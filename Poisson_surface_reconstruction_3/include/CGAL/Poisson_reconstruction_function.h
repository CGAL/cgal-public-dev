// Copyright (c) 2007-09  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Saboret, Pierre Alliez, Tong Zhao, Cédric Portaneri

#ifndef CGAL_IMPLICIT_RECONSTRUCTION_FUNCTION_H
#define CGAL_IMPLICIT_RECONSTRUCTION_FUNCTION_H

#include <CGAL/license/Implicit_surface_reconstruction_3.h>

#include <CGAL/disable_warnings.h>

#ifndef CGAL_DIV_NORMALIZED
#  ifndef CGAL_DIV_NON_NORMALIZED
#    define CGAL_DIV_NON_NORMALIZED 1
#  endif
#endif

#include <vector>
#include <deque>
#include <algorithm>
#include <cmath>
#include <iterator>

#include <CGAL/IO/trace.h>
#include <CGAL/Reconstruction_triangulation_3.h>
#include <CGAL/Covariance_matrix_3.h>
#include <CGAL/spatial_sort.h>
#ifdef CGAL_EIGEN3_ENABLED
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <CGAL/Eigen_solver_traits.h>
#include <Eigen/Eigenvalues>
#else
#endif
#include <CGAL/centroid.h>
#include <CGAL/property_map.h>
#include <CGAL/surface_reconstruction_points_assertions.h>
#include <CGAL/implicit_refine_triangulation.h>
#include <CGAL/Robust_circumcenter_filtered_traits_3.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Timer.h>
#include <CGAL/IO/write_ply_points.h>
#include <CGAL/enum.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Mesh_3/Octree_3.h>

#include <boost/shared_ptr.hpp>
#include <boost/array.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/iterator/indirect_iterator.hpp>
#include <boost/range.hpp>

#include <Spectra/SymEigsSolver.h>
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseCholesky.h>
#include <unsupported/Eigen/SparseExtra>

//Mesh
#include <CGAL/Surface_mesh.h>

/*!
  \file Poisson_reconstruction_function.h
*/

typedef CGAL::cpp11::array<unsigned char, 3>   Color;

namespace CGAL {

  template< class F >
  struct Output_rep< ::Color, F > {
    const ::Color& c;
    static const bool is_specialized = true;
    Output_rep (const ::Color& c) : c(c)
    { }
    std::ostream& operator() (std::ostream& out) const
    {
      if (is_ascii(out))
        out << int(c[0]) << " " << int(c[1]) << " " << int(c[2]);
      else
        out.write(reinterpret_cast<const char*>(&c), sizeof(c));
      return out;
    }
  };

  namespace internal {

    template <class RT>
    bool
    invert(
      const RT& a0,  const RT& a1,  const RT& a2,
      const RT& a3,  const RT& a4,  const RT& a5,
      const RT& a6,  const RT& a7,  const RT& a8,
      RT& i0,   RT& i1,   RT& i2,
      RT& i3,   RT& i4,   RT& i5,
      RT& i6,   RT& i7,   RT& i8)
    {
      // Compute the adjoint.
      i0 = a4*a8 - a5*a7;
      i1 = a2*a7 - a1*a8;
      i2 = a1*a5 - a2*a4;
      i3 = a5*a6 - a3*a8;
      i4 = a0*a8 - a2*a6;
      i5 = a2*a3 - a0*a5;
      i6 = a3*a7 - a4*a6;
      i7 = a1*a6 - a0*a7;
      i8 = a0*a4 - a1*a3;

      RT det = a0*i0 + a1*i3 + a2*i6;

      if(det != 0) {
        RT idet = (RT(1.0))/det;
        i0 *= idet;
        i1 *= idet;
        i2 *= idet;
        i3 *= idet;
        i4 *= idet;
        i5 *= idet;
        i6 *= idet;
        i7 *= idet;
        i8 *= idet;
        return true;
      }

      return false;
    }

  }


/// \cond SKIP_IN_MANUAL
struct Implicit_visitor {
  void before_insertion() const
  {}
};

struct Implicit_skip_vertices {
  double ratio;
  Random& m_random;
  Implicit_skip_vertices(const double ratio, Random& random)
    : ratio(ratio), m_random(random) {}

  template <typename Iterator>
  bool operator()(Iterator) const {
    return m_random.get_double() < ratio;
  }
};

// Given f1 and f2, two sizing fields, that functor wrapper returns
//   max(f1, f2*f2)
// The wrapper stores only pointers to the two functors.
template <typename F1, typename F2>
struct Special_wrapper_of_two_functions_keep_pointers {
  F1 *f1;
  F2 *f2;
  Special_wrapper_of_two_functions_keep_pointers(F1* f1, F2* f2)
    : f1(f1), f2(f2) {}

  template <typename X>
  double operator()(const X& x) const {
    return (std::max)((*f1)(x), CGAL::square((*f2)(x)));
  }

  template <typename X>
  double operator()(const X& x) {
    return (std::max)((*f1)(x), CGAL::square((*f2)(x)));
  }
}; // end struct Special_wrapper_of_two_functions_keep_pointers<F1, F2>

/// \endcond


/*!
\ingroup PkgImplicitSurfaceReconstruction3Ref

\brief Implementation of the Implicit Surface Reconstruction methods.

This class offers 2 algorithms:

1. Poisson Surface Reconstruction

Given a set of 3D points with oriented normals sampled on the boundary
of a 3D solid, the Poisson Surface Reconstruction method \cgalCite{Kazhdan06}
solves for an approximate indicator function of the inferred
solid, whose gradient best matches the input normals. The output
scalar function, represented in an adaptive octree, is then
iso-contoured using an adaptive marching cubes.

We implements a variant of this algorithm which solves for a piecewise
linear function on a 3D Delaunay triangulation instead of an adaptive octree.

2. Spectral Surface Reconstruction

Given a set of 3D points with unoriented normals sampled on the boundary
of a 3D solid, the Spectral Surface Reconstruction Method \cgalCite{cgal:a-vvrup-07}
computes an implicit function by solving a generalized eigenvalue problem
such that its gradient is most aligned with the principal axes of a tensor field.
The principal axes and eccentricities of the tensor field locally represent
respectively the most likely direction of the normal to the surface, and the
confidence in this direction estimation.

The GEP is solved by Spectra library.


\tparam Gt Geometric traits class.
\tparam PointRange Input point range
\tparam NormalMap Normal property map

\cgalModels `ImplicitFunction`

*/
template <class Gt,
          class PointRange,
          class NormalMap>
class Implicit_reconstruction_function
{
// Public types
public:

  /// \name Types
  /// @{
  typedef Gt Geom_traits; ///< Geometric traits class
  /// \cond SKIP_IN_MANUAL
  typedef Reconstruction_triangulation_3<Robust_circumcenter_filtered_traits_3<Gt>, PointRange, NormalMap>
                                                   Triangulation;
  /// \endcond
  typedef typename Triangulation::Cell_handle   Cell_handle;

  // Geometric types
  typedef typename Geom_traits::FT FT; ///< number type.
  typedef typename Geom_traits::Point_3 Point; ///< point type.
  typedef typename Geom_traits::Vector_3 Vector; ///< vector type.
  typedef typename Geom_traits::Sphere_3 Sphere; ///< sphere type

  /// @}

// Private types
private:

  // Internal 3D triangulation, of type Reconstruction_triangulation_3.
  // Note: implicit_refine_triangulation() requires a robust circumcenter computation.

  // Repeat Triangulation types
  typedef typename Triangulation::Triangulation_data_structure Triangulation_data_structure;
  typedef typename Geom_traits::Ray_3 Ray;
  typedef typename Geom_traits::Plane_3 Plane;
  typedef typename Geom_traits::Segment_3 Segment;
  typedef typename Geom_traits::Triangle_3 Triangle;
  typedef typename Geom_traits::Tetrahedron_3 Tetrahedron;
  typedef typename Triangulation::Vertex_handle Vertex_handle;
  typedef typename Triangulation::Cell   Cell;
  typedef typename Triangulation::Vertex Vertex;
  typedef typename Triangulation::Facet  Facet;
  typedef typename Triangulation::Edge   Edge;
  typedef typename Triangulation::Cell_circulator  Cell_circulator;
  typedef typename Triangulation::Facet_circulator Facet_circulator;
  typedef typename Triangulation::Cell_iterator    Cell_iterator;
  typedef typename Triangulation::Facet_iterator   Facet_iterator;
  typedef typename Triangulation::Edge_iterator    Edge_iterator;
  typedef typename Triangulation::Vertex_iterator  Vertex_iterator;
  typedef typename Triangulation::Point_iterator   Point_iterator;
  typedef typename Triangulation::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Triangulation::Finite_cells_iterator    Finite_cells_iterator;
  typedef typename Triangulation::Finite_facets_iterator   Finite_facets_iterator;
  typedef typename Triangulation::Finite_edges_iterator    Finite_edges_iterator;
  typedef typename Triangulation::All_cells_iterator       All_cells_iterator;
  typedef typename Triangulation::Locate_type Locate_type;

  enum Cache_state { UNINITIALIZED, BUSY, INITIALIZED };
  // Thread-safe cache for barycentric coordinates of a cell
  class Cached_bary_coord
  {
  private:
    std::atomic<Cache_state> m_state;
    std::array<double, 9> m_bary;
  public:
    Cached_bary_coord() : m_state (UNINITIALIZED) { }

    // Copy operator to satisfy vector, shouldn't be used
    Cached_bary_coord(const Cached_bary_coord&)
    {
      CGAL_error();
    }

    bool is_initialized()
    {
      Cache_state s = m_state;
      if (s == UNINITIALIZED)
      {
        // If the following line successfully replaces UNINITIALIZED
        // by BUSY, then the current thread in charge of initialization
        if (m_state.compare_exchange_weak(s, BUSY))
          return false;
      }
      // Otherwise, either the thread is BUSY by another thread, or
      // it's already INITIALIZED. Either way, we way until it's INITIALIZED
      else
        while (m_state != INITIALIZED) { }

      // At this point, it's always INITIALIZED
      return true;
    }

    void set_initialized() { m_state = INITIALIZED; }

    const double& operator[] (const std::size_t& idx) const { return m_bary[idx]; }
    double& operator[] (const std::size_t& idx) { return m_bary[idx]; }
  };

  // Wrapper for thread safety of maintained cell hint for fast
  // locate, with conversions atomic<Cell*>/Cell_handle
  class Cell_hint
  {
    std::atomic<Cell*> m_cell;
  public:

    Cell_hint() : m_cell(nullptr) { }

    // Poisson_reconstruction_function should be copyable, although we
    // don't need to copy that
    Cell_hint(const Cell_hint&) : m_cell(nullptr) { }

    Cell_handle get() const
    {
      return Triangulation_data_structure::Cell_range::s_iterator_to(*m_cell);
    }
    void set (Cell_handle ch) { m_cell = ch.operator->(); }
  };

  typedef typename CGAL::Eigen_sparse_matrix<FT>            Matrix;
  typedef typename Eigen::SparseMatrix<FT>                  ESMatrix;
  typedef typename Eigen::Matrix<FT, Eigen::Dynamic, Eigen::Dynamic>  EMatrix;
  typedef typename CGAL::Covariance_matrix_3<Geom_traits>   Covariance;

  typedef typename Spectra::SparseSymMatProd<FT>    OpType;
  typedef typename Spectra::SparseCholesky<FT>      BOpType;
  typedef typename Spectra::SparseSymMatProd<FT>    VopType;

  typedef CGAL::cpp11::array<unsigned char, 3>    Color;
  typedef CGAL::cpp11::tuple<Point, Color>        PC;
  typedef CGAL::Nth_of_tuple_property_map<0, PC>  VF_point_map;
  typedef CGAL::Nth_of_tuple_property_map<1, PC>  VF_color_map;
  typedef std::vector<std::pair<Point, double>>   Point_list;
  typedef std::vector<Color>                      Color_list;

  typedef typename PointRange::const_iterator InputIterator;

  //Mesh
  typedef CGAL::Surface_mesh<Point>                           Mesh;

// Data members.
// Warning: the Surface Mesh Generation package makes copies of implicit functions,
// thus this class must be lightweight and stateless.
private:

  // operator() is pre-computed on vertices of *m_tr by solving
  // ...
  boost::shared_ptr<Triangulation> m_tr;
  mutable std::shared_ptr<std::vector<Cached_bary_coord> > m_bary;

  PointRange m_octree_pwn;
  std::vector<Point> m_octree_steiner;

  mutable std::vector<Point> Dual;
  mutable std::vector<Vector> Normal;

  // contouring and meshing
  Point m_sink; // Point with the minimum value of operator()
  mutable Cell_hint m_hint; // last cell found = hint for next search

  FT m_average_spacing;
  FT m_enlarge_ratio;

  /// function to be used for the different constructors available that are
  /// doing the same thing but with default template parameters
  template <typename PointMap,
            typename Visitor>
  void forward_constructor(
    PointRange& points,
    PointMap point_map,
    NormalMap normal_map,
    Visitor visitor)
  {
    CGAL::Timer task_timer; task_timer.start();
    CGAL_TRACE_STREAM << "Creates Implicit triangulation...\n";

    // Inserts points in triangulation
    m_tr->insert(
      points,
      point_map,
      visitor);

    m_tr->intialize_normal(normal_map);

    // Prints status
    CGAL_TRACE_STREAM << "Creates Implicit triangulation: " << task_timer.time() << " seconds, "
                                                            << std::endl;
  }

// Public methods
public:

  /// \name Creation
  /// @{


  /*!
    Creates a Implicit function from the input range of points.

    \tparam PointMap is a model of `ReadablePropertyMap` with
      a `value_type = Point`.  It can be omitted if `InputIterator`
      `value_type` is convertible to `Point`.

    \param points input point range.
    \param point_map property map: value_type of `InputIterator` -> Point_3.
    \param normal_map property map: value_type of `InputIterator` -> Vector_3.
  */
  Implicit_reconstruction_function()
    : m_tr(new Triangulation)
    , m_bary(new std::vector<Cached_bary_coord>)
  {
  }

  // TOOD: destructors -> delete triangulation and barycenters
  ~Implicit_reconstruction_function()
  {
    m_tr.reset();
    m_bary.reset();
  }


  void reset()
  {
    m_tr.reset();
    m_bary.reset();

    // m_tr = new Triangulation;
    // m_Bary = new std::vector<boost::array<double,9> >;

    Dual.clear();
    Normal.clear();
  }

  template <typename PointMap>
  void initialize_point_map(PointRange& points,
							PointMap point_map,
							NormalMap normal_map,
							bool use_octree,
							bool octree_debug_visu)
  {
	m_average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>
						(points, 6, CGAL::parameters::point_map(point_map));
	Implicit_visitor visitor = Implicit_visitor();
    if(use_octree)
    {
      CGAL_TRACE_STREAM << "init octree...\n";
	  typedef typename OCTREE::Octree<Geom_traits, PointRange, PointMap, NormalMap> Octree;
      Octree octree(points, point_map, normal_map);
	  if(octree_debug_visu)
	  {
		octree.dump_bbox("bbox_scaled_isotrope");
	  }

      CGAL_TRACE_STREAM << "refine octree...\n";
      octree.refine(9, 1);
	  if(octree_debug_visu)
	  {
		// drawing all leafs is same as all nodes but cheaper
        octree.dump_octree("octree_all_nodes", OCTREE::SHOW_ALL_LEAFS);
        octree.dump_octree("octree_non_empty_nodes", OCTREE::SHOW_NON_EMPTY_NODES);
        octree.dump_octree("octree_non_empty_leafs", OCTREE::SHOW_NON_EMPTY_LEAFS);
	  }

      CGAL_TRACE_STREAM << "2:1 octree grading...\n";
      octree.grade();
	  if(octree_debug_visu)
	  {
        octree.dump_octree("balanced_octree_all_nodes", OCTREE::SHOW_ALL_LEAFS);
        octree.dump_octree("balanced_octree_non_empty_nodes", OCTREE::SHOW_NON_EMPTY_NODES);
        octree.dump_octree("balanced_octree_non_empty_leafs", OCTREE::SHOW_NON_EMPTY_LEAFS);
        if(octree.debug_grading())
          CGAL_TRACE_STREAM << " octree correctly balanced!\n";
				else
          CGAL_TRACE_STREAM << " Error: octree not correctly balanced!\n";
	  }

      CGAL_TRACE_STREAM << "generate octree new points with normal...\n";
      octree.generate_points(std::back_inserter(m_octree_pwn), std::back_inserter(m_octree_steiner));
	  if(octree_debug_visu)
	  {
		octree.dump_octree_pwn("octree_pwn", m_octree_pwn);
		octree.dump_octree_point("octree_steiner", m_octree_steiner);
	  }

      CGAL_TRACE_STREAM << "creates implicit triangulation...\n";
      m_tr->insert(m_octree_pwn, point_map, visitor);
      m_tr->intialize_normal(normal_map);
	  for (auto steiner = m_octree_steiner.begin(); steiner != m_octree_steiner.end(); steiner++)
		m_tr->insert(*steiner, Triangulation::Point_type::STEINER, Cell_handle(), visitor);
    }
	else
	{
	  forward_constructor(points, PointMap(), NormalMap(), visitor);
    first_delaunay_refinement(visitor);
	}


  }

  /// \endcond

  /// @}

  /// \name Operations
  /// @{

  /// Returns a sphere bounding the inferred surface.
  Sphere bounding_sphere() const
  {
    return m_tr->bounding_sphere();
  }

  /// \cond SKIP_IN_MANUAL
  const Triangulation& tr() const {
    return *m_tr;
  }

  void initialize_insides() const
  {
    Vertex_handle v, e;

    std::list<Cell_handle> cells;
    typename std::list<Cell_handle>::iterator it;

    for (v = m_tr->finite_vertices_begin(), e = m_tr->finite_vertices_end(); v!= e; ++v){
      cells.clear();
      m_tr->incident_cells(v, std::back_inserter(cells));

      bool flag = true;
      for(it = cells.begin(); it != cells.end(); it++)
      {
        Cell_handle cell = *it;
        if(m_tr->is_infinite(cell)){
          flag = false; break;
        }
      }

      v->position() = flag? static_cast<unsigned char>(Triangulation::INSIDE):
                            static_cast<unsigned char>(Triangulation::BOUNDARY);
    }
  }


  template <class Visitor>
  bool first_delaunay_refinement(Visitor visitor,
                                 const unsigned int max_vertices = (unsigned int)1000000000,
                                 const FT enlarge_ratio = 1.5)
  {
    CGAL::Timer refine_timer;
    CGAL_TRACE_STREAM << "Delaunay refinement...\n";

    // Delaunay refinement
    //const FT radius_edge_ratio_bound = std::max(1.1, -3.32317 * m_average_spacing + 1.4128);
    const FT radius_edge_ratio_bound = 1.3;
    const FT radius = sqrt(bounding_sphere().squared_radius()); // get triangulation's radius
    const FT cell_radius_bound = radius / 5.; // large
    m_enlarge_ratio = enlarge_ratio;

    internal::Implicit::Constant_sizing_field<Triangulation> sizing_field(CGAL::square(cell_radius_bound));
    std::vector<int> NB;

    NB.push_back( delaunay_refinement(radius_edge_ratio_bound, sizing_field, max_vertices, enlarge_ratio));

    while(m_tr->insert_fraction(visitor)) {
      NB.push_back( delaunay_refinement(radius_edge_ratio_bound, sizing_field, max_vertices, enlarge_ratio));
    }

    // Prints status
    CGAL_TRACE_STREAM << "Delaunay refinement: " << "added ";
    for(std::size_t i = 0; i < NB.size()-1; i++){
      CGAL_TRACE_STREAM << NB[i] << " + ";
    }
    CGAL_TRACE_STREAM << NB.back() << " Steiner points, "
                      << refine_timer.time() << " seconds, "
                      << std::endl;

    return true;
  }


  /// split bad tetrahedrons into 2 in order to calculate the normal derivative properly
  unsigned int break_bad_tets_on_boundary()
  {
    int counter = 0;
    std::vector<Cell_handle> cells;

    for(Finite_vertices_iterator v = m_tr->finite_vertices_begin(); v != m_tr->finite_vertices_end(); ++v)
    {
      if(v->position() == Triangulation::BOUNDARY){
        cells.clear();
        m_tr->incident_cells(v, std::back_inserter(cells));

        // iterate over all incident cells
        typename std::vector<Cell_handle>::iterator cellit;
        for(cellit = cells.begin(); cellit != cells.end(); cellit++){
          Cell_handle cell = *cellit;
          // skip infinite cells
          if(!m_tr->is_cell(cell) || m_tr->is_infinite(cell))
            continue;
          // find cells whose 4 vertices are on boundary
          bool flag = true;
          for(int i = 0; i < 4; i++)
            if(cell->vertex(i)->position() != Triangulation::BOUNDARY)
            {
              flag = false; break;
            }
          if(!flag) continue;

          // find the edge to break
          int index_i = -1;
          int index_j = -1;
          for(int j = 0; j < 4; j++)
          {
            Facet facet_j = std::make_pair(cell, j);
            Facet facet_jo = m_tr->mirror_facet(facet_j);

            if(m_tr->is_infinite((facet_jo.first)->vertex(facet_jo.second))){
              if(index_i == -1) index_i = j;
              else index_j = j;
            }
          }

          // break edge
          if(index_i != -1 && index_j != -1 && index_i != index_j){
            m_tr->insert_in_edge(cell, index_i, index_j);
            counter++;
          }
        }
      }
    }

    return counter;
  }


  // Poisson surface reconstruction
  // This variant requires all parameters.
  template <class SparseLinearAlgebraTraits_d,
            class Visitor>
  bool compute_poisson_implicit_function(
                                 SparseLinearAlgebraTraits_d solver, // = SparseLinearAlgebraTraits_d(),
                                 Visitor visitor,
                                 double average_spacing_ratio = 5.0) // this parameter should be passed to second delaunay refinement
  {
    CGAL::Timer task_timer; task_timer.start();

#ifdef CGAL_DIV_NON_NORMALIZED
    CGAL_TRACE_STREAM << "Solve Poisson equation with non-normalized divergence...\n";
#else
    CGAL_TRACE_STREAM << "Solve Poisson equation with normalized divergence...\n";
#endif

    // Computes the Poisson indicator function operator()
    // at each vertex of the triangulation.
    double lambda = 0.1;
    if ( ! solve_poisson(solver, lambda) )
    {
      std::cerr << "Error: cannot solve Poisson equation" << std::endl;
      return false;
    }

    // Shift and orient operator() such that:
    // - operator() = 0 on the input points,
    // - operator() < 0 inside the surface.
    set_contouring_value(median_value_at_input_vertices());

    // Prints status
    CGAL_TRACE_STREAM << "Solve Poisson equation: " << task_timer.time() << " seconds, "
                                                    << std::endl;
    task_timer.reset();

    return true;
  }
  /// \endcond

  /*!
    This function must be called after the
    insertion of oriented points. It computes the piecewise linear scalar
    function operator() by: applying Delaunay refinement, solving for
    operator() at each vertex of the triangulation with a sparse linear
    solver, and shifting and orienting operator() such that it is 0 at all
    input points and negative inside the inferred surface.

    \tparam SparseLinearAlgebraTraits_d Symmetric definite positive sparse linear solver.
    If \ref thirdpartyEigen "Eigen" 3.1 (or greater) is available and `CGAL_EIGEN3_ENABLED`
    is defined, an overload with \link Eigen_solver_traits <tt>Eigen_solver_traits<Eigen::ConjugateGradient<Eigen_sparse_symmetric_matrix<double>::EigenType> ></tt> \endlink
    as default solver is provided.

    \param solver sparse linear solver.
    \param smoother_hole_filling controls if the Delaunay refinement is done for the input points, or for an approximation of the surface obtained from a first pass of the algorithm on a sample of the points.

    \return `false` if the linear solver fails.
  */
  template <class SparseLinearAlgebraTraits_d>
  bool compute_poisson_implicit_function(SparseLinearAlgebraTraits_d solver, bool smoother_hole_filling = false)
  {
    if (smoother_hole_filling)
      return compute_poisson_implicit_function<SparseLinearAlgebraTraits_d,Implicit_visitor>(solver, Implicit_visitor(), 5.);
    else
      return compute_poisson_implicit_function<SparseLinearAlgebraTraits_d,Implicit_visitor>(solver, Implicit_visitor());
  }

  /// \cond SKIP_IN_MANUAL
#ifdef CGAL_EIGEN3_ENABLED
  // This variant provides the default sparse linear traits class = Eigen_solver_traits.
  bool compute_poisson_implicit_function(bool smoother_hole_filling = false)
  {
    typedef Eigen_solver_traits<Eigen::ConjugateGradient<Eigen_sparse_symmetric_matrix<double>::EigenType> > Solver;
    return compute_poisson_implicit_function<Solver>(Solver(), smoother_hole_filling);
  }
#endif


  // Spectral Surface Reconstruction with K = L + N
  // This variant requires all parameters.
  template <class Visitor,
            class ReliabilityMap,
            class ConfidenceMap>
  bool compute_spectral_implicit_function(
                                 ReliabilityMap  reliability_map,
                                 ConfidenceMap  confidence_map,
                                 Visitor    visitor,
                                 double bilaplacian = 1.0,
                                 double laplacian = 0.0, // this parameter is dangerous
                                 double average_spacing_ratio = 5.0) // pass to second_delaunay_refinement
  {
    if(laplacian < 0.0)
	{
      initialize_insides();
      unsigned int bad_tets = break_bad_tets_on_boundary();
      CGAL_TRACE_STREAM << "Break " << bad_tets << " bad tets!" << std::endl;
    }

    CGAL::Timer task_timer; task_timer.start();

    // Computes the Implicit indicator function operator()
    // at each vertex of the triangulation.
    if ( ! solve_spectral(reliability_map, confidence_map, bilaplacian, laplacian) )
    {
      std::cerr << "Error: cannot solve Implicit equation" << std::endl;
      return false;
    }

    // Shift and orient operator() such that:
    // - operator() = 0 on the input points,
    // - operator() < 0 inside the surface.
    set_contouring_value(median_value_at_input_vertices());

    // Prints status
    CGAL_TRACE_STREAM << "Solve Spectral equation: " << task_timer.time() << " seconds, "
                                                    << std::endl;
    task_timer.reset();

    return true;
  }


  template <class ReliabilityMap>
  bool compute_spectral_implicit_function(ReliabilityMap reliability_map,
                                          FT confidence = 15.,
                                          double bilaplacian = 1,
                                          double laplacian = 0,
                                          bool smoother_hole_filling = false)
  {
    typedef typename CGAL::Constant_property_map<InputIterator, FT> CoefficientMap;
    CoefficientMap confidence_map = CGAL::Constant_property_map<InputIterator, FT>(FT(confidence));

    if (smoother_hole_filling)
      return compute_spectral_implicit_function<Implicit_visitor>(reliability_map, confidence_map, Implicit_visitor(), bilaplacian, laplacian, 5.0);
    else
      return compute_spectral_implicit_function<Implicit_visitor>(reliability_map, confidence_map, Implicit_visitor(), bilaplacian, laplacian);
  }


  template <class ConfidenceMap>
  bool compute_spectral_implicit_function(FT reliability,
                                          ConfidenceMap confidence_map,
                                          double bilaplacian = 1.0,
                                          double laplacian = 1.0,
                                          bool smoother_hole_filling = false)
  {
    typedef typename CGAL::Constant_property_map<InputIterator, FT> CoefficientMap;
    CoefficientMap reliability_map = CGAL::Constant_property_map<InputIterator, FT>(FT(reliability));

    if (smoother_hole_filling)
      return compute_spectral_implicit_function<Implicit_visitor>(reliability_map, confidence_map, Implicit_visitor(), bilaplacian, laplacian, 5);
    else
      return compute_spectral_implicit_function<Implicit_visitor>(reliability_map, confidence_map, Implicit_visitor(), bilaplacian, laplacian);
  }
  /// \endcond

  /*!
    This function must be called after the
    insertion of unoriented points. It computes the piecewise linear scalar
    function operator() by: applying Delaunay refinement, solving for
    operator() at each vertex of the triangulation with a sparse linear
    solver, and shifting and orienting operator() such that it is 0 at all
    input points and negative inside the inferred surface.

    \tparam ReliabilityMap property map: value_type of `InputIterator` -> FT.
    \tparam ConfidenceMap property map: value_type of `InputIterator` -> FT.

    \param reliability_map a property map of confidence coefficients for point coordinates
    \param confidence_map a property map of confidence coefficients for point normals
    \param bilaplacian bilaplacian term weight
    \param laplacian laplacian term weight
    \param smoother_hole_filling controls if the Delaunay refinement is done for the input points, or for an approximation of the surface obtained from a first pass of the algorithm on a sample of the points.

    \return `false` if the solver fails.
  */
  template <class ReliabilityMap,
            class ConfidenceMap>
  bool compute_spectral_implicit_function(ReliabilityMap  reliability_map,
                                          ConfidenceMap  confidence_map,
                                          double bilaplacian = 1.0,
                                          double laplacian = 1.0,
                                          bool smoother_hole_filling = false)
  {
    if (smoother_hole_filling)
      return compute_spectral_implicit_function<Implicit_visitor>(reliability_map, confidence_map, Implicit_visitor(), bilaplacian, laplacian, 5);
    else
      return compute_spectral_implicit_function<Implicit_visitor>(reliability_map, confidence_map, Implicit_visitor(), bilaplacian, laplacian);
  }

  /*!
    This function provides an overloaded function so that users can give a global reliability coefficient
    and a global confidence coefficient. Users can also provide a property map for one of them and a global
    coefficient for the other.

    \param reliability the confidence coefficient for point coordinates
    \param confidence the confidence coefficient for point normals
    \param bilaplacian bilaplacian term weight
    \param laplacian laplacian term weight
    \param smoother_hole_filling controls if the Delaunay refinement is done for the input points, or for an approximation of the surface obtained from a first pass of the algorithm on a sample of the points.

    \return `false` if the solver fails.
  */
  bool compute_spectral_implicit_function(FT reliability = 10.,
                                          FT confidence = 10.,
                                          double bilaplacian = 1,
                                          double laplacian = 0.1,
                                          bool smoother_hole_filling = false)
  {
    typedef typename CGAL::Constant_property_map<InputIterator, FT> CoefficientMap;
    CoefficientMap reliability_map = CGAL::Constant_property_map<InputIterator, FT>(FT(reliability));
    CoefficientMap confidence_map = CGAL::Constant_property_map<InputIterator, FT>(FT(confidence));

    if (smoother_hole_filling)
      return compute_spectral_implicit_function<Implicit_visitor>(reliability_map, confidence_map, Implicit_visitor(), bilaplacian, laplacian, 5.0);
    else
      return compute_spectral_implicit_function<Implicit_visitor>(reliability_map, confidence_map, Implicit_visitor(), bilaplacian, laplacian);
  }



  /// \endcond

  /*!
    `ImplicitFunction` interface: evaluates the implicit function at a
    given 3D query point. The function `compute_implicit_function()` must be
    called before the first call to `operator()`.
  */
  FT operator()(const Point& p) const
  {
    Cell_handle hint = m_hint.get();
    hint = m_tr->locate(p, hint);
    m_hint.set(hint);

    if(m_tr->is_infinite(hint)) {
      int i = hint->index(m_tr->infinite_vertex());
      return hint->vertex((i+1)&3)->f();
    }

    FT a,b,c,d;
    barycentric_coordinates(p,hint,a,b,c,d);
    return a * hint->vertex(0)->f() +
           b * hint->vertex(1)->f() +
           c * hint->vertex(2)->f() +
           d * hint->vertex(3)->f();
  }

  /// \cond SKIP_IN_MANUAL
  boost::tuple<FT, Cell_handle, bool> special_func(const Point& p) const
  {
    m_hint = m_tr->locate(p, m_hint); // no hint when we use hierarchy

    if(m_tr->is_infinite(m_hint)) {
      int i = m_hint->index(m_tr->infinite_vertex());
      return boost::make_tuple(m_hint->vertex((i + 1) & 3)->f(),
                               m_hint, true);
    }

    FT a, b, c, d;
    barycentric_coordinates(p, m_hint, a, b, c, d);
    return boost::make_tuple(a * m_hint->vertex(0)->f() +
                             b * m_hint->vertex(1)->f() +
                             c * m_hint->vertex(2)->f() +
                             d * m_hint->vertex(3)->f(),
                             m_hint, false);
  }

  void initialize_cell_indices()
  {
    int i = 0;
    for(Finite_cells_iterator fcit = m_tr->finite_cells_begin();
        fcit != m_tr->finite_cells_end();
        ++fcit){
      fcit->info() = i++;
    }
  }

  void initialize_barycenters() const
  {
    m_bary->clear();
    m_bary->resize(m_tr->number_of_cells());
  }

  void initialize_cell_normals() const
  {
    Normal.resize(m_tr->number_of_finite_cells());
    int i = 0;
    int N = 0;
    for(Finite_cells_iterator fcit = m_tr->finite_cells_begin();
        fcit != m_tr->finite_cells_end();
        ++fcit){
      Normal[i] = cell_normal(fcit);
      if(Normal[i] == NULL_VECTOR){
        N++;
      }
      ++i;
    }
    std::cerr << N << " out of " << i << " cells have NULL_VECTOR as normal" << std::endl;
  }

  void initialize_duals() const
  {
    Dual.resize(m_tr->number_of_finite_cells());
    int i = 0;
    for(Finite_cells_iterator fcit = m_tr->finite_cells_begin();
        fcit != m_tr->finite_cells_end();
        ++fcit){
      Dual[i++] = m_tr->dual(fcit);
    }
  }

  void clear_duals() const
  {
    Dual.clear();
  }

  void clear_normals() const
  {
    Normal.clear();
  }

  void initialize_matrix_entry(Cell_handle ch) const
  {
    Cached_bary_coord& bary = (*m_bary)[ch->info()];

    if (bary.is_initialized())
      return;

    // If the cache was uninitialized, this thread is in charge of
    // initializing it
    const Point& pa = ch->vertex(0)->point();
    const Point& pb = ch->vertex(1)->point();
    const Point& pc = ch->vertex(2)->point();
    const Point& pd = ch->vertex(3)->point();

    Vector va = pa - pd;
    Vector vb = pb - pd;
    Vector vc = pc - pd;

    internal::invert(va.x(), va.y(), va.z(),
                     vb.x(), vb.y(), vb.z(),
                     vc.x(), vc.y(), vc.z(),
                     bary[0],bary[1],bary[2],
                     bary[3],bary[4],bary[5],
                     bary[6],bary[7],bary[8]);

    bary.set_initialized();
  }
  /// \endcond

  /// Returns a point located inside the inferred surface.
  Point get_inner_point() const
  {
    // Gets point / the implicit function is minimum
    return m_sink;
  }

  /// @}

// Private methods:
private:

  /// Delaunay refinement (break bad tetrahedra, where
  /// bad means badly shaped or too big). The normal of
  /// Steiner points is set to zero.
  /// Returns the number of vertices inserted.

  template <typename Sizing_field>
  unsigned int delaunay_refinement(FT radius_edge_ratio_bound, ///< radius edge ratio bound (ignored if zero)
                                   Sizing_field sizing_field, ///< cell radius bound (ignored if zero)
                                   unsigned int max_vertices, ///< number of vertices bound
                                   FT enlarge_ratio) ///< bounding box enlarge ratio
  {
    return delaunay_refinement(radius_edge_ratio_bound,
                               sizing_field,
                               max_vertices,
                               enlarge_ratio,
                               internal::Implicit::Constant_sizing_field<Triangulation>());
  }

  template <typename Sizing_field,
            typename Second_sizing_field>
  unsigned int delaunay_refinement(FT radius_edge_ratio_bound, ///< radius edge ratio bound (ignored if zero)
                                   Sizing_field sizing_field, ///< cell radius bound (ignored if zero)
                                   unsigned int max_vertices, ///< number of vertices bound
                                   FT enlarge_ratio, ///< bounding box enlarge ratio
                                   Second_sizing_field second_sizing_field)
  {
    Sphere elarged_bsphere = enlarged_bounding_sphere(enlarge_ratio);
    unsigned int nb_vertices_added = implicit_refine_triangulation(*m_tr,radius_edge_ratio_bound,sizing_field,second_sizing_field,max_vertices,elarged_bsphere);
    return nb_vertices_added;
  }

  /// Computes enlarged geometric bounding sphere of the embedded triangulation.
  Sphere enlarged_bounding_sphere(FT ratio) const
  {
    Sphere bsphere = bounding_sphere(); // triangulation's bounding sphere
    return Sphere(bsphere.center(), bsphere.squared_radius() * ratio * ratio);
  }


  /// Poisson Surface Reconstruction.
  /// Returns false on error.
  ///
  /// @commentheading Template parameters:
  /// @param SparseLinearAlgebraTraits_d Symmetric definite positive sparse linear solver.
  template <class SparseLinearAlgebraTraits_d>
  bool solve_poisson(
    SparseLinearAlgebraTraits_d solver, ///< sparse linear solver
    double lambda)
  {
    CGAL_TRACE("Calls solve_poisson()\n");

    double time_init = clock();

    double duration_assembly = 0.0;
    double duration_solve = 0.0;

    initialize_cell_indices();
    initialize_barycenters();

    // get #variables
    constrain_one_vertex_on_convex_hull();
    m_tr->index_all_vertices(true); // index all unconstrained vertices
    unsigned int nb_variables = static_cast<unsigned int>(m_tr->number_of_vertices()-1);

    CGAL_TRACE("  Number of variables: %ld\n", (long)(nb_variables));

    // Assemble linear system A*X=B
    typename SparseLinearAlgebraTraits_d::Matrix A(nb_variables); // matrix is symmetric definite positive
    typename SparseLinearAlgebraTraits_d::Vector X(nb_variables), B(nb_variables);

    initialize_duals();
#ifndef CGAL_DIV_NON_NORMALIZED
    initialize_cell_normals();
#endif
    Finite_vertices_iterator v, e;
    for(v = m_tr->finite_vertices_begin(),
        e = m_tr->finite_vertices_end();
        v != e;
        ++v)
    {
      if(!m_tr->is_constrained(v)) {
#ifdef CGAL_DIV_NON_NORMALIZED
        B[v->index()] = div(v); // rhs -> divergent
#else // not defined(CGAL_DIV_NORMALIZED)
        B[v->index()] = div_normalized(v); // rhs -> divergent
#endif // not defined(CGAL_DIV_NORMALIZED)
        assemble_poisson_row<SparseLinearAlgebraTraits_d>(A, v, B, lambda);
      }
    }

    clear_duals();
    clear_normals();
    duration_assembly = (clock() - time_init)/CLOCKS_PER_SEC;
    CGAL_TRACE("  Creates matrix: done (%.2lf s)\n", duration_assembly);

    CGAL_TRACE("  Solve sparse linear system...\n");

    // Solve "A*X = B". On success, solution is (1/D) * X.
    time_init = clock();
    double D;
    if(!solver.linear_solver(A, B, X, D))
      return false;
    CGAL_surface_reconstruction_points_assertion(D == 1.0);
    duration_solve = (clock() - time_init)/CLOCKS_PER_SEC;

    CGAL_TRACE("  Solve sparse linear system: done (%.2lf s)\n", duration_solve);

    // copy function's values to vertices
    unsigned int index = 0;
    for (v = m_tr->finite_vertices_begin(), e = m_tr->finite_vertices_end(); v!= e; ++v)
      if(!m_tr->is_constrained(v))
        v->f() = X[index++];

    CGAL_TRACE("End of solve_poisson()\n");

    return true;
  }

  /// Spectral Surface reconstruction.
  /// Returns false on error.
  ///
  /// @commentheading Template parameters:
  template <class ReliabilityMap,
            class ConfidenceMap>
  bool solve_spectral(ReliabilityMap reliability_map,
                      ConfidenceMap confidence_map,
                      double bilaplacian,
                      double laplacian)
  {
    CGAL_TRACE("Calls solve_spectral()\n");

    double time_init = clock();

    double duration_assembly = 0.0;
    double duration_solve = 0.0;

    initialize_cell_indices();
    initialize_barycenters();
    initialize_insides();

    // get #variables
    m_tr->index_all_vertices();
    const int nb_variables = static_cast<int>(m_tr->number_of_vertices());
    const int nb_input_vertices = m_tr->nb_input_vertices();
    //const int nb_insides = static_cast<int>(m_tr->nb_inside_vertices());
  	CGAL_TRACE("  %d input vertices out of %d\n", nb_input_vertices, nb_variables);

    //std::ofstream oFileT("triangulation.off", std::ios::out);
    //oFileT << *m_tr;

    //check_ratio_radius_edge();

    // Assemble matrices
    Matrix AA(nb_variables), L(nb_variables), F(nb_variables); // matrix is symmetric definite positive
    Matrix V_inv(nb_variables), N(nb_variables), V(nb_variables);

    ESMatrix B(nb_variables, nb_variables);
    EMatrix X(nb_variables, 1);

    bool flag_boundary = (laplacian < 0) ? true : false;

    initialize_duals();

    CGAL_TRACE("  Begin calculation: (%.2lf s)\n", (clock() - time_init) / CLOCKS_PER_SEC);
    Finite_vertices_iterator v, e;
    double duration_cal = 0., duration_assign= 0.;
    for(v = m_tr->finite_vertices_begin(), e = m_tr->finite_vertices_end();
        v != e;
        ++v)
    {
        FT fitting = (v->type() == Triangulation::INPUT) ? get(reliability_map, v->input_iterator()) : 0;
        assemble_spectral_row(v, AA, L, F, V_inv, V, N, fitting, confidence_map, duration_assign, duration_cal, flag_boundary);
    }

    CGAL_TRACE("  Calculate elem: total (%.2lf s)\n", duration_cal / CLOCKS_PER_SEC);
    CGAL_TRACE("  Assign: total (%.2lf s)\n", duration_assign / CLOCKS_PER_SEC);

    //double time_b = clock();

    ESMatrix EL = L.eigen_object(), EA = AA.eigen_object();
    ESMatrix EV_inv = V_inv.eigen_object(), EV = V.eigen_object();

    const FT radius = sqrt(bounding_sphere().squared_radius()); // get triangulation's radius

    EV_inv = EV_inv * ::pow(radius, 3);
    EV = EV / ::pow(radius, 3);

    if(flag_boundary){
      EL = EL + N.eigen_object();
      EL = EL / radius;
      B = EL.transpose() * EV_inv * EL * bilaplacian + EV * F.eigen_object();
    }
    else{
      EL = EL / radius;
      //B = EL * EV_inv * EL * bilaplacian + EL * laplacian + F.eigen_object();
      //B = EL * EV_inv * EL * bilaplacian + EL * laplacian + EV * F.eigen_object();
      //B = EL.transpose() * EL * bilaplacian + EL * laplacian + F.eigen_object();
      B = EL * EL * bilaplacian + EL * laplacian + F.eigen_object();
    }
/*
    dihedral_angle_per_cell();
    save_steiner_point();

    std::cerr << "Check A...." << std::endl;
    check_spd(EA);

    std::cerr << "Check L...." << std::endl;
    check_spd(EL);

    std::cerr << "write l..." << std::endl;
    std::string lvec_outfile("lvec.ply");
    save_smallest_eigvec(EL, lvec_outfile, nb_variables, 1e-7);

    std::cerr << "write a..." << std::endl;
    std::string avec_outfile("avec.ply");
    save_smallest_eigvec(EA, avec_outfile, nb_variables, 1e-8);

    ESMatrix identity(nb_variables, nb_variables);
    identity.setIdentity();
    identity = identity * 1e-8;

    std::cerr << "write l*l..." << std::endl;
    ESMatrix ELL = EL * EL + identity;
    std::string llvec_outfile("llvec.ply");
    save_smallest_eigvec(ELL, llvec_outfile, nb_variables, 1e-4);

    std::cerr << "write l*v*l..." << std::endl;
    ESMatrix ELVL = EL * EV_inv * EL + identity;
    std::string lvlvec_outfile("lvlvec.ply");
    save_smallest_eigvec(ELVL, lvlvec_outfile, nb_variables, 1e-4);
*/
    clear_duals();
    duration_assembly = (clock() - time_init) / CLOCKS_PER_SEC;
    CGAL_TRACE("  Creates matrix: done (%.2lf s)\n", duration_assembly);
    CGAL_TRACE("  Solve generalized eigenvalue problem...\n");

    // Solve generalized eigenvalue problem
    time_init = clock();
    spectral_solver<ESMatrix, EMatrix, Spectra::LARGEST_ALGE>(EA, B, EL, X);

    duration_solve = (clock() - time_init) / CLOCKS_PER_SEC;

    CGAL_TRACE("  Solve generalized eigenvalue problem: done (%.2lf s)\n", duration_solve);

    /*
    EMatrix AX = EA * X;
    EMatrix BX = B * X;
    EMatrix LX = EL * X;*/

    // copy function's values to vertices
    unsigned int index = 0;
    for (v = m_tr->finite_vertices_begin(), e = m_tr->finite_vertices_end(); v!= e; ++v){
      v->f() = X(index, 0);
      //v->af() = AX(index, 0);
      //v->bf() = BX(index, 0);
      //v->lf() = LX(index, 0);
      index += 1;
    }

    CGAL_TRACE("End of solve_spectral()\n");

    return true;
  }


  /// @commentheading Template parameters:
  /// @param MatType The name of the matrix operation class for A and B
  /// @param RMatType The name of the matrix operation class for X
  /// @param SelectionRule An enumeration value indicating the selection rule of the requested eigenvalues
  template <typename MatType, typename RMatType, int SelectionRule>
  void spectral_solver(const MatType& A, const MatType& B, const MatType& L, RMatType& X, int k = 1, int m = 37)
  {
      CGAL_TRACE("Begin solving spectra...\n");
      OpType op(A);
      BOpType Bop(B);
      // Make sure B is positive definite and the decompoition is successful
      assert(Bop.info() == Spectra::SUCCESSFUL);

      Spectra::SymGEigsSolver<FT, SelectionRule, OpType, BOpType, Spectra::GEIGS_CHOLESKY> eigs(&op, &Bop, k, m);
      eigs.init();
      //int nconv = eigs.compute();
	  eigs.compute();

      CGAL_TRACE("Problem solved!\n");

      if(eigs.info() != Spectra::SUCCESSFUL)
        CGAL_TRACE("  Spectra failed! %d\n", eigs.info());

      X = eigs.eigenvectors();
  }


  template <typename MatType>
  void save_smallest_eigvec(const MatType& A, const std::string outfile, const int nb_variables, const double tol = 1e-8, int k = 5, int m = 37, int n = 5)
  {
      CGAL_TRACE("Begin finding smallest eigenvector...\n");
      VopType op(A);

      Spectra::SymEigsSolver<FT, Spectra::SMALLEST_ALGE, VopType> eigs(&op, k, m);
      eigs.init();
      int nconv = eigs.compute(1000, tol);

      Eigen::VectorXcd evalues;
      //EMatrix evecs(nb_variables, k);

      if(eigs.info() == Spectra::SUCCESSFUL)
      {
        evalues = eigs.eigenvalues();
        auto evecs = eigs.eigenvectors();
        for(int i = 0; i < k; i++){
          if(std::abs(evalues[k - i - 1]) > 1e-6){
            std::cerr << evalues[k - i - 1] << std::endl;
            unsigned int index = 0;
            for (Finite_vertices_iterator v = m_tr->finite_vertices_begin(); v!= m_tr->finite_vertices_end(); ++v)
              //std::cerr << evecs(index++, k - i - 1) << std::endl;
              v->check() = evecs(index++, k - i - 1);

            Point my_pmax = draw_xslice_function(500, 0, 4, outfile);
            my_pmax = draw_xslice_function(500, 1, 4, "1_" + outfile);
            my_pmax = draw_xslice_function(500, 2, 4, "2_" + outfile);
            my_pmax = draw_xslice_function(500, 3, 4, "3_" + outfile);

            /*
            typename Triangulation::Locate_type lt;
            int li, lj;
            Cell_handle my_cmax = m_tr -> locate(my_pmax, lt, li, lj);

            std::ofstream out("max_tet.off");
            out << "OFF" << std::endl;
            out << "4 4 0" << std::endl;
            for(int m = 0; m < 4; m++){
              Point pt = my_cmax->vertex(m)->point();
              out << pt.x() << " " << pt.y() << " " << pt.z() << std::endl;
            }
            out << "3 0 1 2" << std::endl;
            out << "3 0 1 3" << std::endl;
            out << "3 0 2 3" << std::endl;
            out << "3 1 2 3" << std::endl;
            out.close();
            */


            /*
            auto chosen_evecs = evecs.col(k - i - 1);
            double min_value = chosen_evecs.minCoeff();

            for(int m = 0; m < n; m++){
              int idx;
              double max_value = chosen_evecs.maxCoeff(&idx);
              Finite_vertices_iterator max_v = m_tr->finite_vertices_begin();

              for(int j = 0; j < idx; j++) max_v++;

              std::vector<Cell_handle> max_cells;
              m_tr->incident_cells(max_v, std::back_inserter(max_cells));
              int num_cells = 0;

              for(auto m_cell = max_cells.begin(); m_cell != max_cells.end(); m_cell++){
                if(!m_tr->is_infinite(*m_cell))
                  num_cells++;
              }

              std::ofstream out(std::to_string(m) + "_tet.off");
              out << "OFF" << std::endl;
              out << std::to_string(num_cells * 4) << " " << std::to_string(num_cells) << " " << std::to_string(0) << std::endl << std::endl;

              for(auto m_cell = max_cells.begin(); m_cell != max_cells.end(); m_cell++){
                if(!m_tr->is_infinite(*m_cell))
                  for(int j = 0; j < 4; j++){
                    Point curr_pt = (*m_cell)->vertex(j)->point();
                    out << std::to_string(curr_pt.x()) << " " << std::to_string(curr_pt.y()) << " " << std::to_string(curr_pt.z()) << std::endl;
                  }
              }

              int my_idx = 0;

              for(int j = 0; j < max_cells.size(); j++){
                if(!m_tr->is_infinite(max_cells[j])){
                  out << std::to_string(4) << " " << std::to_string(my_idx * 4) << " " << std::to_string(my_idx * 4 + 1) << " " << std::to_string(my_idx * 4 + 2) << " " << std::to_string(my_idx * 4 + 3) << std::endl;
                  my_idx++;
                }
              }
              out.close();
              evecs(idx) = min_value;
            }
            */
            break;
          }
        }

        // int list_idx[n];


      }
      else
        std::cerr << "Check failed! " << eigs.info() << std::endl;

  }

  void dihedral_angle_per_cell()
  {
    std::ofstream out("cotan.txt");

    for(Finite_cells_iterator cb = m_tr->finite_cells_begin(); cb != m_tr->finite_cells_end(); cb++)
    {
      //double min_cotan = 1e7;
      Point center = bounding_sphere().center();

      for(int i = 0; i < 3; i++)
        for(int j = i + 1; j < 4; j++)
        {
          double cotan = cotan_per_edge(cb, i, j);
          //if(cotan < min_cotan) min_cotan = cotan;
          out << std::to_string(cotan) << " ";
        }

      Vector point_from_center = cb->vertex(0)->point() - center;
      FT length_pc = std::sqrt(point_from_center * point_from_center);
      out << std::to_string(length_pc) << std::endl;
    }

    out.close();
  }

  void check_ratio_radius_edge()
  {
    std::ofstream out("ratio.txt");
    Point center = bounding_sphere().center();

    for(Finite_cells_iterator cb = m_tr->finite_cells_begin(); cb != m_tr->finite_cells_end(); cb++){
      Point circum = CGAL::circumcenter(cb->vertex(0)->point(),
                                        cb->vertex(1)->point(),
                                        cb->vertex(2)->point(),
                                        cb->vertex(3)->point());
      std::cerr << circum << std::endl;
      double radius = std::sqrt((cb->vertex(0)->point() - circum) * (cb->vertex(0)->point() - circum));

      double short_edge = std::sqrt((cb->vertex(0)->point() - cb->vertex(1)->point()) * (cb->vertex(0)->point() - cb->vertex(1)->point()));
      double min_cotan = 1e7;

      for(int i = 0; i < 3; i++)
        for(int j = i + 1; j < 4; j++)
        {
          double cotan = cotan_per_edge(cb, i, j);
          if(cotan < min_cotan) min_cotan = cotan;

          Vector edge = cb->vertex(i)->point() - cb->vertex(j)->point();
          double len_edge = std::sqrt(edge * edge);
          if(len_edge < short_edge) short_edge = len_edge;
        }

      out << std::to_string(radius / short_edge) << " " << std::to_string(std::sqrt((cb->vertex(0)->point() - center) * (cb->vertex(0)->point() - center)));
      out << " " << min_cotan;
      out << " " << cb->vertex(0)->point().x() << " " << cb->vertex(0)->point().y() << " " << cb->vertex(0)->point().z();
      out << " " << cb->vertex(1)->point().x() << " " << cb->vertex(1)->point().y() << " " << cb->vertex(1)->point().z();
      out << " " << cb->vertex(2)->point().x() << " " << cb->vertex(2)->point().y() << " " << cb->vertex(2)->point().z();
      out << " " << cb->vertex(3)->point().x() << " " << cb->vertex(3)->point().y() << " " << cb->vertex(3)->point().z();
      out << std::endl;
    }

    out.close();
  }

  void save_steiner_point()
  {
    std::ofstream out("steiner.xyz");

    for(Finite_vertices_iterator vb = m_tr->finite_vertices_begin(); vb != m_tr->finite_vertices_end(); vb++)
    {
      if(vb->type() != Triangulation::INPUT){
        Point pb = vb->point();
        out << std::to_string(pb.x()) << " " << std::to_string(pb.y()) << " " << std::to_string(pb.z()) << std::endl;
      }
    }

    out.close();
  }

  FT cotan_per_edge(Cell_handle cell, int i, int j)
  {
    Vertex_handle vi = cell->vertex(i);
    Vertex_handle vj = cell->vertex(j);

    Point pi = vi->point();
    Point pj = vj->point();

    std::vector<Point> vpq;

    for(int i = 0; i < 4; i++)
      if(cell->vertex(i)->index() != vi->index() && cell->vertex(i)->index() != vj->index())
        vpq.push_back(cell->vertex(i)->point());

    Vector ni = CGAL::cross_product(pi - vpq[0], pi - vpq[1]);
    Vector nj = CGAL::cross_product(pj - vpq[0], pj - vpq[1]);

    ni = ni / std::sqrt(ni * ni);
    nj = nj / std::sqrt(nj * nj);

    Vector nij = CGAL::cross_product(ni, nj);
    FT cotan = (ni * nj) / std::sqrt(nij * nij);

    return cotan;
  }

  template <typename MatType>
  void check_spd(const MatType& A, int k = 5, int m = 37)
  {
      CGAL_TRACE("Begin checking eigenvalue...\n");
      VopType op(A);
      // Make sure B is positive definite and the decompoition is successful
      assert(op.info() == Spectra::SUCCESSFUL);

      Spectra::SymEigsSolver<FT, Spectra::SMALLEST_ALGE, VopType> eigs(&op, k, m);
      eigs.init();
      int nconv = eigs.compute(1000, 1e-8);

      Eigen::VectorXcd evalues;
      if(eigs.info() == Spectra::SUCCESSFUL)
      {
        evalues = eigs.eigenvalues();
        for(int i = 0; i < k; i++){
        std::cerr << i << "th value: " << evalues[i] << std::endl;
        }
      }
      else
        std::cerr << "Check failed! " << eigs.info() << std::endl;

  }


  /// Helping functions to assemble matrices

  /// Shift and orient the implicit function such that:
  /// - the implicit function = 0 for points / f() = contouring_value,
  /// - the implicit function < 0 inside the surface.
  ///
  /// Returns the minimum value of the implicit function.
  FT set_contouring_value(FT contouring_value)
  {
    // median value set to 0.0
    shift_f(-contouring_value);

    // Check value on convex hull (should be positive): if more than
    // half the vertices of the convex hull are negative, we flip the
    // sign (this is particularly useful if the surface is open, then
    // it is closed using the smallest part of the sphere).
    std::vector<Vertex_handle> convex_hull;
    m_tr->adjacent_vertices(m_tr->infinite_vertex(), std::back_inserter(convex_hull));
    unsigned int nb_negative = 0;
    for (std::size_t i = 0; i < convex_hull.size (); ++ i)
      if (convex_hull[i]->f() < 0.0)
        ++nb_negative;

    if(nb_negative > convex_hull.size() / 2)
      flip_f();

    // Update m_sink
    FT sink_value = find_sink();
    return sink_value;
  }

  template <class MatrixType>
  FT median_value_at_diagonal(const MatrixType& matrix, const int nb_variables)
  {
    Eigen::VectorXd diag_matrix = matrix.diagonal();
    std::sort(diag_matrix.data(), diag_matrix.data() + nb_variables);

    int mid = nb_variables / 2;

    if(nb_variables % 2 == 0)
      return 0.5 * (diag_matrix(mid) + diag_matrix(mid - 1));
    else
      return diag_matrix(mid);
  }


/// Gets median value of the implicit function over input vertices.
  FT median_value_at_input_vertices() const
  {
    std::deque<FT> values;
    Finite_vertices_iterator v, e;
    for(v = m_tr->finite_vertices_begin(),
        e = m_tr->finite_vertices_end();
        v != e;
        v++)
      if(v->type() == Triangulation::INPUT)
        values.push_back(v->f());

    std::size_t size = values.size();
    if(size == 0)
    {
      std::cerr << "Contouring: no input points\n";
      return 0.0;
    }

    std::sort(values.begin(),values.end());
    std::size_t index = size / 2;
    // return values[size/2];
    return 0.5 * (values[index] + values[index+1]); // avoids singular cases
  }

  void barycentric_coordinates(const Point& p,
                               Cell_handle cell,
                               FT& a,
                               FT& b,
                               FT& c,
                               FT& d) const
  {

    // const Point& pa = cell->vertex(0)->point();
    // const Point& pb = cell->vertex(1)->point();
    // const Point& pc = cell->vertex(2)->point();
    const Point& pd = cell->vertex(3)->point();
#if 1
    //Vector va = pa - pd;
    //Vector vb = pb - pd;
    //Vector vc = pc - pd;
    Vector vp = p - pd;

    //FT i00, i01, i02, i10, i11, i12, i20, i21, i22;
    //internal::invert(va.x(), va.y(), va.z(),
    //       vb.x(), vb.y(), vb.z(),
    //       vc.x(), vc.y(), vc.z(),
    //       i00, i01, i02, i10, i11, i12, i20, i21, i22);
    initialize_matrix_entry(cell);
    const Cached_bary_coord& i = (*m_bary)[cell->info()];

    //    UsedBary[cell->info()] = true;
    a = i[0] * vp.x() + i[3] * vp.y() + i[6] * vp.z();
    b = i[1] * vp.x() + i[4] * vp.y() + i[7] * vp.z();
    c = i[2] * vp.x() + i[5] * vp.y() + i[8] * vp.z();
    d = 1 - ( a + b + c);
#else
    FT v = volume(pa,pb,pc,pd);
    a = std::fabs(volume(pb,pc,pd,p) / v);
    b = std::fabs(volume(pa,pc,pd,p) / v);
    c = std::fabs(volume(pb,pa,pd,p) / v);
    d = std::fabs(volume(pb,pc,pa,p) / v);

    std::cerr << "_________________________________\n";
    std::cerr << aa << "  " << bb << "  " << cc << "  " << dd << std::endl;
    std::cerr << a << "  " << b << "  " << c << "  " << d << std::endl;

#endif
  }

  FT find_sink()
  {
    m_sink = CGAL::ORIGIN;
    FT min_f = 1e38;
    Finite_vertices_iterator v, e;
    for(v = m_tr->finite_vertices_begin(),
        e = m_tr->finite_vertices_end();
        v != e;
        v++)
    {
      if(v->f() < min_f)
      {
        m_sink = v->point();
        min_f = v->f();
      }
    }
    return min_f;
  }

  void shift_f(const FT shift)
  {
    Finite_vertices_iterator v, e;
    for(v = m_tr->finite_vertices_begin(),
        e = m_tr->finite_vertices_end();
        v != e;
        v++)
      v->f() += shift;
  }

  void flip_f()
  {
    Finite_vertices_iterator v, e;
    for(v = m_tr->finite_vertices_begin(),
        e = m_tr->finite_vertices_end();
        v != e;
        v++)
      v->f() = -v->f();
  }

  Vertex_handle any_vertex_on_convex_hull()
  {
    Cell_handle ch = m_tr->infinite_vertex()->cell();
    return  ch->vertex( (ch->index( m_tr->infinite_vertex())+1)%4);
  }

  void constrain_one_vertex_on_convex_hull(const FT value = 0.0)
  {
    Vertex_handle v = any_vertex_on_convex_hull();
    m_tr->constrain(v);
    v->f() = value;
  }

  // TODO: Some entities are computed too often
  // - nn and area should not be computed for the face and its opposite face
  //
  // divergent
  FT div_normalized(Vertex_handle v)
  {
    std::vector<Cell_handle> cells;
    cells.reserve(32);
    m_tr->incident_cells(v,std::back_inserter(cells));

    FT div = 0;
    typename std::vector<Cell_handle>::iterator it;
    for(it = cells.begin(); it != cells.end(); it++)
    {
      Cell_handle cell = *it;
      if(m_tr->is_infinite(cell))
        continue;

      // compute average normal per cell
      Vector n = get_cell_normal(cell);

      // zero normal - no need to compute anything else
      if(n == CGAL::NULL_VECTOR)
        continue;

      // compute n'
      int index = cell->index(v);
      const Point& x = cell->vertex(index)->point();
      const Point& a = cell->vertex((index+1)%4)->point();
      const Point& b = cell->vertex((index+2)%4)->point();
      const Point& c = cell->vertex((index+3)%4)->point();
      Vector nn = (index % 2 == 0) ? CGAL::cross_product(b - a, c - a) : CGAL::cross_product(c - a, b - a);
      nn = nn / std::sqrt(nn * nn); // normalize
      Vector p = a - x;
      Vector q = b - x;
      Vector r = c - x;
      FT p_n = std::sqrt(p * p);
      FT q_n = std::sqrt(q * q);
      FT r_n = std::sqrt(r * r);
      FT solid_angle = p * (CGAL::cross_product(q, r));
      solid_angle = std::abs(solid_angle / (p_n * q_n * r_n + (p * q) * r_n + (q * r) * p_n + (r * p) * q_n));

      FT area = std::sqrt(squared_area(a, b, c));
      FT length = p_n + q_n + r_n;
      div += n * nn * area / length;
    }
    return div * FT(3.0);
  }

  FT div(Vertex_handle v)
  {
    std::vector<Cell_handle> cells;
    cells.reserve(32);
    m_tr->incident_cells(v,std::back_inserter(cells));

    FT div = 0.0;
    typename std::vector<Cell_handle>::iterator it;
    for(it = cells.begin(); it != cells.end(); it++)
    {
      Cell_handle cell = *it;
      if(m_tr->is_infinite(cell))
        continue;

      const int index = cell->index(v);
      const Point& a = cell->vertex(m_tr->vertex_triple_index(index, 0))->point();
      const Point& b = cell->vertex(m_tr->vertex_triple_index(index, 1))->point();
      const Point& c = cell->vertex(m_tr->vertex_triple_index(index, 2))->point();
      const Vector nn = CGAL::cross_product(b - a, c - a);

      div+= nn * (//v->normal() +
                  m_tr->normal(cell->vertex((index + 1) % 4)) +
                  m_tr->normal(cell->vertex((index + 2) % 4)) +
                  m_tr->normal(cell->vertex((index + 3) % 4)));
    }
    return div;
  }

  Vector get_cell_normal(Cell_handle cell)
  {
    return Normal[cell->info()];
  }

  Vector cell_normal(Cell_handle cell) const
  {
    const Vector& n0 = m_tr->normal(cell->vertex(0));
    const Vector& n1 = m_tr->normal(cell->vertex(1));
    const Vector& n2 = m_tr->normal(cell->vertex(2));
    const Vector& n3 = m_tr->normal(cell->vertex(3));
    Vector n = n0 + n1 + n2 + n3;
    if(n != NULL_VECTOR){
      FT sq_norm = n * n;
      if(sq_norm != 0.0){
        return n / std::sqrt(sq_norm); // normalize
      }
    }
    return NULL_VECTOR;
  }

  // cotan formula as area(voronoi face) / len(primal edge)
  FT cotan_geometric(Edge& edge)
  {
    Cell_handle cell = edge.first;
    Vertex_handle vi = cell->vertex(edge.second);
    Vertex_handle vj = cell->vertex(edge.third);

    // primal edge
    const Point& pi = vi->point();
    const Point& pj = vj->point();
    Vector primal = pj - pi;
    FT len_primal = std::sqrt(primal * primal);

    return area_voronoi_face(edge) / len_primal;
  }

  // anisotropic Laplace coefficient (formula in paper)
  FT mcotan_geometric(Edge& edge, const FT cij, const FT ri, const FT rj, const bool convert)
  {
    Cell_handle cell = edge.first;
    Vertex_handle vi = cell->vertex(edge.second);
    Vertex_handle vj = cell->vertex(edge.third);

    // primal edge
    const Point& pi = vi->point();
    const Point& pj = vj->point();
    Vector primal = pj - pi;
    primal = primal / std::sqrt(primal * primal);

    // find normals
    Vector ni = m_tr->normal(vi);
		Vector nj = m_tr->normal(vj);

    // should use covariance to check isotropic
    Covariance cov_i(pi, ni, ri), cov_j(pj, nj, rj);
    Covariance cov_ij(cov_i, cov_j, convert);

    if(cov_ij.isotropic()) return cij;

    // average normals
    FT dot = cov_ij.ut_c_v(primal, primal);

		return cij * dot;
  }

  // anisotropic Laplace coefficient (formula derived from cotan_geometric)
  FT mcotan_geometric_in_metric(Edge& edge, const FT cij, const FT ri, const FT rj, const bool convert)
  {
    Cell_handle cell = edge.first;
    Vertex_handle vi = cell->vertex(edge.second);
    Vertex_handle vj = cell->vertex(edge.third);

    // primal edge
    const Point& pi = vi->point();
    const Point& pj = vj->point();
    Vector primal = pj - pi;

    // find normals
    Vector ni = m_tr->normal(vi);
		Vector nj = m_tr->normal(vj);

    // should use covariance to check isotropic
    Covariance cov_i(pi, ni, ri), cov_j(pj, nj, rj);
    Covariance cov_ij(cov_i, cov_j, convert);
    if(cov_ij.isotropic()) return cij;

    // calculate the voronoi area in a metric
    FT cell_area = area_voronoi_face_in_metric(edge, cov_ij);
    FT len_primal = std::sqrt(cov_ij.ut_c_v(primal, primal));

		return cell_area / len_primal;
  }

  // discrete Laplacian
  // defined as the product of the dihedral angle and the length of the corresponding edge
  FT cotan_laplacian(Edge& edge)
  {
    Cell_handle cell = edge.first;
    Vertex_handle vi = cell->vertex(edge.second);
    Vertex_handle vj = cell->vertex(edge.third);

    // circulate around edge
    Cell_circulator circ = m_tr->incident_cells(edge);
    Cell_circulator done = circ;
    FT cotan = 0;
    do
    {
      cell = circ;
      if(!m_tr->is_infinite(cell))
        cotan += cotan_dihedral_per_cell(cell, vi, vj);

      circ++;
    }
    while(circ != done);

    return cotan / 6;
  }

  // normal derivative
  FT cotan_normal_derivative(Cell_handle cell, int j, int f)
  {
      Vertex_handle vj = cell->vertex(j);
      Vertex_handle vf = cell->vertex(f);

      FT cotan = cotan_dihedral_per_cell(cell, vj, vf);
      return cotan / 6.;
  }

  // Given an edge ij in one cell, calculate the contangent of the dihedral angle and the edge length for its opposite edge
  FT cotan_dihedral_per_cell(Cell_handle cell, Vertex_handle vi, Vertex_handle vj)
  {
    Point pi = vi->point();
    Point pj = vj->point();

    std::vector<Point> vpq;

    for(int i = 0; i < 4; i++)
      if(cell->vertex(i)->index() != vi->index() && cell->vertex(i)->index() != vj->index())
        vpq.push_back(cell->vertex(i)->point());

    Vector ni = CGAL::cross_product(pi - vpq[0], pi - vpq[1]);
    Vector nj = CGAL::cross_product(pj - vpq[0], pj - vpq[1]);

    ni = ni / std::sqrt(ni * ni);
    nj = nj / std::sqrt(nj * nj);

    Vector lpq = vpq[0] - vpq[1];
    FT length_lpq = std::sqrt(lpq * lpq);

    Vector nij = CGAL::cross_product(ni, nj);
    FT cotan = (ni * nj) * length_lpq / std::sqrt(nij * nij);

    return cotan;
  }

  // anisotropic Laplace coefficient
  FT mcotan_laplacian(Edge& edge, const FT cij, const FT ri, const FT rj, const bool convert)
  {
    Cell_handle cell = edge.first;
    Vertex_handle vi = cell->vertex(edge.second);
    Vertex_handle vj = cell->vertex(edge.third);

    // primal edge
    const Point& pi = vi->point();
    const Point& pj = vj->point();

    // find normals
    Vector ni = m_tr->normal(vi);
		Vector nj = m_tr->normal(vj);

    // should use covariance to check isotropic
    Covariance cov_i(pi, ni, ri), cov_j(pj, nj, rj);
    Covariance cov_ij(cov_i, cov_j, convert);
    if(cov_ij.isotropic()) return cij;

    // circulate around edge
    Cell_circulator circ = m_tr->incident_cells(edge);
    Cell_circulator done = circ;
    FT mcotan = 0.;

    do
    {
      cell = circ;
      if(!m_tr->is_infinite(cell))
        mcotan += mcotan_dihedral_per_cell(cell, vi, vj, cov_ij, ri);

      circ++;
    }
    while(circ != done);

    //return 1.0;
    return mcotan / 6.;
  }

  FT mcotan_dihedral_per_cell(Cell_handle cell, Vertex_handle vi, Vertex_handle vj, Covariance& cov_ij, const FT ri)
  {
    Point pi = vi->point();
    Point pj = vj->point();

    std::vector<Point> vpq;
    std::vector<Vector> npq;

    for(int i = 0; i < 4; i++)
      if(cell->vertex(i)->index() != vi->index() && cell->vertex(i)->index() != vj->index()){
        vpq.push_back(cell->vertex(i)->point());
        npq.push_back(m_tr->normal(cell->vertex(i)));
      }

    Covariance cov_p(vpq[0], npq[0], ri), cov_q(vpq[1], npq[1], ri);
    Covariance cov_pq(cov_p, cov_q, true);
    Covariance ctet(cov_ij, cov_pq, true);

    Vector ni = CGAL::cross_product(pi - vpq[0], pi - vpq[1]);
    Vector nj = CGAL::cross_product(pj - vpq[0], pj - vpq[1]);

    ni = ni / std::sqrt(ni * ni);
    nj = nj / std::sqrt(nj * nj);

    Vector lpq = vpq[0] - vpq[1];
    FT length_lpq = std::sqrt(ctet.ut_c_v(lpq, lpq));

    FT dot_pq = ctet.ut_c_v(ni, nj);
    FT cross_pq = std::sqrt(squared_area_in_metric(ni, nj, ctet));

    FT mcotan = dot_pq * length_lpq / cross_pq;

    return mcotan;
  }

  FT squared_area_in_metric(const Point& a, const Point& b, const Point& c, Covariance& cov)
  {
    Vector u = b - a;
    Vector v = c - a;
    FT ut_cov_u = cov.ut_c_v(u, u);
    FT vt_cov_v = cov.ut_c_v(v, v);
    FT ut_cov_v = cov.ut_c_v(u, v);

    return ut_cov_u * vt_cov_v - ut_cov_v * ut_cov_v;
  }

  FT squared_area_in_metric(const Vector& u, const Vector& v, Covariance& cov)
  {
    FT ut_cov_u = cov.ut_c_v(u, u);
    FT vt_cov_v = cov.ut_c_v(v, v);
    FT ut_cov_v = cov.ut_c_v(u, v);

    return ut_cov_u * vt_cov_v - ut_cov_v * ut_cov_v;
  }

  // spin around edge
  // return area(voronoi face)
  FT area_voronoi_face(Edge& edge)
  {
    // circulate around edge
    Cell_circulator circ = m_tr->incident_cells(edge);
    Cell_circulator done = circ;
    std::vector<Point> voronoi_points;
    voronoi_points.reserve(9);
    do
    {
      Cell_handle cell = circ;
      if(!m_tr->is_infinite(cell))
        voronoi_points.push_back(Dual[cell->info()]);
      else // one infinite tet, switch to another calculation
        return area_voronoi_face_boundary(edge);
      circ++;
    }
    while(circ != done);

    if(voronoi_points.size() < 3)
    {
      CGAL_surface_reconstruction_points_assertion(false);
      return 0.0;
    }

    // sum up areas
    FT area = 0.0;
    const Point& a = voronoi_points[0];
    std::size_t nb_triangles = voronoi_points.size() - 1;
    for(std::size_t i = 1; i < nb_triangles; i++)
    {
      const Point& b = voronoi_points[i];
      const Point& c = voronoi_points[i + 1];
      area += std::sqrt(squared_area(a, b, c));
    }
    return area;
  }

  // spin around edge
  // return area(voronoi face) in a specific metric
  FT area_voronoi_face_in_metric(Edge& edge, Covariance& cov_ij)
  {
    // circulate around edge
    Cell_circulator circ = m_tr->incident_cells(edge);
    Cell_circulator done = circ;
    std::vector<Point> voronoi_points;
    voronoi_points.reserve(9);
    do
    {
      Cell_handle cell = circ;
      if(!m_tr->is_infinite(cell))
        voronoi_points.push_back(Dual[cell->info()]);
      else // one infinite tet, switch to another calculation
        return area_voronoi_face_boundary_in_metric(edge, cov_ij);
      circ++;
    }
    while(circ != done);

    if(voronoi_points.size() < 3)
    {
      CGAL_surface_reconstruction_points_assertion(false);
      return 0.0;
    }

    // sum up areas
    FT area = 0.0;
    const Point& a = voronoi_points[0];
    std::size_t nb_triangles = voronoi_points.size() - 1;
    for(std::size_t i=1;i<nb_triangles;i++)
    {
      const Point& b = voronoi_points[i];
      const Point& c = voronoi_points[i+1];
      area += std::sqrt(squared_area_in_metric(a, b, c, cov_ij)) / 2.;
    }
    return area;
  }

  // approximate area when a cell is infinite
  FT area_voronoi_face_boundary(Edge& edge)
  {
    FT area = 0.0;
    Vertex_handle vi = edge.first->vertex(edge.second);
    Vertex_handle vj = edge.first->vertex(edge.third);

    const Point& pi = vi->point();
    const Point& pj = vj->point();
    Point m = CGAL::midpoint(pi, pj);

    // circulate around each incident cell
    Cell_circulator circ = m_tr->incident_cells(edge);
    Cell_circulator done = circ;
    do
    {
      Cell_handle cell = circ;
      if(!m_tr->is_infinite(cell))
      {
        // circumcenter of cell
        Point c = Dual[cell->info()];
        Tetrahedron tet = m_tr->tetrahedron(cell);

        int i = cell->index(vi);
        int j = cell->index(vj);
        int k = Triangulation_utils_3::next_around_edge(i, j);
        int l = Triangulation_utils_3::next_around_edge(j, i);

        Vertex_handle vk = cell->vertex(k);
        Vertex_handle vl = cell->vertex(l);

        const Point& pk = vk->point();
        const Point& pl = vl->point();

        // if circumcenter is outside tet
        // pick barycenter instead
        if(tet.has_on_unbounded_side(c))
        {
          Point cell_points[4] = {pi, pj, pk, pl};
          c = CGAL::centroid(cell_points, cell_points+4);
        }

        Point ck = CGAL::circumcenter(pi, pj, pk);
        Point cl = CGAL::circumcenter(pi, pj, pl);

        area += std::sqrt(squared_area(m, c, ck));
        area += std::sqrt(squared_area(m, c, cl));
      }
      circ++;
    }
    while(circ != done);
    return area;
  }

  // approximate area when a cell is infinite
  FT area_voronoi_face_boundary_in_metric(Edge& edge, Covariance& cov_ij)
  {
    FT area = 0.0;
    Vertex_handle vi = edge.first->vertex(edge.second);
    Vertex_handle vj = edge.first->vertex(edge.third);

    const Point& pi = vi->point();
    const Point& pj = vj->point();
    Point m = CGAL::midpoint(pi,pj);

    // circulate around each incident cell
    Cell_circulator circ = m_tr->incident_cells(edge);
    Cell_circulator done = circ;
    do
    {
      Cell_handle cell = circ;
      if(!m_tr->is_infinite(cell))
      {
        // circumcenter of cell
        Point c = Dual[cell->info()];
        Tetrahedron tet = m_tr->tetrahedron(cell);

        int i = cell->index(vi);
        int j = cell->index(vj);
        int k = Triangulation_utils_3::next_around_edge(i,j);
        int l = Triangulation_utils_3::next_around_edge(j,i);

        Vertex_handle vk = cell->vertex(k);
        Vertex_handle vl = cell->vertex(l);

        const Point& pk = vk->point();
        const Point& pl = vl->point();

        // if circumcenter is outside tet
        // pick barycenter instead
        if(tet.has_on_unbounded_side(c))
        {
          Point cell_points[4] = {pi, pj, pk, pl};
          c = CGAL::centroid(cell_points, cell_points + 4);
        }

        Point ck = CGAL::circumcenter(pi, pj, pk);
        Point cl = CGAL::circumcenter(pi, pj, pl);

        area += std::sqrt(squared_area_in_metric(m, c, ck, cov_ij));
        area += std::sqrt(squared_area_in_metric(m, c, cl, cov_ij));
      }
      circ++;
    }
    while(circ != done);
    return area;
  }

  FT volume_voronoi_cell(Vertex_handle v)
  {
    if(!has_finite_voronoi_cell(v))
      return approx_volume_voronoi_cell(v);

    std::list<Tetrahedron> tetrahedra;
    tessellate_voronoi_cell(v, tetrahedra);
    return volume(tetrahedra);
  }

  bool has_finite_voronoi_cell(Vertex_handle v)
  {
    std::list<Cell_handle> cells;
    m_tr->incident_cells(v, std::back_inserter(cells));

    if(cells.size() == 0)
      return false;

    typename std::list<Cell_handle>::iterator it;
    for(it = cells.begin(); it != cells.end(); it++)
    {
      if(m_tr->is_infinite(*it))
        return false;
    }
    return true;
  }

  FT approx_volume_voronoi_cell(Vertex_handle v)
  {
    FT total_volume = 0.0;

    // get all cells incident to v
    std::list<Cell_handle> cells;
    m_tr->incident_cells(v, std::back_inserter(cells));
    typename std::list<Cell_handle>::iterator it;
    for(it = cells.begin(); it != cells.end(); it++)
    {
      Cell_handle cell = *it;

      if(m_tr->is_infinite(cell))
        continue;

      Tetrahedron tet = m_tr->tetrahedron(cell);
      total_volume += std::abs(tet.volume());
    }
    return 0.25 * total_volume; // approximation! Could use circumenter insted, as this one uses implicitly the
  }

  bool tessellate_voronoi_cell(Vertex_handle v, std::list<Tetrahedron>& tetrahedra,
                                                const bool add_to_vertex = false)
  {
    Point a = v->point();

    // get all vertices incident to v
    std::list<Vertex_handle> vertices;
    m_tr->incident_vertices(v, std::back_inserter(vertices));
    typename std::list<Vertex_handle>::iterator it;
    for(it = vertices.begin(); it != vertices.end(); it++)
    {
      // build edge from two vertices
      Vertex_handle v2 = *it;
      Cell_handle cell;
      int i1, i2;
      if(!m_tr->is_edge(v, v2, cell, i1, i2))
        return false;
      Edge edge(cell, i1, i2);

      // spin around edge to get incident cells
      Cell_circulator c = m_tr->incident_cells(edge);
      Cell_circulator done = c;
      unsigned int degree = 0;
      do
        degree++;
      while(++c != done);
      assert(degree >= 3);

      // choose first as pivot
      Point b = m_tr->dual(c);
      Cell_circulator curr = m_tr->incident_cells(edge);
      curr++;
      Cell_circulator next = m_tr->incident_cells(edge);
      next++;
      next++;
      unsigned int nb_tets = degree - 2;
      for(unsigned int i = 0; i < nb_tets; i++)
      {
        Point c = m_tr->dual(curr);
        Point d = m_tr->dual(next);
        Tetrahedron tet(a, b, c, d);
        //if(add_to_vertex)
          //v->add(tet);
        //else
        tetrahedra.push_back(tet);
        curr++;
        next++;
      }
    }
    return true;
  }

  FT volume(std::list<Tetrahedron>& tetrahedra)
  {
    FT total_volume = 0.0;
    typename std::list<Tetrahedron>::iterator it;
    for(it = tetrahedra.begin(); it != tetrahedra.end(); it++)
    {
      Tetrahedron& tetrahedron = *it;
      total_volume += std::abs(tetrahedron.volume());
    }
    return total_volume;
  }

  FT volume(Cell_handle cell)
  {
    Point a = cell->vertex(0)->point();
    Point b = cell->vertex(1)->point();
    Point c = cell->vertex(2)->point();
    Point d = cell->vertex(3)->point();

    Tetrahedron tet(a, b, c, d);

    return std::abs(tet.volume());
  }


  /// Assemble vi's row of the linear system A*X=B
  ///
  /// @commentheading Template parameters:
  /// @param SparseLinearAlgebraTraits_d Symmetric definite positive sparse linear solver.
  template <class SparseLinearAlgebraTraits_d>
  void assemble_poisson_row(typename SparseLinearAlgebraTraits_d::Matrix& A,
                            Vertex_handle vi,
                            typename SparseLinearAlgebraTraits_d::Vector& B,
                            double lambda)
  {
    // for each vertex vj neighbor of vi
    std::vector<Edge> edges;
    m_tr->incident_edges(vi,std::back_inserter(edges));

    double diagonal = 0.0;

    for(typename std::vector<Edge>::iterator it = edges.begin();
        it != edges.end();
        it++)
      {
        Vertex_handle vj = it->first->vertex(it->third);
        if(vj == vi){
          vj = it->first->vertex(it->second);
        }
        if(m_tr->is_infinite(vj))
          continue;

        // get corresponding edge
        Edge edge(it->first, it->first->index(vi), it->first->index(vj));
        if(vi->index() < vj->index()){
          std::swap(edge.second, edge.third);
        }

        double cij = cotan_geometric(edge);

        if(m_tr->is_constrained(vj)){
          if(! is_valid(vj->f())){
            std::cerr << "vj->f() = " << vj->f() << " is not valid" << std::endl;
          }
          B[vi->index()] -= cij * vj->f(); // change rhs
          if(! is_valid( B[vi->index()])){
            std::cerr << " B[vi->index()] = " <<  B[vi->index()] << " is not valid" << std::endl;
          }

        } else {
          if(! is_valid(cij)){
            std::cerr << "cij = " << cij << " is not valid" << std::endl;
          }
          A.set_coef(vi->index(),vj->index(), -cij, true /*new*/); // off-diagonal coefficient
        }

        diagonal += cij;
      }
    // diagonal coefficient
    if (vi->type() == Triangulation::INPUT){
      A.set_coef(vi->index(),vi->index(), diagonal + lambda, true /*new*/) ;
    } else{
      A.set_coef(vi->index(),vi->index(), diagonal, true /*new*/);
    }
  }

  /// Assemble vi's row of the GEV system
  ///
  /// @commentheading Template parameters:
  template <class ConfidenceMap>
  void assemble_spectral_row(Vertex_handle vi,
                             Matrix& AA,
                             Matrix& L,
                             Matrix& F,
                             Matrix& V_inv,
                             Matrix& V,
                             Matrix& N,
                             const FT fitting,
                             ConfidenceMap confidence_map,
                             FT& duration_assign,
                             FT& duration_cal,
                             bool flag_boundary)
  {
    // for each vertex vj neighbor of vi
    std::vector<Edge> edges;
    m_tr->incident_edges(vi, std::back_inserter(edges));

    double diagonal = 0.0;
    double mdiagonal = 0.0;
    double time_init = 0.0;

    for(typename std::vector<Edge>::iterator it = edges.begin();
        it != edges.end();
        it++)
      {
        Vertex_handle vj = it->first->vertex(it->third);
        if(vj == vi){
          vj = it->first->vertex(it->second);
        }
        if(m_tr->is_infinite(vj))
          continue;

        // get corresponding edge
        Edge edge(it->first, it->first->index(vi), it->first->index(vj));
        if(vi->index() < vj->index()){
          std::swap(edge.second, edge.third);
        }

        time_init = clock();

        FT cij = cotan_laplacian(edge);
        //FT cij = cotan_geometric(edge);

        bool convert = true;
        FT ri = get(confidence_map, it->first->vertex(edge.second)->input_iterator());
        FT rj = get(confidence_map, it->first->vertex(edge.third)->input_iterator());
        FT mcij = mcotan_laplacian(edge, cij, ri, rj, convert);
        //FT mcij = mcotan_geometric(edge, cij, ri, rj, convert);

        duration_cal += clock() - time_init; time_init = clock();

        AA.set_coef(vi->index(), vj->index(), -mcij, true);
        L.set_coef(vi->index(), vj->index(), -cij, true);

        duration_assign += clock() - time_init; time_init = clock();

        diagonal += cij;
        mdiagonal += mcij;
      }
    // diagonal coefficient

    FT vol = volume_voronoi_cell(vi);
    duration_cal += clock() - time_init; time_init = clock();

    time_init = clock();
    AA.set_coef(vi->index(), vi->index(), mdiagonal, true);
    L.set_coef(vi->index(), vi->index(), diagonal, true);
    V_inv.set_coef(vi->index(), vi->index(), std::min(1.0 / vol, 1e7), true);
    V.set_coef(vi->index(), vi->index(), vol, true);

    if(vi->type() == Triangulation::INPUT)
      F.set_coef(vi->index(), vi->index(), fitting, true);

    duration_assign += clock() - time_init;

    // normal derivative for boundary points
    if(flag_boundary && (vi->position() == Triangulation::BOUNDARY))
    {
      std::vector<Facet> facets;
      m_tr->incident_facets(vi, std::back_inserter(facets));

      // get mirror facets
      size_t number_facets = facets.size();
      for(size_t i = 0; i < number_facets; i++)
        facets.push_back(m_tr->mirror_facet(facets[i]));

      for(typename std::vector<Facet>::iterator facet = facets.begin(); facet != facets.end(); facet++){
        Cell_handle cell = facet->first;
        int index_f = facet->second;

        if(m_tr->is_infinite(cell))
          continue;

        bool flag = true;
        for(int i = 0; i < 4; i++)
          if((i != index_f) && (cell->vertex(i)->position() != Triangulation::BOUNDARY)){
            flag = false;
            break;
          }

        if(!flag) continue;

        for(int j = 0; j < 4; j++)
          if(j != index_f){
            FT njf = cotan_normal_derivative(cell, j, index_f);
            N.add_coef(vi->index(), cell->vertex(index_f)->index(), njf);
            N.add_coef(vi->index(), cell->vertex(j)->index(), -njf);
          }

      }
    }
  }

public:

  /// Marching Tetrahedra
  unsigned int marching_tetrahedra(const FT value, Mesh &mesh)
  {
    std::vector<Point> points;
    std::vector< std::vector<std::size_t> > polygons;

    return m_tr->marching_tets(value, mesh, points, polygons);
  }

  Point draw_xslice_function(
		const unsigned int size,
		const double x,
    const int mode,
    const std::string outfile)
  {
    Point_list point_xslice;
    Color_list rgb_xslice;

    Point center = bounding_sphere().center();
    double radius = sqrt(bounding_sphere().squared_radius()) * m_enlarge_ratio;

    double ymin = center.y() - radius, ymax = center.y() + radius;
    double zmin = center.z() - radius, zmax = center.z() + radius;

    const double yincr = (ymax - ymin) / size;
		const double zincr = (zmax - zmin) / size;

    double my_fmin = 1e10;
    double my_fmax = -1e10;

    Point my_pmax;

    Cell_handle hint;
    double y = ymin;
    for(unsigned int i = 0; i < size; i++)
    {
      double z = zmin;

      for(unsigned int j = 0; j < size; j++)
      {
        Point a(x, y ,z);
        double va = 0.;
        //bool ba = locate_and_evaluate_function(a, hint, va);
        bool ba = locate_and_evaluate_function(a, hint, va, mode);

        if(ba)
        {
          if(va < my_fmin) my_fmin = va;
          else if(va > my_fmax){
            my_pmax = a;
            my_fmax = va;
          }

          point_xslice.push_back(std::make_pair(a, va));
        }

        z += zincr;
      }
      y += yincr;
    }

    std::cerr << "fmin: " << my_fmin << std::endl;
    std::cerr << "fmax: " << my_fmax << std::endl;

    for(typename Point_list::iterator e = point_xslice.begin(); e != point_xslice.end(); e++){
    //for(const auto &e : point_xslice){
      Color my_color;
      color_and_vertex_function(e->second, my_color, my_fmin, my_fmax);
      rgb_xslice.push_back(my_color);
    }
    save_slice(point_xslice, rgb_xslice, outfile);

    return my_pmax;
  }

private:

  bool save_slice(Point_list& point_xslice, Color_list& rgb_xslice, const std::string outfile)
  {
    if(rgb_xslice.size() == 0) return false;

    std::vector<PC> pc_xslice;
    std::ofstream out("value_" + outfile);
    CGAL::set_binary_mode(out);

    for(int i = 0; i < (int) rgb_xslice.size(); i++)
      pc_xslice.push_back(CGAL::cpp11::make_tuple(point_xslice[i].first, rgb_xslice[i]));

    point_xslice.clear();
    rgb_xslice.clear();

    CGAL::write_ply_points_with_properties(out, pc_xslice, CGAL::make_ply_point_writer(VF_point_map()),
                                          std::make_tuple(VF_color_map(),
                                          CGAL::PLY_property<unsigned char>("red"),
                                          CGAL::PLY_property<unsigned char>("green"),
                                          CGAL::PLY_property<unsigned char>("blue")));

    return true;
  }

  bool locate_and_evaluate_function(const Point& query, Cell_handle hint, double& value, const int mode)
  {
    typename Triangulation::Locate_type lt;
    int li, lj;
    Cell_handle cell = m_tr -> locate(query, lt, li, lj, hint);
    if(lt == Triangulation::CELL)
    {
      hint = cell;
      FT a, b, c, d;
      barycentric_coordinates(query, cell, a, b, c, d);
      if(mode == 0)
        value =  a * cell->vertex(0)->f() +
          b * cell->vertex(1)->f() +
          c * cell->vertex(2)->f() +
          d * cell->vertex(3)->f();
      else if(mode == 1)
        value =  a * cell->vertex(0)->lf() +
          b * cell->vertex(1)->lf() +
          c * cell->vertex(2)->lf() +
          d * cell->vertex(3)->lf();
      else if(mode == 2)
        value =  a * cell->vertex(0)->bf() +
          b * cell->vertex(1)->bf() +
          c * cell->vertex(2)->bf() +
          d * cell->vertex(3)->bf();
      else if(mode == 3)
        value =  a * cell->vertex(0)->af() +
          b * cell->vertex(1)->af() +
          c * cell->vertex(2)->af() +
          d * cell->vertex(3)->af();
      else if(mode == 4)
        value =  a * cell->vertex(0)->check() +
          b * cell->vertex(1)->check() +
          c * cell->vertex(2)->check() +
          d * cell->vertex(3)->check();
      return true;
    }
    return false;
  }

/*
  bool locate_and_evaluate_function(const Point& query, Cell_handle hint, double& value)
  {
    typename Triangulation::Locate_type lt;
    int li, lj;
    Cell_handle cell = m_tr -> locate(query, lt, li, lj, hint);
    if(lt == Triangulation::CELL)
    {
      hint = cell;
      FT a, b, c, d;
      barycentric_coordinates(query, cell, a, b, c, d);

      value =  a * cell->vertex(0)->f() +
          b * cell->vertex(1)->f() +
          c * cell->vertex(2)->f() +
          d * cell->vertex(3)->f();
      return true;
    }

    return false;
  }*/

  void color_and_vertex_function(const double value, Color& color, const double min_value, const double max_value)
  {
    double ratio = (value - min_value) / (max_value - min_value);
    get_rainbow_color(ratio, color);
  }

  void get_rainbow_color(const double ratio, Color& color)
  {
    int h = int(ratio * 256 * 6);
    int x = h % 256;

    switch(h / 256)
    {
      case 0: color[0] = 255;     color[1] = x;       color[2] = 0; break;
      case 1: color[0] = 255 - x; color[1] = 255;     color[2] = 0; break;
      case 2: color[0] = 0;       color[1] = 255;     color[2] = x; break;
      case 3: color[0] = 0;       color[1] = 255 - x; color[2] = 255; break;
      case 4: color[0] = x;       color[1] = 0;       color[2] = 255; break;
      case 5: color[0] = 255;     color[1] = 0;       color[2] = 255 - x; break;
    }
  }

}; // end of Implicit_reconstruction_function

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_IMPLICIT_RECONSTRUCTION_FUNCTION_H
