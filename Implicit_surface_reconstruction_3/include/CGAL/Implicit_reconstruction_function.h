// Copyright (c) 2007-09  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Laurent Saboret, Pierre Alliez, Tong Zhao

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

#include <CGAL/trace.h>
#include <CGAL/Reconstruction_triangulation_3.h>
#include <CGAL/Covariance_matrix_3.h>
#include <CGAL/spatial_sort.h>
#ifdef CGAL_EIGEN3_ENABLED
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <CGAL/Eigen_solver_traits.h>
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

#include <boost/shared_ptr.hpp>
#include <boost/array.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/range.hpp>

#include <SymGEigsSolver.h>
#include <MatOp/SparseSymMatProd.h>
#include <MatOp/SparseCholesky.h>
#include <unsupported/Eigen/SparseExtra>


/*! 
  \file Implicit_reconstruction_function.h
*/

namespace CGAL {

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
\ingroup PkgImplicitSurfaceReconstruction

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
computes a n implicit function by solving a generalized eigenvalue problem
such that its gradient is most aligned with the principal axes of a tensor field.
The principal axes and eccentricities of the tensor field locally represent
respectively the most likely direction of the normal to the surface, and the 
confidence in this direction estimation.

The GEP is solved by Spectra library.


\tparam Gt Geometric traits class. 

\cgalModels `ImplicitFunction`

*/
template <class Gt>
class Implicit_reconstruction_function
{
// Public types
public:

  /// \name Types 
  /// @{

  typedef Gt Geom_traits; ///< Geometric traits class
  /// \cond SKIP_IN_MANUAL
  typedef Reconstruction_triangulation_3<Robust_circumcenter_filtered_traits_3<Gt> >
                                                   Triangulation;
  /// \endcond
  typedef typename Triangulation::Cell_handle   Cell_handle;

  // Geometric types
  typedef typename Geom_traits::FT FT; ///< number type.
  typedef typename Geom_traits::Point_3 Point; ///< point type.
  typedef typename Geom_traits::Vector_3 Vector; ///< vector type.
  typedef typename Geom_traits::Sphere_3 Sphere; 

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

  typedef typename CGAL::Eigen_sparse_matrix<FT>            Matrix;
  typedef typename Eigen::SparseMatrix<FT>                  ESMatrix;
  typedef typename Eigen::Matrix<FT, Eigen::Dynamic, Eigen::Dynamic>  EMatrix;
  typedef typename CGAL::Covariance_matrix_3<Geom_traits>             Covariance;

  typedef typename Spectra::SparseSymMatProd<FT>   OpType;
  typedef typename Spectra::SparseCholesky<FT>     BOpType;

  typedef CGAL::cpp11::tuple<Point, FT> Point_with_property;
  typedef CGAL::Nth_of_tuple_property_map<0, Point_with_property> PP_point_map;
  typedef CGAL::Nth_of_tuple_property_map<1, Point_with_property> PP_func_map;


// Data members.
// Warning: the Surface Mesh Generation package makes copies of implicit functions,
// thus this class must be lightweight and stateless.
private:

  // operator() is pre-computed on vertices of *m_tr by solving
  // ...
  boost::shared_ptr<Triangulation> m_tr;

  mutable boost::shared_ptr<std::vector<boost::array<double,9> > > m_Bary;
  mutable std::vector<Point> Dual;
  mutable std::vector<Vector> Normal;
  mutable std::vector<FT> Reliability;

  // contouring and meshing
  Point m_sink; // Point with the minimum value of operator()
  mutable Cell_handle m_hint; // last cell found = hint for next search

  FT average_spacing;


  /// function to be used for the different constructors available that are
  /// doing the same thing but with default template parameters
  template <typename PointRange,
            typename PointMap,
            typename NormalPMap,
            typename Visitor
  >
  void forward_constructor(
    PointRange& points,
    PointMap point_map,
    NormalPMap normal_pmap,
    Visitor visitor)
  {
    CGAL::Timer task_timer; task_timer.start();
    CGAL_TRACE_STREAM << "Creates Implicit triangulation...\n";

    // Inserts points in triangulation
    m_tr->insert(
      points,
      point_map,
      normal_pmap,
      visitor);

    // Prints status
    CGAL_TRACE_STREAM << "Creates Implicit triangulation: " << task_timer.time() << " seconds, "
                                                           << std::endl;
  }


// Public methods
public:

  /// \name Creation 
  /// @{


  /*! 
    Creates a Implicit function from the  range of points `[first, beyond)`. 

    \tparam PointRange is a model of `Range`. The value type of
   its iterator is the key type of the named parameter `point_map`.

    \tparam PointMap is a model of `ReadablePropertyMap` with
      a `value_type = Point`.  It can be omitted if `InputIterator`
      `value_type` is convertible to `Point`. 
    
    \tparam NormalPMap is a model of `ReadablePropertyMap`
      with a `value_type = Vector`.
  */ 
  template <typename PointRange,
            typename PointMap,
            typename NormalPMap
  >
  Implicit_reconstruction_function(
    PointRange& points, ///< input point range
    PointMap point_map, ///< property map: `value_type of InputIterator` -> `Point` (the position of an input point).
    NormalPMap normal_pmap ///< property map: `value_type of InputIterator` -> `Vector` (the *oriented* normal of an input point).
)
    : m_tr(new Triangulation), m_Bary(new std::vector<boost::array<double,9> > ),
    average_spacing(CGAL::compute_average_spacing<CGAL::Sequential_tag>
                      (points, 6,
                       CGAL::parameters::point_map(point_map)))
  {
    forward_constructor(points, point_map, normal_pmap, Implicit_visitor());
  }

  /// \cond SKIP_IN_MANUAL
  template <typename PointRange,
            typename PointMap,
            typename NormalPMap,
            typename Visitor
  >
  Implicit_reconstruction_function(
    PointRange& points, ///< input point range
    PointMap point_map, ///< property map: `value_type of InputIterator` -> `Point` (the position of an input point).
    NormalPMap normal_pmap, ///< property map: `value_type of InputIterator` -> `Vector` (the *oriented* normal of an input point).
    Visitor visitor)
    : m_tr(new Triangulation), m_Bary(new std::vector<boost::array<double,9> > ),
    average_spacing(CGAL::compute_average_spacing<CGAL::Sequential_tag>(points, 6,
                                                CGAL::parameters::point_map(point_map)))
  {
    forward_constructor(points, point_map, normal_pmap, visitor);
  }

  // This variant creates a default point property map = Identity_property_map and Visitor=Implicit_visitor
  template <typename PointRange,
            typename NormalPMap
  >
  Implicit_reconstruction_function(
    PointRange& points, ///< input point range
    NormalPMap normal_pmap, ///< property map: `value_type of InputIterator` -> `Vector` (the *oriented* normal of an input point).
    typename boost::enable_if<
      boost::is_convertible<typename std::iterator_traits<typename PointRange::iterator>::value_type, Point>
    >::type* = 0
  )
  : m_tr(new Triangulation), m_Bary(new std::vector<boost::array<double,9> > )
  , average_spacing(CGAL::compute_average_spacing<CGAL::Sequential_tag>(points, 6))
  {
    forward_constructor(points, 
      make_identity_property_map(
      typename std::iterator_traits<typename PointRange::iterator>::value_type()),
      normal_pmap, Implicit_visitor());
    CGAL::Timer task_timer; task_timer.start();
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


  template <class Visitor>
  bool first_delaunay_refinement(Visitor visitor,
                                 const FT approximation_ratio = 0.,
                                 const FT radius_edge_ratio_bound = 2.5,
                                 const unsigned int max_vertices = (unsigned int)1e7,
                                 const FT enlarge_ratio = 1.5)
  {
    CGAL::Timer refine_timer;
    CGAL_TRACE_STREAM << "Delaunay refinement...\n";

    // Delaunay refinement
    const FT radius = sqrt(bounding_sphere().squared_radius()); // get triangulation's radius
    const FT cell_radius_bound = radius/5.; // large

    internal::Implicit::Constant_sizing_field<Triangulation> sizing_field(CGAL::square(cell_radius_bound));

    std::vector<int> NB; 

    NB.push_back( delaunay_refinement(radius_edge_ratio_bound,sizing_field,max_vertices,enlarge_ratio));

    while(m_tr->insert_fraction(visitor))
      NB.push_back( delaunay_refinement(radius_edge_ratio_bound,sizing_field,max_vertices,enlarge_ratio));
    
    /*
    if(approximation_ratio > 0. && 
       approximation_ratio * std::distance(m_tr->input_points_begin(),
                                           m_tr->input_points_end()) > 20)
      second_delaunay_refinement(visitor, approximation_ratio, NB);*/

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

  /*
  template <class Visitor>
  bool second_delaunay_refinement(Visitor visitor,
                                  const FT approximation_ratio,
                                  std::vector<int>& NB,
                                  bool flag_spectral = false)
  {
      // Add a pass of Delaunay refinement.
      //
      // In that pass, the sizing field, of the refinement process of the
      // triangulation, is based on the result of a poisson function with a
      // sample of the input points. The ratio is 'approximation_ratio'.
      //
      // For optimization reasons, the cell criteria of the refinement
      // process uses two sizing fields:
      //
      //   - the minimum of the square of 'coarse_implicit_function' and the
      // square of the constant field equal to 'average_spacing',
      //
      //   - a second sizing field that is constant, and equal to:
      //
      //         average_spacing*average_spacing_ratio
      //
      // If a given cell is smaller than the constant second sizing field,
      // then the cell is considered as small enough, and the first sizing
      // field, more costly, is not evaluated.

      typedef Filter_iterator<typename Triangulation::Input_point_iterator,
                              Implicit_skip_vertices> Some_points_iterator;
      //make it deterministic
      Random random(0);
      Implicit_skip_vertices skip(1.-approximation_ratio,random);
      
      CGAL_TRACE_STREAM << "SPECIAL PASS that uses an approximation of the result (approximation ratio: "
                << approximation_ratio << ")" << std::endl;
      CGAL::Timer approximation_timer; approximation_timer.start();

      CGAL::Timer sizing_field_timer; sizing_field_timer.start();

      typedef boost::iterator_range<Some_points_iterator> Some_points_range;
      Some_points_range some_points = boost::make_iterator_range(Some_points_iterator(m_tr->input_points_end(),
                                                     skip,
                                                     m_tr->input_points_begin()),
                                          Some_points_iterator(m_tr->input_points_end(),
                                                     skip));

      Implicit_reconstruction_function<Geom_traits> 
        coarse_implicit_function(some_points,
                                Normal_of_point_with_normal_map<Geom_traits>() );
      
      if(flag_spectral)
        coarse_implicit_function.compute_spectral_implicit_function(Implicit_visitor(), bilaplacian, laplacian, fitting, ratio, mode, 0.);
      else
        coarse_implicit_function.compute_poisson_implicit_function(solver, Implicit_visitor(), 0.);

      internal::Implicit::Constant_sizing_field<Triangulation> 
        min_sizing_field(CGAL::square(average_spacing));
      internal::Implicit::Constant_sizing_field<Triangulation> 
        sizing_field_ok(CGAL::square(average_spacing*average_spacing_ratio));

      Special_wrapper_of_two_functions_keep_pointers<
        internal::Implicit::Constant_sizing_field<Triangulation>,
        Implicit_reconstruction_function<Geom_traits> > sizing_field2(&min_sizing_field,
                                                            &coarse_implicit_function);
        
      sizing_field_timer.stop();
      std::cerr << "Construction time of the sizing field: " << sizing_field_timer.time() 
                << " seconds" << std::endl;

      NB.push_back( delaunay_refinement(radius_edge_ratio_bound,
                                        sizing_field2,
                                        max_vertices,
                                        enlarge_ratio,
                                        sizing_field_ok) );
      approximation_timer.stop();
      CGAL_TRACE_STREAM << "SPECIAL PASS END (" << approximation_timer.time() <<  " seconds)" << std::endl;

      return true;
  }*/


  // Poisson surface reconstruction
  // This variant requires all parameters.
  template <class SparseLinearAlgebraTraits_d,
            class Visitor>
  bool compute_poisson_implicit_function(
                                 SparseLinearAlgebraTraits_d solver,// = SparseLinearAlgebraTraits_d(),
                                 Visitor visitor,
                                 double approximation_ratio = 0,
                                 double average_spacing_ratio = 5) 
  {
    
    first_delaunay_refinement(visitor);
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
      return compute_poisson_implicit_function<SparseLinearAlgebraTraits_d,Implicit_visitor>(solver,Implicit_visitor(),0.02,5);
    else
      return compute_poisson_implicit_function<SparseLinearAlgebraTraits_d,Implicit_visitor>(solver,Implicit_visitor());
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

   
  // Spectral Surface Reconstruction
  // This variant requires all parameters.
  template <class Visitor>
  bool compute_spectral_implicit_function(
                                 Visitor visitor,
                                 double bilaplacian = 1,
                                 double laplacian = 0.1,
                                 double fitting = 1.,
                                 double ratio = 10., 
                                 int mode = 0,
                                 double approximation_ratio = 0,
                                 double average_spacing_ratio = 5) 
  {
    
    first_delaunay_refinement(visitor);
    CGAL::Timer task_timer; task_timer.start();

    // Computes the Implicit indicator function operator()
    // at each vertex of the triangulation.
    if ( ! solve_spectral(bilaplacian, laplacian, fitting, ratio, mode) )
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
  /// \endcond

  /*!
    This function must be called after the
    insertion of oriented points. It computes the piecewise linear scalar
    function operator() by: applying Delaunay refinement, solving for
    operator() at each vertex of the triangulation with a sparse linear
    solver, and shifting and orienting operator() such that it is 0 at all
    input points and negative inside the inferred surface.

    \param bilaplacian bilaplacian term weight
    \param laplacian laplacian term weight
    \param fitting data fitting coefficient
    \param ratio reliability coefficient
    \param mode choose the using formulation
    \param smoother_hole_filling controls if the Delaunay refinement is done for the input points, or for an approximation of the surface obtained from a first pass of the algorithm on a sample of the points.

    \return `false` if the solver fails. 
  */ 
  bool compute_spectral_implicit_function(double bilaplacian = 1, double laplacian = 0.1, 
                                 double fitting = 1., double ratio = 3.,
                                 int mode = 0, bool smoother_hole_filling = false)
  {
    if (smoother_hole_filling)
      return compute_spectral_implicit_function<Implicit_visitor>(Implicit_visitor(), bilaplacian, laplacian, fitting, ratio, mode, 0.02, 5);
    else
      return compute_spectral_implicit_function<Implicit_visitor>(Implicit_visitor(), bilaplacian, laplacian, fitting, ratio, mode);
  }

  /// \endcond

  /*! 
    `ImplicitFunction` interface: evaluates the implicit function at a 
    given 3D query point. The function `compute_implicit_function()` must be 
    called before the first call to `operator()`. 
  */ 
  FT operator()(const Point& p) const
  {
    m_hint = m_tr->locate(p ,m_hint); 

    if(m_tr->is_infinite(m_hint)) {
      int i = m_hint->index(m_tr->infinite_vertex());
      return m_hint->vertex((i+1)&3)->f();
    }

    FT a,b,c,d;
    barycentric_coordinates(p,m_hint,a,b,c,d);
    return a * m_hint->vertex(0)->f() +
           b * m_hint->vertex(1)->f() +
           c * m_hint->vertex(2)->f() +
           d * m_hint->vertex(3)->f();
  }

  boost::tuple<FT, Cell_handle, bool> special_func(const Point& p) const
  {
    m_hint = m_tr->locate(p  ,m_hint  ); // no hint when we use hierarchy

    if(m_tr->is_infinite(m_hint)) {
      int i = m_hint->index(m_tr->infinite_vertex());
      return boost::make_tuple(m_hint->vertex((i+1)&3)->f(),
                               m_hint, true);
    }

    FT a,b,c,d;
    barycentric_coordinates(p,m_hint,a,b,c,d);
    return boost::make_tuple(a * m_hint->vertex(0)->f() +
                             b * m_hint->vertex(1)->f() +
                             c * m_hint->vertex(2)->f() +
                             d * m_hint->vertex(3)->f(),
                             m_hint, false);
  }
  
  /// \cond SKIP_IN_MANUAL
  void initialize_cell_indices()
  {
    int i=0;
    for(Finite_cells_iterator fcit = m_tr->finite_cells_begin();
        fcit != m_tr->finite_cells_end();
        ++fcit){
      fcit->info()= i++;
    }
  }

  void initialize_barycenters() const
  {
    m_Bary->resize(m_tr->number_of_finite_cells());

    for(std::size_t i=0; i< m_Bary->size();i++){
      (*m_Bary)[i][0]=-1;
    }
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
    boost::array<double,9> & entry = (*m_Bary)[ch->info()];
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
           entry[0],entry[1],entry[2],entry[3],entry[4],entry[5],entry[6],entry[7],entry[8]);
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
    m_tr->index_unconstrained_vertices();
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
        assemble_poisson_row<SparseLinearAlgebraTraits_d>(A,v,B,lambda);
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
  bool solve_spectral(
    double bilaplacian, double laplacian,
    double fitting, double ratio, int mode)
  {
    CGAL_TRACE("Calls solve_spectral()\n");

    double time_init = clock();

    double duration_assembly = 0.0;
    double duration_solve = 0.0;

    initialize_cell_indices();
    initialize_barycenters();

    // get #variables
    //constrain_one_vertex_on_convex_hull();
    m_tr->index_all_vertices();
    const int nb_variables = static_cast<int>(m_tr->number_of_vertices());
    const int nb_input_vertices = m_tr->nb_input_vertices();
  	CGAL_TRACE("  %d input vertices out of %d\n", nb_input_vertices, nb_variables);

    // Assemble isotropic laplacian matrix A
    Matrix AA(nb_variables), L(nb_variables), F(nb_variables); // matrix is symmetric definite positive
    Matrix V(nb_variables), V_inv(nb_variables);
    ESMatrix B(nb_variables, nb_variables);
    EMatrix X(nb_variables, 1);

    initialize_duals();

    CGAL_TRACE("  Begin calculation: (%.2lf s)\n", (clock() - time_init)/CLOCKS_PER_SEC);
    Finite_vertices_iterator v, e; 
    double duration_cal = 0., duration_assign= 0.; 
    for(v = m_tr->finite_vertices_begin(), e = m_tr->finite_vertices_end();
        v != e;
        ++v)
        assemble_spectral_row(v, AA, L, F, V, V_inv, duration_assign, duration_cal, fitting, ratio, mode);

    CGAL_TRACE("  Calculate elem: total (%.2lf s)\n", duration_cal/CLOCKS_PER_SEC);
    CGAL_TRACE("  Assign: total (%.2lf s)\n", duration_assign/CLOCKS_PER_SEC);

    double time_b = clock();

    ESMatrix EL = L.eigen_object(), EA = AA.eigen_object();
    ESMatrix EV = V.eigen_object(), EV_inv = V_inv.eigen_object();
    
    const FT radius = sqrt(bounding_sphere().squared_radius()); // get triangulation's radius

    EL = EL / radius;
    EA = EA / radius;
    EV = EV / ::pow(radius, 3);
    EV_inv = EV_inv * ::pow(radius, 3);

    B = EL * EV_inv * EL * bilaplacian + EL * laplacian + EV * F.eigen_object();
    
    clear_duals();
    duration_assembly = (clock() - time_init)/CLOCKS_PER_SEC;
    CGAL_TRACE("  Creates matrix: done (%.2lf s)\n", duration_assembly);

    CGAL_TRACE("  Solve generalized eigenvalue problem...\n");

    // Solve generalized eigenvalue problem
    time_init = clock();
    spectral_solver<ESMatrix, EMatrix, Spectra::LARGEST_ALGE>(EA, B, EL, X);
    duration_solve = (clock() - time_init)/CLOCKS_PER_SEC;

    CGAL_TRACE("  Solve generalized eigenvalue problem: done (%.2lf s)\n", duration_solve);

    // copy function's values to vertices
    unsigned int index = 0;
    for (v = m_tr->finite_vertices_begin(), e = m_tr->finite_vertices_end(); v!= e; ++v)
      v->f() = X(index++, 0);

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
      OpType op(A);
      BOpType Bop(B);
      // Make sure B is positive definite and the decompoition is successful
      assert(Bop.info() == Spectra::SUCCESSFUL);

      Spectra::SymGEigsSolver<FT, SelectionRule, OpType, BOpType, Spectra::GEIGS_CHOLESKY> eigs(&op, &Bop, k, m);
      eigs.init();
      int nconv = eigs.compute(200); // maxit = 200 to reduce running time for failed cases 
      X = eigs.eigenvectors();
      auto lambd = eigs.eigenvalues()[0];
      auto xtax = X.transpose() * A * X;
      auto xtlx = X.transpose() * L * X;
      auto xtltlx = X.transpose() * L.transpose() * L * X;

      if(eigs.info() != Spectra::SUCCESSFUL)
        CGAL_TRACE("  Spectra failed! %d", eigs.info());

      /*
      // test
      std::cerr << "    lambda:" << lambd << std::endl;
      std::cerr << "    xtax  :" << xtax << std::endl;
      std::cerr << "    xtlx  :" << xtlx << std::endl;
      std::cerr << "    xtltlx:" << xtltlx << std::endl;
      */
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
    m_tr ->adjacent_vertices (m_tr->infinite_vertex (),
			     std::back_inserter (convex_hull));
    unsigned int nb_negative = 0;
    for (std::size_t i = 0; i < convex_hull.size (); ++ i)
      if (convex_hull[i]->f() < 0.0)
        ++ nb_negative;
    
    if(nb_negative > convex_hull.size () / 2)
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
        e= m_tr->finite_vertices_end();
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
    std::size_t index = size/2;
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
    const boost::array<double,9> & i = (*m_Bary)[cell->info()];
    if(i[0]==-1){
      initialize_matrix_entry(cell);
    }
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
        e= m_tr->finite_vertices_end();
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
        v!= e;
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
      Vector nn = (index%2==0) ? CGAL::cross_product(b-a,c-a) : CGAL::cross_product(c-a,b-a);
      nn = nn / std::sqrt(nn*nn); // normalize
      Vector p = a - x;
      Vector q = b - x;
      Vector r = c - x;
      FT p_n = std::sqrt(p*p);
      FT q_n = std::sqrt(q*q);
      FT r_n = std::sqrt(r*r);
      FT solid_angle = p*(CGAL::cross_product(q,r));
      solid_angle = std::abs(solid_angle / (p_n*q_n*r_n + (p*q)*r_n + (q*r)*p_n + (r*p)*q_n));

      FT area = std::sqrt(squared_area(a,b,c));
      FT length = p_n + q_n + r_n;
      div += n * nn * area / length ;
    }
    return div * FT(3.0);
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
      const Vector nn = CGAL::cross_product(b-a,c-a);

      div+= nn * (//v->normal() + 
                  cell->vertex((index+1)%4)->normal() +
                  cell->vertex((index+2)%4)->normal() +
                  cell->vertex((index+3)%4)->normal());
    }
    return div;
  }

  Vector get_cell_normal(Cell_handle cell)
  {
    return Normal[cell->info()];
  }

  Vector cell_normal(Cell_handle cell) const
  {
    const Vector& n0 = cell->vertex(0)->normal();
    const Vector& n1 = cell->vertex(1)->normal();
    const Vector& n2 = cell->vertex(2)->normal();
    const Vector& n3 = cell->vertex(3)->normal();
    Vector n = n0 + n1 + n2 + n3;
    if(n != NULL_VECTOR){
      FT sq_norm = n*n;
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

  // anisotropic Laplace coefficient
  FT mcotan_dot_new(Edge& edge, const FT cij, const FT ratio, const bool convert = true, const bool inverse = false)
  {
    Cell_handle cell = edge.first;
    Vertex_handle vi = cell->vertex(edge.second);
    Vertex_handle vj = cell->vertex(edge.third);

    // primal edge
    const Point& pi = vi->point();
    const Point& pj = vj->point();
    Vector primal = pj - pi;
    FT len_primal = std::sqrt(primal * primal);

    // find normals
    Vector na = vi->normal();
		Vector nb = vj->normal();
    if(na * nb < 0.0)
			na = -na;

    // should use covariance to check isotropic
    Covariance ca(pi, na, ratio), cb(pj, nb, ratio);
    Covariance cab(ca, cb, convert);
    if(cab.isotropic()) return cij;
    
    // average normals
		FT dot = cab.ut_c_v(primal, primal);

		//return cij * len_primal / dot;
    return cij / dot;
  }

  // anisotropic Laplace coefficient
  FT mcotan_dot(Edge& edge, const FT cij, const FT ratio, const bool convert, const bool inverse = false)
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
    Vector na = vi->normal();
		Vector nb = vj->normal();
    if(na * nb < 0.0)
			na = -na;

    // should use covariance to check isotropic
    Covariance ca(pi, na, ratio), cb(pj, nb, ratio);
    //Covariance ca(na, ratio), cb(nb, ratio);
    Covariance cab(ca, cb, convert);
    if(cab.isotropic()) return cij;
    
    // average normals
    FT dot = cab.ut_c_v(primal, primal);

		if(inverse)
			dot = 1.0 - dot;

		return cij * dot;
  }

  // anisotropic Laplace coefficient
  FT mcotan_dot_in_metric(Edge& edge, const FT cij, const FT ratio, const bool convert, const bool inverse = false)
  {
    Cell_handle cell = edge.first;
    Vertex_handle vi = cell->vertex(edge.second);
    Vertex_handle vj = cell->vertex(edge.third);

    // primal edge
    const Point& pi = vi->point();
    const Point& pj = vj->point();
    Vector primal = pj - pi;

    // find normals
    Vector na = vi->normal();
		Vector nb = vj->normal();
    if(na * nb < 0.0)
			na = -na;

    // should use covariance to check isotropic
    Covariance ca(pi, na, ratio), cb(pj, nb, ratio);
    Covariance cab(ca, cb, convert);
    if(cab.isotropic()) return cij;
    
    // calculate the voronoi area in a metric
    FT cell_area = area_voronoi_face_in_metric(edge, cab);
    FT len_primal = std::sqrt(cab.ut_c_v(primal, primal));

		return cell_area / len_primal;
  }

  FT mcotan_dot_2007(Edge& edge, const FT cij, const FT ratio, const bool convert = true, const bool inverse = false)
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
    Vector na = vi->normal();
		Vector nb = vj->normal();
    if(na * nb < 0.0)
			na = -na;
    Vector n = na + nb;
    FT sqnorm = n * n;

    if(sqnorm == 0.0)
      return cij;

    // should use covariance to check isotropic
    Covariance ca(pi, na, ratio), cb(pj, nb, ratio);
    //Covariance ca(na, ratio), cb(nb, ratio);
    Covariance cab(ca, cb, convert);
    if(cab.isotropic()) return cij;

    n = n / std::sqrt(sqnorm);

    double dot = std::fabs(n * primal);
    return ratio * cij * ::pow(dot, 2);
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
    for(std::size_t i=1;i<nb_triangles;i++)
    {
      const Point& b = voronoi_points[i];
      const Point& c = voronoi_points[i+1];
      area += std::sqrt(squared_area(a,b,c));
    }
    return area;
  }

  // spin around edge
  // return area(voronoi face) in a specific metric
  FT area_voronoi_face_in_metric(Edge& edge, Covariance& cab)
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
        return area_voronoi_face_boundary_in_metric(edge, cab);
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
      area += std::sqrt(squared_area_in_metric(a,b,c,cab));
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
        int k =  Triangulation_utils_3::next_around_edge(i,j);
        int l =  Triangulation_utils_3::next_around_edge(j,i);

        Vertex_handle vk = cell->vertex(k);
        Vertex_handle vl = cell->vertex(l);

        const Point& pk = vk->point();
        const Point& pl = vl->point();

        // if circumcenter is outside tet
        // pick barycenter instead
        if(tet.has_on_unbounded_side(c))
        {
          Point cell_points[4] = {pi,pj,pk,pl};
          c = CGAL::centroid(cell_points, cell_points+4);
        }

        Point ck = CGAL::circumcenter(pi,pj,pk);
        Point cl = CGAL::circumcenter(pi,pj,pl);

        area += std::sqrt(squared_area(m,c,ck));
        area += std::sqrt(squared_area(m,c,cl));
      }
      circ++;
    }
    while(circ != done);
    return area;
  }

  // approximate area when a cell is infinite
  FT area_voronoi_face_boundary_in_metric(Edge& edge, Covariance& cab)
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
        int k =  Triangulation_utils_3::next_around_edge(i,j);
        int l =  Triangulation_utils_3::next_around_edge(j,i);

        Vertex_handle vk = cell->vertex(k);
        Vertex_handle vl = cell->vertex(l);

        const Point& pk = vk->point();
        const Point& pl = vl->point();

        // if circumcenter is outside tet
        // pick barycenter instead
        if(tet.has_on_unbounded_side(c))
        {
          Point cell_points[4] = {pi,pj,pk,pl};
          c = CGAL::centroid(cell_points, cell_points+4);
        }

        Point ck = CGAL::circumcenter(pi,pj,pk);
        Point cl = CGAL::circumcenter(pi,pj,pl);

        area += std::sqrt(squared_area_in_metric(m,c,ck,cab));
        area += std::sqrt(squared_area_in_metric(m,c,cl,cab));
      }
      circ++;
    }
    while(circ != done);
    return area;
  }


  /// Computes enlarged geometric bounding sphere of the embedded triangulation.
  Sphere enlarged_bounding_sphere(FT ratio) const
  {
    Sphere bsphere = bounding_sphere(); // triangulation's bounding sphere
    return Sphere(bsphere.center(), bsphere.squared_radius() * ratio*ratio);
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
      Cell_handle cell = *it;
      if(m_tr->is_infinite(cell))
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
    return 0.25 * total_volume;
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
        Edge edge( it->first, it->first->index(vi), it->first->index(vj));
        if(vi->index() < vj->index()){
          std::swap(edge.second,  edge.third);
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
  void assemble_spectral_row(Vertex_handle vi, Matrix& AA, 
                             Matrix& L, Matrix& F, Matrix& V, Matrix& V_inv,
                             FT& duration_assign, FT& duration_cal,
                             const FT fitting, 
                             const FT ratio, 
                             const int mode)
  {
    // for each vertex vj neighbor of vi
    std::vector<Edge> edges;
    m_tr->incident_edges(vi,std::back_inserter(edges));

    double diagonal = 0.0;
    double mdiagonal = 0.0;
    double time_init;

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
        Edge edge( it->first, it->first->index(vi), it->first->index(vj));

        time_init = clock();

        if(vi->index() < vj->index()){
          std::swap(edge.second,  edge.third);
        }

        FT cij = cotan_geometric(edge);

        bool convert = true;
        FT mcij;
        switch(mode) {
          case 0: mcij = mcotan_dot(edge, cij, ratio, convert);
          case 1: mcij = mcotan_dot_new(edge, cij, ratio, convert);
          default: mcij = mcotan_dot_in_metric(edge, cij, ratio, convert);
        }
       
        duration_cal += clock() - time_init; time_init = clock();


        AA.set_coef(vi->index(), vj->index(), -mcij, true);
        L.set_coef(vi->index(), vj->index(), -cij, true);

        duration_assign += clock() - time_init;

        diagonal += cij;
        mdiagonal += mcij;
      }
    // diagonal coefficient

    const FT vol = volume_voronoi_cell(vi);

    time_init = clock();
    AA.set_coef(vi->index(),vi->index(), mdiagonal, true);
    L.set_coef(vi->index(),vi->index(), diagonal, true);
    V.set_coef(vi->index(), vi->index(), vol, true);
    V_inv.set_coef(vi->index(), vi->index(), std::min(1.0 / vol, 1e7), true);
    
    if (vi->type() == Triangulation::INPUT)
      F.set_coef(vi->index(),vi->index(), fitting, true);
      
     duration_assign += clock() - time_init;
  }


public:

  // Write function value to ply file (for testing the algorithm)
  bool write_func_to_ply(const std::string outfile){
    std::vector<Point_with_property> my_pts;

    Finite_vertices_iterator v, e;
    for(v = m_tr->finite_vertices_begin(),
        e= m_tr->finite_vertices_end();
        v != e; 
        v++)
        my_pts.push_back(CGAL::cpp11::make_tuple(v->point(), v->f()));

    std::string outname = outfile.substr(0, outfile.find_last_of('.'));

    std::ofstream f("func_" + outname + ".ply");
    CGAL::write_ply_points_with_properties(f, my_pts, CGAL::make_ply_point_writer(PP_point_map()),
            std::make_pair(PP_func_map(), CGAL::PLY_property<FT>("function_value")));

    return true;
  }

  /// Marching Tetrahedra
  unsigned int marching_tetrahedra(const FT value, const std::string outfile)
  {
    std::vector<Point> points;
    std::vector< std::vector<std::size_t> > polygons;
    std::ofstream out("iso_facet_" + outfile);

    return m_tr->marching_tets(value, out, points, polygons);
  }

}; // end of Implicit_reconstruction_function


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_IMPLICIT_RECONSTRUCTION_FUNCTION_H
