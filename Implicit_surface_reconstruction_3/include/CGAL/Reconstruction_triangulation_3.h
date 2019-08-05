// Copyright (c) 2007  INRIA (France).
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
//
// Author(s)     : Laurent Saboret, Pierre Alliez, Tong Zhao, Cédric Portaneri


#ifndef CGAL_IMPLICIT_FCT_DELAUNAY_TRIANGULATION_H
#define CGAL_IMPLICIT_FCT_DELAUNAY_TRIANGULATION_H

#include <CGAL/license/Implicit_surface_reconstruction_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Point_with_normal_3.h>
#include <CGAL/Lightweight_vector_3.h>
#include <CGAL/property_map.h>
#include <CGAL/surface_reconstruction_points_assertions.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_3.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <CGAL/algorithm.h>
#include <CGAL/bounding_box.h>
#include <boost/random/random_number_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>

//Mesh
#include <CGAL/Surface_mesh.h>

#include <vector>
#include <iterator>
#include <unordered_map>

namespace CGAL {

/// \internal
/// The Reconstruction_vertex_base_3 class is the default
/// vertex class of the Reconstruction_triangulation_3 class.
///
/// It provides the interface requested by the Implicit_reconstruction_function class:
/// - Each vertex stores a normal vector.
/// - A vertex is either an input point or a Steiner point added by Delaunay refinement.
/// - In order to solve a linear system over the triangulation, a vertex may be constrained
///   or not (i.e. may contribute to the right or left member of the linear system),
///   and has a unique index.
///
/// @param Gt   Geometric traits class / Point_3 is a typedef to Point_with_normal_3.
/// @param Cb   Vertex base class, model of TriangulationVertexBase_3.

template < class Gt,
           class PointRange,
           class Vb = Triangulation_vertex_base_3<Gt> >
class Reconstruction_vertex_base_3 : public Vb
{
// Public types
public:

  /// Geometric traits class / Point_3 is a typedef to Point_with_normal_3.
  typedef Gt Geom_traits;

  // Repeat Triangulation_vertex_base_3 public types
  /// \cond SKIP_IN_MANUAL
  typedef typename Vb::Cell_handle Cell_handle;
  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other                       Vb2;
    typedef Reconstruction_vertex_base_3<Geom_traits, PointRange, Vb2> Other;
  };
  /// \endcond

  // Geometric types
  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Vector_3 Vector;           ///< typedef to Vector_3
  typedef typename Geom_traits::Point_3 Point;             ///< typedef to Point_with_normal_3
  typedef typename PointRange::const_iterator InputIterator;

// data members
private:

  // TODO: reduce memory footprint
  FT m_f; // value of the implicit function // float precise enough?
  //bool m_constrained; // is vertex constrained? // combine constrained and type
  unsigned char m_type; // INPUT or STEINER
  unsigned char m_position; // INSIDE or BOUNDARY
  unsigned int m_index; // index in matrix (to be stored outside)
  InputIterator m_iter; // the associated PointRange::const_iterator
  FT m_lf;
  FT m_af;
  FT m_bf;
  FT m_check;

// Public methods
public:

  Reconstruction_vertex_base_3()
    : Vb(), m_f(FT(0.0)), m_type(0), m_index(0)
  {}

  Reconstruction_vertex_base_3(const Point& p)
    : Vb(p), m_f(FT(0.0)), m_type(0), m_index(0)
  {}

  Reconstruction_vertex_base_3(const Point& p, Cell_handle c)
    : Vb(p, c), m_f(FT(0.0)), m_type(0), m_index(0)
  {}

  Reconstruction_vertex_base_3(Cell_handle c)
    : Vb(c), m_f(FT(0.0)), m_type(0), m_index(0)
  {}


  /// Gets/sets the value of the implicit function.
  /// Default value is 0.0.
  FT  f() const { return m_f; }
  FT& f()       { return m_f; }

  FT  lf() const { return m_lf; }
  FT& lf()       { return m_lf; }

  FT  bf() const { return m_bf; }
  FT& bf()       { return m_bf; }

  FT  af() const { return m_af; }
  FT& af()       { return m_af; }

  FT  check() const { return m_check; }
  FT& check()       { return m_check; }

  /// Gets/sets the type = INPUT or STEINER.
  unsigned char  type() const { return m_type; }
  unsigned char& type()       { return m_type; }

  /// Gets/sets the type = INSIDE or BOUNDARY.
  unsigned char  position() const { return m_position; }
  unsigned char& position()       { return m_position; }

  /// Gets/sets the index in matrix.
  unsigned int  index() const { return m_index; }
  unsigned int& index()       { return m_index; }

  InputIterator input_iterator() const { return m_iter; }
  InputIterator& input_iterator()      { return m_iter; }

  /*
  /// Gets/sets normal vector.
  /// Default value is null vector.
  const Vector& normal() const { return this->point().normal(); }
  Vector&       normal()       { return this->point().normal(); }
  */

// Private methods
private:

    /// Copy constructor and operator =() are not implemented.
    Reconstruction_vertex_base_3(const Reconstruction_vertex_base_3& toCopy);
    Reconstruction_vertex_base_3& operator =(const Reconstruction_vertex_base_3& toCopy);

}; // end of Reconstruction_vertex_base_3


/// \internal
/// Same as Triangulation_cell_base_with_info_3 
/// but with circumcenter() from Delaunay_triangulation_cell_base_3 
///
/// Added because of the circumcenter() removal from Triangulation_cell_base_3
///
/// @param Info_ Type of the reference to the object stored in the cell
/// @param Gt    Geometric traits class / Point_3 is a typedef to Point_with_normal_3.
/// @param Cb    Vertex base class, model of TriangulationVertexBase_3.

template < typename Info_, typename GT,
           typename Cb = Delaunay_triangulation_cell_base_3<GT> >
class Delaunay_triangulation_cell_base_with_info_3
  : public Cb
{
  Info_ _info;
public:
  typedef typename Cb::Vertex_handle                   Vertex_handle;
  typedef typename Cb::Cell_handle                     Cell_handle;
  typedef Info_                                        Info;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Cb::template Rebind_TDS<TDS2>::Other       Cb2;
    typedef Delaunay_triangulation_cell_base_with_info_3<Info, GT, Cb2>  Other;
  };

  Delaunay_triangulation_cell_base_with_info_3()
    : Cb() {}

  Delaunay_triangulation_cell_base_with_info_3(Vertex_handle v0, Vertex_handle v1,
                                      Vertex_handle v2, Vertex_handle v3)
    : Cb(v0, v1, v2, v3) {}

  Delaunay_triangulation_cell_base_with_info_3(Vertex_handle v0, Vertex_handle v1,
                                      Vertex_handle v2, Vertex_handle v3,
                                      Cell_handle   n0, Cell_handle   n1,
                                      Cell_handle   n2, Cell_handle   n3)
    : Cb(v0, v1, v2, v3, n0, n1, n2, n3) {}

  const Info& info() const { return _info; }
  Info&       info()       { return _info; }
};

/*
/// \internal
/// Helper class:
/// Reconstruction_triangulation_default_geom_traits_3
/// changes in a geometric traits class the Point_3 type to
/// Point_with_normal_3<BaseGt>.
///
/// @param BaseGt   Geometric traits class.
template <class BaseGt>
struct Reconstruction_triangulation_default_geom_traits_3 : public BaseGt
{
  typedef Point_with_normal_3<BaseGt> Point_3;
};*/


/// \internal
/// The Reconstruction_triangulation_3 class
/// provides the interface requested by the Implicit_reconstruction_function class:
/// - A vertex is either an input point or a Steiner point added by Delaunay refinement.
/// - In order to solve a linear system over the triangulation, a vertex may be constrained
///   or not (i.e. may contribute to the right or left member of the linear system),
///   and has a unique index.
/// The vertex class must derive from Reconstruction_vertex_base_3.
///
/// @param BaseGt   Geometric traits class.
/// @param Gt       Geometric traits class / Point_3 is a typedef to Point_with_normal_3<BaseGt>.
/// @param Tds      Model of TriangulationDataStructure_3. The vertex class
///                 must derive from Reconstruction_vertex_base_3.

template <class Gt,
          class PointRange,
          class NormalMap,
          class Tds_ = Triangulation_data_structure_3<Reconstruction_vertex_base_3<Gt, PointRange>, Delaunay_triangulation_cell_base_with_info_3<int, Gt> > >
class Reconstruction_triangulation_3 : public Delaunay_triangulation_3<Gt, Tds_>
{
// Private types
private:

  // Base class
  typedef Delaunay_triangulation_3<Gt, Tds_>  Base;

  // Auxiliary class to build an iterator over input points.
  class Is_steiner_point
  {
  public:
      typedef typename Base::Finite_vertices_iterator Finite_vertices_iterator;

      bool operator()(const Finite_vertices_iterator& v) const
      {
        return (v->type() == Reconstruction_triangulation_3::STEINER);
      }
  };

// Public types
public:

  /// Geometric traits class / Point_3 is a typedef to Point_with_normal_3<BaseGt>.
  typedef Gt  Geom_traits;

  // Repeat base class' types
  /// \cond SKIP_IN_MANUAL
  typedef Tds_ Triangulation_data_structure;
  typedef typename Base::Segment      Segment;
  typedef typename Base::Triangle     Triangle;
  typedef typename Base::Tetrahedron  Tetrahedron;
  typedef typename Base::Line         Line;
  typedef typename Base::Ray          Ray;
  typedef typename Base::Object       Object;
  typedef typename Base::Cell_handle   Cell_handle;
  typedef typename Base::Vertex_handle Vertex_handle;
  typedef typename Base::Cell   Cell;
  typedef typename Base::Vertex Vertex;
  typedef typename Base::Facet  Facet;
  typedef typename Base::Edge   Edge;
  typedef typename Base::Cell_circulator  Cell_circulator;
  typedef typename Base::Facet_circulator Facet_circulator;
  typedef typename Base::Cell_iterator    Cell_iterator;
  typedef typename Base::Facet_iterator   Facet_iterator;
  typedef typename Base::Edge_iterator    Edge_iterator;
  typedef typename Base::Vertex_iterator  Vertex_iterator;
  //typedef typename Base::Point_iterator Point_iterator;
  typedef typename Base::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Base::Finite_cells_iterator    Finite_cells_iterator;
  typedef typename Base::Finite_facets_iterator   Finite_facets_iterator;
  typedef typename Base::Finite_edges_iterator    Finite_edges_iterator;
  typedef typename Base::All_cells_iterator       All_cells_iterator;
  //typedef typename Base::All_vertices_iterator       All_vertices_iterator;
  typedef typename Base::Locate_type Locate_type;
  /// \endcond

  // Geometric types
  typedef typename PointRange::const_iterator InputIterator;
  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Vector_3 Vector; ///< typedef to Vector_3<BaseGt>
  typedef typename Geom_traits::Point_3 Point;  ///< typedef to Point_with_normal_3<BaseGt>
  //typedef typename Geom_traits::Point_3 Point_with_normal; ///< Point_with_normal_3<BaseGt>
  typedef typename std::pair<Point, InputIterator> Point_with_iterator;
  typedef typename Geom_traits::Sphere_3 Sphere;
  typedef typename Geom_traits::Iso_cuboid_3 Iso_cuboid;

  typedef CGAL::Polyhedron_3<Geom_traits> Polyhedron;
  typedef CGAL::Spatial_sort_traits_adapter_3<Geom_traits, CGAL::First_of_pair_property_map<Point_with_iterator> > Search_traits_3;

  //Mesh
  typedef CGAL::Surface_mesh<Point>                           Mesh;

  // Marching Tet
  struct HashEdgePair {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1, T2> &p) const {
        size_t hash = 0;
        boost::hash_combine(hash, p.first);
        boost::hash_combine(hash, p.second);

        return hash;  
    }
  };

  typedef std::pair<unsigned int, unsigned int>               Edge_pair;
  typedef std::unordered_map<Edge_pair, size_t, HashEdgePair> Edge_hash_map;

  /// Point type
  enum Point_type {
    INPUT=0,    ///< Input point.
    STEINER=1   ///< Steiner point created by Delaunay refinement.
  };

  enum Position_type{
    INSIDE=0,
    BOUNDARY=1
  };

  /// Iterator over input vertices.
  typedef Filter_iterator<Finite_vertices_iterator, Is_steiner_point>
                                                    Input_vertices_iterator;

  /// Iterator over input points.
  typedef Iterator_project<Input_vertices_iterator,
                           Project_point<Vertex> >  Input_point_iterator;

  mutable Sphere sphere;
  std::vector<Point_with_iterator> points;
  std::size_t fraction;
  std::list<double> fractions;
  Vertex_handle constrained_vertex;
  NormalMap normals;


public:

  /// Default constructor.
  Reconstruction_triangulation_3()
  {}

  ~Reconstruction_triangulation_3()
  {}

  // Default copy constructor and operator =() are fine.

  // Repeat base class' public methods used below
  /// \cond SKIP_IN_MANUAL
  using Base::points_begin;
  using Base::points_end;
  using Base::number_of_vertices;
  using Base::number_of_finite_cells;
  using Base::number_of_facets;
  using Base::finite_vertices_begin;
  using Base::finite_vertices_end;
  using Base::finite_edges_begin;
  using Base::finite_edges_end;
  using Base::all_vertices_begin;
  using Base::all_vertices_end;

  using Base::geom_traits;
  /// \endcond

  /// Gets first iterator over input vertices.
  Input_vertices_iterator input_vertices_begin() const
  {
      return Input_vertices_iterator(finite_vertices_end(), Is_steiner_point(),
                                     finite_vertices_begin());
  }
  /// Gets past-the-end iterator over input vertices.
  Input_vertices_iterator input_vertices_end() const
  {
      return Input_vertices_iterator(finite_vertices_end(), Is_steiner_point());
  }

  /// Gets iterator over the first input point.
  Input_point_iterator input_points_begin() const
  {
      return Input_point_iterator(input_vertices_begin());
  }
  /// Gets past-the-end iterator over the input points.
  Input_point_iterator input_points_end() const
  {
      return Input_point_iterator(input_vertices_end());
  }

  /// Gets the bounding sphere of input points.


  Sphere bounding_sphere() const
  {
    return sphere;
  }

  void initialize_bounding_sphere() const
  {
    boost::function<Point(Point_with_iterator)> f = boost::bind(&Point_with_iterator::first, _1);

    Iso_cuboid ic = bounding_box(boost::make_transform_iterator(points.begin(), f), 
                                 boost::make_transform_iterator(points.end(), f));
    Point center = midpoint((ic.min)(), (ic.max)());
    sphere = Sphere(center, squared_distance(center, (ic.max)()));
  }

  /// Insert point in the triangulation.
  /// Default type is STEINER.
  template <typename Visitor>
  Vertex_handle insert(const Point& p,
                       Point_type type,// = INPUT,
                       Cell_handle start,// = Cell_handle(),
                       Visitor visitor)
  {

    if(type == INPUT){
      visitor.before_insertion();
    }
    if(this->dimension() < 3){
      Vertex_handle v = Base::insert(p, start);
      v->type() = static_cast<unsigned char>(type);
      return v;
    }
    typename Base::Locate_type lt;
    int li, lj;
    Cell_handle ch = Base::locate(p, lt, li, lj, start);

    Vertex_handle v = Base::insert(p, lt, ch, li, lj);
    v->type() = static_cast<unsigned char>(type);

    return v;
    
  }

  /// Insert point in the triangulation.
  /// Default type is INPUT.
  template <typename Visitor>
  Vertex_handle insert(const Point_with_iterator& p,
                       Point_type type,// = INPUT,
                       Cell_handle start,// = Cell_handle(),
                       Visitor visitor)
  {

    if(type == INPUT){
      visitor.before_insertion();
    }
    if(this->dimension() < 3){
      Vertex_handle v = Base::insert(p.first, start);
      v->type() = static_cast<unsigned char>(type);
      v->input_iterator() = p.second;
      return v;
    }
    typename Base::Locate_type lt;
    int li, lj;
    Cell_handle ch = Base::locate(p.first, lt, li, lj, start);

    Vertex_handle v = Base::insert(p.first, lt, ch, li, lj);
    v->type() = static_cast<unsigned char>(type);
    v->input_iterator() = p.second;

    return v;
    
  }

  /// Insert the [first, beyond) range of points in the triangulation using a spatial sort.
  /// Default type is INPUT.
  ///
  /// @commentheading Template Parameters:
  /// @param PointRange is a model of `Range`. The value type of
  ///        its iterator is the key type of the named parameter `point_map`.
  /// @param PointMap is a model of `ReadablePropertyMap` with a value_type = Point_3.
  ///        It can be omitted if InputIterator value_type is convertible to Point_3.
  ///
  /// @return the number of inserted points.

  // This variant requires all parameters.
  template <typename PointMap,
            typename Visitor
  >
  int insert(
    PointRange& pts, ///< input point range
    PointMap point_map, ///< property map: `value_type of InputIterator` -> `Point_3` (the position of an input point).
    Visitor visitor)
  {
    if(! points.empty()){
      std::cerr << "WARNING: not all points inserted yet" << std::endl;
    }
    // Convert input points to Point_with_normal_3
    for (InputIterator it = pts.begin(); it != pts.end(); ++it)
      points.push_back(std::make_pair(get(point_map, *it), it));
    
    std::size_t n = points.size();

    initialize_bounding_sphere();

    // typedef typename PointRange::iterator Iterator_traits;
    typedef typename PointRange::difference_type Diff_t;
    boost::rand48 random;
    boost::random_number_generator<boost::rand48, Diff_t> rng(random);
    CGAL::cpp98::random_shuffle (points.begin(), points.end(), rng);
    fraction = 0;

    fractions.clear();
    fractions.push_back(1.0);
    
    double m = static_cast<double>(n);
    
    while(m > 500){
      m /= 2;
      fractions.push_front(m / n);
    }
    
    insert_fraction(visitor);
    return 0;
  }

  template <typename Visitor>
  bool insert_fraction(Visitor visitor)
  {
    if(fractions.empty()){
      points.clear();
      return false;
    }
    double frac = fractions.front();
    fractions.pop_front();
    std::size_t more = (std::size_t)(points.size() * frac) - fraction;
    if((fraction + more) > points.size()){
      more = points.size() - fraction;
    }
    Cell_handle hint;

    // spatial sort
    Search_traits_3 traits;
    spatial_sort (points.begin() + fraction, points.begin() + fraction + more, traits);

    for (typename std::vector<Point_with_iterator>::const_iterator p = points.begin()+fraction;
         p != points.begin()+fraction+more; ++p)
    {
      Vertex_handle v = insert(*p, INPUT, hint, visitor);
      hint = v->cell();
    }
    fraction += more;
    return true;
  }

  /// \cond SKIP_IN_MANUAL
  // This variant creates a default point property map = Identity_property_map.
  template <typename Visitor>
  int insert(
    PointRange& pts, ///< input point range
    Visitor visitor)
  {
    return insert(
      pts,
      make_identity_property_map(
      typename std::iterator_traits<typename PointRange::iterator>::value_type()),
      visitor);
  }
  /// \endcond

  /// Delaunay refinement callback:
  /// insert STEINER point in the triangulation.
  template <class CellIt>
  Vertex_handle
  insert_in_hole(const Point& p, CellIt cell_begin, CellIt cell_end,
	         Cell_handle begin, int i, Point_type type = STEINER)
  {
      Vertex_handle v = Base::insert_in_hole(p, cell_begin, cell_end, begin, i);
      v->type() = static_cast<unsigned char>(type);
      return v;
  }

  Vertex_handle
  insert_in_edge(Cell_handle cell, int i, int j, Point_type type = STEINER)
  {
      Point pi = cell->vertex(i)->point();
      Point pj = cell->vertex(j)->point();
      Point mid = CGAL::midpoint(pi, pj);

      Vertex_handle v = Base::insert_in_edge(mid, cell, i, j);
      
      v->type() = static_cast<unsigned char>(type);
      v->position() = INSIDE;
      return v;
  }

  void intialize_normal(NormalMap normal_map) { normals = normal_map; }

  Vector normal(Vertex_handle v) const
  {
    if(v->type() == INPUT)
      return get(normals, *(v->input_iterator()));
    else
      return CGAL::NULL_VECTOR;
  }

  /*
  /// Index unconstrained vertices following the order of Finite_vertices_iterator.
  /// @return the number of unconstrained vertices.
  unsigned int index_unconstrained_vertices()
  {
    unsigned int index = 0;
    for (Finite_vertices_iterator v = finite_vertices_begin(),
         e = finite_vertices_end();
         v!= e;
         ++v)
    {
      if(! is_constrained(v))
        v->index() = index++;
    }
    return index;
  }*/

  unsigned int index_all_vertices(bool flag_constrained = false)
  {
    unsigned int index = 0;

    // index all unconstrained vertices
    for (Finite_vertices_iterator v = finite_vertices_begin(),
         e = finite_vertices_end();
         v!= e;
         ++v)
    {
      if(!flag_constrained)
        v->index() = index++;
      else if(!is_constrained(v))
        v->index() = index++;
    }

    // index the constrained vertices
    if(flag_constrained)
    {
      for (Finite_vertices_iterator v = finite_vertices_begin(),
         e = finite_vertices_end();
         v!= e;
         ++v)
         if(is_constrained(v))
          v->index() == index++;
    }
    return index;
  }

  /*
  unsigned int index_all_inside_vertices()
  {
    unsigned int index = 0;
    unsigned int iindex = 0;
    for (Finite_vertices_iterator v = finite_vertices_begin(),
         e = finite_vertices_end();
         v!= e;
         ++v)
    {
      v->index() = index++;
      if(v->position() == INSIDE)
        v->iindex() = iindex++;
    }
    return index;
  }*/

  unsigned int nb_input_vertices()
  {
	  unsigned int count = 0;
	  for (Finite_vertices_iterator v = finite_vertices_begin(),
		  e = finite_vertices_end();
		  v != e; v++)
		  if (v->type() == INPUT)
			  count++;
	
	  return count;
  }

  unsigned int nb_inside_vertices()
  {
	  unsigned int count = 0;
	  for (Finite_vertices_iterator v = finite_vertices_begin(),
		  e = finite_vertices_end();
		  v != e; v++)
		  if (v->position() == INSIDE)
			  count++;
	
	  return count;
  }

  /// Is vertex constrained, i.e.
  /// does it contribute to the right or left member of the linear system?

  bool is_constrained(Vertex_handle v) const
  {
    return v == constrained_vertex;
  }

  void constrain(Vertex_handle v)
  {
    constrained_vertex = v;
  }

  /// Marching Tets

  template <class Point_3, class Polygon_3>
  unsigned int marching_tets(const FT value,
                             Mesh&    mesh,
                             std::vector< Point_3 >&   m_contour_points,
                             std::vector< Polygon_3 >& m_contour_polygons)
  {
    Edge_hash_map m_edge_map;
    std::vector<Vector> m_directions;
    size_t num_points = find_level_point_set(value, m_edge_map, m_contour_points, m_directions);
    size_t num_faces  = find_level_polygon_set(m_edge_map, m_contour_points, m_directions, m_contour_polygons);
    m_edge_map.clear();
    m_directions.clear();

    bool flag_manifold = CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(m_contour_polygons);

    if(!flag_manifold) {
      std::cerr << "Marching Tetrahedron failed!" << std::endl;
      return 0;
    }
    else {
      CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(m_contour_points, m_contour_polygons, mesh);
      return num_faces;
    }  
  }

  template <class Point_3>
  size_t find_level_point_set(const FT              value,
                              Edge_hash_map&        m_edge_map,
                              std::vector<Point_3>& m_pts,
                              std::vector<Vector>&  m_dirs)
  {
    size_t index = 0;

    for(Finite_edges_iterator edge = this->finite_edges_begin();
        edge != this->finite_edges_end();
        ++edge)
    {
      Cell_handle cell = edge->first;     
      double v1 = cell->vertex(edge->second)->f() - value;
      double v2 = cell->vertex(edge->third)->f() - value;
      
      if((v1 * v2 <= 0) && !((std::abs(v1) < 1e-8) && (std::abs(v2) < 1e-8)))
      {
        const Point& p1 = cell->vertex(edge->second)->point();
        const Point& p2 = cell->vertex(edge->third)->point();

        unsigned int i1 = std::min(cell->vertex(edge->second)->index(), cell->vertex(edge->third)->index());
        unsigned int i2 = std::max(cell->vertex(edge->second)->index(), cell->vertex(edge->third)->index());

        if(std::abs(v1) < 1e-8) 
          m_pts.push_back(p1);
        else if(std::abs(v2) < 1e-8) 
          m_pts.push_back(p2);
        else if(v1 > 0 && v2 < 0) {
          double ratio = (0. - v1) / (v2 - v1);
          Point_3 p = p1 + ratio * (p2 - p1);
          m_pts.push_back(p);
        }
        else {
          double ratio = (0. - v2) / (v1 - v2);
          Point_3 p = p2 + ratio * (p1 - p2);
          m_pts.push_back(p);
        }
        // save hash pair
        m_edge_map.insert({std::make_pair(i1, i2), index++});
        // save direction
        Vector direction = (v1 > v2) ? p1 - p2 : p2 - p1;
        m_dirs.push_back(direction);
      }
    }

    return index;
  }

  template <class Point_3, class Polygon_3>
  size_t find_level_polygon_set(Edge_hash_map&          m_edge_map,
                                std::vector<Point_3>&   m_pts,
                                std::vector<Vector>&    m_dirs,
                                std::vector<Polygon_3>& m_polys)
  {
    size_t num_faces = 0;

    for(Finite_cells_iterator v = this->finite_cells_begin();
        v != this->finite_cells_end();
        ++v)
    {
      std::vector<size_t> cell_points;
      Vector direction;
      for(size_t i = 0; i < 3; i++)
        for(size_t j = i + 1; j < 4; j++)
        {
          size_t i1 = v->vertex(i)->index();
          size_t i2 = v->vertex(j)->index();

          typename Edge_hash_map::const_iterator got = m_edge_map.find(std::make_pair(std::min(i1, i2), std::max(i1, i2)));

          if (got != m_edge_map.end()) {
            cell_points.push_back(got->second);
            direction = m_dirs[got->second];
          } 
        }

      if(cell_points.size() == 3)
      {
        Vector u = m_pts[cell_points[1]] - m_pts[cell_points[0]];
        Vector v = m_pts[cell_points[2]] - m_pts[cell_points[0]];

        Vector n = CGAL::cross_product(u, v);

        if(n * direction >= 0) {
          Polygon_3 m_idx{cell_points[0], cell_points[1], cell_points[2]};
          m_polys.push_back(m_idx);
        } 
        else {
          Polygon_3 m_idx{cell_points[0], cell_points[2], cell_points[1]};
          m_polys.push_back(m_idx);        
        }

        num_faces += 1;
      }
      else if(cell_points.size() == 4)
      {
        Vector u = m_pts[cell_points[1]] - m_pts[cell_points[0]];
        Vector v = m_pts[cell_points[2]] - m_pts[cell_points[0]];

        Vector n = CGAL::cross_product(u, v);

        if(n * direction <= 0){
          Polygon_3 m_idx_1{cell_points[0], cell_points[2], cell_points[3]},
                    m_idx_2{cell_points[0], cell_points[3], cell_points[1]};
          m_polys.push_back(m_idx_1);
          m_polys.push_back(m_idx_2);
        }
        else{
          Polygon_3 m_idx_1{cell_points[0], cell_points[1], cell_points[3]},
                    m_idx_2{cell_points[0], cell_points[3], cell_points[2]};
          m_polys.push_back(m_idx_1);
          m_polys.push_back(m_idx_2);
        }

        num_faces += 2;
      } 
    }

    return num_faces;
  }


  // old version
  template <class Point_3, class Polygon_3>
  unsigned int marching_tets_old(const FT value, 
                             //const std::string outfile,
                             Mesh &mesh,
                             std::vector< Point_3 >& m_contour_points,
                             std::vector< Polygon_3 >& m_contour_polygons)
  {
    unsigned int nb_tri = 0;
    Finite_cells_iterator v, e; 

    for(v = this->finite_cells_begin(),
        e = this->finite_cells_end();
        v != e;
        ++v)
        nb_tri += contour(v, value, m_contour_points, m_contour_polygons);

    
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(m_contour_points, m_contour_polygons, mesh);
    //if (CGAL::is_closed(mesh) && (!CGAL::Polygon_mesh_processing::is_outward_oriented(mesh)))
    //  CGAL::Polygon_mesh_processing::reverse_face_orientations(mesh);

    return nb_tri;
  }

  template <class Point_3, class Polygon_3>
  unsigned int contour(Cell_handle cell, const FT value, 
                        std::vector< Point_3 >& m_pts, 
                        std::vector< Polygon_3 >& m_polys)
  {
    std::list<Point> cell_points;
    Vector direction;

    if(!extract_level_set_points(cell, value, cell_points, direction))
      return 0;

    if(cell_points.size() == 3)
    {
      typename std::list<Point>::iterator it = cell_points.begin();
      const Point& a = (*it); it++;
      const Point& b = (*it); it++;
      const Point& c = (*it);

      Vector n = CGAL::cross_product((b - a), (c - a));

      m_pts.push_back(a);
      m_pts.push_back(b);
      m_pts.push_back(c);

      if(n * direction >= 0){
        std::vector<std::size_t> m_idx{m_pts.size() - 3, m_pts.size() - 2, m_pts.size() - 1};
        m_polys.push_back(m_idx);
      } 
      else{
        std::vector<std::size_t> m_idx{m_pts.size() - 3, m_pts.size() - 1, m_pts.size() - 2};
        m_polys.push_back(m_idx);        
      }
      return 1;
    }
    else if(cell_points.size() == 4)
    {
      typename std::list<Point>::iterator it = cell_points.begin();
      std::vector<Point> p(4);
      for(int i = 0; i < 4; i++)
      {
        p[i] = (*it);
        it++;
      }
      // compute normal
      Vector u = p[1] - p[0];
      Vector v = p[2] - p[0];
      Vector n = CGAL::cross_product(u, v);

      m_pts.push_back(p[0]);
      m_pts.push_back(p[1]);
      m_pts.push_back(p[2]);
      m_pts.push_back(p[3]);

      if(n * direction <= 0){
        std::vector<std::size_t> m_idx_1{m_pts.size() - 4, m_pts.size() - 2, m_pts.size() - 1},
                                 m_idx_2{m_pts.size() - 4, m_pts.size() - 1, m_pts.size() - 3};
        m_polys.push_back(m_idx_1);
        m_polys.push_back(m_idx_2);
      }
      else{
        std::vector<std::size_t> m_idx_1{m_pts.size() - 4, m_pts.size() - 3, m_pts.size() - 1},
                                 m_idx_2{m_pts.size() - 4, m_pts.size() - 1, m_pts.size() - 2};
        m_polys.push_back(m_idx_1);
        m_polys.push_back(m_idx_2);
      }
      return 2;
    }
    return 0;   
  }

  bool extract_level_set_points(Cell_handle cell, const FT value, std::list<Point>& points, Vector& direction)
  {
    Point point;
    if(level_set(cell, value, 0, 1, point, direction)) points.push_back(point);
		if(level_set(cell, value, 0, 2, point, direction)) points.push_back(point);
		if(level_set(cell, value, 0, 3, point, direction)) points.push_back(point);
		if(level_set(cell, value, 1, 2, point, direction)) points.push_back(point);
		if(level_set(cell, value, 1, 3, point, direction)) points.push_back(point);
		if(level_set(cell, value, 2, 3, point, direction)) points.push_back(point);
		return points.size() != 0;
  }

  bool level_set(Cell_handle cell, const FT value, const int i1, const int i2, Point& p, Vector& direction)
  {
    const Point& p1 = cell->vertex(i1)->point();
    const Point& p2 = cell->vertex(i2)->point();
    double v1 = cell->vertex(i1)->f();
    double v2 = cell->vertex(i2)->f();

    if(v1 <= value && v2 >= value)
    {
      double ratio = (value - v1) / (v2 - v1);
      p = p1 + ratio * (p2 - p1);
      direction = p2 - p1;
      return true;
    }
    else if(v2 <= value && v1 >= value)
    {
      double ratio = (value - v2) / (v1 - v2);
      p = p2 + ratio * (p1 - p2);
      direction = p1 - p2;
      return true;
    }
    return false;
  }

}; // end of Reconstruction_triangulation_3

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_IMPLICIT_FCT_DELAUNAY_TRIANGULATION_H
