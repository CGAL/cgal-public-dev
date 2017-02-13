// Copyright (c) 1999-2004   INRIA Sophia-Antipolis (France).
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

#ifndef DUAL_H // TODO: give a better name
#define DUAL_H

#include <utility>
#include <iterator>
#include <vector>
#include <map>

#include <boost/bind.hpp>
#include <boost/function_output_iterator.hpp>
#include <boost/range/join.hpp>

#include <CGAL/Bbox_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>

namespace CGAL {

// v must be an interior triangulation vertex.
// Construct a cell using only low-cost combinatorial operations
template < class Polyhedron_3, class DelaunayTriangulation_3 >
Polyhedron_3
dual(
  const DelaunayTriangulation_3& t,
  typename DelaunayTriangulation_3::Vertex_handle v
) {
  typedef typename DelaunayTriangulation_3::Point Point;
  typedef typename DelaunayTriangulation_3::Edge Edge;
  typedef typename DelaunayTriangulation_3::Vertex_handle Vertex_handle;
  typedef typename DelaunayTriangulation_3::Cell_handle Cell_handle;
  typedef typename DelaunayTriangulation_3::Cell_circulator Cell_circulator;

  CGAL_triangulation_precondition(v != Vertex_handle());
  CGAL_triangulation_expensive_precondition(tds().is_vertex(v));
  CGAL_triangulation_precondition(!t.is_infinite(v));
  CGAL_triangulation_precondition(t.dimension() == 3);

  std::vector<Edge> finite_inc_edges;
  finite_inc_edges.reserve(16); // < 16 in average?
  t.finite_incident_edges(v, std::back_inserter(finite_inc_edges));
  // Vertex is interior
  CGAL_triangulation_expensive_precondition(finite_inc_edges.size() == t.degree(v));

  typedef std::map<Cell_handle, size_t> Cell_id_map;
  // Indices correspond to cells incident to v and to its dual vertices.
  Cell_id_map cell_ids;
  {
    typedef std::pair<typename Cell_id_map::iterator, bool> Iter_bool_pair;
    typedef typename Cell_id_map::value_type Value_type;
    typedef Iter_bool_pair(Cell_id_map::*Insert)(const Value_type&);
    t.incident_cells(v, boost::make_function_output_iterator(
      boost::bind(static_cast<Insert>(&Cell_id_map::insert), &cell_ids,
        // Const references are necessary for C++11 compability.
        boost::bind(&std::make_pair<const Cell_handle&, const size_t&>,
          _1, boost::bind(&Cell_id_map::size, &cell_ids)))));
  }

  // Dual points of finite cells incident to v.
  std::vector<Point> points(cell_ids.size());
  for(typename Cell_id_map::iterator cit = cell_ids.begin(), end = cell_ids.end(); cit != end; ++cit) {
    points[cit->second] = t.dual(cit->first);
  }

  Polyhedron_3 result;
  Polyhedron_incremental_builder_3<typename Polyhedron_3::HalfedgeDS> builder(result.hds()/*TODO: keep it?*/, true);
  builder.begin_surface(cell_ids.size(), finite_inc_edges.size(), 3 * cell_ids.size());
  for(typename std::vector<Point>::iterator pit = points.begin(), pend = points.end(); pit != pend; ++pit) {
    builder.add_vertex(*pit);
  }

  for(typename std::vector<Edge>::iterator eit = finite_inc_edges.begin(), eend = finite_inc_edges.end(); eit != eend; ++eit) {
    Cell_circulator ccir = t.incident_cells(*eit);
    Cell_circulator cend = ccir;
    builder.begin_facet();
    do {
      builder.add_vertex_to_facet(cell_ids[ccir]);
      ++ccir;
    } while ( ccir != cend );
    builder.end_facet();
  }
  builder.end_surface();

  CGAL_triangulation_expensive_postcondition(result.is_valid());
  CGAL_triangulation_expensive_postcondition(result.is_closed());

  return result;
}

// Construct a cell as intersection of natural voronoi planes and user-defined bounding planes
// Uses an expesive halfspace_intersection_3 function
// The user-defined planes should be correctly defined (no checks) in order to obtain a valid result
template < class Polyhedron_3, class DelaunayTriangulation_3, class PlaneRange >
Polyhedron_3
dual(
  const DelaunayTriangulation_3& t,
  typename DelaunayTriangulation_3::Vertex_handle v,
  const PlaneRange& planes
) {
  typedef typename DelaunayTriangulation_3::Geom_traits::Plane_3 Plane;
  typedef typename DelaunayTriangulation_3::Vertex_handle Vertex_handle;

  CGAL_triangulation_precondition(v != Vertex_handle());
  CGAL_triangulation_expensive_precondition(tds().is_vertex(v));
  CGAL_triangulation_precondition(!t.is_infinite(v));
  CGAL_triangulation_precondition(t.dimension() == 3);

  std::vector<Vertex_handle> finite_adj_vertices;
  finite_adj_vertices.reserve(16); // TODO: < 16 in average?
  t.finite_adjacent_vertices(v, std::back_inserter(finite_adj_vertices));
  CGAL_assertion(!finite_adj_vertices.empty());

  std::vector<Plane> voronoi_planes;
  voronoi_planes.reserve(finite_adj_vertices.size());
  // Construct dual planes from incident edges
  for(typename std::vector<Vertex_handle>::iterator vit = finite_adj_vertices.begin(), vend = finite_adj_vertices.end(); vit != vend; ++vit) {
    voronoi_planes.push_back(t.geom_traits().construct_bisector_3_object()((*vit)->point(), v->point()));
  }

  boost::range::joined_range<std::vector<Plane>, const PlaneRange>
  all_planes =  boost::join(voronoi_planes, planes);
  Polyhedron_3 result;
  halfspace_intersection_3(boost::begin(all_planes), boost::end(all_planes), result);

  // TODO : should it be included if it's not possible to verify
  // the corresponding preconditions?
//   CGAL_triangulation_expensive_postcondition(result.is_valid());
//   CGAL_triangulation_expensive_postcondition(result.is_closed());

  return result;
}

template < class Polyhedron_3, class DelaunayTriangulation_3 >
Polyhedron_3
dual(
  const DelaunayTriangulation_3& t,
  typename DelaunayTriangulation_3::Vertex_handle v,
  const Bbox_3& bbox
) {
  typedef typename DelaunayTriangulation_3::Geom_traits::Plane_3 Plane;

  std::vector<Plane> bbox_planes;
  bbox_planes.reserve(6);
  bbox_planes.push_back(Plane(-1, 0, 0, bbox.xmin()));
  bbox_planes.push_back(Plane(1, 0, 0, -bbox.xmax()));
  bbox_planes.push_back(Plane(0, -1, 0, bbox.ymin()));
  bbox_planes.push_back(Plane(0, 1, 0, -bbox.ymax()));
  bbox_planes.push_back(Plane(0, 0, -1, bbox.zmin()));
  bbox_planes.push_back(Plane(0, 0, 1, -bbox.zmax()));

  return dual<Polyhedron_3>(t, v, bbox_planes);
}

} //namespace CGAL

#endif // DUAL_H