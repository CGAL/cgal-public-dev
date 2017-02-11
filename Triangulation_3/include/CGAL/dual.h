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

template < class Polyhedron_3, class DelaunayTriangulation_3, class PlaneIterator >
Polyhedron_3
dual(const DelaunayTriangulation_3& t, typename DelaunayTriangulation_3::Vertex_handle v, PlaneIterator planes_begin, PlaneIterator planes_end)
{
  typedef typename DelaunayTriangulation_3::Point Point;
  typedef typename DelaunayTriangulation_3::Geom_traits::Plane_3 Plane;
  typedef typename DelaunayTriangulation_3::Edge Edge;
  typedef typename DelaunayTriangulation_3::Vertex_handle Vertex_handle;
  typedef typename DelaunayTriangulation_3::Cell_handle Cell_handle;
  typedef typename DelaunayTriangulation_3::Cell_circulator Cell_circulator;

  CGAL_triangulation_precondition(v != Vertex_handle());
  CGAL_triangulation_precondition(!t.is_infinite(v));
  CGAL_triangulation_precondition(t.dimension() == 3);
  CGAL_triangulation_expensive_precondition(tds().is_vertex(v));

  Polyhedron_3 result;
  {
    std::vector<Edge> finite_inc_edges;
    finite_inc_edges.reserve(16); // < 16 in average?
    t.finite_incident_edges(v, std::back_inserter(finite_inc_edges));
    CGAL_assertion(!finite_inc_edges.empty());

    if(finite_inc_edges.size() == t.degree(v) && planes_begin == planes_end) {
      // Vertex v is NOT on the triangulation convex hull and there are no user-defined boundary planes.
      // Dual cell can be constructed directly from dual points of the triangulation cells incident to v
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
      CGAL_triangulation_postcondition(result.is_closed());

      return result;
    }
  }

  // Vertex v is on the convex hull or user-defined boundary planes exist.
  // Dual cell will be constructed through the intersection of the voronoi dual planes and the input planes.
  std::vector<Vertex_handle> finite_adj_vertices;
  finite_adj_vertices.reserve(16); // TODO: < 16 in average?
  t.finite_adjacent_vertices(v, std::back_inserter(finite_adj_vertices));
  CGAL_assertion(!finite_adj_vertices.empty());

  std::vector<Plane> voronoi_planes;
  voronoi_planes.reserve(finite_adj_vertices.size());
  for(typename std::vector<Vertex_handle>::iterator vit = finite_adj_vertices.begin(), vend = finite_adj_vertices.end(); vit != vend; ++vit) {
    voronoi_planes.push_back(t.geom_traits().construct_bisector_3_object()((*vit)->point(), v->point()));
  }

  boost::iterator_range<PlaneIterator>
  input_planes = boost::make_iterator_range(planes_begin, planes_end);
  boost::range::joined_range<std::vector<Plane>, boost::iterator_range<PlaneIterator> >
  plane_range =  boost::join(voronoi_planes, input_planes);
  halfspace_intersection_3(boost::begin(plane_range), boost::end(plane_range), result);

  return result;
}

template < class Polyhedron_3, class DelaunayTriangulation_3 >
Polyhedron_3
dual(const DelaunayTriangulation_3& t, typename DelaunayTriangulation_3::Vertex_handle v, const Bbox_3& bbox)
{
  typedef typename DelaunayTriangulation_3::Geom_traits::Plane_3 Plane;

  std::vector<Plane> bbox_planes;
  bbox_planes.reserve(6);
  bbox_planes.push_back(Plane(-1, 0, 0, bbox.xmin()));
  bbox_planes.push_back(Plane(1, 0, 0, -bbox.xmax()));
  bbox_planes.push_back(Plane(0, -1, 0, bbox.ymin()));
  bbox_planes.push_back(Plane(0, 1, 0, -bbox.ymax()));
  bbox_planes.push_back(Plane(0, 0, -1, bbox.zmin()));
  bbox_planes.push_back(Plane(0, 0, 1, -bbox.zmax()));

  return dual<Polyhedron_3>(t, v, bbox_planes.begin(), bbox_planes.end());
}

} //namespace CGAL

#endif // DUAL_H