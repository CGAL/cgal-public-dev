// Copyright (c) 2017, 2018 GeometryFactory (France).
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
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_POLYLINE_TRACING_POINT_SET_TRACER_H
#define CGAL_POLYLINE_TRACING_POINT_SET_TRACER_H

#include <CGAL/Polyline_tracing/Motorcycle.h>
#include <CGAL/Polyline_tracing/Motorcycle_graph_node_dictionary.h>
#include <CGAL/Polyline_tracing/internal/robust_intersections.h>

#include <CGAL/Polygon_mesh_processing/locate.h>

#include <CGAL/assertions.h>
#include <CGAL/use.h>
#include <CGAL/Origin.h>
#include <CGAL/boost/graph/helpers.h>

#include <boost/optional.hpp>
#include <boost/variant.hpp>

#include <utility>
#include <vector>

namespace CGAL {

namespace Polyline_tracing {

template<typename MotorcycleGraph,
         typename DestinationContainer =
           std::vector<typename MotorcycleGraph::Geom_traits::Face_location> >
class Point_set_tracer
{
  typedef Point_set_tracer<MotorcycleGraph, DestinationContainer>           Self;

public:
  typedef typename MotorcycleGraph::Geom_traits                             Geom_traits;
  typedef typename Geom_traits::Triangle_mesh                               Triangle_mesh;

  typedef DestinationContainer                                              Destination_container;
  typedef typename Destination_container::value_type                        Destination;

  typedef typename Geom_traits::FT                                          FT;
  typedef typename Geom_traits::Point_d                                     Point;
  typedef typename Geom_traits::Segment_d                                   Segment;
  typedef typename Geom_traits::Vector_d                                    Vector;
  typedef typename Geom_traits::Ray_d                                       Ray;

  typedef typename MotorcycleGraph::Nodes                                   Nodes;
  typedef typename Nodes::Node_ptr                                          Node_ptr;

  typedef typename Geom_traits::Face_location                               Face_location;

  typedef typename MotorcycleGraph::Motorcycle                              Motorcycle;

  typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor      face_descriptor;

  // - bool: whether we have found a destination or not
  // - Node_ptr: the origin of the path (might be different from mc.current_position()
  //             if on the border)
  // - Node_ptr: the destination
  // - FT: the time at the destination
  // - bool: is the destination final
  typedef boost::tuple<bool, Node_ptr, Node_ptr, FT, bool>                  result_type;

  // Access
  // @todo add destination with a point ? dangerous if the point is on an edge
  void add_destination(const Destination& dest) { dests.push_back(dest); }
  Destination_container& destinations() { return dests; }
  const Destination_container& destinations() const { return dests; }

  // Constructor
  Point_set_tracer() : dests(), pos(0) { }
  Point_set_tracer(const Destination_container& destinations) : dests(destinations), pos(0) { }

  // Functions
  result_type set_next_destination(const Motorcycle& mc, Nodes& points, const Triangle_mesh& mesh) const;

  result_type operator()(const vertex_descriptor vd, const Motorcycle& mc,
                         Nodes& points, const Triangle_mesh& mesh) const;
  result_type operator()(const halfedge_descriptor hd, const Motorcycle& mc,
                         Nodes& points, const Triangle_mesh& mesh) const;
  result_type operator()(const face_descriptor fd, const Motorcycle& mc,
                         Nodes& points, const Triangle_mesh& mesh) const;

private:
  // A vector of destinations, with the conditions that two consecutive destinations
  // are on the same face of the mesh
  Destination_container dests;
  mutable std::size_t pos;
};

// -----------------------------------------------------------------------------

template<typename MotorcycleGraphTraits, typename DestinationContainer>
typename Point_set_tracer<MotorcycleGraphTraits, DestinationContainer>::result_type
Point_set_tracer<MotorcycleGraphTraits, DestinationContainer>::
set_next_destination(const Motorcycle& mc, Nodes& points, const Triangle_mesh& mesh) const
{
  CGAL_precondition(!destinations().empty());

  if(destinations().empty() || pos >= destinations().size())
  {
    std::cerr << "Warning: No destination or all destinations are already reached" << std::endl;
    return boost::make_tuple(false, mc.current_position(), mc.current_position(),
                             mc.current_time(), true /*final destination*/);
  }

  // intentional copies of both locations
  Face_location mc_loc = mc.current_location();
  Face_location dest_loc = destinations()[pos];
  Polygon_mesh_processing::locate_in_common_face(mc_loc, dest_loc, mesh);

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
  CGAL::Polygon_mesh_processing::internal::snap_location_to_border<Triangle_mesh>(dest_loc);
  CGAL_postcondition(CGAL::Polygon_mesh_processing::is_on_face_border(dest_loc, mesh));
#endif

  face_descriptor next_fd = mc_loc.first;
  Node_ptr origin_in_next_face = points.get_sibling(mc.current_position(), next_fd);
  std::pair<Node_ptr, bool> destination = points.insert(dest_loc, mesh);

  const Point& destination_point = destination.first->point();
  FT time_at_destination = mc.current_time() +
    CGAL::sqrt(CGAL::squared_distance(origin_in_next_face->point(),
                                      destination_point)) / mc.speed();

  // the last destination is marked as final
  const bool is_final_destination = (pos == destinations().size() - 1);

  ++pos;

  return boost::make_tuple(true, origin_in_next_face, destination.first,
                           time_at_destination, is_final_destination);
}

template<typename MotorcycleGraphTraits, typename DestinationContainer>
typename Point_set_tracer<MotorcycleGraphTraits, DestinationContainer>::result_type
Point_set_tracer<MotorcycleGraphTraits, DestinationContainer>::
operator()(vertex_descriptor vd, const Motorcycle& mc,
           Nodes& points, const Triangle_mesh& mesh) const
{
  CGAL_USE(vd);
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << " Point set tracing from a point on the vertex " << vd << std::endl;
#endif
  return set_next_destination(mc, points, mesh);
}

template<typename MotorcycleGraphTraits, typename DestinationContainer>
typename Point_set_tracer<MotorcycleGraphTraits, DestinationContainer>::result_type
Point_set_tracer<MotorcycleGraphTraits, DestinationContainer>::
operator()(halfedge_descriptor hd, const Motorcycle& mc,
           Nodes& points, const Triangle_mesh& mesh) const
{
  CGAL_USE(hd);
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << " Point set tracing from a point on the halfedge " << hd << std::endl;
#endif
  return set_next_destination(mc, points, mesh);
}

template<typename MotorcycleGraphTraits, typename DestinationContainer>
typename Point_set_tracer<MotorcycleGraphTraits, DestinationContainer>::result_type
Point_set_tracer<MotorcycleGraphTraits, DestinationContainer>::
operator()(face_descriptor fd, const Motorcycle& mc,
           Nodes& points, const Triangle_mesh& mesh) const
{
  CGAL_USE(fd);
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << " Point set tracing from a point in the face " << fd << std::endl;
#endif
  return set_next_destination(mc, points, mesh);
}

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_POINT_SET_TRACER_H
