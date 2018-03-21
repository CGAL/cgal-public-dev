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
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_POLYLINE_TRACING_POINT_SET_TRACER_H
#define CGAL_POLYLINE_TRACING_POINT_SET_TRACER_H

#include <CGAL/Polyline_tracing/Dictionary.h>
#include <CGAL/Polyline_tracing/Motorcycle.h>
#include <CGAL/Polyline_tracing/internal/robust_intersections.h>

#include <CGAL/Origin.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polygon_mesh_processing/locate.h>

#include <boost/optional.hpp>
#include <boost/variant.hpp>

#include <utility>
#include <vector>

namespace CGAL {

namespace Polyline_tracing {

template<typename MotorcycleGraphTraits, typename Tracer>
class Motorcycle;

template<typename MotorcycleGraphTraits>
class Point_set_tracer
{
  typedef Point_set_tracer<MotorcycleGraphTraits>            Self;

public:
  typedef MotorcycleGraphTraits                               Geom_traits;
  typedef typename Geom_traits::Triangle_mesh                 Triangle_mesh;

  typedef typename Geom_traits::FT                            FT;
  typedef typename Geom_traits::Point_d                       Point;
  typedef typename Geom_traits::Segment_d                     Segment;
  typedef typename Geom_traits::Vector_d                      Vector;
  typedef typename Geom_traits::Ray_d                         Ray;

  typedef Dictionary<Geom_traits>                             Dictionary;
  typedef typename Dictionary::DEC_it                         DEC_it;

  typedef typename Geom_traits::Face_location                 Face_location;

  typedef Motorcycle<Geom_traits, Self>                       Motorcycle;

  typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor      face_descriptor;

  // - bool: whether we have found a destination or not
  // - DEC_it: the source of the path (might be different from mc.current_position()
  //           if on the border)
  // - DEC_it: the destination
  // - FT: the time at the destination
  // - bool: is the destination final
  typedef boost::tuple<bool, DEC_it, DEC_it, FT, bool>        result_type;

  // Access
  void set_destinations(const std::vector<Face_location>& dests) { destinations = dests; }

  // Constructor
  Point_set_tracer() : destinations(), pos(-1) { }
  Point_set_tracer(const std::vector<Face_location>& destinations) : destinations(destinations), pos(-1) { }

  // Functions
  result_type operator()(vertex_descriptor vd, const Motorcycle& mc,
                         Dictionary& points, const Triangle_mesh& mesh) const;
  result_type operator()(halfedge_descriptor hd, const Motorcycle& mc,
                         Dictionary& points, const Triangle_mesh& mesh) const;
  result_type operator()(face_descriptor fd, const Motorcycle& mc,
                         Dictionary& points, const Triangle_mesh& mesh) const;

private:
  // A vector of destination, with the conditions that two consecutive destinations
  // are on the same face of the mesh
  std::vector<Face_location> destinations; //@todo make it variant of point/face_location
  mutable std::size_t pos;
};

// -----------------------------------------------------------------------------

template<typename MotorcycleGraphTraits>
typename Point_set_tracer<MotorcycleGraphTraits>::result_type
Point_set_tracer<MotorcycleGraphTraits>::
operator()(vertex_descriptor vd, const Motorcycle& mc,
           Dictionary& /*points*/, const Triangle_mesh& /*mesh*/) const
{
  // check which face we should move into, find it, compute.
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << " Point set tracing from a point on the vertex " << vd << std::endl;
#endif
  CGAL_precondition(!destinations.empty());

  // @todo
  CGAL_assertion(false);

  ++pos;
  if(pos >= destinations.size())
  {
    std::cerr << "Warning: tried to get a destination but we have already reached all destinations" << std::endl;
    return boost::make_tuple(true, mc.current_position(), mc.current_position(),
                             mc.current_time(), true /*final destination*/);
  }

  return boost::make_tuple(true, mc.current_position(), mc.current_position(),
                           mc.current_time(), true /*final destination*/);
}

template<typename MotorcycleGraphTraits>
typename Point_set_tracer<MotorcycleGraphTraits>::result_type
Point_set_tracer<MotorcycleGraphTraits>::
operator()(halfedge_descriptor hd, const Motorcycle& mc,
           Dictionary& points, const Triangle_mesh& mesh) const
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << " Point set tracing from a point on the halfedge " << hd << std::endl;
#endif
  CGAL_precondition(!destinations.empty());

  ++pos;
  if(pos >= destinations.size())
  {
    std::cerr << "Warning: tried to get a destination but we have already reached all destinations" << std::endl;
    return boost::make_tuple(true, mc.current_position(), mc.current_position(),
                             mc.current_time(), true /*final destination*/);
  }

  // Compute the position of the motorcycle in the opposite face
  halfedge_descriptor opp_hd = opposite(hd, mesh);
  face_descriptor opp_fd = face(opp_hd, mesh);
  CGAL_assertion(opp_fd != boost::graph_traits<Triangle_mesh>::null_face());
  Face_location opp_loc =
    CGAL::Polygon_mesh_processing::locate_in_adjacent_face(mc.current_position()->location(), opp_fd, mesh);

  // Insert the source seen from the opposite face in the dictionary
  std::pair<DEC_it, bool> source_in_next_face = points.insert(opp_loc,
                                                              mc.current_position()->point(),
                                                              mesh);

  // Now deal with the destination
  const Face_location& loc = destinations[pos];
  std::pair<DEC_it, bool> destination = points.insert(loc, mesh);
  const Point& destination_point = destination.first->point();

  FT time_at_destination = mc.current_time() +
    CGAL::sqrt(CGAL::squared_distance(source_in_next_face.first->point(),
                                      destination_point)) / mc.speed();

  // last destination is marked as final
  bool is_final_destination = (pos == destinations.size() - 1);

  return boost::make_tuple(true, source_in_next_face.first, destination.first,
                           time_at_destination, is_final_destination);
}

template<typename MotorcycleGraphTraits>
typename Point_set_tracer<MotorcycleGraphTraits>::result_type
Point_set_tracer<MotorcycleGraphTraits>::
operator()(face_descriptor fd, const Motorcycle& mc,
           Dictionary& points, const Triangle_mesh& mesh) const
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << " Point set tracing from a point in the face " << fd << std::endl;
#endif
  CGAL_precondition(!destinations.empty());

  ++pos;
  if(pos >= destinations.size())
  {
    std::cerr << "Warning: tried to get a destination but we have already reached all destinations" << std::endl;
    return boost::make_tuple(true, mc.current_position(), mc.current_position(),
                             mc.current_time(), true /*final destination*/);
  }

  const Face_location& loc = destinations[pos];

  std::cout << "Destination n° " << pos << " location: " << loc.first
            << " bar: " << loc.second[0] << " " << loc.second[1] << " " << loc.second[2] << std::endl;

  CGAL_assertion(loc.first == fd);

  std::pair<DEC_it, bool> destination = points.insert(loc, mesh);
  const Point& destination_point = destination.first->point();

  FT time_at_destination = mc.current_time() +
    CGAL::sqrt(CGAL::squared_distance(mc.source()->point(), destination_point)) / mc.speed();

  // last destination is marked as final
  bool is_final_destination = (pos == destinations.size() - 1);

  return boost::make_tuple(true, mc.current_position(), destination.first,
                           time_at_destination, is_final_destination);
}

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_POINT_SET_TRACER_H
