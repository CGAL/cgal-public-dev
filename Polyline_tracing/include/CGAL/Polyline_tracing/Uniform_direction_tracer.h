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

#ifndef CGAL_POLYLINE_TRACING_UNIFORM_TRACER_H
#define CGAL_POLYLINE_TRACING_UNIFORM_TRACER_H

#include <CGAL/Polyline_tracing/Motorcycle.h>
#include <CGAL/Polyline_tracing/Motorcycle_graph_node_dictionary.h>
#include <CGAL/Polyline_tracing/internal/robust_intersections.h>

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Origin.h>
#include <CGAL/Polygon_mesh_processing/locate.h>

#include <boost/optional.hpp>
#include <boost/tuple/tuple.hpp>

#include <utility>

namespace CGAL {

namespace Polyline_tracing {

template<typename MotorcycleGraph>
class Uniform_direction_tracer_visitor
{
  typedef Uniform_direction_tracer_visitor<MotorcycleGraph>                 Self;

public:
  typedef typename MotorcycleGraph::Geom_traits                             Geom_traits;
  typedef typename Geom_traits::Triangle_mesh                               Triangle_mesh;

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
  typedef boost::variant<vertex_descriptor,
                         halfedge_descriptor,
                         face_descriptor>                                   descriptor_variant;

  // - bool: whether we have found a destination or not
  // - Node_ptr: the origin of the path (might be different from mc.current_position()
  //             when mc.current_position() is on a face border)
  // - Node_ptr: the destination
  // - FT: the time at the destination
  // - bool: is the destination final
  typedef boost::tuple<bool, Node_ptr, Node_ptr, FT, bool>                  result_type;

  // Constructors
  Uniform_direction_tracer_visitor(const Vector& direction = CGAL::NULL_VECTOR,
                                   const Geom_traits& gt = Geom_traits())
    :
      dir(direction),
      gt(gt)
  {
    // @tmp this tracer only makes sense in 2D for now
    CGAL_static_assertion((Geom_traits::dim == 2));
  }

  // Access
  Vector& direction() { return dir; }
  const Vector& direction() const { return dir; }

  // Functions
  bool is_direction_null(const Motorcycle& mc) const;
  result_type compute_next_destination(const Node_ptr start_point,
                                       const face_descriptor fd, const Motorcycle& mc,
                                       Nodes& points, const Triangle_mesh& mesh) const;

  result_type operator()(const vertex_descriptor vd, const Motorcycle& mc,
                         Nodes& points, const Triangle_mesh& mesh)  const;
  result_type operator()(const halfedge_descriptor hd, const Motorcycle& mc,
                         Nodes& points, const Triangle_mesh& mesh) const;
  result_type operator()(const face_descriptor fd, const Motorcycle& mc,
                         Nodes& points, const Triangle_mesh& mesh) const;

private:
  mutable Vector dir;
  Geom_traits gt;
};

// -----------------------------------------------------------------------------

template<typename MotorcycleGraph>
bool
Uniform_direction_tracer_visitor<MotorcycleGraph>::
is_direction_null(const Motorcycle& mc) const
{
  // If the direction is null, try to compute it from the origin and the destination,
  // if both were provided.
  if(CGAL::NULL_VECTOR == dir && bool(mc.input_destination()))
  {
    CGAL_assertion(mc.origin() != Node_ptr() && mc.destination() != Node_ptr());
    Vector actual_dir = Vector(mc.origin()->point(), mc.destination()->point());
    dir = actual_dir;
    return (actual_dir == CGAL::NULL_VECTOR);
  }

  return false;
}

template<typename MotorcycleGraph>
typename Uniform_direction_tracer_visitor<MotorcycleGraph>::result_type
Uniform_direction_tracer_visitor<MotorcycleGraph>::
compute_next_destination(const Node_ptr start_point, const face_descriptor fd,
                         const Motorcycle& mc, Nodes& points, const Triangle_mesh& mesh) const
{
  CGAL_precondition(start_point->face() == fd);
  CGAL_precondition(num_vertices(mesh) != 0);

  typedef typename property_map_selector<Triangle_mesh, boost::vertex_point_t>::const_type VertexPointMap;
  VertexPointMap vpmap = get_const_property_map(vertex_point, mesh);

  // @todo do I really need to look at all the halfedges? why not simply break
  // as soon as it's a valid destination with time different from current_time?
  Point farthest_destination;
  FT time_at_farthest_destination = mc.current_time(); // minimum allowed time value

  Ray r(start_point->point(), direction());

  typedef CGAL::Halfedge_around_face_circulator<Triangle_mesh>  Halfedge_around_facet_circulator;
  Halfedge_around_facet_circulator hcir_begin(halfedge(fd, mesh), mesh);
  Halfedge_around_facet_circulator hcir = hcir_begin;
  do
  {
    halfedge_descriptor hd = *hcir++;
    Segment s(get(vpmap, source(hd, mesh)), get(vpmap, target(hd, mesh)));

    if(gt.do_intersect_2_object()(r, s))
    {
      // returns a point because we ignore the degenerate configuration of the ray
      // and segment being aligned (the next halfedge will give us an intersection
      // at a vertex descriptor, which is the point we need)

      boost::optional<Point> res = internal::robust_intersection<Geom_traits>(r, s, gt);
      if(!res)
        continue;

      const Point new_destination = *res;

      // compute the time at destination
      FT time_at_new_destination = mc.current_time() +
        CGAL::sqrt(CGAL::squared_distance(start_point->point(), new_destination)) / mc.speed();

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
      std::cout << "new potential destination: " << new_destination
                << " a time: " << time_at_new_destination << std::endl;
#endif

      if(time_at_new_destination > time_at_farthest_destination)
      {
        farthest_destination = new_destination;
        time_at_farthest_destination = time_at_new_destination;
      }
    }
  } while (hcir != hcir_begin);

  CGAL_postcondition(time_at_farthest_destination >= mc.current_time());
  if(time_at_farthest_destination == mc.current_time())
  {
    std::cerr << "Warning: motorcycle has no interesting intersection "
              << "with the border of the face: " << fd << std::endl;

    // Since we have already dealt with the case of a null direction, if the current
    // position is the new destination, then the direction is pointing outside
    // of this face.
    CGAL_assertion(direction() != CGAL::NULL_VECTOR);
    return boost::make_tuple(false, Node_ptr(), Node_ptr(), 0., false /*not final*/);
  }

  Face_location loc = CGAL::Polygon_mesh_processing::locate_in_face(
                        farthest_destination, fd, mesh);

  // A uniform tracer will trace until it reaches a boundary. It is important
  // that the location of this new destination reflects that it is on the boundary
  // (that is, one of its barycentric coordinates should be 0).
  //
  // @todo a proper snap to closer halfedge / vertex. detect when we are walking an
  // edge, which should yield a target on a vertex (pass the bool in parameter).
  CGAL::Polygon_mesh_processing::internal::snap_location_to_border<Triangle_mesh>(loc);
  CGAL_postcondition(CGAL::Polygon_mesh_processing::is_on_face_border(loc, mesh));

  std::pair<Node_ptr, bool> destination = points.insert(loc, farthest_destination, mesh);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "new origin at      : " << &*start_point
            << " p: " << start_point->point()
            << " loc: " << start_point->face()
            << " bc: [" << start_point->location().second[0] << " "
                        << start_point->location().second[1] << " "
                        << start_point->location().second[2] << "] "
            << " time: " << mc.current_time() << std::endl;
  std::cout << "new destination at: " << &*destination.first
            << " p: " << destination.first->point()
            << " loc: " << destination.first->face()
            << " bc: [" << destination.first->location().second[0] << " "
                        << destination.first->location().second[1] << " "
                        << destination.first->location().second[2] << "] "
            << " time: " << time_at_farthest_destination << std::endl;
#endif

  return boost::make_tuple(true, start_point, destination.first,
                           time_at_farthest_destination, false /*not final*/);
}

template<typename MotorcycleGraph>
typename Uniform_direction_tracer_visitor<MotorcycleGraph>::result_type
Uniform_direction_tracer_visitor<MotorcycleGraph>::
operator()(vertex_descriptor vd, const Motorcycle& mc,
           Nodes& points, const Triangle_mesh& mesh) const
{
  // check which face we should move into, find it, compute.
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << " Uniform tracing from a point on the vertex " << vd;
  std::cout << " with direction: " << direction() << std::endl;
#endif

  if(is_direction_null(mc))
  {
    return boost::make_tuple(true, mc.current_position(), mc.current_position(),
                             mc.current_time(), true /*final destination*/);
  }

  halfedge_descriptor hd = halfedge(vd, mesh);

  // loop the incident faces of vd
  typedef CGAL::Face_around_target_circulator<Triangle_mesh> face_around_target_circulator;
  face_around_target_circulator fatc(hd, mesh), done(fatc);
  do
  {
    face_descriptor fd = *fatc;
    if(fd == boost::graph_traits<Triangle_mesh>::null_face())
    {
      ++fatc;
      continue;
    }

    // Compute the position of the motorcycle in the current face
    Node_ptr origin_in_fd = points.get_sibling(mc.current_position(), fd);
    result_type res = compute_next_destination(origin_in_fd, fd, mc, points, mesh);

    // Since direction == NULL_VECTOR has been filtered in Tracer.h, the destination
    // should not be equal to the origin
    // @todo This check would fail if one is manipulating a mesh with a completely
    // degenerate face
    if(res.template get<0>() && res.template get<2>() != origin_in_fd)
      return res;

    ++fatc;
  }
  while(fatc != done);

  // If we couldn't find a destination, then we must be on the border of the mesh
  // with a direction pointing out.
  CGAL_assertion(bool(is_border(vd, mesh)));
  return boost::make_tuple(false, mc.current_position(), mc.current_position(),
                           mc.current_time(), true /*final destination*/);
}

template<typename MotorcycleGraph>
typename Uniform_direction_tracer_visitor<MotorcycleGraph>::result_type
Uniform_direction_tracer_visitor<MotorcycleGraph>::
operator()(halfedge_descriptor hd, const Motorcycle& mc,
           Nodes& points, const Triangle_mesh& mesh) const
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << " Uniform tracing from a point on the halfedge " << hd;
  std::cout << " with direction: " << direction() << std::endl;
#endif
  namespace PMP = CGAL::Polygon_mesh_processing;

  if(is_direction_null(mc))
  {
    return boost::make_tuple(true, mc.current_position(), mc.current_position(),
                             mc.current_time(), true /*final destination*/);
  }

  // When we reach the border at the interior of a halfedge, the path continues
  // on the adjacent face.

  // Exception cases:
  // - We are computing the first destination
  // - The motorcycle is walking on the halfedge
  // In those cases, first try to find a valid destination on face(hd, mesh)
  const Face_location& origin_loc = mc.origin()->location();
  bool is_motorcycle_walking_hd = PMP::is_on_halfedge(origin_loc, hd, mesh);

  if(mc.input_destination() == boost::none || is_motorcycle_walking_hd)
  {
    face_descriptor fd = face(hd, mesh);
    result_type res = compute_next_destination(mc.current_position(), fd, mc, points, mesh);

    // Since (direction == NULL_VECTOR) has been filtered in Tracer.h, the destination
    // should not be equal to the origin
    // @todo This check would fail if one is manipulating a mesh with a degenerate (flat) face
    if(res.template get<0>() &&
       res.template get<2>() != mc.current_position())
      return res;
  }

  // Check the other face incident to hd
  halfedge_descriptor opp_hd = opposite(hd, mesh);

  if(is_border(opp_hd, mesh))
  {
    // Origin is on the border of the mesh and the direction is pointing out,
    return boost::make_tuple(false, mc.current_position(), mc.current_position(),
                             mc.current_time(), true /*final destination*/);
  }

  // Compute the position of the motorcycle in the opposite face
  face_descriptor opp_fd = face(opp_hd, mesh);
  CGAL_assertion(opp_fd != boost::graph_traits<Triangle_mesh>::null_face());

  // Insert the origin seen from the opposite face in the dictionary
  Node_ptr origin_in_next_face = points.get_sibling(mc.current_position(), opp_fd);
  result_type opp_res = compute_next_destination(origin_in_next_face, opp_fd, mc, points, mesh);
  CGAL_assertion(opp_res.template get<0>());

  return opp_res;
}

template<typename MotorcycleGraph>
typename Uniform_direction_tracer_visitor<MotorcycleGraph>::result_type
Uniform_direction_tracer_visitor<MotorcycleGraph>::
operator()(face_descriptor fd, const Motorcycle& mc,
           Nodes& points, const Triangle_mesh& mesh) const
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << " Uniform tracing from a point in the face " << fd;
  std::cout << " with direction: " << direction() << std::endl;
#endif

  if(is_direction_null(mc))
  {
    return boost::make_tuple(true, mc.current_position(), mc.current_position(),
                             mc.current_time(), true /*final destination*/);
  }

  return compute_next_destination(mc.current_position(), fd, mc, points, mesh);
}

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_UNIFORM_TRACER_H
