// Copyright (c) 2017 GeometryFactory (France).
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

#ifndef CGAL_POLYLINE_TRACING_UNIFORM_TRACER_H
#define CGAL_POLYLINE_TRACING_UNIFORM_TRACER_H

#include <CGAL/Polyline_tracing/Dictionary.h>
#include <CGAL/Polyline_tracing/Motorcycle.h>
#include <CGAL/Polyline_tracing/internal/robust_intersections.h>

#include <CGAL/Origin.h>
#include <CGAL/boost/graph/helpers.h>

#include <boost/variant.hpp>

#include <utility>

namespace CGAL {

namespace Polyline_tracing {

template<typename K, typename PolygonMesh>
class Motorcycle;

template<typename K, typename PolygonMesh>
class Uniform_direction_tracer_visitor
  : public boost::static_visitor<
             boost::tuple<typename Dictionary_entry<K, PolygonMesh>::Face_location,
                          typename Dictionary_entry<K, PolygonMesh>::Face_location,
                          typename K::FT> >
{
public:
  typedef typename K::FT                                            FT;
  typedef typename K::Point_2                                       Point;
  typedef typename K::Segment_2                                     Segment;
  typedef typename K::Vector_2                                      Vector;
  typedef typename K::Ray_2                                         Ray;

  typedef Dictionary<K, PolygonMesh>                                Dictionary;
  typedef typename Dictionary::DEC_it                               DEC_it;
  typedef Dictionary_entry<K, PolygonMesh>                          Dictionary_entry;
  typedef typename Dictionary_entry::Face_location                  Face_location;

  typedef Motorcycle<K, PolygonMesh>                                Motorcycle;

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor      face_descriptor;

  typedef boost::tuple<DEC_it, DEC_it, FT>                          result_type;

  Uniform_direction_tracer_visitor(const Motorcycle* mc,
                                   Dictionary& points,
                                   const PolygonMesh& mesh);

  result_type compute_next_destination(const DEC_it start_point,
                                       const face_descriptor fd) const;

  result_type operator()(vertex_descriptor vd) const;
  result_type operator()(halfedge_descriptor hd) const;
  result_type operator()(face_descriptor fd) const;

private:
  const Motorcycle* mc;
  Dictionary& points;
  const PolygonMesh& mesh;
};

// -----------------------------------------------------------------------------

template<typename K, typename PolygonMesh>
Uniform_direction_tracer_visitor<K, PolygonMesh>::
Uniform_direction_tracer_visitor(const Motorcycle* mc,
                                 Dictionary& points,
                                 const PolygonMesh& mesh)
  :
    mc(mc),
    points(points),
    mesh(mesh)
{ }

template<typename K, typename PolygonMesh>
typename Uniform_direction_tracer_visitor<K, PolygonMesh>::result_type
Uniform_direction_tracer_visitor<K, PolygonMesh>::
compute_next_destination(const DEC_it start_point,
                         const face_descriptor fd) const
{
  CGAL_precondition(start_point->location().first == fd);
  CGAL_precondition(mc->direction()); // direction must be known

  Vector mc_dir = *(mc->direction());

  // degenerate case
  if(mc_dir == CGAL::NULL_VECTOR)
  {
    std::cerr << "Warning: direction is null, thus the next destination is the current position" << std::endl;
    return boost::make_tuple(mc->position(), mc->position(), mc->current_time());
  }

  CGAL_precondition(num_vertices(mesh) != 0);

  typedef CGAL::Halfedge_around_face_circulator<PolygonMesh>  Halfedge_around_facet_circulator;

  Point farthest_destination;
  FT time_at_farthest_destination = std::numeric_limits<FT>::min();

  Ray r(mc->position()->point(), mc_dir);

  // @todo should be a function that computes the face

  Halfedge_around_facet_circulator hcir_begin(halfedge(fd, mesh), mesh);
  Halfedge_around_facet_circulator hcir = hcir_begin;

  do
  {
    halfedge_descriptor hd = *hcir++;

    // @BGL point
    Segment s(mesh.point(source(hd, mesh)), mesh.point(target(hd, mesh)));
    std::cout << "ray; segment r:" << r << " s: " << s << std::endl;

    if(K().do_intersect_2_object()(r, s))
    {
      const Point new_destination = internal::robust_intersection<K>(r, s);
      std::cout << "new potential destination: " << new_destination << std::endl;

      // compute time at destination
      FT time_at_new_destination = mc->current_time() + // @todo tracer
        CGAL::sqrt(CGAL::squared_distance(mc->position()->point(), new_destination)) / mc->speed();

      if(time_at_new_destination > time_at_farthest_destination)
      {
        farthest_destination = new_destination;
        time_at_farthest_destination = time_at_new_destination;
      }
    }
  } while (hcir != hcir_begin);

  // no intersection with the border... Is the point not in the face?
  CGAL_assertion(time_at_farthest_destination != std::numeric_limits<FT>::min());

  // @todo handle the case where the new destination is the current_point
  // (I guess the motorcycle crashes in that case... ?)

  std::cout << "new destination at : " << farthest_destination
            << " time: " << time_at_farthest_destination << std::endl;

  Face_location loc = internal::compute_location<Point, Face_location, PolygonMesh>(farthest_destination, fd, mesh);
  DEC_it new_destination_it = points.insert(farthest_destination, loc);

  return boost::make_tuple(start_point,
                           new_destination_it,
                           time_at_farthest_destination);
}

template<typename K, typename PolygonMesh>
typename Uniform_direction_tracer_visitor<K, PolygonMesh>::result_type
Uniform_direction_tracer_visitor<K, PolygonMesh>::
operator()(vertex_descriptor vd) const
{
  // check which face we should move into, find it, compute.
  // be careful of null_face
  CGAL_assertion(false); // todo
  return result_type();
}

template<typename K, typename PolygonMesh>
typename Uniform_direction_tracer_visitor<K, PolygonMesh>::result_type
Uniform_direction_tracer_visitor<K, PolygonMesh>::
operator()(halfedge_descriptor hd) const
{
  // find out which face we should move into, move into it !
  // be careful of null_face
  CGAL_assertion(false);
  return result_type();
}

template<typename K, typename PolygonMesh>
typename Uniform_direction_tracer_visitor<K, PolygonMesh>::result_type
Uniform_direction_tracer_visitor<K, PolygonMesh>::
operator()(face_descriptor fd) const
{
  // @todo
  CGAL_assertion(false);
  return compute_next_destination(DEC_it(), fd);
}

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_UNIFORM_TRACER_H
