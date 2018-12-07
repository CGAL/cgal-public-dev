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

#ifndef CGAL_POLYLINE_TRACING_MOTORCYCLE_H
#define CGAL_POLYLINE_TRACING_MOTORCYCLE_H

#include <CGAL/Polyline_tracing/Motorcycle_graph_node_dictionary.h>
#include <CGAL/Polyline_tracing/Track.h>

#include <CGAL/Polygon_mesh_processing/locate.h>

#include <CGAL/function.h>

#include <boost/optional.hpp>
#include <boost/parameter.hpp>
#include <boost/variant.hpp>

#include <fstream>
#include <iostream>
#include <list>
#include <set>
#include <sstream>

namespace CGAL {

namespace Polyline_tracing {

namespace internal {

template<typename MotorcycleGraphTraits>
class Target_point_set_comparer
{
  typedef MotorcycleGraphTraits                                    Geom_traits;
  typedef typename Geom_traits::FT                                 FT;
  typedef Motorcycle_graph_node_dictionary<Geom_traits>            Nodes;
  typedef typename Nodes::Node_ptr                                 Node_ptr;

public:
  // Given two pairs <point, visiting time>, compare the visiting times
  bool operator()(const std::pair<Node_ptr, FT>& lhs, const std::pair<Node_ptr, FT>& rhs) const {
    return lhs.second < rhs.second;
  }
};

} // namespace internal

// -----------------------------------------------------------------------------

template<typename MotorcycleGraphTraits>
class Motorcycle
{
  typedef Motorcycle<MotorcycleGraphTraits>                                 Self;

public:
  typedef MotorcycleGraphTraits                                             Geom_traits;
  typedef typename Geom_traits::Triangle_mesh                               Triangle_mesh;

  typedef typename Geom_traits::FT                                          FT;
  typedef typename Geom_traits::Point_d                                     Point;
  typedef typename Geom_traits::Vector_d                                    Vector;

  typedef Motorcycle_graph_node_dictionary<Geom_traits>                     Nodes;
  typedef typename Nodes::Node                                              Node;
  typedef typename Nodes::Node_ptr                                          Node_ptr;
  typedef typename Node::Face_location                                      Face_location;

  typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor      face_descriptor;

  typedef boost::variant<vertex_descriptor,
                         halfedge_descriptor,
                         face_descriptor>                                   descriptor_variant;

  typedef boost::variant<Point, Face_location>                              Point_or_location;
  typedef boost::optional<Point_or_location>                                Optional_point_or_location;

  // Target points are put in a set sorted by time. It is not a multiset because
  // same times should be equal to same points.
  typedef std::pair<Node_ptr, FT>                                           Track_point;
  typedef std::set<Track_point,
                   internal::Target_point_set_comparer<Geom_traits> >       Target_point_container;
  typedef typename Target_point_container::iterator                         TPC_iterator;

  typedef Motorcycle_track<Geom_traits>                                     Track;
  typedef typename Track::Track_segment                                     Track_segment;
  typedef typename Track::iterator                                          Track_iterator;

  // - bool: whether we have found a destination or not
  // - Node_ptr: the origin of the path
  // - Node_ptr: the destination
  // - FT: the time at the destination
  // - bool: is the destination final
  typedef boost::tuple<bool, Node_ptr, Node_ptr, FT, bool>                  result_type;

  // Using type erasure to avoid templating the motorcycle with the tracer type
  typedef CGAL::cpp11::function<result_type(const vertex_descriptor,
                                            const Self&, Nodes&,
                                            const Triangle_mesh&)>          Vertex_tracer;
  typedef CGAL::cpp11::function<result_type(const halfedge_descriptor,
                                            const Self&, Nodes&,
                                            const Triangle_mesh&)>          Halfedge_tracer;
  typedef CGAL::cpp11::function<result_type(const face_descriptor,
                                            const Self&, Nodes&,
                                            const Triangle_mesh&)>          Face_tracer;

  enum Motorcycle_status
  {
    IN_MOTION = 0,
    CRASHED
  };

  enum Motorcycle_nature
  {
    FREE_MOTORCYCLE = 0,
    CONSTRAINT_MOTORCYCLE,
    BORDER_MOTORCYCLE
  };

  // Access
  int id() const { return id_; }
  void set_id(int id) { id_ = id; }
  Motorcycle_status& status() { return status_; }
  const Motorcycle_status& status() const { return status_; }
  const Motorcycle_nature& nature() const { return nature_; }

  const Point_or_location& input_origin() const { return input_orig; }
  Optional_point_or_location& input_destination() { return input_dest; }
  const Optional_point_or_location& input_destination() const { return input_dest; }

  Node_ptr& origin() { return orig; }
  const Node_ptr& origin() const { CGAL_precondition(orig != Node_ptr()); return orig; }
  FT& time_at_origin() { return origin_time; }
  const FT& time_at_origin() const { return origin_time; }

  bool is_initialized() const { return conf != Node_ptr() && dest != Node_ptr(); }
  Node_ptr& current_position() { return conf; }
  const Node_ptr& current_position() const { CGAL_precondition(conf != Node_ptr()); return conf; }
  const Face_location& current_location() const { return current_position()->location(); }
  face_descriptor current_face() const { return current_position()->face(); }
  FT& current_time() { return current_time_; }
  const FT& current_time() const { return current_time_; }

  Node_ptr closest_target() const;
  FT time_at_closest_target() const;

  Node_ptr& destination() { return dest; }
  const Node_ptr destination() const { CGAL_precondition(dest != Node_ptr()); return dest; }
  FT& time_at_destination() { return destination_time; }
  const FT& time_at_destination() const { return destination_time; }

  bool& is_destination_final() { return is_dest_final; }
  bool is_destination_final() const { return is_dest_final; }

  FT speed() const { return speed_; }

  Target_point_container& targets() { return target_points; }
  const Target_point_container& targets() const { return target_points; }

  Track& track() { return track_; }
  const Track& track() const { return track_; }

  Track_iterator newest_track_segment() { CGAL_precondition(!track_.empty()); return --(track_.end()); }

  // Constructor
  template<typename Tracer,
           typename NamedParameters = cgal_bgl_named_params<bool, internal_np::all_default_t> >
  Motorcycle(const Point_or_location& origin, const Tracer& tracer,
             const NamedParameters& np = CGAL::parameters::all_default());

  // Functions
  void add_target(const Node_ptr target_point, const FT time_at_target);
  void clear_targets();

  result_type compute_next_destination(Nodes& points, const Triangle_mesh& mesh);

  std::pair<TPC_iterator, bool> has_target(const Node_ptr e) const;
  std::pair<TPC_iterator, bool> has_target(const Face_location loc) const;
  std::pair<TPC_iterator, bool> has_target_at_time(const FT visiting_time) const;
  std::pair<TPC_iterator, bool> has_target_at_time(const FT min_visiting_time, const FT max_visiting_time) const;
  bool has_left_starting_position() const;
  bool has_target_at_time(const Node_ptr e, const FT visiting_time) const;
  bool has_reached_blocked_point() const;
  bool is_tentative_track_degenerate() const;

  void remove_closest_target_from_targets();

  void crash();
  bool drive_to_closest_target();

  // Output
  friend std::ostream& operator<<(std::ostream& out, const Self& mc)
  {
    out << "Motorcycle #" << mc.id() << " (status? " << mc.status() << ") ";
    if(!mc.is_initialized())
      return out;

    out << "going from (current) origin: (" << mc.origin()->point() << ")"
        << " to destination: (" << mc.destination()->point() << ")" << std::endl
        << "  currently at position: " << &*(mc.current_position()) << " (" << mc.current_position()->point() << ")"
        << " [L: " << mc.current_face() << "]"
        << " at time: " << mc.current_time() << std::endl
        << "  with targets:" << std::endl;
    typename Target_point_container::const_iterator tpc_it = mc.targets().begin();
    typename Target_point_container::const_iterator end = mc.targets().end();
    for(; tpc_it!=end; ++tpc_it)
    {
      out << "\t " << &*(tpc_it->first)
          << " P: (" << tpc_it->first->point() << ") T: " << tpc_it->second
          << " B: " << tpc_it->first->is_blocked()
          << " L: " << tpc_it->first->face()
          << " bc: [" << tpc_it->first->location().second[0] << " "
                      << tpc_it->first->location().second[1] << " "
                      << tpc_it->first->location().second[2] << "]" << std::endl;
    }

    return out;
  }

  void output_origin_and_destination() const;
  void output_track() const;

protected:
  // ID and status
  int id_;
  Motorcycle_status status_;

  // @fixme 'nature' shouldn't exist, info() should be used, defined in the user code
  const Motorcycle_nature nature_; // a motorcycle is like a scorpion

  // The very first origin and destination points, before insertion in the dictionary
  const Point_or_location input_orig;
  boost::optional<Point_or_location> input_dest;

  // Below might change when we move in the mesh
  Node_ptr orig; // origin
  Node_ptr dest; // destination
  Node_ptr conf; // current position (last confirmed position)

  // indicates whether we should stop at the destination or try to trace
  bool is_dest_final;

  const FT speed_; // speed of the motorcycle
  FT origin_time; // time at the current origin
  FT current_time_; // time at the current position
  FT destination_time; // time at the destination

  // The tentative targets, ordered by increasing distance from 'conf'
  Target_point_container target_points;

  Track track_;

  // Tracers
  Vertex_tracer vertex_tracer;
  Halfedge_tracer halfedge_tracer;
  Face_tracer face_tracer;

private:
  // Explanation about disallowing all copy/moves operators:
  // - minor: This class is heavy
  // - main: The three tracers above are constructed in a very particular way,
  //         due to the usage of type erasure for the tracer (which avoids having to
  //         template this Motorcycle class with a Tracer class): when we construct
  //         'vertex_tracer', we intentionally copy the tracer, because we don't want
  //         to assume that the tracer will have a long-enough lifetime. However,
  //         we do not want to copy it for 'halfedge/face_tracers' because tracers
  //         usually have states. Thus, we initialize them using the tracer contained
  //         within 'vertex_tracer'.
  //         Problem: if you copy, the default operator is not a deep copy
  //         and the references in 'halfedge/face_tracers' do not point
  //         to the new copied vertex_tracer... and you run in trouble.
  //         Solution: disable all copy/move operations.

  // disable copy operators
  Motorcycle& operator=(const Motorcycle& other);
  Motorcycle(const Motorcycle& other);

  // disable move operators
  Motorcycle (Motorcycle&& other);
  Motorcycle& operator=(Motorcycle&& other);
};

template<typename MotorcycleGraphTraits>
template <typename Tracer, typename NamedParameters>
Motorcycle<MotorcycleGraphTraits>::
Motorcycle(const Point_or_location& origin, const Tracer& tracer, const NamedParameters& np)
  :
    id_(-1),
    status_(IN_MOTION),
    nature_(boost::choose_param(boost::get_param(np, internal_np::nature), FREE_MOTORCYCLE)),
    input_orig(origin),
    input_dest(boost::choose_param(boost::get_param(np, internal_np::destination), boost::none)),
    orig(), dest(), conf(),
    is_dest_final(false),
    speed_(boost::choose_param(boost::get_param(np, internal_np::speed), 1.)),
    origin_time(boost::choose_param(boost::get_param(np, internal_np::initial_time), 0.)),
    current_time_(origin_time),
    target_points(internal::Target_point_set_comparer<MotorcycleGraphTraits>()),
    track_(),
    vertex_tracer(tracer),
    // about below, see remark above
    halfedge_tracer(std::ref(*(vertex_tracer.template target<Tracer>()))),
    face_tracer(std::ref(*(vertex_tracer.template target<Tracer>())))
{
  CGAL_precondition(speed() > 0.);
}

template<typename MotorcycleGraphTraits>
void
Motorcycle<MotorcycleGraphTraits>::
add_target(const Node_ptr target_point, const FT time_at_target)
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << " > Adding target: " << &*target_point
            << " at time: " << time_at_target
            << " to motorcycle #" << id_ << std::endl;
#endif

  CGAL_precondition(target_point != Node_ptr());

  // Don't add targets to a crashed motorcycle
  CGAL_precondition(status() != CRASHED);

  // Don't want to insert targets in other faces
  CGAL_precondition(target_point->face() == current_face());

  // Don't want to insert the same point twice...
  CGAL_expensive_precondition(!has_target(target_point).second);
  // ... or accidentally ignore adding a point due to another point having the same time
  CGAL_expensive_precondition(!has_target_at_time(time_at_target).second);

  // No target should be inserted before the current point
  // Note: equality is important because we might re-insert the current point
  //       with the current time to bump the priority to the top of the queue
  //       (e.g. after computing a new destination)
  CGAL_precondition(time_at_target >= current_time());

  targets().insert(std::make_pair(target_point, time_at_target));
}

template<typename MotorcycleGraphTraits>
void
Motorcycle<MotorcycleGraphTraits>::
crash()
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~X" << std::endl;
  std::cout << "Crashing " << *this << std::endl;
#endif

  CGAL_precondition(status() != CRASHED);
  status() = CRASHED;

  clear_targets();
}

template<typename MotorcycleGraphTraits>
void
Motorcycle<MotorcycleGraphTraits>::
clear_targets()
{
  // Go through the next targets of the motorcycle and remove it from the list
  // of motorcycles that might reach the target
  TPC_iterator it = targets().begin(), end = targets().end();
  for(; it!=end; ++it)
  {
    Node_ptr target_point = it->first;
    CGAL_assertion(it->second > current_time());

    target_point->remove_motorcycle(id());
    // @todo if 'target_point' does not contain any motorcycle after this remove, delete it ?
    // keeping it can maybe help with snapping?
  }

  targets().clear();
}

template<typename MotorcycleGraphTraits>
typename Motorcycle<MotorcycleGraphTraits>::Node_ptr
Motorcycle<MotorcycleGraphTraits>::
closest_target() const
{
  CGAL_precondition(!targets().empty());
  return targets().begin()->first;
}

template<typename MotorcycleGraphTraits>
void
Motorcycle<MotorcycleGraphTraits>::
remove_closest_target_from_targets()
{
  CGAL_precondition(!targets().empty());
  targets().erase(targets().begin());
}

template<typename MotorcycleGraphTraits>
std::pair<typename Motorcycle<MotorcycleGraphTraits>::TPC_iterator, bool>
Motorcycle<MotorcycleGraphTraits>::
has_target(const Node_ptr e) const
{
  // Note that since the set is sorted on the visting time, we have no choice
  // but to loop linearly thus this runs in O(n). It's bad, but usually the list
  // of targets is small. @todo ?
  TPC_iterator tpit = targets().begin(), end = targets().end();
  for(; tpit!=end; ++tpit)
    if(tpit->first == e)
      return std::make_pair(tpit, true);

  return std::make_pair(end, false);
}

template<typename MotorcycleGraphTraits>
std::pair<typename Motorcycle<MotorcycleGraphTraits>::TPC_iterator, bool>
Motorcycle<MotorcycleGraphTraits>::
has_target(const Face_location loc) const
{
  // Same remark as the 'has_target' function above.
  TPC_iterator tpit = targets().begin(), end = targets().end();
  for(; tpit!=end; ++tpit)
    if(tpit->first->location() == loc)
      return std::make_pair(tpit, true);

  return std::make_pair(end, false);
}

template<typename MotorcycleGraphTraits>
std::pair<typename Motorcycle<MotorcycleGraphTraits>::TPC_iterator, bool>
Motorcycle<MotorcycleGraphTraits>::
has_target_at_time(const FT visiting_time) const
{
  TPC_iterator res = targets().find(std::make_pair(Node_ptr(), visiting_time));
  return std::make_pair(res, (res != targets().end()));
}

template<typename MotorcycleGraphTraits>
std::pair<typename Motorcycle<MotorcycleGraphTraits>::TPC_iterator, bool>
Motorcycle<MotorcycleGraphTraits>::
has_target_at_time(const FT min_visiting_time, const FT max_visiting_time) const
{
  std::cout << "checking for target in interval: " << min_visiting_time << " || " << max_visiting_time << std::endl;
  CGAL_precondition(min_visiting_time <= max_visiting_time);

  TPC_iterator tit = targets().lower_bound(std::make_pair(Node_ptr(), min_visiting_time));
  bool is_valid_iterator = (tit != targets().end() && tit->second <= max_visiting_time);

  // "while" but actually the function just returns the first it finds (there
  // should only be one if max-min is small - as it should be).
  while(is_valid_iterator)
  {
    Node_ptr target = tit->first;
    const FT visiting_time = tit->second;
    CGAL_assertion(target != Node_ptr() && visiting_time >= min_visiting_time);

    if(visiting_time <= max_visiting_time)
      return std::make_pair(tit, true);

    ++tit;
    is_valid_iterator = (tit != targets().end() && tit->second <= max_visiting_time);
  }

  return std::make_pair(targets().end(), false);
}

template<typename MotorcycleGraphTraits>
bool
Motorcycle<MotorcycleGraphTraits>::
has_left_starting_position() const
{
  return (current_position() != origin());
}

template<typename MotorcycleGraphTraits>
bool
Motorcycle<MotorcycleGraphTraits>::
has_target_at_time(const Node_ptr e, const FT visiting_time) const
{
  TPC_iterator res = targets().find(std::make_pair(e, visiting_time));
  return (res != targets().end() && res->first == e);
}

template<typename MotorcycleGraphTraits>
bool
Motorcycle<MotorcycleGraphTraits>::
has_reached_blocked_point() const
{
  CGAL_precondition(is_initialized());
  return current_position()->is_blocked();
}

template<typename MotorcycleGraphTraits>
bool
Motorcycle<MotorcycleGraphTraits>::
is_tentative_track_degenerate() const
{
  CGAL_precondition(is_initialized());
  return (targets().empty() ||
          current_position() == closest_target());
}

template<typename MotorcycleGraphTraits>
typename Motorcycle<MotorcycleGraphTraits>::FT
Motorcycle<MotorcycleGraphTraits>::
time_at_closest_target() const
{
  // Motorcycles are in the priority queue before they are properly initialized
  // (i.e. before their source/destinations are put in the dictionary and before
  // the list of initial targets is set up).
  if(!is_initialized())
    return current_time();

  CGAL_precondition(!targets().empty());
  return targets().begin()->second;
}

template<typename MotorcycleGraphTraits>
bool
Motorcycle<MotorcycleGraphTraits>::
drive_to_closest_target()
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>" << std::endl;
  std::cout << "Driving " << *this << std::endl;
#endif

  CGAL_precondition(status() == IN_MOTION);
  CGAL_precondition(!targets().empty());

  // Don't create degenerate track segments if they are not required (otherwise
  // it increases the cost of computing intersections between tracks)
  // @todo don't insert the first one that is degenerate (see also mg::drive_to_closest_target())
  bool created_new_track_segment = false;
  if(track().empty() || current_position() != closest_target())
  {
    Track_segment ts(id(), current_position(), current_time(), closest_target(), time_at_closest_target());
    track().push_back(ts);
    created_new_track_segment = true;
  }

  current_time() = time_at_closest_target();
  current_position() = closest_target();
  remove_closest_target_from_targets();

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "Reached point:" << std::endl << *(current_position()) << std::endl;
#endif

  return created_new_track_segment;
}

template<typename MotorcycleGraphTraits>
typename Motorcycle<MotorcycleGraphTraits>::result_type
Motorcycle<MotorcycleGraphTraits>::
compute_next_destination(Nodes& points, const Triangle_mesh& mesh)
{
  CGAL_precondition(targets().empty());

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*" << std::endl;
  std::cout << "Computing the next path for motorcycle #" << id() << std::endl;
  std::cout << "Current position: " << conf->point() << std::endl
            << "Location: " << current_face() << " b: "
            << current_location().second[0] << " "
            << current_location().second[1] << " "
            << current_location().second[2] << std::endl;
#endif

  const Face_location& loc = current_location();
  descriptor_variant dv = CGAL::Polygon_mesh_processing::get_descriptor_from_location(loc, mesh);

  if(const vertex_descriptor* v = boost::get<vertex_descriptor>(&dv))
  {
    return vertex_tracer(*v, *this, points, mesh);
  }
  else if(const halfedge_descriptor* h = boost::get<halfedge_descriptor>(&dv))
  {
    return halfedge_tracer(*h, *this, points, mesh);
  }
  else
  {
    const face_descriptor* f = boost::get<face_descriptor>(&dv);
    CGAL_assertion(f);
    return face_tracer(*f, *this, points, mesh);
  }
}

template<typename MotorcycleGraphTraits>
void
Motorcycle<MotorcycleGraphTraits>::
output_origin_and_destination() const
{
  std::stringstream oss_orig, oss_dest;
  oss_orig << "results_" << Geom_traits::dimension() << "/motorcycles_origins.xyz" << std::ends;
  oss_dest << "results_" << Geom_traits::dimension() << "/motorcycles_destinations.xyz" << std::ends;
  std::ofstream oof, odf;

  oof.open(oss_orig.str().c_str(), std::ios::app);
  odf.open(oss_dest.str().c_str(), std::ios::app);

  oof.precision(17);
  odf.precision(17);

  oof << origin()->point();
  odf << destination()->point();

  if(Geom_traits::dimension() == 2) // The '.xyz' format expects 3D points
  {
    oof << " 0";
    odf << " 0";
  }

  oof << '\n';
  odf << '\n';
}

template<typename MotorcycleGraphTraits>
void
Motorcycle<MotorcycleGraphTraits>::
output_track() const
{
  std::ostringstream out_filename;
  out_filename << "results_" << Geom_traits::dimension() << "/track_" << id() << ".off" << std::ends;
  std::ofstream out(out_filename.str().c_str());
  out.precision(17);

  return track().write_off(out);
}

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_MOTORCYCLE_H
