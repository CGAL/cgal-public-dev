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

#ifndef CGAL_POLYLINE_TRACING_MOTORCYCLE_GRAPH_H
#define CGAL_POLYLINE_TRACING_MOTORCYCLE_GRAPH_H

#include <CGAL/Polyline_tracing/Dictionary.h>
#include <CGAL/Polyline_tracing/Motorcycle.h>
#include <CGAL/Polyline_tracing/Motorcycle_priority_queue.h>
#include <CGAL/Polyline_tracing/internal/robust_intersections.h>
#include <CGAL/Polyline_tracing/internal/VPM_selector.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/array.h>
#include <CGAL/assertions.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/enum.h>
#include <CGAL/number_utils.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/result_of.h>

#include <boost/optional.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/unordered_map.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>

namespace CGAL {

namespace Polyline_tracing {

namespace internal {

template<typename Motorcycle_graph>
struct Collision_information
{
  typedef typename Motorcycle_graph::FT                      FT;
  typedef typename Motorcycle_graph::DEC_it                  DEC_it;
  typedef typename Motorcycle_graph::Face_location           Face_location;

  Collision_information()
    :
      closest_collision_location(),
      closest_collision(),
      is_closest_collision_location_already_in_dictionary(false),
      foreign_mc_id(-1),
      time_at_closest_collision(std::numeric_limits<FT>::max()),
      foreign_time_at_closest_collision(std::numeric_limits<FT>::max())
  { }

  Face_location closest_collision_location;
  DEC_it closest_collision;
  bool is_closest_collision_location_already_in_dictionary;
  std::size_t foreign_mc_id;
  FT time_at_closest_collision;
  FT foreign_time_at_closest_collision;
};

} // namespace internal

template<typename MotorcycleGraphTraits>
class Motorcycle_graph
{
  typedef Motorcycle_graph<MotorcycleGraphTraits>             Self;

  typedef internal::Collision_information<Self>               Collision_information;
  typedef typename MotorcycleGraphTraits::Kernel              K;

public:
  typedef MotorcycleGraphTraits                               Geom_traits;
  typedef typename Geom_traits::Triangle_mesh                 Triangle_mesh;

  typedef typename Geom_traits::FT                            FT;

  typedef typename Geom_traits::Point_2                       Point_2;
  typedef typename Geom_traits::Segment_2                     Segment_2;
  typedef typename Geom_traits::Vector_2                      Vector_2;
  typedef typename Geom_traits::Ray_2                         Ray_2;

  typedef typename Geom_traits::Point_d                       Point;
  typedef typename Geom_traits::Segment_d                     Segment;
  typedef typename Geom_traits::Vector_d                      Vector;
  typedef typename Geom_traits::Ray_d                         Ray;

  typedef typename Geom_traits::Bbox_d                        Bbox;

  typedef Dictionary<Geom_traits>                             Dictionary;
  typedef Dictionary_entry<Geom_traits>                       Dictionary_entry;
  typedef typename Dictionary::DEC_it                         DEC_it;

  typedef typename Geom_traits::Barycentric_coordinates       Barycentric_coordinates;
  typedef typename Geom_traits::Face_location                 Face_location;

  typedef Motorcycle_impl_base<Geom_traits>                   Motorcycle;
  typedef boost::shared_ptr<Motorcycle>                       Motorcycle_ptr;
  typedef std::vector<Motorcycle_ptr>                         Motorcycle_container;

  typedef typename Motorcycle::Track                          Track;
  typedef typename Motorcycle::TPC_iterator                   TPC_iterator;

  typedef Motorcycle_priority_queue<Geom_traits>              Motorcycle_PQ;
  typedef Motorcycle_priority_queue_entry<Geom_traits>        Motorcycle_PQE;

  typedef typename Geom_traits::halfedge_descriptor           halfedge_descriptor;
  typedef typename Geom_traits::face_descriptor               face_descriptor;

  typedef boost::tuple<std::size_t, DEC_it, FT, DEC_it, FT>   Track_segment;
  typedef std::list<Track_segment>                            Track_segment_container;
  typedef boost::unordered_map<face_descriptor,
                               Track_segment_container>       Track_face_map;
  typedef typename Track_face_map::iterator                   TFM_iterator;

  // Access
  Motorcycle& motorcycle(const std::size_t id)
  {
    CGAL_precondition(id >= 0 && id < motorcycles.size());
    return *(motorcycles[id]);
  }
  const Motorcycle& motorcycle(const std::size_t id) const
  {
    CGAL_precondition(id >= 0 && id < motorcycles.size());
    return *(motorcycles[id]);
  }

  const Geom_traits& geom_traits() const { return gt; }

  // Constructor
  Motorcycle_graph(Triangle_mesh& mesh, const Geom_traits& gt = Geom_traits());

  // Functions
  void add_motorcycle(Motorcycle_ptr mc);
  void add_motorcycle(Motorcycle_ptr mc, std::size_t new_id);

  template<typename MotorcycleContainerIterator>
  void add_motorcycles(MotorcycleContainerIterator mit, MotorcycleContainerIterator last);

  /// \param fd face in which the segment belongs
  /// \param id the id of the motorcycle
  /// \param s the source of the oriented segment
  /// \param t the target of the oriented segment
  ///
  /// \return iterator in the tracking map
  TFM_iterator add_track_segment_to_map(face_descriptor fd, std::size_t id,
                                        DEC_it s, const FT time_at_s,
                                        DEC_it t, const FT time_at_t);

  /// \param p first point
  /// \param p_time time at first point
  /// \param q second point
  /// \param q_time time at second point
  ///
  /// \return new point and time at the new point
  std::pair<DEC_it, FT> compute_halving_point(const Motorcycle& mc,
                                              DEC_it p, const FT p_time,
                                              DEC_it q, const FT q_time);

  /// \param p first point
  /// \param p_time time at first point
  /// \param q second point
  /// \param q_time time at second point
  ///
  /// \return new point and time at the new point
  std::pair<DEC_it, FT> compute_middle_point(DEC_it p, const FT p_time,
                                             DEC_it q, const FT q_time);

  bool compute_motorcycle_next_path(Motorcycle& mc);
  void crash_motorcycle(Motorcycle& mc);
  void crash_motorcycles_with_same_source_and_direction();
  void drive_to_closest_target(Motorcycle& mc);

  void find_collision_at_tentative_track_destination(const Motorcycle& mc, const Motorcycle& fmc,
                                                     const FT fmc_visiting_time, Collision_information& tc) const;
  void find_collision_between_collinear_tracks(const Motorcycle& mc, const Segment_2& mcs,
                                               const Motorcycle& fmc, const Track_segment& fmc_track, const Segment_2& fmcs,
                                               const bool is_fmc_moving_on_track, Collision_information& tc) const;
  void find_collision_between_tracks(const Motorcycle& mc, const Segment_2& mcs,
                                     const Motorcycle& fmc, const Track_segment& fmc_track,
                                     const bool is_fmc_moving_on_track, Collision_information& tc) const;
  void find_collision_with_complete_track(Motorcycle& mc, const Segment_2& mcs,
                                          const Track_segment& fmc_track, Collision_information& tc);
  void find_collision_with_live_motorcycle(Motorcycle& mc, const Segment_2& mcs,
                                           const Motorcycle& fmc, Collision_information& tc);

  /// \return collision point (if any), time at the collision for `mc`, id of
  ///         the foreign intersecting motorcycle, time at the collision for
  ///         the foreign motorcycle.
  boost::tuple<DEC_it, FT, std::size_t, FT> find_collision(Motorcycle& mc);
  void generate_enclosing_face();
  bool has_motorcycle_reached_final_destination(const Motorcycle& mc) const;
  void initialize_motorcycles();
  bool has_motorcycle_reached_crashing_point(const Motorcycle& mc) const;

  template<typename MotorcycleContainerIterator>
  void trace_graph(MotorcycleContainerIterator mit, MotorcycleContainerIterator last);

  // Post-tracing checks
  bool is_valid() const;

  // Output
  void output_all_dictionary_points() const;
  void output_motorcycles_sources_and_destinations() const;

private:
  Geom_traits gt;

  Dictionary points; // All the points that will be used throughout the algorithm
  Motorcycle_container motorcycles;
  Motorcycle_PQ motorcycle_pq; // motorcycle priority queue

  bool using_enclosing_bbox; // indicates whether a mesh is passed input
  Triangle_mesh& mesh; // not 'const' in case we need to create it

  // map to keep in memory the completed tracks of the motorcycles for each face
  Track_face_map track_face_map;
};

// -----------------------------------------------------------------------------
template<typename MotorcycleGraphTraits>
Motorcycle_graph<MotorcycleGraphTraits>::
Motorcycle_graph(Triangle_mesh& mesh, const Geom_traits& gt)
  :
    gt(gt),
    points(),
    motorcycles(),
    motorcycle_pq(),
    using_enclosing_bbox(false),
    mesh(mesh),
    track_face_map()
{
  if(num_vertices(mesh) == 0)
  {
    std::cerr << " Warning: empty mesh in input" << std::endl;
    using_enclosing_bbox = true;
  }
  else
  {
    // Input must be a mesh with triangle faces
    CGAL_precondition(CGAL::is_triangle_mesh(mesh));
  }

  //@tmp disabled while I find out what to do with the "no mesh provided option"
  // The issue is that the points are identified by a location described with barycentric
  // coordinates. I guess, I could generate a bbox, then a triangle that includes
  // the box ? Pretty ugly, though...
  CGAL_assertion(!using_enclosing_bbox);
}

template<typename MotorcycleGraphTraits>
void
Motorcycle_graph<MotorcycleGraphTraits>::
add_motorcycle(Motorcycle_ptr mc)
{
  return add_motorcycle(mc, motorcycles.size());
}

template<typename MotorcycleGraphTraits>
void
Motorcycle_graph<MotorcycleGraphTraits>::
add_motorcycle(Motorcycle_ptr mc, std::size_t new_id)
{
  mc->set_id(new_id);

  boost::optional<Point>& destination_point = mc->initial_destination_point();
  boost::optional<Vector>& direction = mc->direction();

  if(destination_point == boost::none && direction == boost::none)
    std::cerr << "Warning: Neither destination nor direction are provided." << std::endl;

  motorcycles.push_back(mc);
}

template<typename MotorcycleGraphTraits>
template<typename MotorcycleContainerIterator>
void
Motorcycle_graph<MotorcycleGraphTraits>::
add_motorcycles(MotorcycleContainerIterator mit, MotorcycleContainerIterator last)
{
  if(!motorcycles.empty())
    std::cerr << "Warning: motorcycle container was already not empty when calling add_motorcycles()" << std::endl;

  motorcycles.reserve(motorcycles.size() + std::distance(mit, last));

  // unique motorcycle id, starting at motorcycles.size() in case we have
  // already added some motorcycles
  std::size_t counter = motorcycles.size();

  while(mit != last)
    add_motorcycle(*mit++, counter++);
}

template<typename MotorcycleGraphTraits>
typename Motorcycle_graph<MotorcycleGraphTraits>::TFM_iterator
Motorcycle_graph<MotorcycleGraphTraits>::
add_track_segment_to_map(face_descriptor fd, std::size_t id,
                         DEC_it s, const FT time_at_s,
                         DEC_it t, const FT time_at_t)
{
  CGAL_precondition(s->location().first == fd);
  CGAL_precondition(t->location().first == fd);
  CGAL_precondition(id >=0 && id < motorcycles.size());
  CGAL_precondition(time_at_s <= time_at_t);

  Track_segment tr = boost::make_tuple(id, s, time_at_s, t, time_at_t);
  Track_segment_container l;
  l.push_back(tr);

  std::pair<typename Motorcycle_graph<MotorcycleGraphTraits>::TFM_iterator, bool>
    is_insert_success = track_face_map.insert(std::make_pair(fd, l));

  if(!is_insert_success.second)
    is_insert_success.first->second.push_back(tr);

  return is_insert_success.first;
}

template<typename MotorcycleGraphTraits>
std::pair<typename Motorcycle_graph<MotorcycleGraphTraits>::DEC_it,
          typename Motorcycle_graph<MotorcycleGraphTraits>::FT>
Motorcycle_graph<MotorcycleGraphTraits>::
compute_halving_point(const Motorcycle& m, DEC_it p, const FT p_time,
                                           DEC_it q, const FT q_time)
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "*** " << std::endl;
  std::cout << " Computing halving point on motorcycle #" << m.id() << "'s track"
            << "points are: " << std::endl << "  " << *p << std::endl
                                           << "  " << *q << std::endl;
#endif
  CGAL_precondition(p != q);
  CGAL_precondition(p->location().first == q->location().first);

#ifdef CGAL_MOTORCYCLE_GRAPH_USE_ADVANCED_HALVING
  // interface with the halving data structure @todo
#else
  return compute_middle_point(p, p_time, q, q_time);
#endif
}

template<typename MotorcycleGraphTraits>
std::pair<typename Motorcycle_graph<MotorcycleGraphTraits>::DEC_it,
          typename Motorcycle_graph<MotorcycleGraphTraits>::FT>
Motorcycle_graph<MotorcycleGraphTraits>::
compute_middle_point(DEC_it p, const FT p_time, DEC_it q, const FT q_time)
{
  if(p->location().first != q->location().first)
  {
    std::cerr << "Error: middle point computation with different faces" << std::endl;
    // asserting because using p.loc().first is too dangerous if r is not guaranteed
    // to be on p's face
    CGAL_assertion(false);
  }

  const Barycentric_coordinates& p_coords = p->location().second;
  const Barycentric_coordinates& q_coords = q->location().second;

  Barycentric_coordinates middle_coords = CGAL::make_array(0.5*(p_coords[0] + q_coords[0]),
                                                           0.5*(p_coords[1] + q_coords[1]),
                                                           0.5*(p_coords[2] + q_coords[2]));
  Face_location middle_loc = std::make_pair(p->location().first, middle_coords);
  const FT time_at_r = 0.5 * (p_time + q_time);
  std::pair<DEC_it, bool> entry = points.insert(middle_loc, mesh);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "  New middle point: (" << entry.first->point()
                                       << ") at time: " << time_at_r << std::endl;
  std::cout << "    Location: " << p->location().first
            << " bc: " << middle_coords[0] << " "
                       << middle_coords[1] << " "
                       << middle_coords[2] << std::endl;
#endif

  return std::make_pair(entry.first, time_at_r);
}

template<typename MotorcycleGraphTraits>
bool
Motorcycle_graph<MotorcycleGraphTraits>::
compute_motorcycle_next_path(Motorcycle& mc)
{
  boost::tuple<bool, DEC_it, DEC_it, FT, bool> next_path =
    mc.compute_next_destination(points, mesh);

  if(!next_path.template get<0>()) // couldn't find a next path
    return false;

  const DEC_it& next_source = next_path.template get<1>();
  const DEC_it& next_destination = next_path.template get<2>();
  const FT time_at_next_destination = next_path.template get<3>();
  const bool is_destination_final = next_path.template get<4>();

  if(next_source != mc.position())
  {
    // if 'next source' is different from the current position, it should only
    // be a location change, not a position change
    CGAL_assertion(next_source->point() == mc.position()->point());

    // block the previous destination
    mc.position()->block();

    // If the next source is the current position, the point already knows the motorcycle
    next_source->add_motorcycle(mc.id(), mc.current_time());
  }

  // Add the next source as target, even if it is equal to the current position
  mc.add_target(next_source, mc.current_time());

  if(next_source != next_destination)
  {
    // No need to add the same information twice
    mc.add_target(next_destination, time_at_next_destination);
    next_destination->add_motorcycle(mc.id(), time_at_next_destination);
  }

  mc.source() = next_source;
  mc.time_at_source() = mc.current_time();
  mc.destination() = next_destination;
  mc.position() = mc.source();

  mc.set_destination_finality(is_destination_final);

  return true;
}

template<typename MotorcycleGraphTraits>
void
Motorcycle_graph<MotorcycleGraphTraits>::
crash_motorcycle(Motorcycle& mc)
{
  if(mc.is_crashed())
    return;

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~X" << std::endl;
  std::cout << "Crashing " << mc;
#endif

  mc.crash();

  // go through the next targets of the motorcycle and remove it from the list
  // of motorcycles that might reach the target
  typename Motorcycle::Target_point_container::iterator it = mc.targets().begin();
  typename Motorcycle::Target_point_container::iterator end = mc.targets().end();
  for(; it!=end; ++it)
  {
    DEC_it target_point = it->first;
    target_point->remove_motorcycle(mc.id());
  }
  mc.targets().clear(); // might be unnecessary @todo

  motorcycle_pq.erase(mc);
}

template<typename MotorcycleGraphTraits>
void
Motorcycle_graph<MotorcycleGraphTraits>::
crash_motorcycles_with_same_source_and_direction()
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "Checking for motorcycles with same sources and directions" << std::endl;
#endif

  // brute force, for now
  // A smarter version is to sort motorcycles by direction (slope),
  // and check for consecutive entries @todo
  std::size_t number_of_motorcycles = motorcycles.size();
  for(std::size_t mc_id=0; mc_id<number_of_motorcycles; ++mc_id)
  {
    Motorcycle& mc = motorcycle(mc_id);

    if(mc.source() == mc.destination() || mc.is_crashed())
      continue;

    for(std::size_t fmc_id=0; fmc_id<number_of_motorcycles; ++fmc_id)
    {
      Motorcycle& fmc = motorcycle(fmc_id);

      // Note: not ignoring crashed motorcycles in case of > 2 motorcycles with
      // same source and destination

      if(fmc.id() == mc.id() ||
         fmc.source() == fmc.destination() || // a degenerate track does not block anything
         mc.source() != fmc.source()) // must have identical sources
        continue;

      CGAL_assertion(mc.current_location().first == fmc.current_location().first);

      Point_2 bcs_mc_s(mc.source()->location().second[0], mc.source()->location().second[1]);
      Point_2 bcs_mc_d(mc.destination()->location().second[0], mc.destination()->location().second[1]);
      Point_2 bcs_fmc_d(fmc.destination()->location().second[0], fmc.destination()->location().second[1]);
      Point_2 bcs_fmc_s(fmc.source()->location().second[0], fmc.source()->location().second[1]);

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
      Vector_2 bcs_mc_v(bcs_mc_s, bcs_mc_d);
      Vector_2 bcs_fmc_v(bcs_fmc_s, bcs_fmc_d);

      FT mc_v_n = bcs_mc_v * bcs_mc_v;
      FT fmc_v_n = bcs_fmc_v * bcs_fmc_v;

      FT sp = gt.compute_scalar_product_2_object()(bcs_mc_v, bcs_fmc_v);

      std::cout << "SProduct: " << sp << std::endl;
      std::cout << "SProduct normalized " << sp * sp / (fmc_v_n * mc_v_n ) << std::endl;

      // @fixme hardcoded value, but numeric_limits::eps not small enough
      // due to the multiple intermediary computations
      if( CGAL::abs( 1 - sp * sp / (fmc_v_n * mc_v_n) ) < 1e-15)
      {
        std::cout << "Crashing degenerate motorcycles: "
                  << mc.id() << " and " << fmc.id() << std::endl;
        crash_motorcycle(mc);
        crash_motorcycle(fmc);
        break;
      }
#else
      // only aligned tracks block one another
      if(!gt.collinear_2_object()(bcs_mc_s, // == fmc.source()->point()
                                  bcs_mc_d, bcs_fmc_d))
        continue;

      std::cout << "Collinear tracks with the same source" << std::endl;

      // Moving away from each other from the same point is allowed.
      // use two calls to ordered_along_line() instead? @todo
      if(gt.angle_2_object()(bcs_mc_s, bcs_mc_d, bcs_fmc_s, bcs_fmc_d) == CGAL::ACUTE)
      {
        std::cout << "Crashing degenerate motorcycles: "
                  << mc.id() << " and " << fmc.id() << std::endl;
        crash_motorcycle(mc);
        crash_motorcycle(fmc);
        break;
      }
#endif
    }
  }
}

template<typename MotorcycleGraphTraits>
void
Motorcycle_graph<MotorcycleGraphTraits>::
drive_to_closest_target(Motorcycle& mc)
{
  CGAL_assertion(!mc.is_crashed());
  CGAL_assertion(!mc.targets().empty());

  DEC_it closest_target = mc.closest_target();

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>" << std::endl;
    std::cout << "Driving " << mc;
#endif

  mc.position() = closest_target;
  mc.current_time() = mc.targets().begin()->second;
  mc.track().insert(std::make_pair(closest_target, mc.current_time()));
  mc.remove_closest_target_from_targets();

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "  now at: (" << mc.position()->point() << ")" << std::endl;
#endif
}

template<typename MotorcycleGraphTraits>
void
Motorcycle_graph<MotorcycleGraphTraits>::
find_collision_at_tentative_track_destination(const Motorcycle& mc,
                                              const Motorcycle& fmc,
                                              const FT fmc_visiting_time,
                                              Collision_information& tc) const
{
  FT time_at_collision = mc.time_at_closest_target();
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "  /!\\ Tentative path collides with track: " << fmc.id()
            << " at its end. Time: " << time_at_collision << std::endl;
#endif

  // Compare with other possible collisions to keep the closest
  if(time_at_collision < tc.time_at_closest_collision)
  {
    tc.is_closest_collision_location_already_in_dictionary = true;

    tc.closest_collision = mc.closest_target();
    tc.time_at_closest_collision = time_at_collision;

    tc.foreign_mc_id = fmc.id();
    tc.foreign_time_at_closest_collision = fmc_visiting_time;
  }
}

template<typename MotorcycleGraphTraits>
void
Motorcycle_graph<MotorcycleGraphTraits>::
find_collision_between_collinear_tracks(const Motorcycle& mc,
                                        const Segment_2& mcs,
                                        const Motorcycle& fmc,
                                        const Track_segment& fmc_track,
                                        const Segment_2& fmcs,
                                        const bool is_fmc_moving_on_track,
                                        Collision_information& tc) const
{
  bool is_fmcs_degenerate = gt.is_degenerate_2_object()(fmcs);

  // Check if the tentative tracks actually intersect
  if((is_fmcs_degenerate &&
      !gt.collinear_are_ordered_along_line_2_object()(mcs.source(),
                                                      fmcs.source(),
                                                      mcs.target()) &&
       !gt.collinear_are_ordered_along_line_2_object()(mcs.source(),
                                                       fmcs.target(),
                                                       mcs.target())) ||
    !gt.do_intersect_2_object()(mcs, fmcs))
  {
    // No intersection
    std::cout << "  No intersection" << std::endl;
    return;
  }

  // Many different configurations exist, e.g. (_S is for source, _T for target):
  //  MC_S  ---- FMC_S ---- FMC_T ---- MC_T
  //  FMC_T ---- MC_S  ---- FMC_S ---- MC_T
  // etc.
  // If, on the ray MC_S->MC_T,
  // - if fmc_S is "before" MC_S, then it doesn't matter for MC
  // - if fmc_S is "after" MC_S, then it depends on the motorcycles' directions

  // We need a ray because fmcs' source can be outside of mc's tentative track
  Ray_2 mcr(mcs.source(), mcs.target());
  if(!mcr.has_on(fmcs.source()))
    return;

  // Compute the respective direction of the two motorcycles:
  // use two calls to collinear_are_aligned_along_line() instead of angle ? @todo
  bool are_motorcycles_moving_in_the_same_direction =
    (is_fmcs_degenerate ||
     gt.angle_2_object()(mcs.source(), mcs.target(),
                         fmcs.source(), fmcs.target()) == CGAL::ACUTE);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "  is degen: " << is_fmcs_degenerate << std::endl;
  std::cout << "  angle: " << gt.angle_2_object()(mcs.source(), mcs.target(),
                                                  fmcs.source(), fmcs.target()) << std::endl;
  std::cout << "  are motorcycles moving in the same direction: "
            << are_motorcycles_moving_in_the_same_direction << std::endl;
#endif

  FT time_at_collision = 0.;
  const DEC_it fmc_track_source = fmc_track.template get<1>();
  const FT time_at_fmc_track_source = fmc_track.template get<2>();
  const DEC_it fmc_track_destination = fmc_track.template get<3>();
  const FT time_at_fmc_track_destination = fmc_track.template get<4>();

  // The motorcycles move in the same direction ==> mc will impact fmcs' source
  if(are_motorcycles_moving_in_the_same_direction)
  {
    // The weird configuration of both motorcycles moving in the same direction
    // AND with the same source is checked at the very start, see function
    // 'crash_motorcycles_with_same_source_and_direction()'
    CGAL_assertion(is_fmcs_degenerate || mcs.source() != fmcs.source());

    // @todo find something nicer than computing time like that
    time_at_collision = mc.current_time() +
      CGAL::sqrt(CGAL::squared_distance(mc.position()->point(), fmc_track_source->point())) / mc.speed();
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "  Motorcycles #" << mc.id() << " crashes into the source of Motorcycle #"
              << fmc.id() << " at time: " << time_at_collision << std::endl;
#endif

    // Compare with other possible collisions to keep the closest
    if(time_at_collision < tc.time_at_closest_collision)
    {
      tc.is_closest_collision_location_already_in_dictionary = true;
      tc.time_at_closest_collision = time_at_collision;
      tc.closest_collision = fmc_track_source;
      tc.foreign_mc_id = fmc.id();
      tc.foreign_time_at_closest_collision = time_at_fmc_track_source;
    }
  }
  else // Motorcycles are moving towards each other
  {
    // Ignore the case where the motorcycles are moving away from each other
    // from the same point (not a configuration where they should crash)
    if(mcs.source() == fmcs.source())
      return;

    // If the foreign motorcycle is (also) moving, we return the middle point
    // that both motorcycles would reach at the same time.
    // Note that this point might not actually be reached by either motorcycle,
    // e.g. if each motorcycle reaches its destination first.
    if(is_fmc_moving_on_track)
    {
      // @todo, if speeds are ever allowed to change, the speed of fmc here
      // must be changed to the speed on the track segment 'fmc_track'
      time_at_collision = mc.current_time() +
        (CGAL::sqrt(CGAL::squared_distance(mc.position()->point(), fmc_track_source->point())) -
           fmc.speed() * (mc.current_time() - time_at_fmc_track_source)) /
             (mc.speed() + fmc.speed());
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
      std::cout << "  mc and fmc would meet at time: " << time_at_collision << std::endl;
#endif
      CGAL_postcondition(time_at_collision > mc.current_time());

      // Compare with other possible collisions to keep the closest
      if(time_at_collision < tc.time_at_closest_collision)
      {
        tc.is_closest_collision_location_already_in_dictionary = false;

        // @todo do something nicer to compute time ?
        const Vector_2 mcv(mcs);
        const FT ratio = (time_at_collision - mc.current_time()) /
                           (mc.time_at_closest_target() - mc.current_time());
        const Point_2 collision = mcs.source() + ratio * mcv;

        tc.closest_collision_location = std::make_pair(
          mc.current_location().first, CGAL::make_array(collision[0],
                                                        collision[1],
                                                        1. - collision[0] - collision[1]));
        tc.time_at_closest_collision = time_at_collision;

        tc.foreign_mc_id = fmc.id();
        tc.foreign_time_at_closest_collision = time_at_collision;
      }
    }
    // If fmc is not moving, then mc crashes into the final position of the track
    else
    {
      time_at_collision = mc.current_time() +
        CGAL::sqrt(CGAL::squared_distance(mc.position()->point(), fmc_track_destination->point())) / mc.speed();
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
      std::cout << "  mc crashes into fmc's final position at: " << time_at_collision << std::endl;
#endif

      // Compare with other possible collisions to keep the closest
      if(time_at_collision < tc.time_at_closest_collision)
      {
        tc.is_closest_collision_location_already_in_dictionary = true;

        tc.closest_collision = fmc_track_destination;
        tc.time_at_closest_collision = time_at_collision;

        tc.foreign_mc_id = fmc.id();
        tc.foreign_time_at_closest_collision = time_at_fmc_track_destination;
      }
    }
  }
}

template<typename MotorcycleGraphTraits>
void
Motorcycle_graph<MotorcycleGraphTraits>::
find_collision_between_tracks(const Motorcycle& mc,
                              const Segment_2& mcs,
                              const Motorcycle& fmc,
                              const Track_segment& fmc_track,
                              const bool is_fmc_moving_on_track,
                              // below are out parameters
                              Collision_information& tc) const
{
  // Non degenerate mc segment
  CGAL_precondition(mcs.source() != mcs.target());

  const DEC_it fmc_track_source = fmc_track.template get<1>();
  const FT time_at_fmc_track_source = fmc_track.template get<2>();
  const DEC_it fmc_track_destination = fmc_track.template get<3>();
  const FT time_at_fmc_track_destination = fmc_track.template get<4>();

  // Both tracks must be on the same face
  CGAL_precondition(fmc_track_source->location().first == fmc_track_destination->location().first);

  FT time_at_collision = 0.;

  const Point_2 s = gt.construct_point_2_object()(fmc_track_source->location().second[0],
                                                  fmc_track_source->location().second[1]);
  const Point_2 t = gt.construct_point_2_object()(fmc_track_destination->location().second[0],
                                                  fmc_track_destination->location().second[1]);
  const Segment_2 fmcs = gt.construct_segment_2_object()(s, t);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "¤¤ Checking collision with track of motorcycle #" << fmc.id() << std::endl;
  std::cout << " + source: " << std::endl << *fmc_track_source << std::endl;
  std::cout << " + target: " << std::endl << *fmc_track_destination << std::endl;
  std::cout << " BCS Segment: " << s << " || " << t << std::endl;
#endif

  // Ignore the degenerate case of a degenerate fmc track starting at the same source
  bool is_fmcs_degenerate = gt.is_degenerate_2_object()(fmcs);
  if(is_fmcs_degenerate && mcs.source() == fmcs.source())
  {
    std::cout << "degenerate fmc and mcs.source() == fmcs.source()" << std::endl;
    return;
  }

  // Detect whether the motorcycles share the same supporting line.
  // Note that we know that 'mcs' is not degenerate.
  // @todo should collinearity checks handle a tolerance ?
  if(gt.collinear_2_object()(mcs.source(), mcs.target(), fmcs.source()) &&
     gt.collinear_2_object()(mcs.source(), mcs.target(), fmcs.target()))
  {
    std::cout << "  /!\\ Tracks are aligned" << std::endl;
    return find_collision_between_collinear_tracks(mc, mcs, fmc, fmc_track, fmcs, is_fmc_moving_on_track, tc);
  }

  // --- From here on, the tracks are not collinear ---

  // The next two "if" are two checks to avoid having to compute the intersection
  // - #1: Check if the current position of mc is a known intersection with the foreign track
  // - #2: Check if the closest target of mc is a known intersection with the foreign track

  // Check #1
  if(mc.position()->has_motorcycle(fmc.id(), time_at_fmc_track_source, time_at_fmc_track_destination))
  {
    // Ignore this intersection: since we are seeking collision in the tentative track,
    // it means that the position was not blocked
    return;
  }

  // Check #2
  FT fmc_visiting_time;
  if(mc.closest_target()->has_motorcycle(fmc.id(), time_at_fmc_track_source,
                                         time_at_fmc_track_destination, fmc_visiting_time))
  {
    return find_collision_at_tentative_track_destination(mc, fmc, fmc_visiting_time, tc);
  }

  // --- The general case: the intersection must be computed ---

  if(is_fmcs_degenerate)
  {
    CGAL_assertion(!gt.do_intersect_2_object()(mcs, fmcs));

    // If there is an intersection, it should have been caught by the first part
    // of that function, branching: "collinear > moving in the same direction"
    std::cout << "  No intersection" << std::endl;
    return;
  }

  if(!gt.do_intersect_2_object()(mcs, fmcs))
  {
    // No intersection, move to the next motorcycle
    std::cout << "  No intersection" << std::endl;
    return;
  }

  // Below computes the intersection in the barycentric coordinates system
  Point_2 collision = internal::robust_intersection<Geom_traits>(mcs, fmcs, gt);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "  /!\\ collision between motorcycles #" << mc.id() << " and #" << fmc.id() << std::endl;
#endif

  // Convert it to a location in the ambiant dimension
  Barycentric_coordinates coords = CGAL::make_array(collision[0], collision[1],
                                                    1. - collision[0] - collision[1]);
  Face_location collision_location = std::make_pair(mc.current_location().first, coords);

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
  // @todo snap here to any existing point instead of border ?
  // Harder than it seems: how to snap close points that are on different sides
  // of an edge, etc. ?
  CGAL::Polygon_mesh_processing::internal::snap_location_to_border<Triangle_mesh>(collision_location);
#endif

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "Collision location: " << collision_location.first << " bc: "
            << collision_location.second[0] << " " << collision_location.second[1] << " " << collision_location.second[2] << std::endl;
#endif

  // Although we might not have known that these two tracks do intersect,
  // their intersection might be a point that has already been used
  std::pair<DEC_it, bool> is_already_in_dictionary = points.find(collision_location);
  if(is_already_in_dictionary.second)
  {
    DEC_it collision_point = is_already_in_dictionary.first;

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "Already in the dictionary at: " << &*collision_point << std::endl << *collision_point << std::endl;
#endif

    // The point already exists, check if either 'mc' or 'fmc' are already visiting it
    // so that we can extract a coherent time
    if(collision_point->has_motorcycle(mc.id(), mc.current_time(), mc.time_at_closest_target()))
    {
      if(collision_point == mc.closest_target())
      {
        time_at_collision = mc.time_at_closest_target();
      }
      else
      {
        // The tentative track of 'mc' can only be intersected at a known point that has 'mc'
        // if that known point is the closest target:
        // - it can't be 'mc.source()' otherwise we sjpimd have found the intersection
        //   on the previous call to find_collision() when the current source
        //   was the closest target
        // - it can't be another point otherwise closest_target is not actually the closest
        CGAL_assertion(false);
      }
    }
    else // collision_point is a known point but has not (yet) visited by 'mc'
    {
      // No choice but to compute the visiting time
      time_at_collision = mc.current_time() +
        CGAL::sqrt(CGAL::squared_distance(mc.position()->point(),
                                          collision_point->point())) / mc.speed();

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
      // Although we have found an existing point at the location of the intersection,
      // this point was neither the source or the closest target of 'mc'.

      // We check that the times are not identical to an existing point on the
      // tentative track, otherwise we have a problem because we have two points
      // that already exist and are absurdly close...
      //
      // Asserting it for now: need a global snapping method to handle this type
      // of problem...
      CGAL_assertion(time_at_collision != mc.current_time());
      CGAL_assertion(time_at_collision != mc.time_at_closest_target());
#endif
    }

    // Compare with other collisions to keep the closest
    if(time_at_collision < tc.time_at_closest_collision)
    {
      tc.is_closest_collision_location_already_in_dictionary = true;

      tc.closest_collision = collision_point;
      tc.time_at_closest_collision = time_at_collision;

      tc.foreign_mc_id = fmc.id();
      FT fmc_visiting_time;
      if(collision_point->has_motorcycle(fmc.id(), time_at_fmc_track_source,
                                         time_at_fmc_track_destination, fmc_visiting_time))
      {
        // The collision point is visited by 'fmc' at time 'fmc_visiting_time'
        tc.foreign_time_at_closest_collision = fmc_visiting_time;
      }
      else  // collision_point is a known point but has not (yet) visited by 'fmc'
      {
        // No choice but to compute the time
        tc.foreign_time_at_closest_collision = time_at_fmc_track_source +
          CGAL::sqrt(CGAL::squared_distance(fmc_track_source->point(),
                                            mc.closest_target()->point())) / fmc.speed();

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
        // Although we have found an existing point at the location of the intersection,
        // this point was not known on the foreign track

        // We check that the times are not identical to an existing point on the
        // foreign track, otherwise we have a problem because we have two points
        // that already exist and are absurdly close...
        //
        // Asserting it for now: need a global snapping method to handle this type
        // of problem...
        CGAL_assertion(!fmc.has_target_at_time(tc.foreign_time_at_closest_collision).second);
#endif
      }
    }
  }
  else // the collision location has never been seen before!
  {
    Point collision_point = CGAL::Polygon_mesh_processing::internal::loc_to_point(collision_location, mesh);
    time_at_collision = mc.current_time() +
      CGAL::sqrt(CGAL::squared_distance(mc.position()->point(), collision_point)) / mc.speed();

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "Location never seen before, corresponds to point "
              << collision_point << " at time: " << time_at_collision << std::endl;
#endif

    // If the collision time is greater than the time at the destination, then
    // the collision is not interesting
    CGAL_assertion(mc.targets().rbegin()->first == mc.destination());
    if(time_at_collision > mc.targets().rbegin()->second)
      return;

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
    // Different locations but same times... to be handled @todo
    if(time_at_collision == mc.time_at_source())
    {
      CGAL_assertion(false); // @todo
    }
    else if(time_at_collision == mc.time_at_closest_target())
    {
      // 'time_at_collision' says that we are intersecting at the closest target,
      // but it's a new location. Ignore the new location and use the closest
      // target as location for the collision

      FT ftime = time_at_fmc_track_source +
        CGAL::sqrt(CGAL::squared_distance(fmc_track_source->point(),
                                          mc.closest_target()->point())) / fmc.speed();
      CGAL_assertion(!fmc.has_target_at_time(ftime).second);

      return find_collision_at_tentative_track_destination(mc, fmc, ftime, tc);
    }
#endif

    // Compare with other collisions to keep the closest
    if(time_at_collision < tc.time_at_closest_collision)
    {
      tc.is_closest_collision_location_already_in_dictionary = false;

      tc.closest_collision_location = collision_location;
      tc.time_at_closest_collision = time_at_collision;

      tc.foreign_mc_id = fmc.id();
      tc.foreign_time_at_closest_collision = time_at_fmc_track_source +
        CGAL::sqrt(CGAL::squared_distance(fmc_track_source->point(),
                                          collision_point)) / fmc.speed();

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
      // Although we haven't found an existing point at the location of the intersection
      // this point might be at the same time from the source of the track
      std::pair<TPC_iterator, bool> res =
        fmc.has_target_at_time(tc.foreign_time_at_closest_collision);

      if(res.second)
      {
        // for the times to be equal, the points should be very close
        TPC_iterator target_point = res.first;
        CGAL_assertion(target_point->second == tc.foreign_time_at_closest_collision);
        DEC_it alternate_collision = target_point->first;
        CGAL_assertion(CGAL::squared_distance(alternate_collision->point(), collision_point)
                         < std::numeric_limits<FT>::epsilon());

        // @todo should I recompute time_at_collision, re-check vs time_at_closest_collision
        // and all that jazz... ? Assert for now...
        // '<=' since 'time_at_closest_collision' is now 'time_at_collision'
        CGAL_assertion(mc.current_time() + CGAL::sqrt(CGAL::squared_distance(
          mc.position()->point(), alternate_collision->point())) / mc.speed()
            <= tc.time_at_closest_collision);

        // for now, simply saying that the collision is now that point
        tc.is_closest_collision_location_already_in_dictionary = true;
        tc.closest_collision = alternate_collision;
      }
#endif
    }
  }
}

template<typename MotorcycleGraphTraits>
void
Motorcycle_graph<MotorcycleGraphTraits>::
find_collision_with_complete_track(Motorcycle& mc, const Segment_2& mcs,
                                   const Track_segment& fmc_track,
                                   // below are out parameters
                                   Collision_information& tc)
{
  const std::size_t fmc_id = fmc_track.template get<0>();
  const Motorcycle& fmc = motorcycle(fmc_id);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
std::cout << "¤ Checking for intersection with the complete track of motorcycle " << fmc.id() << std::endl;
#endif

  // 'false' because the motorcycle is not moving on that track
  return find_collision_between_tracks(mc, mcs, fmc, fmc_track, false, tc);
}

template<typename MotorcycleGraphTraits>
void
Motorcycle_graph<MotorcycleGraphTraits>::
find_collision_with_live_motorcycle(Motorcycle& mc, const Segment_2& mcs,
                                    const Motorcycle& fmc,
                                    // below are out parameters
                                    Collision_information& tc)
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
std::cout << "¤ Checking for intersection with live motorcycle #" << fmc.id() << std::endl;
#endif

  if(// the motorcycles must be different
     mc.id() == fmc.id() ||
     // the motorcycles must be in the same face
     mc.current_location().first != fmc.current_location().first ||
     // the foreign motorcycle must be in motion
     fmc.is_crashed())
  {
    std::cout << " ignoring fmc... " << std::endl;
    std::cout << "  > motorcycles #" << mc.id() << " and #" << fmc.id() << std::endl;
    std::cout << "  > faces: " << mc.current_location().first << " and " << fmc.current_location().first << std::endl;
    std::cout << "  > crashed status: " << fmc.is_crashed() << std::endl;
    return;
  }

  Track_segment fmc_track = boost::make_tuple(fmc.id(), fmc.source(), fmc.time_at_source(),
                                              fmc.closest_target(), fmc.time_at_closest_target());

  // 'true' because fmc is currently moving on that track
  return find_collision_between_tracks(mc, mcs, fmc, fmc_track, true, tc);
}

// search for a possible collision with another motorcycle between the current
// position of mc and the next target
template<typename MotorcycleGraphTraits>
boost::tuple<typename Motorcycle_graph<MotorcycleGraphTraits>::DEC_it,
             typename Motorcycle_graph<MotorcycleGraphTraits>::FT,
             std::size_t,
             typename Motorcycle_graph<MotorcycleGraphTraits>::FT>
Motorcycle_graph<MotorcycleGraphTraits>::
find_collision(Motorcycle& mc)
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "Checking for collisions on motorcycle #" << mc.id() << "'s track" << std::endl
            << "Currently on face: " << mc.current_location().first << std::endl;
#endif

  CGAL_precondition(!mc.is_crashed());
  CGAL_precondition(!mc.targets().empty());

  // A bunch of output parameters are regrouped into the 'Collision_information' struct,
  // which describes the best (closest to mc.position()) tentative collision.
  Collision_information tc;

  // The motorcycles must be on the same face
  CGAL_precondition(mc.current_location().first == mc.closest_target()->location().first);

  // Use the barycentric coordinate systems to compute intersections
  const Point_2 s = gt.construct_point_2_object()(mc.current_location().second[0],
                                                  mc.current_location().second[1]);
  const Point_2 t = gt.construct_point_2_object()(mc.closest_target()->location().second[0],
                                                  mc.closest_target()->location().second[1]);
  const Segment_2 mc_tentative_track = gt.construct_segment_2_object()(s, t);

  std::cout << "MC tentative track: " << std::endl << *(mc.position()) << std::endl
                                                   << *(mc.closest_target()) << std::endl;
  std::cout << "BCS Segment: " << std::endl << s << std::endl << t << std::endl;

  // A degenerate tentative track has no interesting collisions
  if(mc_tentative_track.is_degenerate())
    return boost::make_tuple(tc.closest_collision, tc.time_at_closest_collision,
                             tc.foreign_mc_id, tc.foreign_time_at_closest_collision);

  // Checking for intersection is done in two steps:
  // - Check with complete tracks in the face
  // - Check the motorcycles that are currently traveling in the face

  // Step 1: check complete tracks
  const face_descriptor mc_face = mc.current_location().first;
  TFM_iterator it = track_face_map.find(mc_face);
  if(it != track_face_map.end())
  {
    const Track_segment_container& face_tracks = it->second;

    typename Track_segment_container::const_iterator tl_it = face_tracks.begin();
    typename Track_segment_container::const_iterator tl_end = face_tracks.end();
    for(; tl_it!=tl_end; ++tl_it)
    {
      const Track_segment& track = *tl_it;
      find_collision_with_complete_track(mc, mc_tentative_track, track, tc);
    }
  }

  // Step 2: check incomplete tracks (path of a motorcycle currently moving in the same face)
  std::size_t number_of_motorcycles = motorcycles.size();
  for(std::size_t fmc_id = 0; fmc_id<number_of_motorcycles; ++fmc_id)
  {
    Motorcycle& fmc = motorcycle(fmc_id);
    find_collision_with_live_motorcycle(mc, mc_tentative_track, fmc, tc);
  }

  // Now, the closest collision (if any) is known
  if(tc.foreign_mc_id != std::size_t(-1) && // there is a collision
     !tc.is_closest_collision_location_already_in_dictionary)
  {
    // Insert the collision point in the dictionary. Motorcycle info will be added later.
    std::pair<DEC_it, bool> entry = points.insert(tc.closest_collision_location, mesh);
    tc.closest_collision = entry.first;
  }

  return boost::make_tuple(tc.closest_collision, tc.time_at_closest_collision,
                           tc.foreign_mc_id, tc.foreign_time_at_closest_collision);
}

template<typename MotorcycleGraphTraits>
void
Motorcycle_graph<MotorcycleGraphTraits>::
generate_enclosing_face()
{
  // generate a bbox that includes all known positions and all crash points
  // 2D only for now @todo
  CGAL_precondition(Geom_traits::dimension == 2);
  Bbox bbox;

  std::size_t number_of_motorcycles = motorcycles.size();
  for(std::size_t mc_id = 0; mc_id<number_of_motorcycles; ++mc_id)
  {
    Motorcycle& mc = motorcycle(mc_id);
    bbox += mc.initial_source_point().bbox();

    if(mc.initial_destination_point() != boost::none)
      bbox += mc.initial_destination_point()->bbox();

    // this part is brute force for now, but can be done in O(nlogn) by sorting
    // according to the slopes (farthest intersections happen when the slopes
    // of the motorcycles are close)
    for(std::size_t fmc_id = 0; fmc_id<number_of_motorcycles; ++fmc_id)
    {
      // segment - segment, segment - ray, or ray-ray intersections @todo
    }
  }

  // Slightly increase the size of the bbox to strictly contain all points

  // Manually create the mesh with Euler operations
}

template<typename MotorcycleGraphTraits>
bool
Motorcycle_graph<MotorcycleGraphTraits>::
has_motorcycle_reached_final_destination(const Motorcycle& mc) const
{
  return mc.is_motorcycle_destination_final();
}

template<typename MotorcycleGraphTraits>
void
Motorcycle_graph<MotorcycleGraphTraits>::
initialize_motorcycles()
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "Initialize motorcycles" << std::endl;
#endif

  // if no mesh has been given in input, generate a mesh made of a single quad face
  // that contains all the interesting motorcycle interactions (crashes)
  if(using_enclosing_bbox)
    generate_enclosing_face();

  typedef CGAL::internal::P2_or_P3_to_P3<Triangle_mesh>                 P2_or_P3_to_P3;
  typedef CGAL::P2_to_P3_VPM<Triangle_mesh>                             VPM;
  VPM vpm(mesh);

  typedef CGAL::AABB_face_graph_triangle_primitive<Triangle_mesh, VPM>  AABB_face_graph_primitive;
  typedef CGAL::AABB_traits<K, AABB_face_graph_primitive>               AABB_face_graph_traits;
  CGAL::AABB_tree<AABB_face_graph_traits> tree;

  CGAL::Polygon_mesh_processing::build_aabb_tree(
    mesh, tree, CGAL::Polygon_mesh_processing::parameters::vertex_point_map(vpm));

  std::size_t number_of_motorcycles = motorcycles.size();
  for(std::size_t mc_id = 0; mc_id<number_of_motorcycles; ++mc_id)
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "      _" << std::endl;
    std::cout << "    D/_" << std::endl;
    std::cout << "    /(__`=-/" << std::endl;
    std::cout << "  (o)     (o)" << std::endl;
    std::cout << "Initializing motorcycle #" << mc_id << std::endl;
#endif

    Motorcycle& mc = motorcycle(mc_id);
    const FT speed = mc.speed();
    boost::optional<Vector>& direction = mc.direction();

    // Add the source to the dictionary
    const Point& ini_source_point = mc.initial_source_point();
    const FT time_at_source = mc.current_time();

    // An AABB tree is a 3D structure, so we need to convert the point to a Point_3.
    // If the point is already a Point_3, this doesn't do anything.
    // @todo handle weird point types
    P2_or_P3_to_P3 to_p3;
    const typename P2_or_P3_to_P3::Point_3& source_point = to_p3(ini_source_point);

    const Face_location source_location =
      CGAL::Polygon_mesh_processing::locate(source_point, tree, mesh,
        CGAL::Polygon_mesh_processing::parameters::vertex_point_map(vpm));

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
    // @todo snap to border here ?
#endif

    std::pair<DEC_it, bool> source = points.insert(source_location,
                                                   ini_source_point,
                                                   mc_id, time_at_source);
    mc.source() = source.first;
    mc.position() = source.first;

    // @todo if source or destination is on a border of the mesh, we must find
    // a common face to both...

    // add the target to the dictionary
    boost::optional<Point>& opt_destination_point = mc.initial_destination_point();
    DEC_it destination;
    FT time_at_destination;

    if(opt_destination_point == boost::none) // destination was not provided
    {
      boost::tuple<bool, DEC_it, DEC_it, FT, bool> res = mc.compute_next_destination(points, mesh);
      if(!res.template get<0>())
      {
        // Couldn't find an initial destination ==> the motorcycle instantly crashes
        destination = mc.source();
        time_at_destination = mc.current_time();
      }
      else // found a destination
      {
        // the location algorithm might change the source to make sure that source
        // and destination are on the same face
        if(mc.source() != res.template get<1>())
        {
          std::cerr << "Warning: source has changed!" << std::endl
                    << "Previously: " << std::endl << *(mc.source()) << std::endl
                    << "Now: " << std::endl << *(res.template get<1>()) << std::endl;

          // The source change must only be a change of face descriptor, not of position
          CGAL_assertion(mc.source()->point() == res.template get<1>()->point());

          // if the source point wasn't previously in the dictionary, it can be cleaned off
          if(source.second)
          {
            // make sure that the only motorcycle visiting that point is 'mc'
            CGAL_assertion(source.first->visiting_motorcycles().size() == 1);
            points.erase(source.first);
          }
          // --- WARNING: 'source' SHOULD NOT BE USED ANYMORE FROM HERE ON ---

          mc.source() = res.template get<1>();
          mc.position() = mc.source();
          mc.source()->add_motorcycle(mc_id, time_at_source);
        }

        destination = res.template get<2>();
        time_at_destination = res.template get<3>();

        mc.set_destination_finality(res.template get<4>());
      }

      destination->add_motorcycle(mc_id, time_at_destination);
      opt_destination_point = destination->point();
    }
    else // destination is known, only need to compute the time of arrival
    {
      const Point& destination_point = *opt_destination_point;

      // source and destination should be on the same face, so we use the overload
      // where a face is provided @fixme find common face
      const Face_location destination_location =
        CGAL::Polygon_mesh_processing::locate(source_location.first,
                                              destination_point, mesh);

      // @todo should this  be computed by the tracer (?)
      time_at_destination = time_at_source +
        CGAL::sqrt(CGAL::squared_distance(ini_source_point, destination_point)) / speed;

      std::pair<DEC_it, bool> destination_entry = points.insert(destination_location,
                                                                destination_point,
                                                                mc_id, time_at_destination);
      destination = destination_entry.first;
    }

    mc.destination() = destination;

    // Sanity checks:
    // - source and destination must be on the same face
    CGAL_postcondition(mc.source()->location().first == mc.destination()->location().first);
    CGAL_postcondition(mc.destination()->location().second[0] >= 0. && mc.destination()->location().second[0] <= 1. &&
                       mc.destination()->location().second[1] >= 0. && mc.destination()->location().second[1] <= 1. &&
                       mc.destination()->location().second[2] >= 0. && mc.destination()->location().second[2] <= 1.);

    // Initialize the motorcycle target queue
    mc.targets().insert(std::make_pair(mc.source(), time_at_source));

    if(mc.source() != mc.destination())
      mc.targets().insert(std::make_pair(destination, time_at_destination));

    // this is useful to not get an empty track when sour=dest
    // but it creates duplicates @fixme
    mc.track().insert(std::make_pair(mc.source(), mc.current_time()));

    // Fill direction if needed
    if(direction == boost::none)
    {
      mc.direction() = Vector(mc.source()->point(), mc.destination()->point());
      std::cout << "Computing direction from destination: " << *(mc.direction()) << std::endl;
    }

    // sanity check: (destination - source) is collinear with the direction
    Ray r(mc.source()->point(), *(mc.direction()));
    if(!r.has_on(mc.destination()->point()))
    {
      std::cerr << "Error: Incompatible destination and direction: " << std::endl
                << "- destination: " << mc.destination()->point() << std::endl
                << "- direction: " << *(mc.direction()) << std::endl;
//      CGAL_assertion(false); // @tmp (this is usually wrong due to numerical errors, need an "almost_has_on")
    }
  }
}

template<typename MotorcycleGraphTraits>
bool
Motorcycle_graph<MotorcycleGraphTraits>::
has_motorcycle_reached_crashing_point(const Motorcycle& mc) const
{
  return  // mc has reached the track of a foreign motorcycle
         (mc.has_reached_blocked_point() ||
          // multiple motorcycles will reach mc's current position at the same time
          mc.has_reached_simultaneous_collision_point());
}

template<typename MotorcycleGraphTraits>
template<typename MotorcycleContainerIterator>
void
Motorcycle_graph<MotorcycleGraphTraits>::
trace_graph(MotorcycleContainerIterator mit, MotorcycleContainerIterator last)
{
  add_motorcycles(mit, last);
  initialize_motorcycles();
  motorcycle_pq.initialize(motorcycles);

#ifdef CGAL_MOTORCYCLE_GRAPH_OUTPUT
  output_motorcycles_sources_and_destinations();
#endif

  // this can only happen at the beginning, simpler to get it out the way immediately
  crash_motorcycles_with_same_source_and_direction();

  while(!motorcycle_pq.empty())
  {
    // get the earliest available event
    Motorcycle_PQE pqe = motorcycle_pq.top();
    Motorcycle& mc = pqe.motorcycle();

    // move the motorcycle to the target (which becomes the confirmed position)
    drive_to_closest_target(mc);

    if(mc.position() == mc.destination())
    {
      // Add the track source -- destination to the track map
      add_track_segment_to_map(mc.current_location().first, mc.id(),
                               mc.source(), mc.time_at_source(),
                               mc.destination(), mc.current_time());

      if(has_motorcycle_reached_final_destination(mc) ||
         has_motorcycle_reached_crashing_point(mc))
      {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
        std::cout << "Reached motorcycle's crashing point: " << std::endl
                  << " - final destination: " << has_motorcycle_reached_final_destination(mc) << std::endl
                  << " - blocked: " << mc.has_reached_blocked_point() << std::endl
                  << " - simultaneous collision: " << mc.has_reached_simultaneous_collision_point() << std::endl;
#endif
        crash_motorcycle(mc);
      }
      // not crashing, try to compute the next path
      else
      {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
        std::cout << "Reached destination: " << mc.destination()->point();
        std::cout << " Now computing motorcycle's next path..." << std::endl;
#endif
        if(compute_motorcycle_next_path(mc))
        {
          // a new path was found and set up, update the queue and continue
          motorcycle_pq.update(mc);

          // The purpose of this 'goto' is to avoid blocking the current position
          // when moving to the next path: the idea is to have the next target
          // be exactly the point at which we are. This is done to then analyze
          // the new path for potential collisions, similarly to how the first target
          // of a new motorcycle is its source.
          //
          // Since the next item is exactly the point at which we are, it'll be
          // first in queue and it is not dangerous to not block that point.
          // (Note that we could imagine having at the same time multiple motorcycles
          // whose closest target is their current position, but not blocking here
          // does not create issues either.)
          //
          // (It might as well be a "break;", but I want to print the priority queue.)
          goto next_item;
        }
        else
        {
          // couldn't find a next destination, crash the motorcycle
          crash_motorcycle(mc);
        }
      }
    }
    // motorcycle has not reached its destination, but still might be crashing
    else if(has_motorcycle_reached_crashing_point(mc) &&
            // only to prevent multiple motorcycles starting from the same source
            // (but with different directions) from blocking each other
            mc.current_time() != 0.)
    {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
      std::cout << "Reached crashing point: " << std::endl
                << " - blocked: " << mc.has_reached_blocked_point() << std::endl
                << " - simultaneous collision: " << mc.has_reached_simultaneous_collision_point() << std::endl;
#endif
      // Add the track source -- crash position to the track map
      add_track_segment_to_map(mc.current_location().first, mc.id(),
                               mc.source(), mc.time_at_source(),
                               mc.position(), mc.current_time());

      crash_motorcycle(mc);
    }
    // the motorcycle can continue without issue towards its destination
    else
    {
      // check for potential foreign tracks intersecting the next move of mc
      boost::tuple<DEC_it, FT, int, FT> res = find_collision(mc);
      DEC_it collision_point = res.template get<0>();
      const FT time_at_collision_point = res.template get<1>();
      const std::size_t foreign_motorcycle_id = res.template get<2>();
      const FT foreign_time_at_collision_point = res.template get<3>();

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
      std::cout << "--------" << std::endl << "Find collision results:" << std::endl << "--------" << std::endl;
      if(foreign_motorcycle_id != std::size_t(-1))
      {
        std::cout << " - collision_point: " << std::endl << *(collision_point) << std::endl
                  << " - motorcycle: #" << mc.id() << std::endl
                  << " - foreign_motorcycle: #" << foreign_motorcycle_id << std::endl
                  << " - time_at_collision_point: " << time_at_collision_point << std::endl
                  << " - foreign_time_at_collision_point: " << foreign_time_at_collision_point << std::endl;
      }
      else
      {
        std::cout << " No collision! " << std::endl;
      }
#endif
      // Treat the potential intersection

      if(// there is intersection
         foreign_motorcycle_id != std::size_t(-1) &&
         // the impact is closer than the next target
         time_at_collision_point <= mc.time_at_closest_target() &&
         // the collision is not the next target of 'mc' or the foreign track
         // does not know this collision point yet
         (collision_point != mc.closest_target() ||
          !collision_point->has_motorcycle(foreign_motorcycle_id,
                                           foreign_time_at_collision_point)))
      {
        if(!collision_point->has_motorcycle(mc.id(), time_at_collision_point))
        {
          // Call the halving structure to create a new point
          std::pair<DEC_it, FT> halving_entity =
            compute_halving_point(mc, mc.position(), mc.current_time(),
                                  collision_point, time_at_collision_point);
          DEC_it halving_point = halving_entity.first;
          const FT time_at_halving_point = halving_entity.second;

          // Degeneracies should have been caught before
          CGAL_postcondition(halving_point != mc.position() &&
                             halving_point != collision_point);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
          std::cout << "Adding collision: " << &*collision_point
                    << " and halving: " << &*halving_point
                    << " to motorcycle #" << mc.id() << std::endl;
#endif
          mc.add_target(collision_point, time_at_collision_point);
          mc.add_target(halving_point, time_at_halving_point);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
          std::cout << "Adding motorcycle #" << mc.id()
                    << " to collision: " << &*collision_point
                    << " and havling: " << &*halving_point << std::endl;
#endif
          collision_point->add_motorcycle(mc.id(), time_at_collision_point);
          halving_point->add_motorcycle(mc.id(), time_at_halving_point);

          CGAL_postcondition(mc.has_target_at_time(collision_point, time_at_collision_point));
        }

        // @todo factorize this a bit so we don't have multiple calls to "has_motorcycle"
        // followed by "add_motorcycle": this is all in log(n) complexity (admittedly,
        // log(n) on something that is unlikely to contain more than 2 elements, but still)
        Motorcycle& foreign_mc = motorcycle(foreign_motorcycle_id);
        if(!collision_point->has_motorcycle(foreign_motorcycle_id,
                                            foreign_time_at_collision_point))
        {
          // it is useful to know that the collision point is on the foreign track,
          // even if the collision point is on the confirmed part of the track
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
          std::cout << "Adding #" << foreign_motorcycle_id << " to collision point: " << &*collision_point << std::endl;
#endif
          collision_point->add_motorcycle(foreign_motorcycle_id, foreign_time_at_collision_point);

          if(// the collision point is not on the confirmed track for the foreign mc
             foreign_time_at_collision_point > foreign_mc.current_time())
          {
            // Call the halving structure to create a new point
            std::pair<DEC_it, FT> foreign_halving_entity =
              compute_halving_point(foreign_mc,
                                    foreign_mc.position(), foreign_mc.current_time(),
                                    collision_point, foreign_time_at_collision_point);
            DEC_it foreign_halving_point = foreign_halving_entity.first;
            const FT foreign_time_at_halving_point = foreign_halving_entity.second;

            // Degeneracies should have been caught before
            CGAL_postcondition(foreign_halving_point != foreign_mc.position() &&
                               foreign_halving_point != collision_point);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
            std::cout << "adding collision: " << &*collision_point
                      << " and halving: " << &*foreign_halving_point
                      << " to motorcycle #" << foreign_mc.id() << std::endl;
#endif
            foreign_mc.add_target(collision_point, foreign_time_at_collision_point);
            foreign_mc.add_target(foreign_halving_point, foreign_time_at_halving_point);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
            std::cout << "adding motorcycle #" << foreign_mc.id()
                      << " to halving point: " << &*foreign_halving_point << std::endl;
#endif
            foreign_halving_point->add_motorcycle(foreign_motorcycle_id, foreign_time_at_halving_point);

            // The target list of the foreign motorcycle was modified and the queue must be updated
            motorcycle_pq.update(foreign_mc);

            CGAL_postcondition(foreign_mc.has_target_at_time(collision_point, foreign_time_at_collision_point));
          }
          else
          {
            // this is a new point for the foreign motorcycle, but it belongs to
            // its confirmed track, and must therefore be blocked
            collision_point->block();

            // Add it to the track of the foreign motorcycle (useful to check
            // the validity of the final graph)
            foreign_mc.track().insert(std::make_pair(collision_point, foreign_time_at_collision_point));
          }
        }
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
        if(foreign_motorcycle_id != std::size_t(-1))
        {
          std::cout << "Post-treatment collision point:" << std::endl << *collision_point << std::endl;
          std::cout << "Motorcycles involved: " << std::endl << mc << std::endl
                    << motorcycle(foreign_motorcycle_id) << std::endl;
        }
#endif
      }

      // The target list of mc was modified and the PQ must be updated
      CGAL_assertion(!mc.is_crashed());
      motorcycle_pq.update(mc);
    }

    // block the newly reached position
    mc.position()->block();

next_item:
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "---" << std::endl;
    std::cout << "Queue after update:" << std::endl << motorcycle_pq << std::endl;
#endif
  }
}

template<typename MotorcycleGraphTraits>
bool
Motorcycle_graph<MotorcycleGraphTraits>::
is_valid() const
{
  // mega brute force validity check
  // @todo do something nice

  std::size_t number_of_motorcycles = motorcycles.size();
  for(std::size_t mc_id = 0; mc_id<number_of_motorcycles; ++mc_id)
  {
    const Motorcycle& mc = motorcycle(mc_id);
    const Track& mc_track = mc.track();
    CGAL_assertion(mc_track.size() > 0);
    if(mc_track.size() <= 2) // ignore degenerate tracks
      continue;

    typename Track::const_iterator lit = mc_track.begin(), end = mc_track.end();
    DEC_it current = lit->first, next;

    while(++lit != end)
    {
      next = lit->first;
      if(current->location().first != next->location().first)
      {
        std::cout << "Should be equal: " << current->point() << " and " << next->point() << std::endl;
        std::cout << "id: " << mc_id << std::endl;
        CGAL_assertion(CGAL::squared_distance(current->point(), next->point()) < std::numeric_limits<FT>::epsilon());
        current = next;
        continue;
      }

      if(current == next)
      {
        // current = next; // unneeded but left for clarity
        continue;
      }

      const Point_2 ts = gt.construct_point_2_object()(current->location().second[0],
                                                       current->location().second[1]);
      const Point_2 tt = gt.construct_point_2_object()(next->location().second[0],
                                                       next->location().second[1]);
      Segment_2 s = gt.construct_segment_2_object()(ts, tt);

      for(std::size_t fmc_id = 0; fmc_id<number_of_motorcycles; ++fmc_id)
      {
        if(fmc_id == mc_id)
          continue;

        const Motorcycle& fmc = motorcycle(fmc_id);
        const Track& fmc_track = fmc.track();
        CGAL_assertion(fmc_track.size() > 0);

        typename Track::const_iterator flit = fmc_track.begin(), fend = fmc_track.end();
        DEC_it fcurrent = flit->first, fnext = fcurrent;

        // degenerate fmc track
        if(fmc_track.size() == 1)
        {
          const Point_2 fts = gt.construct_point_2_object()(fcurrent->location().second[0],
                                                            fcurrent->location().second[1]);
          const Point_2 ftt = gt.construct_point_2_object()(fnext->location().second[0],
                                                            fnext->location().second[1]);
          Segment_2 fs = gt.construct_segment_2_object()(fts, ftt);

          if(gt.do_intersect_2_object()(s, fs))
          {
            std::cout << "Intersection ¤~~~~~~~~~~~~~~~~~¤ " << std::endl;
            std::cout << "motorcycle #" << mc_id << " (track size: " << mc_track.size();
            std::cout << ") with motorcycle #" << fmc_id << " (track size: " << fmc_track.size() << ")" << std::endl;
            std::cout << "cu/ne: " << std::endl << current->point() << " ## " << next->point() << std::endl;
            std::cout << "fcu/fne: " << std::endl << fcurrent->point() << " ## " << fnext->point() << std::endl;
            std::cout << "DECITs:" << std::endl << &*current << std::endl << &*next << std::endl << &*fcurrent << std::endl << &*fnext << std::endl;
            std::cout << "BCS points: " << std::endl << ts << std::endl << tt << std::endl << fts << std::endl << ftt << std::endl;

            // Xor
            CGAL_assertion((current == fcurrent && next != fcurrent) ||
                           (current != fcurrent && next == fcurrent));
          }
        }

        while(++flit != fend)
        {
          fnext = flit->first;

          // different face locations
          if(current->location().first != fcurrent->location().first)
          {
            // @todo handle intersections at an edge
            fcurrent = fnext;
            continue;
          }

          if(fcurrent->location().first != fnext->location().first)
          {
            std::cout << "Should be equal: " << fcurrent->point() << " and " << fnext->point() << std::endl;
            std::cout << "id: " << fmc_id << std::endl;
            CGAL_assertion(CGAL::squared_distance(fcurrent->point(), fnext->point()) < std::numeric_limits<FT>::epsilon());
            fcurrent = fnext;
            continue;
          }

          const Point_2 fts = gt.construct_point_2_object()(fcurrent->location().second[0],
                                                            fcurrent->location().second[1]);
          const Point_2 ftt = gt.construct_point_2_object()(fnext->location().second[0],
                                                            fnext->location().second[1]);
          Segment_2 fs = gt.construct_segment_2_object()(fts, ftt);

          if(gt.do_intersect_2_object()(s, fs))
          {
            std::cout << "Intersection ¤~~~~~~~~~~~~~~~~~¤ " << std::endl;
            std::cout << "motorcycle #" << mc_id << " (track size: " << mc_track.size();
            std::cout << ") with motorcycle #" << fmc_id << " (track size: " << fmc_track.size() << ")" << std::endl;
            std::cout << "DECITs:" << *current << std::endl << *next << std::endl << *fcurrent << std::endl << *fnext << std::endl;

            if(fcurrent == fnext) // degenerate fmc track
            {
              CGAL_assertion((current == fcurrent && next != fcurrent) ||
                             (current != fcurrent && next == fcurrent));
            }
            else
            {
            CGAL_assertion((current == fcurrent && current != fnext &&
                            next != fcurrent && next != fnext) ||
                           (current != fcurrent && current == fnext &&
                            next != fcurrent && next != fnext) ||
                           (current != fcurrent && current != fnext &&
                            next == fcurrent && next != fnext) ||
                           (current != fcurrent && current != fnext &&
                            next != fcurrent && next == fnext));
            }

            // Any intersection that is not at the source must crash the motorcycle
            // if the motorcycle reaches the collision point at a later time
            // than the other motorcycle
            typename Track::const_iterator flitb = flit;
            if((next == fcurrent && lit->second >= (--flitb)->second) ||
                (next == fnext && lit->second >= flit->second))
            {
              CGAL_assertion(lit == --(mc.track().end()));
            }
          }
          fcurrent = fnext;
        }
      }
      current = next;
    }
  }

  return true;
}

template<typename MotorcycleGraphTraits>
void
Motorcycle_graph<MotorcycleGraphTraits>::
output_all_dictionary_points() const
{
  typename Dictionary::Dictionary_entry_container::const_iterator dit = points.all_entries().begin();
  typename Dictionary::Dictionary_entry_container::const_iterator end = points.all_entries().end();

  std::stringstream oss;
  oss << "results_" << gt.dimension << "/dictionary_points.xyz" << std::ends;
  std::ofstream os(oss.str().c_str());
  oss.precision(17);

  for(; dit!=end; ++dit)
  {
    os << dit->point();
    if(gt.dimension == 2) // The '.xyz' format expects 3D points
      os << " 0";
    os << '\n';
  }
}

template<typename MotorcycleGraphTraits>
void
Motorcycle_graph<MotorcycleGraphTraits>::
output_motorcycles_sources_and_destinations() const
{
  std::stringstream oss_sour, oss_dest;
  oss_sour << "results_" << gt.dimension << "/motorcycles_sources.xyz" << std::ends;
  oss_dest << "results_" << gt.dimension << "/motorcycles_destinations.xyz" << std::ends;
  std::ofstream oss(oss_sour.str().c_str());
  std::ofstream osd(oss_dest.str().c_str());
  oss.precision(17);
  osd.precision(17);

  for(std::size_t i=0; i<motorcycles.size(); ++i)
  {
    oss << motorcycle(i).source()->point();
    osd << motorcycle(i).destination()->point();

    if(gt.dimension == 2) // The '.xyz' format expects 3D points
    {
      oss << " 0";
      osd << " 0";
    }

    oss << '\n';
    osd << '\n';
  }
}

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_MOTORCYCLE_GRAPH_H
