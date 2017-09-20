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
#include <CGAL/Polyline_tracing/Tracer.h>
#include <CGAL/Polyline_tracing/internal/robust_intersections.h>
#include <CGAL/Polyline_tracing/internal/VPM_selector.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/assertions.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/enum.h>
#include <CGAL/number_utils.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/result_of.h>

#include <boost/optional.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/unordered_map.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

namespace CGAL {

namespace Polyline_tracing {

namespace internal {

template<typename Motorcycle_graph>
struct Collision_information
{
  typedef typename Motorcycle_graph::FT                      FT;
  typedef typename Motorcycle_graph::Point                   Point;
  typedef typename Motorcycle_graph::DEC_it                  DEC_it;

  Collision_information()
    :
      closest_collision_point(),
      closest_collision(),
      is_closest_collision_point_already_in_dictionary(false),
      foreign_mc_id(-1),
      time_at_closest_collision(std::numeric_limits<FT>::max()),
      foreign_time_at_closest_collision(std::numeric_limits<FT>::max())
  { }

  Point closest_collision_point;
  DEC_it closest_collision;
  bool is_closest_collision_point_already_in_dictionary;
  std::size_t foreign_mc_id;
  FT time_at_closest_collision;
  FT foreign_time_at_closest_collision;
};

} // namespace internal

template<typename K, typename TriangleMesh>
class Motorcycle_graph
{
  typedef Motorcycle_graph<K, TriangleMesh>                    Self;

  typedef internal::Collision_information<Self>               Collision_information;

public:
  typedef typename K::FT                                      FT;
  typedef typename K::Point_2                                 Point;
  typedef typename K::Segment_2                               Segment;
  typedef typename K::Vector_2                                Vector;
  typedef typename K::Ray_2                                   Ray;

  typedef Dictionary<K, TriangleMesh>                         Dictionary;
  typedef Dictionary_entry<K, TriangleMesh>                   Dictionary_entry;
  typedef typename Dictionary::DEC_it                         DEC_it;
  typedef typename Dictionary_entry::Barycentric_coordinates  Barycentric_coordinates;
  typedef typename Dictionary_entry::Face_location            Face_location;

  typedef Motorcycle<K, TriangleMesh>                         Motorcycle;
  typedef std::vector<Motorcycle*>                            Motorcycle_container;

  typedef Motorcycle_priority_queue<K, TriangleMesh>           Motorcycle_PQ;
  typedef Motorcycle_priority_queue_entry<K, TriangleMesh>     Motorcycle_PQE;

  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor    halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor        face_descriptor;

  typedef boost::tuple<std::size_t, DEC_it, FT, DEC_it, FT>   Track;
  typedef std::list<Track>                                    Track_list;
  typedef boost::unordered_map<face_descriptor , Track_list>  Track_face_map;
  typedef typename Track_face_map::iterator                   TFM_iterator;

  // Access
  Motorcycle& motorcycle(const std::size_t id) { return *(motorcycles[id]); }
  const Motorcycle& motorcycle(const std::size_t id) const { return *(motorcycles[id]); }

  // Constructor
  Motorcycle_graph(TriangleMesh& mesh);

  // Functions

  template<typename MotorcycleContainerIterator>
  void add_motorcycles(MotorcycleContainerIterator mit, MotorcycleContainerIterator last);

  /// \param fd face in which the segment belongs
  /// \param id the id of the motorcycle
  /// \param s the source of the oriented segment
  /// \param t the target of the oriented segment
  ///
  /// \return iterator in the tracking map
  TFM_iterator add_track_to_map(face_descriptor fd, std::size_t id,
                                DEC_it s, const FT time_at_s,
                                DEC_it t, const FT time_at_t);

  /// \param p first point
  /// \param p_time time at first point
  /// \param q second point
  /// \param q_time time at second point
  ///
  /// \return new point and time at the new point
  std::pair<DEC_it, FT> compute_halving_point(const Motorcycle& mc, DEC_it p, const FT p_time, DEC_it q, const FT q_time);

  /// \param p first point
  /// \param p_time time at first point
  /// \param q second point
  /// \param q_time time at second point
  ///
  /// \return new point and time at the new point
  std::pair<DEC_it, FT> compute_middle_point(DEC_it p, const FT p_time, DEC_it q, const FT q_time);
  bool compute_motorcycle_next_path(Motorcycle& mc);
  void crash_motorcycle(Motorcycle& mc);
  void crash_motorcycles_with_same_source_and_direction();
  void drive_to_closest_target(Motorcycle& mc);

  void find_collision_with_track(Motorcycle& mc, const Segment& mcs,
                                 const Motorcycle& fmc, const Track& fmc_track,
                                 bool is_fmc_moving_on_track, Collision_information& tc);
  void find_collision_with_complete_track(Motorcycle& mc, const Segment& mcs,
                                          const Track& fmc_track, Collision_information& tc);
  void find_collision_with_live_motorcycle(Motorcycle& mc, const Segment& mcs,
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
  Dictionary points; // All the points that will be used throughout the algorithm
  Motorcycle_container motorcycles;
  Motorcycle_PQ motorcycle_pq; // motorcycle priority queue

  bool using_enclosing_bbox; // indicates whether a mesh is passed input
  TriangleMesh& mesh; // not 'const' in case we need to create it

  // map to keep in memory the completed tracks of the motorcycles for each face
  Track_face_map track_face_map;
};

// -----------------------------------------------------------------------------
template<typename K, typename TriangleMesh>
Motorcycle_graph<K, TriangleMesh>::
Motorcycle_graph(TriangleMesh& mesh)
  : points(),
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
  // The issue is that the points are based on a location described by barycentric
  // coordinates. I guess, I could generate a bbox, then a triangle that includes
  // the box ? Pretty ugly, though...
  CGAL_assertion(!using_enclosing_bbox);
}

template<typename K, typename TriangleMesh>
template<typename MotorcycleContainerIterator>
void
Motorcycle_graph<K, TriangleMesh>::
add_motorcycles(MotorcycleContainerIterator mit, MotorcycleContainerIterator last)
{
  motorcycles.reserve(std::distance(mit, last));

  std::size_t counter = 0; // unique motorcycle ids
  while(mit != last)
  {
    Motorcycle& mc = *mit++;
    mc.set_id(counter);

    boost::optional<Point>& destination_point = mc.initial_destination_point();
    boost::optional<Vector>& direction = mc.direction();

    CGAL_precondition_msg(destination_point != boost::none || direction != boost::none,
      "A motorcycle must have least a destination or a direction.");

    motorcycles.push_back(&mc);
    ++counter;
  }
}

template<typename K, typename TriangleMesh>
typename Motorcycle_graph<K, TriangleMesh>::TFM_iterator
Motorcycle_graph<K, TriangleMesh>::
add_track_to_map(face_descriptor fd, std::size_t id,
                 DEC_it s, const FT time_at_s,
                 DEC_it t, const FT time_at_t)
{
  Track tr = boost::make_tuple(id, s, time_at_s, t, time_at_t);
  Track_list l;
  l.push_back(tr);

  std::pair<typename Motorcycle_graph<K, TriangleMesh>::TFM_iterator, bool>
    is_insert_success = track_face_map.insert(std::make_pair(fd, l));

  if(!is_insert_success.second)
    is_insert_success.first->second.push_back(tr);

  return is_insert_success.first;
}

template<typename K, typename TriangleMesh>
std::pair<typename Motorcycle_graph<K, TriangleMesh>::DEC_it,
          typename Motorcycle_graph<K, TriangleMesh>::FT>
Motorcycle_graph<K, TriangleMesh>::
compute_halving_point(const Motorcycle& m, DEC_it p, const FT p_time,
                                           DEC_it q, const FT q_time)
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << " *** " << std::endl;
  std::cout << "Computing half on motorcycle track: " << m.id() << std::endl
            << "  " << *p << std::endl << "  " << *q << std::endl;
#endif
  // Assert that both points are in the same face
  CGAL_precondition(p->location().first == q->location().first);

#ifdef CGAL_MOTORCYCLE_GRAPH_USE_ADVANCED_HALVING
  // interface with the halving data structure @todo
#else
  return compute_middle_point(p, p_time, q, q_time);
#endif
}

template<typename K, typename TriangleMesh>
std::pair<typename Motorcycle_graph<K, TriangleMesh>::DEC_it,
          typename Motorcycle_graph<K, TriangleMesh>::FT>
Motorcycle_graph<K, TriangleMesh>::
compute_middle_point(DEC_it p, const FT p_time, DEC_it q, const FT q_time)
{
  if(p->location().first != q->location().first)
  {
    std::cerr << "Warning: middle point computation with different faces" << std::endl;
    // asserting because using p.loc().first is too dangerous if r is not guaranteed
    // to be on p's face
    CGAL_assertion(false);
  }

  Point r = K().construct_midpoint_2_object()(p->point(), q->point());
  const FT time_at_r = 0.5 * (p_time + q_time);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "  New middle point: (" << r  << ") at time: " << time_at_r << std::endl;
#endif

  // As long as the faces are equal, the barycenter coordinates, the coordinates
  // of the middle point can be deduced from the (known) coordinates of 'p' and 'q'
  // but is it better (accuracy + consistency) ? @todo
  const Face_location middle_point_location =
    CGAL::Polygon_mesh_processing::locate(p->location().first, r, mesh);

  std::pair<DEC_it, bool> entry = points.insert(middle_point_location, r);
  return std::make_pair(entry.first, time_at_r);
}

template<typename K, typename TriangleMesh>
bool
Motorcycle_graph<K, TriangleMesh>::
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

  mc.source() = next_source;
  mc.time_at_source() = mc.current_time();
  mc.destination() = next_destination;
  mc.position() = mc.source();

  next_source->add_motorcycle(mc.id(), mc.current_time());
  next_destination->add_motorcycle(mc.id(), time_at_next_destination);

  mc.add_target(next_source, mc.current_time());

  if(next_source != next_destination)
    mc.add_target(next_destination, time_at_next_destination);

  mc.set_destination_finality(is_destination_final);

  return true;
}

template<typename K, typename TriangleMesh>
void
Motorcycle_graph<K, TriangleMesh>::
crash_motorcycle(Motorcycle& mc)
{
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

template<typename K, typename TriangleMesh>
void
Motorcycle_graph<K, TriangleMesh>::
crash_motorcycles_with_same_source_and_direction()
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "Checking for motorcycles with same sources and directions" << std::endl;
#endif

  // brute force, for now
  // A smarter version is to sort motorcycles by direction (slope),
  // and check for consecutive entries @todo
  std::size_t number_of_motorcycles = motorcycles.size();
  for(std::size_t mc_id = 0; mc_id<number_of_motorcycles; ++mc_id)
  {
    Motorcycle& mc = motorcycle(mc_id);

    if(mc.source() == mc.destination() || mc.is_crashed())
      continue;

    for(std::size_t fmc_id = 0; fmc_id<number_of_motorcycles; ++fmc_id)
    {
      Motorcycle& fmc = motorcycle(fmc_id);

      if(fmc.id() == mc.id() ||
         fmc.source() == fmc.destination() || // a degenerate track does not block anything
         mc.source() != fmc.source()) // must have identical sources
        continue;

      // @todo Snap near degeneracies here ?
//      std::cout << "checking: " << mc_id << " and " << fmc_id << std::endl;
//      typename K::Vector_2 v1(mc.source()->point(), mc.destination()->point());
//      typename K::Vector_2 v2(fmc.source()->point(), fmc.destination()->point());
//      std::cout << *(mc.direction()) << std::endl;
//      std::cout << *(fmc.direction()) << std::endl;
//      std::cout << v1 / CGAL::sqrt(v1 * v1) << std::endl;
//      std::cout << v2 / CGAL::sqrt(v2 * v2) << std::endl;
//      std::cout << "scalar product: " << K().compute_scalar_product_2_object()(v1, v2) << std::endl;

      // only aligned tracks block one another
      if(!K().collinear_2_object()(mc.source()->point(), // == fmc.source()->point()
                                   mc.destination()->point(),
                                   fmc.destination()->point()))
        continue;

      std::cout << "collinear" << std::endl;

      // Moving away from each other from the same point is allowed.
      // use two calls to ordered_along_line() instead? @todo
      if(K().angle_2_object()(mc.source()->point(),
                              mc.destination()->point(),
                              fmc.source()->point(),
                              fmc.destination()->point()) == CGAL::ACUTE)
      {
        std::cout << "Crashing degenerate motorcycles: "
                  << mc.id() << " and " << fmc.id() << std::endl;
        crash_motorcycle(mc);
        crash_motorcycle(fmc);
        break;
      }
    }
  }
}

template<typename K, typename TriangleMesh>
void
Motorcycle_graph<K, TriangleMesh>::
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
  mc.track().push_back(closest_target);
  mc.erase_closest_target();

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "  now at: (" << mc.position()->point() << ")" << std::endl;
#endif
}

template<typename K, typename TriangleMesh>
void
Motorcycle_graph<K, TriangleMesh>::
generate_enclosing_face()
{
  // generate a bbox that includes all known positions and all crash points
  // 2D only for now @todo
  Bbox_2 bbox;

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
template<typename K, typename TriangleMesh>
bool
Motorcycle_graph<K, TriangleMesh>::
has_motorcycle_reached_final_destination(const Motorcycle& mc) const
{
  return mc.is_motorcycle_destination_final();
}

template<typename K, typename TriangleMesh>
void
Motorcycle_graph<K, TriangleMesh>::
find_collision_with_track(Motorcycle& mc, const Segment& mcs,
                          const Motorcycle& fmc, const Track& fmc_track,
                          bool is_fmc_moving_on_track,
                          // below are out parameters
                          Collision_information& tc)
{
  CGAL_precondition(mcs.source() != mcs.target());

  const DEC_it fmc_track_source = fmc_track.template get<1>();
  const FT time_at_fmc_track_source = fmc_track.template get<2>();
  const DEC_it fmc_track_destination = fmc_track.template get<3>();
  const FT time_at_fmc_track_destination = fmc_track.template get<4>();

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "Checking collision with track of motorcycle: " << fmc.id() << std::endl;
  std::cout << " source: " << fmc_track_source->point() << " at t: " << time_at_fmc_track_source << std::endl;
  std::cout << " target: " << fmc_track_destination->point() << " at t: " << time_at_fmc_track_destination << std::endl;
#endif

  FT time_at_collision = 0.;

  const Segment fmcs = K().construct_segment_2_object()(
                         fmc_track_source->point(), fmc_track_destination->point());
  bool is_fmcs_degenerate = fmcs.is_degenerate();

  // Ignore the degenerate case of a degenerate fmc track starting at the same source
  if(is_fmcs_degenerate && mcs.source() == fmcs.source())
  {
    std::cout << "degenerate fmc and mcs.source() == fmcs.source()" << std::endl;
    return;
  }

  // Detect whether the motorcycles share the same supporting line.
  // Note that we know that 'mcs' is not degenerate.
  if(K().collinear_2_object()(mcs.source(), mcs.target(), fmcs.source()) &&
     K().collinear_2_object()(mcs.source(), mcs.target(), fmcs.target()))
  {
    std::cerr << "  /!\\ Tracks are aligned" << std::endl;

    // Check if the tentative tracks actually intersect
    if(!K().do_intersect_2_object()(mcs, fmcs))
    {
      // No intersection, try the next motorcycle
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
    Ray mcr(mcs.source(), mcs.target());
    if(!mcr.has_on(fmcs.source()))
      return;

    // Compute the respective direction of the two motorcycles:
    // use two calls to collinear_are_aligned_along_line() instead of angle ? @todo
    bool are_motorcycles_moving_in_the_same_direction =
      (is_fmcs_degenerate ||
       K().angle_2_object()(mcs.source(), mcs.target(), fmcs.source(), fmcs.target()) == CGAL::ACUTE);

    std::cout << "  is degen: " << is_fmcs_degenerate << std::endl;
    std::cout << "  angle: " << K().angle_2_object()(mcs.source(), mcs.target(),
                                                     fmcs.source(), fmcs.target()) << std::endl;
    std::cout << "  are motorcycles moving in the same direction: "
              << are_motorcycles_moving_in_the_same_direction << std::endl;

    // The motorcycles move in the same direction ==> mc will impact fmcs' source
    if(are_motorcycles_moving_in_the_same_direction)
    {
      // The weird configuration of both motorcycles moving in the same direction
      // AND with the same source is checked at the very start, see function
      // 'crash_motorcycles_with_same_source_and_direction()'
      CGAL_assertion(is_fmcs_degenerate || mcs.source() != fmcs.source());

      time_at_collision = mc.current_time() +
        CGAL::sqrt(CGAL::squared_distance(mcs.source(), fmcs.source())) / mc.speed();
      std::cout << "  mc crashes into fmcs' source at time: " << time_at_collision << std::endl;

      // Compare with other possible collisions to keep the closest
      if(time_at_collision < tc.time_at_closest_collision)
      {
        tc.is_closest_collision_point_already_in_dictionary = true;
        tc.time_at_closest_collision = time_at_collision;
        tc.closest_collision = fmc_track_source;
        tc.foreign_mc_id = fmc.id();
        tc.foreign_time_at_closest_collision = time_at_fmc_track_source;
      }
    }
    else // Motorcycles are moving towards each other
    {
      // Ignore the case where the motorcycles are moving away from each other
      // from the same point (not problematic)
      if(mcs.source() == fmcs.source())
        return;

      // If the foreign motorcycle is (also) moving, we return the middle point
      // that both motorcycles would reach at the same time.
      // Note that this point might not actually be reached by either motorcycle,
      // e.g. if they reach their destination first.
      if(is_fmc_moving_on_track)
      {
        // @fixme, if speeds are ever allowed to change, the speed of fmc here
        // must be changed to the speed on that segment
        time_at_collision = mc.current_time() +
          (CGAL::sqrt(CGAL::squared_distance(mcs.source(), fmcs.source())) -
             fmc.speed() * (mc.current_time() - time_at_fmc_track_source)) /
               (mc.speed() + fmc.speed());
        std::cout << "  mc and fmc would meet at time: " << time_at_collision << std::endl;
        CGAL_postcondition(time_at_collision > mc.current_time());

        // Compare with other possible collisions to keep the closest
        if(time_at_collision < tc.time_at_closest_collision)
        {
          tc.is_closest_collision_point_already_in_dictionary = false;
          tc.time_at_closest_collision = time_at_collision;

          Vector mcv(mcs);
          tc.closest_collision_point = mcs.source() +
            time_at_collision * mc.speed() * mcv / CGAL::sqrt(mcv.squared_length());

          tc.foreign_mc_id = fmc.id();
          tc.foreign_time_at_closest_collision = time_at_collision;
        }
      }
      // If fmc is not moving, then mc crashes into the final position of the track
      else
      {
        time_at_collision = mc.current_time() +
                              CGAL::sqrt(CGAL::squared_distance(
                                mcs.source(), fmcs.target())) / mc.speed();
        std::cout << "  mc crashes into fmc's final position at: " << time_at_collision << std::endl;

        // Compare with other possible collisions to keep the closest
        if(time_at_collision < tc.time_at_closest_collision)
        {
          tc.is_closest_collision_point_already_in_dictionary = true;
          tc.time_at_closest_collision = time_at_collision;
          tc.closest_collision = fmc_track_destination;
          tc.foreign_mc_id = fmc.id();
          tc.foreign_time_at_closest_collision = time_at_fmc_track_destination;
        }
      }
    }
  }
  // --- From here on, the tracks are not collinear ---
  else
  {
    // The next two "if" are two checks for robustness : a non-degenerate
    // intersection of the two tracks is at most one point. To avoid numerical errors,
    // the function 'motorcycle_visiting_time' is used to detect whether an extremity
    // of a track is already that known intersection.
    std::pair<bool, FT> reaching_time = mc.position()->motorcycle_visiting_time(fmc.id());

    // Check if the current position of mc is a known intersection with the foreign track
    if(// the foreign motorcycle passes through the current position of mc
        reaching_time.first &&
        // Check that this point is indeed on the foreign track
        reaching_time.second >= time_at_fmc_track_source &&
        reaching_time.second <= time_at_fmc_track_destination)
    {
      // Ignore this intersection: since we are in this function, it means that
      // the position was not blocked and we can move on.
      return;
    }

    // Check if the closest target of mc is a known intersection with the foreign track
    reaching_time = mc.closest_target()->motorcycle_visiting_time(fmc.id());
    if(// the foreign motorcycle passes through the closest_target of mc
       reaching_time.first &&
       // Check that this point is indeed on the foreign track
       reaching_time.second >= time_at_fmc_track_source &&
       reaching_time.second <= time_at_fmc_track_destination)
    {
      time_at_collision = mc.time_at_closest_target();
      std::cout << "  /!\\ tentative path collides with track: " << fmc.id()
                << " at its end. Time: " << time_at_collision << std::endl;

      // Compare with other possible collisions to keep the closest
      if(time_at_collision < tc.time_at_closest_collision)
      {
        tc.is_closest_collision_point_already_in_dictionary = true;
        tc.closest_collision = mc.closest_target();
        tc.time_at_closest_collision = time_at_collision;
        tc.foreign_mc_id = fmc.id();
        tc.foreign_time_at_closest_collision =
          tc.closest_collision->find_motorcycle(fmc.id())->second;
      }
      return;
    }

    // --- Nothing is easy, the intersection must be computed ---
    if(!K().do_intersect_2_object()(mcs, fmcs))
    {
      // No intersection, move to the next motorcycle
      std::cout << "  No intersection" << std::endl;
      return;
    }

    std::cout << "  /!\\ collision with motorcycle: " << fmc.id()
              << " between tracks ("
              << mcs.source() << ")--(" << mcs.target() << ") and ("
              << fmcs.source() << ")--(" << fmcs.target() << ")" << std::endl;

    const Point& collision_point = internal::robust_intersection<K>(mcs, fmcs);

    time_at_collision = mc.current_time() +
      CGAL::sqrt(CGAL::squared_distance(mcs.source(), collision_point)) / mc.speed();
    std::cout << "  collision at: (" << collision_point
              << ") at time: " << time_at_collision << std::endl;

    // if the collision time is greater than the time at destination, then
    // the collision is not interesting
    CGAL_assertion(mc.targets().rbegin()->first == mc.destination());
    if(time_at_collision > mc.targets().rbegin()->second)
      return;

    // Compare with other collisions to keep the closest
    if(time_at_collision < tc.time_at_closest_collision)
    {
      tc.is_closest_collision_point_already_in_dictionary = false;
      tc.time_at_closest_collision = time_at_collision;
      tc.foreign_mc_id = fmc.id();
      tc.closest_collision_point = collision_point;
      tc.foreign_time_at_closest_collision = time_at_fmc_track_source +
        CGAL::sqrt(CGAL::squared_distance(fmcs.source(), collision_point)) / fmc.speed();
    }
  }
}

template<typename K, typename TriangleMesh>
void
Motorcycle_graph<K, TriangleMesh>::
find_collision_with_complete_track(Motorcycle& mc, const Segment& mcs,
                                   const Track& fmc_track,
                                   // below are out parameters
                                   Collision_information& tc)
{
  const std::size_t fmc_id = fmc_track.template get<0>();
  const Motorcycle& fmc = motorcycle(fmc_id);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
std::cout << "Checking for intersection with the complete track of motorcycle " << fmc.id() << std::endl;
#endif

  return find_collision_with_track(mc, mcs, fmc, fmc_track,
                                   false /*the motorcycle is not moving on that track*/,
                                   tc);
}

template<typename K, typename TriangleMesh>
void
Motorcycle_graph<K, TriangleMesh>::
find_collision_with_live_motorcycle(Motorcycle& mc, const Segment& mcs,
                                    const Motorcycle& fmc,
                                    // below are out parameters
                                    Collision_information& tc)
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
std::cout << "Checking for intersection with the live motorcycle " << fmc.id() << std::endl;
#endif

  if(mc.id() == fmc.id() ||
     mc.current_location().first != fmc.current_location().first || // the motorcycles must be in the same face
     fmc.is_crashed()) // the track of a crashed motorcycle track is complete
  {
    std::cout << " ignoring... " << std::endl;
    std::cout << "  ids: " << mc.id() << " " << fmc.id() << std::endl;
    std::cout << "  faces: " << mc.current_location().first << " and " << fmc.current_location().first << std::endl;
    std::cout << "  crashed status: " << fmc.is_crashed() << std::endl;
    return;
  }

  Track fmc_track = boost::make_tuple(fmc.id(),
                                      fmc.source(), fmc.time_at_source(),
                                      fmc.closest_target(),
                                      fmc.time_at_closest_target());

  return find_collision_with_track(mc, mcs, fmc, fmc_track,
                                   true /*fmc is currently moving on that track*/,
                                   tc);
}

// search for a possible collision with another motorcycle between the current
// position of mc and the next target
template<typename K, typename TriangleMesh>
boost::tuple<typename Motorcycle_graph<K, TriangleMesh>::DEC_it,
             typename Motorcycle_graph<K, TriangleMesh>::FT,
             std::size_t,
             typename Motorcycle_graph<K, TriangleMesh>::FT>
Motorcycle_graph<K, TriangleMesh>::
find_collision(Motorcycle& mc)
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "Checking for collisions on motorcycle " << mc.id() << "'s track" << std::endl
            << "Currently on face: " << mc.current_location().first << std::endl;
#endif

  CGAL_precondition(!mc.is_crashed());
  CGAL_precondition(!mc.targets().empty());

  // A bunch of output parameters are regrouped into the 'Collision_information' struct,
  // which describes the best (closest to mc.position()) tentative collision.
  Collision_information tc;

  Segment mc_tentative_track =
    K().construct_segment_2_object()(mc.position()->point(), mc.closest_target()->point());

  std::cout << "MC tentative track: " << std::endl << *(mc.position()) << std::endl
                                                   << *(mc.closest_target()) << std::endl;

  // A degenerate tentative track has no interesting collisions
  if(mc_tentative_track.is_degenerate())
    return boost::make_tuple(tc.closest_collision, tc.time_at_closest_collision,
                             tc.foreign_mc_id, tc.foreign_time_at_closest_collision);

  // -----------------------------------------------------------------------
  // Checking for intersection is done in two steps:
  // - Check with complete tracks in the face
  // - Check the motorcycles that are currently traveling in the face

  // Step 1: check complete tracks
  const face_descriptor mc_face = mc.current_location().first;
  TFM_iterator it = track_face_map.find(mc_face);
  if(it != track_face_map.end())
  {
    const Track_list& face_tracks = it->second;

    typename Track_list::const_iterator tl_it = face_tracks.begin();
    typename Track_list::const_iterator tl_end = face_tracks.end();
    for(; tl_it!=tl_end; ++tl_it)
    {
      const Track& track = *tl_it;
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
     !tc.is_closest_collision_point_already_in_dictionary)
  {
    // Insert the collision point in the dictionary. Motorcycle info will be added later.
    const Face_location collision_location =
      CGAL::Polygon_mesh_processing::locate(
        mc.current_location().first, tc.closest_collision_point, mesh);

    std::pair<DEC_it, bool> entry = points.insert(collision_location,
                                                  tc.closest_collision_point);
    tc.closest_collision = entry.first;
  }

  return boost::make_tuple(tc.closest_collision, tc.time_at_closest_collision,
                           tc.foreign_mc_id, tc.foreign_time_at_closest_collision);
}

template<typename K, typename TriangleMesh>
void
Motorcycle_graph<K, TriangleMesh>::
initialize_motorcycles()
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "Initialize motorcycles" << std::endl;
#endif

  // if no mesh has been given in input, generate a mesh made of a single quad face
  // that contains all the interesting motorcycle interactions (crashes)
  if(using_enclosing_bbox)
    generate_enclosing_face();

  typedef CGAL::internal::P2_or_P3_to_P3<TriangleMesh>                 P2_or_P3_to_P3;
  typedef CGAL::P2_to_P3_VPM<TriangleMesh>                             VPM;
  VPM vpm(mesh);

  typedef CGAL::AABB_face_graph_triangle_primitive<TriangleMesh, VPM>  AABB_face_graph_primitive;
  typedef CGAL::AABB_traits<K, AABB_face_graph_primitive>             AABB_face_graph_traits;
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
    std::cout << "Initializing motorcycle: " << mc_id << std::endl;
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
      CGAL_assertion(res.template get<0>()); // a destination must be found

      // the location algorithm might change source to make sure source and destination
      // are on the same face
      if(mc.source() != res.template get<1>())
      {
        std::cerr << "Warning: source has changed!" << std::endl
                  << "Previously: " << *(mc.source()) << std::endl
                  << "Now: " << *(res.template get<1>()) << std::endl;

        // The source change must only be a change of face descriptor, not of position
        CGAL_assertion(mc.source()->point() == res.template get<1>()->point());

        // if the source point wasn't previously in the dictionary, clean it off
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

      destination->add_motorcycle(mc_id, time_at_destination);
      opt_destination_point = destination->point();
    }
    else // destination is known, only need to compute the time of arrival
    {
      const Point& destination_point = *opt_destination_point;

      // source and destination should be on the same face, so we use the overload
      // where a face is provided @fixme
      const Face_location destination_location =
        CGAL::Polygon_mesh_processing::locate(source_location.first,
                                              destination_point, mesh);

      // check that the destination is indeed in the same face as the source
      CGAL_precondition(destination_location.second[0] >= 0. &&
                        destination_location.second[0] <= 1. &&
                        destination_location.second[1] >= 0. &&
                        destination_location.second[1] <= 1. &&
                        destination_location.second[2] >= 0. &&
                        destination_location.second[2] <= 1.);

      if(mc.id()%3 == 0)// @tmp, for fun, don't always stop at the destination
        mc.set_destination_finality(true);

      // @todo this should be computed by the tracer (?)
      time_at_destination = time_at_source +
        CGAL::sqrt(CGAL::squared_distance(ini_source_point, destination_point)) / speed;

      std::pair<DEC_it, bool> destination_entry = points.insert(destination_location,
                                                                destination_point,
                                                                mc_id, time_at_destination);
      destination = destination_entry.first;

      if(direction == boost::none)
      {
        mc.direction() = Vector(ini_source_point, destination_point);
      }
      else // both the destination and the direction are known
      {
        // sanity check: (destination - source) is collinear with the direction
        Ray r(ini_source_point, *(mc.direction()));
        if(!r.has_on(destination_point))
        {
          std::cerr << "Error: Incompatible destination and direction: " << std::endl
                    << "- destination: " << destination_point << std::endl
                    << "- direction: " << *(mc.direction()) << std::endl;
          CGAL_assertion(false);
        }
      }
    }

    mc.destination() = destination;

    // Initialize the motorcycle target queue
    mc.targets().insert(std::make_pair(mc.source(), time_at_source));

    if(mc.source() != mc.destination())
      mc.targets().insert(std::make_pair(destination, time_at_destination));

    // this is useful to not get an empty track when sour=dest
    // but it creates duplicates @fixme
    mc.track().push_back(mc.source());
  }
}

template<typename K, typename TriangleMesh>
bool
Motorcycle_graph<K, TriangleMesh>::
has_motorcycle_reached_crashing_point(const Motorcycle& mc) const
{
  return  // mc has reached the track of a foreign motorcycle
         (mc.has_reached_blocked_point() ||
          // multiple motorcycles will reach mc's current position at the same time
          mc.has_reached_simultaneous_collision_point());
}

template<typename K, typename TriangleMesh>
template<typename MotorcycleContainerIterator>
void
Motorcycle_graph<K, TriangleMesh>::
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
      add_track_to_map(mc.current_location().first, mc.id(),
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
        std::cout << "Computing motorcycle's next path" << std::endl;
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
          // couldn't find a next path, crash the motorcycle
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
      add_track_to_map(mc.current_location().first, mc.id(),
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
      std::cout << "Find collision results: " << std::endl;
      if(foreign_motorcycle_id != std::size_t(-1))
      {
        std::cout << " - collision_point: " << *(collision_point) << std::endl
                  << " - foreign_motorcycle_id: " << foreign_motorcycle_id << std::endl
                  << " - time_at_collision_point: " << time_at_collision_point << std::endl
                  << " - foreign_time_at_collision_point: " << foreign_time_at_collision_point << std::endl;
      }
      else
      {
        std::cout << " No collision! " << std::endl;
      }
#endif

      if(// there is intersection
         foreign_motorcycle_id != std::size_t(-1) &&
         // the impact is closer than the next target
         time_at_collision_point <= mc.time_at_closest_target() &&
         // the collision is not the next target of 'mc' or the foreign track
         // does not know this collision point yet
         (collision_point != mc.closest_target() ||
          !collision_point->has_motorcycle(foreign_motorcycle_id)))
      {
        if(!collision_point->has_motorcycle(mc.id()))
        {
          // Call the halving structure to create a new point
          std::pair<DEC_it, FT> halving_entity =
            compute_halving_point(mc,
                                  mc.position(), mc.current_time(),
                                  collision_point, time_at_collision_point);
          DEC_it halving_point = halving_entity.first;
          const FT time_at_halving_point = halving_entity.second;

          // Degeneracies should be caught before
          CGAL_postcondition(halving_point != mc.position() &&
                             halving_point != collision_point);

          mc.add_target(collision_point, time_at_collision_point);
          mc.add_target(halving_point, time_at_halving_point);

          halving_point->add_motorcycle(mc.id(), time_at_halving_point);
          collision_point->add_motorcycle(mc.id(), time_at_collision_point);
        }

        if(!collision_point->has_motorcycle(foreign_motorcycle_id))
        {
          // it is useful to know that the collision point is on the foreign track,
          // even if the collision point is on the confirmed part of the track
          collision_point->add_motorcycle(foreign_motorcycle_id, foreign_time_at_collision_point);

          Motorcycle& foreign_mc = motorcycle(foreign_motorcycle_id);
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

            // Degeneracies should be caught before
            CGAL_postcondition(foreign_halving_point != foreign_mc.position() &&
                               foreign_halving_point != collision_point);

            foreign_mc.add_target(collision_point, foreign_time_at_collision_point);
            foreign_mc.add_target(foreign_halving_point, foreign_time_at_halving_point);

            foreign_halving_point->add_motorcycle(foreign_motorcycle_id, foreign_time_at_halving_point);

            // The target list of the foreign motorcycle was modified and the queue must be updated
            motorcycle_pq.update(foreign_mc);
          }
          else
          {
            // this is a new point for the foreign motorcycle, but it belongs to
            // its confirmed track, and must therefore be blocked
            collision_point->block();
          }
        }

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
        std::cout << "Post-treatment collision point: "  << *collision_point << std::endl;
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

template<typename K, typename TriangleMesh>
bool
Motorcycle_graph<K, TriangleMesh>::
is_valid() const
{
  // @todo
  // check that every point on the track is only:
  // - one of extremities of the track and in that case it can be anywhere
  //   on another path
  // - an extremity of another track

  return true;
}

template<typename K, typename TriangleMesh>
void
Motorcycle_graph<K, TriangleMesh>::
output_all_dictionary_points() const
{
  typename Dictionary::Dictionary_entry_container::const_iterator dit = points.all_entries().begin();
  typename Dictionary::Dictionary_entry_container::const_iterator end = points.all_entries().end();

  std::ofstream os("dictionary_points.xyz");
  for(; dit!=end; ++dit)
    os << dit->point() << " 0" << '\n';
}

template<typename K, typename TriangleMesh>
void
Motorcycle_graph<K, TriangleMesh>::
output_motorcycles_sources_and_destinations() const
{
  // must be adapted to surfaces @todo

  std::ofstream oss("out_motorcycles_sources.xyz");
  std::ofstream osd("out_motorcycles_destinations.xyz");
  for(std::size_t i=0; i<motorcycles.size(); ++i)
  {
    oss << motorcycle(i).source()->point() << " 0" << '\n';
    osd << motorcycle(i).destination()->point() << " 0" << '\n';
  }
}

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_MOTORCYCLE_GRAPH_H
