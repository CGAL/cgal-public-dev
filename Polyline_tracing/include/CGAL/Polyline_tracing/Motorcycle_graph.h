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

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Polyline_tracing/Dictionary.h>
#include <CGAL/Polyline_tracing/Motorcycle.h>
#include <CGAL/Polyline_tracing/Motorcycle_priority_queue.h>
#include <CGAL/Polyline_tracing/Tracer.h>

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/enum.h>
#include <CGAL/intersection_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/result_of.h>

#include <boost/optional.hpp>
#include <boost/tuple/tuple.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

namespace CGAL {

namespace Polyline_tracing {

template<typename K, typename PolygonMesh>
class Motorcycle_graph
{
public:
  typedef typename K::FT                                      FT;
  typedef typename K::Point_2                                 Point;
  typedef typename K::Segment_2                               Segment;
  typedef typename K::Vector_2                                Vector;
  typedef typename K::Ray_2                                   Ray;

  typedef Dictionary<K, PolygonMesh>                          Dictionary;
  typedef Dictionary_entry<K, PolygonMesh>                    Dictionary_entry;
  typedef typename Dictionary::DEC_it                         DEC_it;
  typedef typename Dictionary_entry::Barycentric_coordinates  Barycentric_coordinates;
  typedef typename Dictionary_entry::Face_location            Face_location;

  typedef Motorcycle<K, PolygonMesh>                          Motorcycle;
  typedef std::vector<Motorcycle*>                            Motorcycle_container;

  typedef Motorcycle_priority_queue<K, PolygonMesh>           Motorcycle_PQ;
  typedef Motorcycle_priority_queue_entry<K, PolygonMesh>     Motorcycle_PQE;

  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;



  // Access
  Motorcycle& motorcycle(const std::size_t id) { return *(motorcycles[id]); }
  const Motorcycle& motorcycle(const std::size_t id) const { return *(motorcycles[id]); }

  // Constructor
  Motorcycle_graph(PolygonMesh& mesh);

  // Functions
  template<typename MotorcycleContainerIterator>
  void add_motorcycles(MotorcycleContainerIterator mit, MotorcycleContainerIterator last);
  std::pair<DEC_it, FT> compute_halving_point(const Motorcycle& mc, DEC_it p, const FT p_time, DEC_it q, const FT q_time);
  std::pair<DEC_it, FT> compute_middle_point(DEC_it p, const FT p_time, DEC_it q, const FT q_time);
  std::pair<DEC_it, FT> compute_motorcycle_next_destination(const Motorcycle& mc);
  void crash_motorcycle(Motorcycle& mc);
  void crash_motorcycles_with_same_source_and_direction();
  void drive_to_closest_target(Motorcycle& mc);
  boost::tuple<DEC_it, FT, std::size_t, FT> find_collisions(Motorcycle& mc);
  void generate_enclosing_face();
  bool has_motorcycle_reached_final_destination(const Motorcycle& mc) const;
  void initialize_motorcycles();
  void locate_motorcycle(Motorcycle& mc);
  Point robust_intersection(const Segment& s, const Segment& t) const; // @todo put that in traits
  Point robust_intersection(const Ray& s, const Segment& t) const; // @todo put that in traits

  template<typename MotorcycleContainerIterator>
  void trace_motorcycle_graph(MotorcycleContainerIterator mit, MotorcycleContainerIterator last);

  // Post-tracing checks
  bool is_valid() const;

  // Output
  void output_all_dictionary_points() const;
  void output_motorcycles_sources_and_destinations() const;

private:
  Dictionary points;
  Motorcycle_container motorcycles;
  Motorcycle_PQ motorcycle_pq;
  PolygonMesh& mesh; // not const in case we need to create it
  bool using_enclosing_bbox;
};

// -----------------------------------------------------------------------------
template<typename K, typename PolygonMesh>
Motorcycle_graph<K, PolygonMesh>::
Motorcycle_graph(PolygonMesh& mesh)
  : points(),
    motorcycles(),
    motorcycle_pq(),
    mesh(mesh),
    using_enclosing_bbox(false)
{
  if(num_vertices(mesh))
  {
    std::cerr << " Warning: empty mesh in input" << std::endl;
    using_enclosing_bbox = true;
  }
}

template<typename K, typename PolygonMesh>
template<typename MotorcycleContainerIterator>
void
Motorcycle_graph<K, PolygonMesh>::
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

template<typename K, typename PolygonMesh>
std::pair<typename Motorcycle_graph<K, PolygonMesh>::DEC_it,
          typename Motorcycle_graph<K, PolygonMesh>::FT>
Motorcycle_graph<K, PolygonMesh>::
compute_middle_point(DEC_it p, const FT p_time, DEC_it q, const FT q_time)
{
  Point r = K().construct_midpoint_2_object()(p->point(), q->point());
  const FT time_at_r = 0.5 * (p_time + q_time);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "  New middle point: (" << r  << ") at time: " << time_at_r << std::endl;
#endif

  return std::make_pair(points.insert(r), time_at_r);
}

template<typename K, typename PolygonMesh>
std::pair<typename Motorcycle_graph<K, PolygonMesh>::DEC_it,
          typename Motorcycle_graph<K, PolygonMesh>::FT>
Motorcycle_graph<K, PolygonMesh>::
compute_halving_point(const Motorcycle& m, DEC_it p, const FT p_time,
                                           DEC_it q, const FT q_time)
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << " *** " << std::endl;
  std::cout << "Computing half on motorcycle track: " << m.id() << std::endl
            << "  " << *p << std::endl << "  " << *q << std::endl;
#endif
  // assert that p & q are on the same face @todo

#ifdef CGAL_MOTORCYCLE_GRAPH_USE_ADVANCED_HALVING
  // interface with the tracer data structure @todo
#else
  return compute_middle_point(p, p_time, q, q_time);
#endif
}

template<typename K, typename PolygonMesh>
std::pair<typename Motorcycle_graph<K, PolygonMesh>::DEC_it,
          typename Motorcycle_graph<K, PolygonMesh>::FT>
Motorcycle_graph<K, PolygonMesh>::
compute_motorcycle_next_destination(const Motorcycle& mc)
{
  // interface with the tracer data structure @todo

  CGAL_precondition(mc.direction());
  Vector mc_dir = *(mc.direction());

  // degenerate case
  if(mc_dir == CGAL::NULL_VECTOR)
  {
    std::cerr << "Warning: direction is null, thus the next destination is the current position" << std::endl;
    return std::make_pair(mc.position(), mc.current_time());
  }

  CGAL_precondition(num_vertices(mesh) != 0);

  typedef CGAL::Halfedge_around_face_circulator<PolygonMesh>  Halfedge_around_facet_circulator;

  Point farthest_destination;
  FT time_at_farthest_destination = std::numeric_limits<FT>::min();

  Ray r(mc.position()->point(), mc_dir);

  // @todo should be a function that computes the face
  face_descriptor fd = mc.current_location().first;

  Halfedge_around_facet_circulator hcir_begin(halfedge(fd, mesh), mesh);
  Halfedge_around_facet_circulator hcir = hcir_begin;

  do
  {
    halfedge_descriptor hd = *hcir++;

    Segment s(mesh.point(source(hd, mesh)), mesh.point(target(hd, mesh))); // @BGL
    std::cout << "ray; segment r:" << r << " s: " << s << std::endl;

    if(K().do_intersect_2_object()(r, s))
    {
      const Point new_destination = robust_intersection(r, s);
      std::cout << "new potential destination: " << new_destination << std::endl;

      // compute time at destination
      FT time_at_new_destination = mc.current_time() + // @todo tracer
        CGAL::sqrt(CGAL::squared_distance(mc.position()->point(), new_destination)) / mc.speed();

      if(time_at_new_destination > time_at_farthest_destination)
      {
        farthest_destination = new_destination;
        time_at_farthest_destination = time_at_new_destination;
      }
    }
  } while (hcir != hcir_begin);

  // no intersection with the border... Is the point not in the face?
  CGAL_assertion(time_at_farthest_destination != std::numeric_limits<FT>::min());

  // @todo handle the case where the new point is the current_point
  // (I guess the motorcycle crashes in that case...)

  std::cout << "new destination at : " << farthest_destination
            << " time: " << time_at_farthest_destination << std::endl;
  DEC_it new_destination_it = points.insert(farthest_destination);
  return std::make_pair(new_destination_it, time_at_farthest_destination);
}

template<typename K, typename PolygonMesh>
void
Motorcycle_graph<K, PolygonMesh>::
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

template<typename K, typename PolygonMesh>
void
Motorcycle_graph<K, PolygonMesh>::
crash_motorcycles_with_same_source_and_direction()
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "Detect degenerate motorcycle configurations" << std::endl;
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

    for(std::size_t mc_2_id = 0; mc_2_id<number_of_motorcycles; ++mc_2_id)
    {
      Motorcycle& mc_2 = motorcycle(mc_2_id);

      if(mc_2.id() == mc.id() ||
         mc_2.source() == mc_2.destination() || // degenerate track does not block anything
         mc.source() != mc_2.source()) // must have identical sources
        continue;

      // only aligned tracks block one another
      if(!K().collinear_2_object()(mc.source()->point(), // == mc_2.source()->point()
                                   mc.destination()->point(),
                                   mc_2.destination()->point()))
        continue;

      // Moving away from each other from the same point is allowed.
      // use ordered_along_line? @todo
      if(K().angle_2_object()(mc.source()->point(),
                              mc.destination()->point(),
                              mc_2.source()->point(),
                              mc_2.destination()->point()) == CGAL::ACUTE)
      {
        std::cout << "Crashing degenerate motorcycles: "
                  << mc.id() << " and " << mc_2.id() << std::endl;
        crash_motorcycle(mc);
        crash_motorcycle(mc_2);
        break;
      }
    }
  }
}

template<typename K, typename PolygonMesh>
void
Motorcycle_graph<K, PolygonMesh>::
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

template<typename K, typename PolygonMesh>
void
Motorcycle_graph<K, PolygonMesh>::
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
    for(std::size_t mc_2_id = 0; mc_2_id<number_of_motorcycles; ++mc_2_id)
    {
      // segment - segment, segment - ray, or ray-ray intersections @todo
    }
  }

  // Slightly increase the size of the bbox to strictly contain all points

  // Manually create the mesh with Euler operations

}
template<typename K, typename PolygonMesh>
bool
Motorcycle_graph<K, PolygonMesh>::
has_motorcycle_reached_final_destination(const Motorcycle& mc) const
{
  if(mc.is_motorcycle_destination_final())
    return true;


  // determine if we've reached the border of the mesh
// if(is_border(mc.current_location(), mesh))
//   return true;

  return false;
}

// search for a possible collision with another motorcycle between the current
// position of mc and the next target
template<typename K, typename PolygonMesh>
boost::tuple<typename Motorcycle_graph<K, PolygonMesh>::DEC_it,
             typename Motorcycle_graph<K, PolygonMesh>::FT, std::size_t,
             typename Motorcycle_graph<K, PolygonMesh>::FT>
Motorcycle_graph<K, PolygonMesh>::
find_collisions(Motorcycle& mc)
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "Checking for collisions on motorcycle " << mc.id() << "'s track" << std::endl;
#endif

  CGAL_precondition(!mc.is_crashed());
  CGAL_precondition(!mc.targets().empty());

  // interface with the ray shooting data structure @todo
  // THIS WHOLE FUNCTION IS @TMP

  // brute force for now:
  Point closest_collision;
  std::size_t foreign_mc_id(-1);
  FT time_at_closest_collision = std::numeric_limits<FT>::max();
  FT foreign_time_at_closest_collision = std::numeric_limits<FT>::max();
  DEC_it closest_collision_point;
  bool is_closest_collision_point_already_in_dictionary = false;

  // Only consider the part of mc's tentative track that is unconfirmed
  Segment mc_tentative_track =
    K().construct_segment_2_object()(mc.position()->point(),
                                     mc.closest_target()->point());

  std::cout << "MC tentative track: " << std::endl
            << *(mc.position()) << std::endl
            << *(mc.closest_target()) << std::endl;

  // A degenerate tentative track has no interesting collisions
  if(mc_tentative_track.is_degenerate())
    return boost::make_tuple(closest_collision_point, time_at_closest_collision,
                             foreign_mc_id, foreign_time_at_closest_collision);

  std::size_t number_of_motorcycles = motorcycles.size();
  for(std::size_t mc_2_id = 0; mc_2_id<number_of_motorcycles; ++mc_2_id)
  {
    if(mc.id() == mc_2_id)
      continue;

    Motorcycle& mc_2 = motorcycle(mc_2_id);

    FT time_at_collision = 0.;

    // mc_2's tentative track. If the foreign motorcycle has already crashed,
    // its tentative track stops at its current (and final) position.
    DEC_it mc_2_tentative_track_source = mc_2.source();
    DEC_it mc_2_tentative_track_target = mc_2.is_crashed() ? mc_2.position() :
                                                             mc_2.closest_target();
    Segment mc_2_tentative_track =
      K().construct_segment_2_object()(mc_2_tentative_track_source->point(),
                                       mc_2_tentative_track_target->point());

    std::cout << "  MC2: " << mc_2.id() << " tentative track ("
              << mc_2_tentative_track.source() << ")--("
              << mc_2_tentative_track.target() << ")" << std::endl;

    // Detect whether the motorcycles share the same supporting line
    // note that we know that 'mc_tentative_track' is not degenerate
    if(K().collinear_2_object()(mc_tentative_track.source(),
                                mc_tentative_track.target(),
                                mc_2.source()->point()) &&
       K().collinear_2_object()(mc_tentative_track.source(),
                                mc_tentative_track.target(),
                                mc_2.destination()->point()))
    {
      std::cerr << "  /!\\ Tracks are aligned" << std::endl;

      // Check if the tentative tracks actually intersect
      if(!K().do_intersect_2_object()(mc_tentative_track, mc_2_tentative_track))
      {
        // No intersection, try the next motorcycle
        std::cout << "  No intersection" << std::endl;
        continue;
      }

      // Many different configurations exist, e.g. (_S is for source, _T for target):
      //  MC_S   ---- MC_2_S ---- MC_2_T ---- MC_T
      //  MC_2_T ---- MC_S   ---- MC_2_S ---- MC_T
      // etc.
      // If, on the ray MC_S->MC_T,
      // - if MC_2_S is "before" MC_S, then it doesn't matter for MC
      // - if MC_2_S is "after" MC_S, then it depends on the motorcycles' directions

      // We need a ray because mc_2_s can be outside of mc's tentative track
      Ray mc_r(mc_tentative_track.source(), mc_tentative_track.target());
      if(!mc_r.has_on(mc_2_tentative_track_source->point()))
        continue;

      // Compute the respective direction of the two motorcycles:
      // use collinear_are_aligned_along_line ? @todo
      bool is_mc_2_degenerate = (mc_2.source() == mc_2.destination());
      bool are_motorcycles_moving_in_the_same_direction = (is_mc_2_degenerate ||
        K().angle_2_object()(mc.source()->point(), mc.destination()->point(),
                             mc_2.source()->point(), mc_2.destination()->point())
          == CGAL::ACUTE);

      std::cout << "  is degen: " << is_mc_2_degenerate << std::endl;
      std::cout << "  angle: " << K().angle_2_object()(mc.source()->point(), mc.destination()->point(),
                                                     mc_2.source()->point(), mc_2.destination()->point()) << std::endl;
      std::cout << "  are motorcycles moving in the same direction: "
                << are_motorcycles_moving_in_the_same_direction << std::endl;

      // The motorcycles move in the same direction ==> mc will impact mc_2's source
      if(are_motorcycles_moving_in_the_same_direction)
      {
        // The weird configuration of both motorcycles moving in the same direction
        // AND with the same source is checked at the very start, see function
        // 'crash_motorcycles_with_same_source_and_direction()'
        CGAL_assertion(is_mc_2_degenerate || mc.source() != mc_2.source());

        // @fixme for surfaces
        time_at_collision = CGAL::sqrt(CGAL::squared_distance(
          mc.source()->point(), mc_2.source()->point())) / mc.speed();
        std::cout << "  mc crashes into mc_2's source at time: " << time_at_collision << std::endl;

        // Compare with other possible collisions to keep the closest
        if(time_at_collision < time_at_closest_collision)
        {
          is_closest_collision_point_already_in_dictionary = true;
          time_at_closest_collision = time_at_collision;
          closest_collision_point = mc_2.source();
          foreign_mc_id = mc_2.id();
          foreign_time_at_closest_collision = 0; //@fixme surface
          continue;
        }
      }
      else // Motorcycles are moving towards each other
      {
        // Ignore the case where the motorcycles are moving away from each other
        // from the same point
        if(mc.source() == mc_2.source())
          continue;

        // If mc_2 has already crashed, then mc crashes into the final position of mc_2
        if(mc_2.is_crashed())
        {
          // @fixme for surfaces (see notes)
          time_at_collision = CGAL::sqrt(CGAL::squared_distance(
            mc.source()->point(), mc_2.position()->point())) / mc.speed();
          std::cout << "  mc crashes into mc_2's final position at: " << time_at_collision << std::endl;

          // Compare with other possible collisions to keep the closest
          if(time_at_collision < time_at_closest_collision)
          {
            is_closest_collision_point_already_in_dictionary = true;
            time_at_closest_collision = time_at_collision;
            closest_collision_point = mc_2.position();
            foreign_mc_id = mc_2.id();
            foreign_time_at_closest_collision = mc_2.current_time();
            continue;
          }
        }
        else // mc_2 has not yet crashed
        {
          // We now simply return the middle point that both motorcycles would reach
          // at the same time. Not that this point might not be reached by either
          // motorcycle, e.g. if they reach their destination first.

          // @fixme for surfaces (see notes)
          time_at_collision = CGAL::sqrt(CGAL::squared_distance(
            mc.source()->point(), mc_2.source()->point())) / (mc.speed() + mc_2.speed());
          std::cout << "  mc and mc_2 would meet at time: " << time_at_collision << std::endl;

          // Compare with other possible collisions to keep the closest
          if(time_at_collision < time_at_closest_collision)
          {
            is_closest_collision_point_already_in_dictionary = false;
            time_at_closest_collision = time_at_collision;
            Vector mc_v(mc_tentative_track);
            closest_collision = // @fixme for surfaces
              mc.source()->point() + time_at_collision * mc.speed() * mc_v / CGAL::sqrt(mc_v.squared_length());
            foreign_mc_id = mc_2.id();

            // @fixme for surfaces
            foreign_time_at_closest_collision = time_at_collision;
            continue;
          }
        }
      }
    }
    // --- From here on, the tracks are not collinear ---
    // Ignore the collision if it is the current position of the motorcycle
    else if(mc.position()->has_motorcycle(mc_2.id()))
    {
      continue;
    }
    // Check if the closest target of mc is the collision with mc_2
    else if(mc.closest_target()->has_motorcycle(mc_2.id()))
    {
      time_at_collision = mc.time_at_closest_target();
      std::cout << "  /!\\ tentative path collides with track: " << mc_2.id()
                << " at its end. Time: " << time_at_collision << std::endl;

      // Compare with other possible collisions to keep the closest
      if(time_at_collision < time_at_closest_collision)
      {
        is_closest_collision_point_already_in_dictionary = true;
        closest_collision_point = mc.closest_target();
        time_at_closest_collision = time_at_collision;
        foreign_mc_id = mc_2.id();
        foreign_time_at_closest_collision =
          closest_collision_point->find_motorcycle(mc_2.id())->second;
      }
    }
    // Nothing is easy, the intersection must be computed
    else
    {
      if(!K().do_intersect_2_object()(mc_tentative_track, mc_2_tentative_track))
      {
        // No intersection, try the next motorcycle
        std::cout << "  No intersection" << std::endl;
        continue;
      }

      std::cout << "  /!\\ collision with motorcycle: " << mc_2.id()
                << " between tracks ("
                << mc_tentative_track.source() << ")--("
                << mc_tentative_track.target() << ") and ("
                << mc_2_tentative_track.source() << ")--("
                << mc_2_tentative_track.target() << ")" << std::endl;

      const Point& collision_point = robust_intersection(mc_tentative_track,
                                                         mc_2_tentative_track);

      // Check if the intersection is within the face, ignore otherwise @todo

      // @fixme for surfaces
      time_at_collision =
          CGAL::sqrt(CGAL::squared_distance(mc.source()->point(), collision_point)) / mc.speed();
      std::cout << "  collision at: (" << collision_point << ") at time: " << time_at_collision << std::endl;

      // Compare with other collisions to keep the closest
      if(time_at_collision < time_at_closest_collision)
      {
        is_closest_collision_point_already_in_dictionary = false;
        time_at_closest_collision = time_at_collision;
        foreign_mc_id = mc_2.id();
        closest_collision = collision_point;
        foreign_time_at_closest_collision =
          CGAL::sqrt(CGAL::squared_distance(mc_2.source()->point(), collision_point)) / mc_2.speed();
      }
    }
  }

  // Insert the point in the dictionary. The motorcycle info will be added later.
  if(foreign_mc_id != std::size_t(-1) &&
     !is_closest_collision_point_already_in_dictionary)
    closest_collision_point = points.insert(closest_collision);

  return boost::make_tuple(closest_collision_point, time_at_closest_collision,
                           foreign_mc_id, foreign_time_at_closest_collision);
}

template<typename K, typename PolygonMesh>
void
Motorcycle_graph<K, PolygonMesh>::
initialize_motorcycles()
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "Initialize motorcycles" << std::endl;
#endif

  // if no mesh has been given in input, generate a mesh made of a single quad face
  // that contains all the interesting motorcycle interactions (crashes)
  if(using_enclosing_bbox)
    generate_enclosing_face();

  std::size_t number_of_motorcycles = motorcycles.size();
  for(std::size_t mc_id = 0; mc_id<number_of_motorcycles; ++mc_id)
  {
    Motorcycle& mc = motorcycle(mc_id);

    const Point& source_point = mc.initial_source_point();
    const FT speed = mc.speed();
    const FT time_at_source = mc.current_time();
    boost::optional<Point>& destination_point = mc.initial_destination_point();
    boost::optional<Vector>& direction = mc.direction();

    // add the source to the dictionary
    DEC_it source_entry = points.insert(source_point, mc_id, time_at_source);
    mc.source() = source_entry;
    mc.position() = source_entry;

    // locate the motorcycle in the graph
    locate_motorcycle(mc);

    // add the target to the dictionary
    DEC_it destination_entry;
    FT time_at_destination;
    if(destination_point == boost::none) // destination was not provided
    {
      std::pair<DEC_it, FT> destination = compute_motorcycle_next_destination(mc);
      destination_entry = destination.first;
      time_at_destination = destination.second;

      destination_entry->add_motorcycle(mc_id, destination.second);
      destination_point = destination_entry->point();
    }
    else // destination is known (but still need to compute the time)
    {
      // @tmp this should be computed by the tracer
      time_at_destination = time_at_source
        + CGAL::sqrt(CGAL::squared_distance(source_point, *destination_point)) / speed;
      destination_entry = points.insert(*destination_point, mc_id, time_at_destination);

      if(direction == boost::none)
      {
        mc.direction() = Vector(source_point, *destination_point);
      }
//    else @todo
//      check if (destination - source) is (roughly?) collinear with the direction
    }

    mc.destination() = destination_entry;

    // Initialize the motorcycle target queue
    mc.targets().insert(std::make_pair(source_entry, time_at_source));
    mc.targets().insert(std::make_pair(destination_entry, time_at_destination));

    // this is useful to not have empty track when sour=dest but creates duplicates @fixme
    mc.track().push_back(source_entry);
  }
}

template<typename K, typename PolygonMesh>
void
Motorcycle_graph<K, PolygonMesh>::
locate_motorcycle(Motorcycle& mc)
{
  if(using_enclosing_bbox || true) // @tmp
  {
    face_descriptor fd = *(faces(mesh).begin());
    // when using an enclosing quad face, barycentric coordinates are both unused
    // and meaningless
    Barycentric_coordinates bc;
    mc.set_location(std::make_pair(fd, bc));
  }

  // actual locate @todo
  face_descriptor fd = *(faces(mesh).begin());
  Barycentric_coordinates bc = { {0,0,0} }; // @todo actually compute them
  mc.set_location(std::make_pair(fd, bc));
}

template<typename K, typename PolygonMesh>
typename Motorcycle_graph<K, PolygonMesh>::Point
Motorcycle_graph<K, PolygonMesh>::
robust_intersection(const Segment &s, const Segment &t) const
{
  // this function should only be called in the case of non collinear segments
  // that are known to have a non-empty intersection
  CGAL_precondition(!K().collinear_2_object()(s.source(), s.target(), t.source()) ||
                    !K().collinear_2_object()(s.source(), s.target(), t.target()));
  CGAL_precondition(K().do_intersect_2_object()(s, t));

  const Point* pp;
  Point intersection_point;

  typedef typename CGAL::cpp11::result_of<
      typename K::Intersect_2(Segment, Segment)>::type Intersection_result;
  Intersection_result res = K().intersect_2_object()(s, t);

  bool need_to_use_exact = false;
  if(!res) // There should be an intersection
  {
    need_to_use_exact = true;
  }
  else
  {
    // The intersection cannot be a segment since the tracks are not collinear
    pp = boost::get<Point>(&*res);
    if(!pp || !s.has_on(*pp) || !t.has_on(*pp)) // another type of sanity check ? @todo
      need_to_use_exact = true;
  }

  if(need_to_use_exact)
  {
    typedef CGAL::Exact_predicates_exact_constructions_kernel EPECK;
    typedef EPECK::Point_2                                    Exact_point;
    typedef EPECK::Segment_2                                  Exact_segment;
    typedef typename CGAL::cpp11::result_of<
      typename EPECK::Intersect_2(Exact_segment, Exact_segment)>::type
                                                              Exact_intersection_result;

    EPECK epeck;
    CGAL::Cartesian_converter<K, EPECK> to_exact;
    CGAL::Cartesian_converter<EPECK, K> from_exact;
    Exact_segment exact_s = to_exact(s);
    Exact_segment exact_t = to_exact(t);

    CGAL_precondition(!epeck.collinear_2_object()(exact_s.source(), exact_s.target(), exact_t.source()) ||
                      !epeck.collinear_2_object()(exact_s.source(), exact_s.target(), exact_t.target()));
    CGAL_precondition(epeck.do_intersect_2_object()(exact_s, exact_t));

    Exact_intersection_result exact_res = epeck.intersect_2_object()(exact_s,
                                                                     exact_t);
    CGAL_assertion(exact_res);
    const Exact_point* epp = boost::get<Exact_point>(&*exact_res);
    const Exact_point exact_intersection_point = *epp;
    CGAL_assertion(epp);
    CGAL_postcondition(exact_s.has_on(exact_intersection_point) &&
                       exact_t.has_on(exact_intersection_point));
    intersection_point = from_exact(exact_intersection_point);
  }
  else
  {
    intersection_point = *pp;
  }

  return intersection_point;
}

template<typename K, typename PolygonMesh>
typename Motorcycle_graph<K, PolygonMesh>::Point
Motorcycle_graph<K, PolygonMesh>::
robust_intersection(const Ray &r, const Segment &s) const
{
  // This function should only be called if we are sure of an intersection between
  // the ray and the segment
  CGAL_precondition(K().do_intersect_2_object()(r, s));

  const Point* pp;
  Point intersection_point;

  typedef typename CGAL::cpp11::result_of<
      typename K::Intersect_2(Ray, Segment)>::type Intersection_result;
  Intersection_result res = K().intersect_2_object()(r, s);

  bool need_to_use_exact = false;
  if(!res) // There should be an intersection
  {
    need_to_use_exact = true;
  }
  else
  {
    // try to enforce a point return by switching to exact if it returns a segment
    // @fixme do something clean
    pp = boost::get<Point>(&*res);
    if(!pp || !s.has_on(*pp)) // another type of sanity check ? @todo
      need_to_use_exact = true;
  }

  if(need_to_use_exact)
  {
    typedef CGAL::Exact_predicates_exact_constructions_kernel EPECK;
    typedef EPECK::Point_2                                    Exact_point;
    typedef EPECK::Ray_2                                      Exact_ray;
    typedef EPECK::Segment_2                                  Exact_segment;
    typedef typename CGAL::cpp11::result_of<
      typename EPECK::Intersect_2(Exact_ray, Exact_segment)>::type
                                                              Exact_intersection_result;

    EPECK epeck;
    CGAL::Cartesian_converter<K, EPECK> to_exact;
    CGAL::Cartesian_converter<EPECK, K> from_exact;
    Exact_ray exact_r = to_exact(r);
    Exact_segment exact_s = to_exact(s);

    Exact_intersection_result exact_res = epeck.intersect_2_object()(exact_r,
                                                                     exact_s);
    CGAL_assertion(exact_res);

    const Exact_point* epp = boost::get<Exact_point>(&*exact_res);
    if(!epp)
    {
      std::cerr << "Warning: intersection is not a point" << std::endl;
      return r.source(); // ugly tmp hack @fixme
    }

    const Exact_point exact_intersection_point = *epp;
    CGAL_postcondition(exact_s.has_on(exact_intersection_point));
    intersection_point = from_exact(exact_intersection_point);
  }
  else
  {
    intersection_point = *pp;
  }

  return intersection_point;
}

template<typename K, typename PolygonMesh>
template<typename MotorcycleContainerIterator>
void
Motorcycle_graph<K, PolygonMesh>::
trace_motorcycle_graph(MotorcycleContainerIterator mit, MotorcycleContainerIterator last)
{
  add_motorcycles(mit, last);
  initialize_motorcycles();
  motorcycle_pq.initialize(motorcycles);

  // this can only happen at the beginning, simpler to get it out the way immediately
  crash_motorcycles_with_same_source_and_direction();

  while(!motorcycle_pq.empty())
  {
    // get the earliest available event
    Motorcycle_PQE pqe = motorcycle_pq.top();
    Motorcycle& mc = pqe.motorcycle();

    // move the motorcycle to the target (new confirmed position)
    drive_to_closest_target(mc);

    if(mc.position() == mc.destination())
    {
      if(has_motorcycle_reached_final_destination(mc))
      {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
        std::cout << "Reached final destination" << std::endl;
#endif
        crash_motorcycle(mc);
      }
      else
      {
        // the new destination is a pair: point, time
        std::pair<DEC_it, FT> new_destination = compute_motorcycle_next_destination(mc);
        new_destination.first->add_motorcycle(mc.id(), new_destination.second);
        mc.set_new_destination(new_destination.first, new_destination.second);

        // change the source @todo
        // change the location @todo
      }
    }
    else if(// only to prevent multiple motorcycles starting from the same source
            // (but with different directions) from blocking each other
            mc.current_time() != 0. &&
            // mc has reached the track of a foreign motorcycle
            (mc.has_reached_blocked_point() ||
            // multiple motorcycles will reach this point at the same time
             mc.has_reached_simultaneous_collision_point()))
    {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
      std::cout << "Reached crashing point: " << std::endl
                << " - blocked: " << mc.has_reached_blocked_point() << std::endl
                << " - simultaneous collision: " << mc.has_reached_simultaneous_collision_point() << std::endl;
#endif

      crash_motorcycle(mc);
    }
    else // the motorcycle can continue without issue towards its destination
    {
      // check for potential foreign tracks intersecting the next move of mc
      boost::tuple<DEC_it, FT, int, FT> res = find_collisions(mc);
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
         // does not know that collision point yet
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

          mc.add_new_target(collision_point, time_at_collision_point);
          mc.add_new_target(halving_point, time_at_halving_point);

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

            foreign_mc.add_new_target(collision_point, foreign_time_at_collision_point);
            foreign_mc.add_new_target(foreign_halving_point, foreign_time_at_halving_point);

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

      // The target list of 'mc' was modified and the queue must be updated
      if(!mc.is_crashed())
        motorcycle_pq.update(mc);
    }

    // block the newly reached position
    mc.position()->block();

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "---" << std::endl;
    std::cout << "Queue after update:" << std::endl << motorcycle_pq << std::endl;
#endif
  }
}

template<typename K, typename PolygonMesh>
bool
Motorcycle_graph<K, PolygonMesh>::
is_valid() const
{
  // @todo
  return true;
}

template<typename K, typename PolygonMesh>
void
Motorcycle_graph<K, PolygonMesh>::
output_all_dictionary_points() const
{
  typename Dictionary::Dictionary_entry_container::const_iterator dit = points.all_entries().begin();
  typename Dictionary::Dictionary_entry_container::const_iterator end = points.all_entries().end();

  std::ofstream os("dictionary_points.xyz");
  for(; dit!=end; ++dit)
    os << dit->point() << " 0" << '\n';
}

template<typename K, typename PolygonMesh>
void
Motorcycle_graph<K, PolygonMesh>::
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
