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

#include <CGAL/assertions.h>
#include <CGAL/enum.h>
#include <CGAL/intersection_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/result_of.h>

#include <boost/tuple/tuple.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

namespace CGAL {

namespace Polyline_tracing {

template<typename K>
class Motorcycle_graph
{
public:
  typedef typename K::FT                              FT;
  typedef typename K::Point_2                         Point;
  typedef typename K::Segment_2                       Segment;
  typedef typename K::Vector_2                        Vector;
  typedef typename K::Ray_2                           Ray;

  typedef Dictionary<K>                               Dictionary;
  typedef Dictionary_entry<K>                         Dictionary_entry;
  typedef typename Dictionary::DEC_it                 DEC_it;

  typedef Motorcycle<K>                               Motorcycle;
  typedef std::vector<Motorcycle>                     Motorcycle_container;

  typedef Motorcycle_priority_queue<K>                Motorcycle_PQ;
  typedef Motorcycle_priority_queue_entry<K>          Motorcycle_PQE;

  Motorcycle& motorcycle(const int id) { return motorcycles[id]; }
  const Motorcycle& motorcycle(const int id) const { return motorcycles[id]; }

  std::pair<DEC_it, FT> compute_halving_point(const Motorcycle& m,
                                              DEC_it p, const FT p_time,
                                              DEC_it q, const FT q_time);
  std::pair<DEC_it, FT> compute_middle_point(DEC_it p, const FT p_time,
                                             DEC_it q, const FT q_time);
  std::pair<DEC_it, FT> compute_motorcycle_next_destination(const Motorcycle& m) const;
  void crash_motorcycle(Motorcycle& m);
  void drive_to_closest_target(Motorcycle& m);
  boost::tuple<DEC_it, FT, int, FT> find_collisions(const Motorcycle& m);

  template<typename MotorcycleSourcesIterator, typename MotorcycleDestinationsIterator>
  void initialize_motorcycles(MotorcycleSourcesIterator mit, MotorcycleSourcesIterator lastm,
                              MotorcycleDestinationsIterator dit, MotorcycleDestinationsIterator lastd);

  template<typename MotorcycleSourcesIterator, typename MotorcycleDestinationsIterator>
  void trace_motorcycle_graph(MotorcycleSourcesIterator mit, MotorcycleSourcesIterator lastm,
                              MotorcycleDestinationsIterator dit, MotorcycleDestinationsIterator lastd);

  // Output
  void output_all_dictionary_points() const;
  void output_motorcycles_sources_and_destinations() const;

private:
  Dictionary points;
  Motorcycle_container motorcycles;
  Motorcycle_PQ motorcycle_pq;
};
// -----------------------------------------------------------------------------

template<typename K>
std::pair<typename Motorcycle_graph<K>::DEC_it, typename Motorcycle_graph<K>::FT>
Motorcycle_graph<K>::
compute_middle_point(DEC_it p, const FT p_time, DEC_it q, const FT q_time)
{
  Point r = K().construct_midpoint_2_object()(p->point(), q->point());
  const FT time_at_r = 0.5 * (p_time + q_time);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "  New middle point: (" << r  << ") at time: " << time_at_r << std::endl;
#endif

  return std::make_pair(points.insert(r), time_at_r);
}

template<typename K>
std::pair<typename Motorcycle_graph<K>::DEC_it, typename Motorcycle_graph<K>::FT>
Motorcycle_graph<K>::
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

template<typename K>
std::pair<typename Motorcycle_graph<K>::DEC_it, typename Motorcycle_graph<K>::FT>
Motorcycle_graph<K>::
compute_motorcycle_next_destination(const Motorcycle& /*m*/) const
{
  // interface with the tracer data structure @todo
  return std::make_pair(DEC_it(), FT());
}

template<typename K>
void
Motorcycle_graph<K>::
crash_motorcycle(Motorcycle& m)
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~X" << std::endl;
  std::cout << "Crashing " << m;
#endif

  m.crash();

  // go through the next targets of the motorcycle and remove it from the list
  // of motorcycles that might reach the target
  typename Motorcycle::Target_point_container::iterator it = m.targets().begin();
  typename Motorcycle::Target_point_container::iterator end = m.targets().end();
  for(; it!=end; ++it)
  {
    DEC_it target_point = it->first;
    target_point->remove_motorcycle(m.id());
  }

  motorcycle_pq.erase(m);
}

template<typename K>
void
Motorcycle_graph<K>::
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

template<typename K>
boost::tuple<typename Motorcycle_graph<K>::DEC_it, typename Motorcycle_graph<K>::FT,
             int, typename Motorcycle_graph<K>::FT>
Motorcycle_graph<K>::
find_collisions(const Motorcycle& mc)
{
  // THIS WHOLE FUNCTION IS @TMP
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "Checking for collisions on motorcycle " << mc.id() << "'s track" << std::endl;
#endif

  CGAL_precondition(!mc.is_crashed());
  CGAL_precondition(!mc.targets().empty());
  // search for a possible collision between the current position of mc
  // and the next target

  // interface with the ray shooting data structure @todo

  // brute force for now:
  Point closest_collision_point;
  int foreign_mc_id = -1;
  FT time_at_closest_collision = std::numeric_limits<FT>::max();
  FT foreign_time_at_closest_collision = std::numeric_limits<FT>::max();
  DEC_it new_point;
  bool is_new_point_already_in_dictionary = false;

  // only consider the part of the tentative track that is unconfirmed
  Segment mc_tentative_track =
    K().construct_segment_2_object()(mc.position()->point(),
                                     mc.closest_target()->point());

  std::cout << "MC tentative track: " << std::endl
            << *(mc.position()) << std::endl
            << *(mc.closest_target()) << std::endl;

  // a degenerate tentative track has no intersecting collisions
  if(mc_tentative_track.is_degenerate())
    return boost::make_tuple(new_point, time_at_closest_collision,
                             foreign_mc_id, foreign_time_at_closest_collision);

  typename Motorcycle_container::const_iterator mc_it = motorcycles.begin();
  typename Motorcycle_container::const_iterator end = motorcycles.end();
  for(; mc_it!=end; ++mc_it)
  {
    const Motorcycle& mc_2 = *mc_it;

    if(mc.id() == mc_2.id())
      continue;

    FT time_at_collision = 0.;

    // mc_2's tentative track. If the foreign motorcycle has already crashed,
    // its tentative track stops at the current (and final) position.
    DEC_it mc_2_tentative_track_source = mc_2.source();
    DEC_it mc_2_tentative_track_target = mc_2.is_crashed() ? mc_2.position() :
                                                             mc_2.closest_target();
    Segment mc_2_tentative_track =
        K().construct_segment_2_object()(mc_2_tentative_track_source->point(),
                                         mc_2_tentative_track_target->point());

    std::cout << "  MC2: " << mc_2.id()
              << " tentative track ("
              << mc_2_tentative_track.source() << ")--("
              << mc_2_tentative_track.target() << ")" << std::endl;

    // detect whether the motorcycles share the same supporting line
    // note that we know that 'mc_tentative_track' is not degenerate
    if(K().collinear_2_object()(mc_tentative_track.source(),
                                mc_tentative_track.target(),
                                mc_2.source()->point()) &&
       K().collinear_2_object()(mc_tentative_track.source(),
                                mc_tentative_track.target(),
                                mc_2.destination()->point()))
    {
      std::cerr << "  /!\\ Tracks are aligned" << std::endl;

      // check if the tentative tracks actually intersect
      if(!K().do_intersect_2_object()(mc_tentative_track, mc_2_tentative_track))
      {
        // no intersection, try the next motorcycle
        std::cout << "  No intersection" << std::endl;
        continue;
      }

      // Many different configurations exist, e.g. (_S is for source, _T for target):
      //  MC_S   ---- MC_2_S ---- MC_2_T ---- MC_T
      //  MC_2_T ---- MC_S   ---- MC_2_S ---- MC_T
      // etc.
      // - if MC_2_S is "before" MC_S, then it doesn't matter for MC
      // - if MC_2_S is "after" MC_S, then it depends on the motorcycles' directions

      // We need a ray because mc_2_s can be outside of mc's tentative track
      Ray mc_r(mc_tentative_track.source(), mc_tentative_track.target());
      if(!mc_r.has_on(mc_2_tentative_track_source->point()))
        continue;

      // compute the respective movement of the two motorcycles:
      bool is_mc_2_degenerate = (mc_2.source() == mc_2.destination());

      std::cout << "is degen: " << is_mc_2_degenerate << std::endl;
      std::cout << "angle: " << K().angle_2_object()(mc.source()->point(), mc.destination()->point(),
                                                     mc_2.source()->point(), mc_2.destination()->point()) << std::endl;

      bool are_motorcycles_moving_in_the_same_direction = (is_mc_2_degenerate ||
        K().angle_2_object()(mc.source()->point(), mc.destination()->point(),
                             mc_2.source()->point(), mc_2.destination()->point())
          == CGAL::ACUTE);

      std::cout << "  are motorcycles moving in the same direction: "
                << are_motorcycles_moving_in_the_same_direction << std::endl;

      // motorcycles move in the same direction, mc will impact mc_2's source
      if(are_motorcycles_moving_in_the_same_direction)
      {
        // both motorcycles moving in the same direction AND with the same source
        // is just weird... Disable it for now.
        // @todo handle the case by crashing both at their sources ?
        CGAL_precondition(mc.source() != mc_2.source());

        // @fixme for surfaces
        time_at_collision = CGAL::sqrt(CGAL::squared_distance(
          mc.source()->point(), mc_2.source()->point())) / mc.speed();
        std::cout << "  mc crashes into mc_2's source at time: " << time_at_collision << std::endl;

        if(time_at_collision < time_at_closest_collision)
        {
          is_new_point_already_in_dictionary = true;
          time_at_closest_collision = time_at_collision;
          new_point = mc_2.source();
          foreign_mc_id = mc_2.id();
          foreign_time_at_closest_collision = 0; //@fixme surface
          continue;
        }
      }
      else // motorcycles are moving towards each other
      {
        // Ignore the case where the motorcycles are moving away from each other
        // from the same point
        if(mc.source() == mc_2.source())
          continue;

        // crash into the destination of mc_2
        if(mc_2.is_crashed())
        {
          // @fixme for surfaces (see notes)
          time_at_collision = CGAL::sqrt(CGAL::squared_distance(
            mc.source()->point(), mc_2.position()->point())) / mc.speed();
          std::cout << "  mc crashes into mc_2's final position at: " << time_at_collision << std::endl;

          if(time_at_collision < time_at_closest_collision)
          {
            is_new_point_already_in_dictionary = true;
            time_at_closest_collision = time_at_collision;
            new_point = mc_2.position();
            foreign_mc_id = mc_2.id();
            foreign_time_at_closest_collision = mc_2.current_time();
            continue;
          }
        }
        else // mc_2 is not crashed
        {
          // we now simply return the middle point that both motorcycles reach
          // at the same time. Not that this point might not be reached by either
          // motorcycle, if they reach their destination first.

          // @fixme for surfaces (see notes)
          time_at_collision = CGAL::sqrt(CGAL::squared_distance(
            mc.source()->point(), mc_2.source()->point())) / (mc.speed() + mc_2.speed());
          std::cout << "  mc and mc_2 would meet at time: " << time_at_collision << std::endl;

          if(time_at_collision < time_at_closest_collision)
          {
            is_new_point_already_in_dictionary = false;
            time_at_closest_collision = time_at_collision;
            Vector mc_v(mc_tentative_track);
            closest_collision_point = // @fixme for surfaces
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
    // Ignore the intersection if it is the current position of the motorcycle
    else if(mc.position()->has_motorcycle(mc_2.id()))
    {
      continue;
    }
    // Check if the target of mc is the intersection with mc_2
    else if(mc.closest_target()->has_motorcycle(mc_2.id()))
    {
      time_at_collision = mc.time_at_closest_target();
      std::cout << "  /!\\ tentative path collides with track: " << mc_2.id()
                << " at its end. Time: " << time_at_collision << std::endl;

      // compare with other intersections to keep the closest
      if(time_at_collision < time_at_closest_collision)
      {
        is_new_point_already_in_dictionary = true;
        new_point = mc.closest_target();
        time_at_closest_collision = time_at_collision;
        foreign_mc_id = mc_2.id();
        foreign_time_at_closest_collision =
          new_point->find_motorcycle(mc_2.id())->second;
      }
    }
    // nothing is easy, the intersection must be computed
    else
    {
      if(!K().do_intersect_2_object()(mc_tentative_track, mc_2_tentative_track))
      {
        // no intersection, try the next motorcycle
        std::cout << "  No intersection" << std::endl;
        continue;
      }

      std::cout << "  /!\\ intersection with motorcycle: " << mc_2.id()
                << " between tracks ("
                << mc_tentative_track.source() << ")--("
                << mc_tentative_track.target() << ") and ("
                << mc_2_tentative_track.source() << ")--("
                << mc_2_tentative_track.target() << ")" << std::endl;

      typedef typename CGAL::cpp11::result_of<
          typename K::Intersect_2(Segment, Segment)>::type Intersection_result;
      Intersection_result res = CGAL::intersection(mc_tentative_track,
                                                   mc_2_tentative_track);

      // there must be an intersection, and it cannot be a segment
      // since the tracks are not collinear
      const Point* pp = boost::get<Point>(&*res);
      CGAL_assertion(pp);
      const Point& collision_point = *pp;

      // check that the intersection is in the face @todo

      // @fixme for surfaces
      time_at_collision =
          CGAL::sqrt(CGAL::squared_distance(mc.source()->point(), collision_point)) / mc.speed();
      std::cout << "  intersection at: (" << collision_point << ") at time: " << time_at_collision << std::endl;

      // compare with other intersections to keep the closest
      if(time_at_collision < time_at_closest_collision)
      {
        is_new_point_already_in_dictionary = false;
        time_at_closest_collision = time_at_collision;
        foreign_mc_id = mc_2.id();
        closest_collision_point = collision_point;
        foreign_time_at_closest_collision =
          CGAL::sqrt(CGAL::squared_distance(mc_2.source()->point(), collision_point)) / mc_2.speed();
      }
    }
  }

  // only try to create the new point if there is an intersection, except if it
  // is already known because it was the tentative point of 'mc' (in that case,
  // 'new_point' is already set up)
  if(foreign_mc_id != -1 && !is_new_point_already_in_dictionary)
    new_point = points.insert(closest_collision_point);

  return boost::make_tuple(new_point, time_at_closest_collision,
                           foreign_mc_id, foreign_time_at_closest_collision);
}

template<typename K>
template<typename MotorcycleSourcesIterator, typename MotorcycleDestinationsIterator>
void
Motorcycle_graph<K>::
initialize_motorcycles(MotorcycleSourcesIterator mit, MotorcycleSourcesIterator lastm,
                       MotorcycleDestinationsIterator dit,
                       MotorcycleDestinationsIterator CGAL_precondition_code(lastd))
{
  CGAL_precondition(std::distance(mit, lastm) == std::distance(dit, lastd));
  motorcycles.reserve(std::distance(mit, lastm));

  int counter = 0; // unique motorcycle ids
  while(mit != lastm)
  {
    const Point s = *mit++;
    const FT speed = 1.; // this should be provided in input @todo

    // @tmp if not provided, this should be computed by the tracer
    const Point d = *dit++;

    if(s == d)
    {
      std::cerr << "Warning: source and destination (" << s << ") "
                << "are equal for motorcycle: " << counter << std::endl;
    }

    // @tmp this should be computed by the tracer (to take e.g. speed into account)
    const FT time_at_source = 0.;
    const FT time_at_destination = time_at_source + CGAL::sqrt(CGAL::squared_distance(s, d)) / speed;

    // add the points to the dictionary
    DEC_it source_entry = points.insert(s, counter, time_at_source);
    DEC_it destination_entry = points.insert(d, counter, time_at_destination);

    // create the motorcycle
    Motorcycle m(counter, source_entry, destination_entry, speed,
                 time_at_source, time_at_destination);
    motorcycles.push_back(m);

    ++counter;
  }

  // some safety checks (two motorcycles sharing a supporting line, ...) @todo

}

template<typename K>
template<typename MotorcyclesInputIterator, typename MotorcycleDestinationsIterator>
void
Motorcycle_graph<K>::
trace_motorcycle_graph(MotorcyclesInputIterator mit, MotorcyclesInputIterator lastm,
                       MotorcycleDestinationsIterator dit, MotorcycleDestinationsIterator lastd)
{
  initialize_motorcycles(mit, lastm, dit, lastd);
  motorcycle_pq.initialize(motorcycles);
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "Initial queue: " << std::endl << motorcycle_pq << std::endl;
#endif

  while(!motorcycle_pq.empty())
  {
    // get the earliest available event
    Motorcycle_PQE pqe = motorcycle_pq.top();
    Motorcycle& mc = pqe.motorcycle();

    // move the motorcycle to the target (new confirmed position)
    drive_to_closest_target(mc);

    if(mc.position() == mc.destination())
    {
      if(mc.is_motorcycle_destination_final())
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
        mc.set_new_destination(new_destination.first, new_destination.second);
      }
    }
    else if(// just to avoid multiple motorcycles starting from the same source
            // blocking each other
            mc.current_time() != 0. &&
            // this motorcycle has impacted an existing trace
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
      // check if other traces exist between the confirmed point and the next
      // tentative point
      boost::tuple<DEC_it, FT, int, FT> res = find_collisions(mc);
      DEC_it collision_point = res.template get<0>();
      const FT time_at_collision_point = res.template get<1>();
      const int foreign_motorcycle_id = res.template get<2>();
      const FT foreign_time_at_collision_point = res.template get<3>();

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
      std::cout << "Find collision results: " << std::endl
                << " - foreign_motorcycle_id: " << foreign_motorcycle_id << std::endl;
      if(foreign_motorcycle_id != -1)
      {
        std::cout << " - collision_point: " << *(collision_point) << std::endl
                  << " - time_at_collision_point: " << time_at_collision_point << std::endl
                  << " - foreign_time_at_collision_point: " << foreign_time_at_collision_point << std::endl;
      }
#endif

      if(// there is intersection
         foreign_motorcycle_id != -1 &&
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

          Motorcycle& foreign_mc = motorcycles[foreign_motorcycle_id];
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
      motorcycle_pq.update(mc);
    }

    // block the current position
    mc.position()->block();

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "---" << std::endl;
    std::cout << "Queue after update:" << std::endl << motorcycle_pq << std::endl;
#endif
  }
}

template<typename K>
void
Motorcycle_graph<K>::
output_all_dictionary_points() const
{
  typename Dictionary::Dictionary_entry_container::const_iterator dit = points.all_entries().begin();
  typename Dictionary::Dictionary_entry_container::const_iterator end = points.all_entries().end();

  std::ofstream os("dictionary_points.xyz");
  for(; dit!=end; ++dit)
    os << dit->point() << " 0" << '\n';
}

template<typename K>
void
Motorcycle_graph<K>::
output_motorcycles_sources_and_destinations() const
{
  // must be adapted to surfaces @todo

  std::ofstream oss("out_motorcycles_sources.xyz");
  std::ofstream osd("out_motorcycles_destinations.xyz");
  for(std::size_t i=0; i<motorcycles.size(); ++i)
  {
    oss << motorcycles[i].source()->point() << " 0" << '\n';
    osd << motorcycles[i].destination()->point() << " 0" << '\n';
  }
}

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_MOTORCYCLE_GRAPH_H
