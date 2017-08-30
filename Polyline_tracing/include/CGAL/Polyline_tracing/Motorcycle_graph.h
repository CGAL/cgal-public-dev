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

#include <boost/tuple/tuple.hpp>

#include <algorithm>
#include <vector>

namespace CGAL {

namespace Polyline_tracing {

template<typename K>
class Motorcycle_graph
{
public:
  typedef typename K::FT                              FT;
  typedef typename K::Point_2                         Point;

  typedef Dictionary<K>                               Dictionary;
  typedef Dictionary_entry<K>                         Dictionary_entry;
  typedef typename Dictionary::DEC_it                 DEC_it;

  typedef Motorcycle<K>                               Motorcycle;
  typedef std::vector<Motorcycle>                     Motorcycle_container;

  typedef Motorcycle_priority_queue<K>                Motorcycle_PQ;
  typedef Motorcycle_priority_queue_entry<K>          Motorcycle_PQE;

  std::pair<DEC_it, FT> compute_motorcycle_next_destination(const Motorcycle& m);
  void crash_motorcycle(Motorcycle& m);
  boost::tuple<int, DEC_it, FT> find_collisions(const Motorcycle& m) const;

  template<typename MotorcycleSourcesIterator, typename MotorcycleDestinationsIterator>
  void initialize_motorcycles(MotorcycleSourcesIterator mit, MotorcycleSourcesIterator lastm,
                              MotorcycleDestinationsIterator dit, MotorcycleDestinationsIterator lastd);
  void drive_to_closest_target(Motorcycle& m);

  template<typename MotorcycleSourcesIterator, typename MotorcycleDestinationsIterator>
  void trace_motorcycle_graph(MotorcycleSourcesIterator mit, MotorcycleSourcesIterator lastm,
                              MotorcycleDestinationsIterator dit, MotorcycleDestinationsIterator lastd);

private:
  Dictionary points;
  Motorcycle_container motorcycles;
  Motorcycle_PQ motorcycle_pq;
};
// -----------------------------------------------------------------------------

template<typename K>
void
Motorcycle_graph<K>::
crash_motorcycle(Motorcycle& m)
{
  m.crash();

  // go through the next targets of the motorcycle and remove it from the list
  // of motorcycles that might reach the target
  typename Motorcycle::Target_point_container::iterator it = m.targets().begin();
  typename Motorcycle::Target_point_container::iterator end = m.targets().end();
  for(; it!=end; ++it)
  {
    DEC_it target_point = it->first;
    target_point->second.remove_motorcycle(m.id());
  }

  motorcycle_pq.erase(motorcycle_pq.handle(m));
}

template<typename K>
boost::tuple<int, typename Motorcycle_graph<K>::DEC_it, typename Motorcycle_graph<K>::FT>
Motorcycle_graph<K>::
find_collisions(const Motorcycle& /*m*/) const
{
  // interface with the ray shooting data structure
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

    // @tmp if not provided, this should be computed by the tracer
    const Point d = *dit++;

    // @tmp this should be computed by the tracer (to take e.g. speed into account)
    const FT distance = CGAL::squared_distance(s, d);

    // add the points to the dictionary
    DEC_it source_entry = points.insert(s, counter, 0.);
    DEC_it destination_entry = points.insert(d, counter, distance);

    // create the motorcycle
    Motorcycle m(counter, source_entry, destination_entry, distance);
    motorcycles.push_back(m);

    ++counter;
  }

  // some safety checks (two motorcycles sharing a supporting line, ...) @todo

}

template<typename K>
void
Motorcycle_graph<K>::
drive_to_closest_target(Motorcycle& mc)
{
  CGAL_assertion(!mc.is_crashed());
  CGAL_assertion(!mc.targets().empty());

  DEC_it closest_target = mc.closest_target();
  closest_target->second.block();

  mc.position() = closest_target;
  mc.distance() = mc.targets().begin()->second;
  mc.path().push_back(closest_target);
  mc.erase_closest_target();

  // @fixme update the priority queue already ? Probably can wait till new
  // points have been inserted.

}

template<typename K>
std::pair<typename Motorcycle_graph<K>::DEC_it, typename Motorcycle_graph<K>::FT>
Motorcycle_graph<K>::
compute_motorcycle_next_destination(const Motorcycle& /*m*/)
{
  // interface with the tracer data structure
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

  while(!motorcycle_pq.empty())
  {
    // get the earliest available event
    Motorcycle_PQE pqe = motorcycle_pq.top();
    Motorcycle& mc = pqe.motorcycle();

    // move the motorcycle to the target (new confirmed position)
    drive_to_closest_target(mc);

    if(mc.position() == mc.destination())
    {
      // the destination should be the last target
      CGAL_assertion(mc.targets().empty());

      if(mc.is_motorcycle_destination_final())
      {
        crash_motorcycle(mc);
      }
      else
      {
        // the new destination is a pair: point, distance
        std::pair<DEC_it, FT> new_destination = compute_motorcycle_next_destination(mc);
        mc.set_new_destination(new_destination.first, new_destination.second);
      }
    }
    else if(// this motorcycle has impacted an existing trace
            mc.has_reached_blocked_point() ||
            // multiple motorcycles will reach this point at the same time
            mc.has_reached_simultaneous_collision_point())
    {
      crash_motorcycle(mc);
    }
    else // the motorcycle can continue without issue towards its destination
    {
      // check if other traces exist between the confirmed point and the next
      // tentative point
      boost::tuple<int, DEC_it, FT> res = find_collisions(mc);
      const int foreign_motorcycle_id = res.get<1>();
      DEC_it collision_point = res.get<2>();
      const FT distance_at_collision_point = res.get<3>();

      if(// there is intersection
         foreign_motorcycle_id != -1 &&
         // the impact is closer than the next target
         distance_at_collision_point <= mc.distance_at_closest_target() &&
         // the collision is not the next target that would also be on the track
         (collision_point != mc.closest_target() ||
          collision_point->has_motorcycle(foreign_motorcycle_id)))
      {
        if(!collision_point->has_motorcycle(mc.id()))
        {
          // Call the halving structure to create a new point
          DEC_it middle_point; // @todo
          const FT distance_at_middle_point; // @todo

          // add as targets
          mc.add_new_target(middle_point, distance_at_middle_point);
          mc.add_new_target(collision_point, distance_at_collision_point);

          // add to M(p')
          middle_point->add_motorcycle(mc.id(), distance_at_middle_point);
          collision_point->add_motorcycle(mc.id(), distance_at_collision_point);
        }

        // deal with the foreign motorcyle
        if(true)
        {

        }
      }
    }
  }
}

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_MOTORCYCLE_GRAPH_H
