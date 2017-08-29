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

  void crash_motorcycle(Motorcycle& m);
  void find_intersections();

  template<typename MotorcycleSourcesIterator, typename MotorcycleDestinationsIterator>
  void initialize_motorcycles(MotorcycleSourcesIterator mit, MotorcycleSourcesIterator lastm,
                              MotorcycleDestinationsIterator dit, MotorcycleDestinationsIterator lastd);
  void reach_closest_target(Motorcycle& m);
  void trace_motorcycle();

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
  motorcycle_pq.erase(motorcycle_pq.handle(m));
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

  int counter = 0;
  while(mit != lastm)
  {
    const Point s = *mit++;

    // @tmp if not provided, this must be computed by the tracer
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
find_intersections()
{
  // interface with the ray shooting data structure
}

template<typename K>
void
Motorcycle_graph<K>::
reach_closest_target(Motorcycle& mc)
{
  CGAL_assertion(!mc.targets().empty());
  DEC_it closest_target = mc.closest_target();
  mc.position() = closest_target;

  // remove the point that we have just reached from the list of targets
  mc.erase_closest_target();

  // @fixme update the priority queue already ? Probably can wait till new
  // points have been inserted.
}

template<typename K>
void
Motorcycle_graph<K>::
trace_motorcycle()
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

  while(! motorcycle_pq.empty())
  {
    // get the earliest available event
    Motorcycle_PQE pqe = motorcycle_pq.top();
    Motorcycle& m = pqe.motorcycle();

    // move the motorcycle to the target (new confirmed position)
    reach_closest_target(m);

    if(m.position() == m.destination()) // reached the destination
    {
      // check if we have reached a final destination (border) @todo
      if(true)
      {
        crash_motorcycle(m);
      }
      else
      {
        // check that S_i is empty @todo
        CGAL_assertion(true);

        // otherwise, compute the new destination
        trace_motorcycle();

        // set the new tentative point
      }
    }
    else if(m.position()->second.is_blocked() || // impacted an existing polyline
            false) // @todo reached the point at the same time as another polyline
    {
      crash_motorcycle(m);
    }
    else // the motorcycle can continue without issue towards its destination
    {
      // check for intersections
      find_intersections();

      if(true)
      {

      }
      else
      {

      }
    }
  }
}

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_MOTORCYCLE_GRAPH_H
