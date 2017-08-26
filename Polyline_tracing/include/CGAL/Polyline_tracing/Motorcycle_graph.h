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

#ifndef CGAL_MOTORCYCLE_GRAPH_H
#define CGAL_MOTORCYCLE_GRAPH_H

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
  typedef Motorcycle<K>                               Motorcycle;
  typedef Motorcycle_priority_queue<K>                Motorcycle_PQ;

  void find_intersections();

  template<typename MotorcycleSourcesIterator, typename MotorcycleDestinationsIterator>
  void initialize_motorcycles(MotorcycleSourcesIterator mit, MotorcycleSourcesIterator lastm,
                              MotorcycleDestinationsIterator dit, MotorcycleDestinationsIterator lastd);

  void trace_motorcycle();

  template<typename MotorcycleSourcesIterator, typename MotorcycleDestinationsIterator>
  void trace_motorcycle_graph(MotorcycleSourcesIterator mit, MotorcycleSourcesIterator lastm,
                              MotorcycleDestinationsIterator dit, MotorcycleDestinationsIterator lastd);

private:
  Dictionary target_points;
  std::vector<Motorcycle> motorcycles;
  Motorcycle_PQ motorcycle_pq;
};

template<typename K>
template<typename MotorcycleSourcesIterator, typename MotorcycleDestinationsIterator>
void
Motorcycle_graph<K>::
initialize_motorcycles(MotorcycleSourcesIterator mit, MotorcycleSourcesIterator lastm,
                       MotorcycleDestinationsIterator dit, MotorcycleDestinationsIterator lastd)
{
  CGAL_precondition(std::distance(mit, lastm) == std::distance(dit, lastd));
  motorcycles.reserve(std::distance(mit, lastm));

  int counter = 0;
  while(mit != lastm)
  {
    const Point s = *mit++;
    const Point d = *dit++;
    const FT distance = CGAL::squared_distance(s, d);

    Motorcycle m(counter, s, d);
    motorcycles.push_back(m);
    target_points.insert(m.position(), counter, 0.);
    target_points.insert(m.destination(), counter, distance);
    ++counter;
  }

  // some safety checks

}

template<typename K>
void
Motorcycle_graph<K>::
trace_motorcycle()
{

}

template<typename K>
template<typename MotorcyclesInputIterator, typename MotorcycleDestinationsIterator>
void
Motorcycle_graph<K>::
trace_motorcycle_graph(MotorcyclesInputIterator mit, MotorcyclesInputIterator lastm,
                       MotorcycleDestinationsIterator dit, MotorcycleDestinationsIterator lastd)
{
  initialize_motorcycles(mit, lastm, dit, lastd);
  target_points.initialize();
  motorcycle_pq.initialize();

  while(! motorcycle_pq.empty())
  {
    // get the earliest available event
    motorcycle_pq.top();

    // set confirmed

    motorcycle_pq.pop();

    // check for crashes
    if(true)
    {

    }
    else
    {
      // compute the destination if needed (?)
      trace_motorcycle();

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

#endif // CGAL_MOTORCYCLE_GRAPH_H
