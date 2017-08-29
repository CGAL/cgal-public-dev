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

#ifndef CGAL_POLYLINE_TRACING_MOTORCYCLE_PRIORITY_QUEUE_H
#define CGAL_POLYLINE_TRACING_MOTORCYCLE_PRIORITY_QUEUE_H

#include <CGAL/Polyline_tracing/Dictionary.h>
#include <CGAL/Polyline_tracing/Motorcycle.h>
#include <CGAL/Polyline_tracing/Motorcycle_priority_queue_entry.h>

#include <boost/heap/fibonacci_heap.hpp>

#include <utility>

namespace CGAL {

namespace Polyline_tracing {

template<typename K>
class Motorcycle_priority_queue
{
public:
  typedef Motorcycle<K>                            Motorcycle;
  typedef std::vector<Motorcycle>                  Motorcycle_container;

  // Picked a fibonacci_heap for now. Would it better to simply use std::priority_queue
  // and ignore+pop values that are meaningless ?
  typedef Motorcycle_priority_queue_entry<K>       MPQ_entry;
  typedef boost::heap::fibonacci_heap<MPQ_entry>   MPQ;
  typedef typename MPQ::handle_type                handle_type;

  bool empty() const { return queue.empty(); }
  const MPQ_entry& top() const { return queue.top(); }
  handle_type handle(const Motorcycle& m) const { return handles[m.id()]; }

  handle_type push(const Motorcycle& m) { queue.push(m); }
  void pop() { return queue.pop(); }
  void update(handle_type handle) { return queue.update(handle); }
  void erase(handle_type handle) { return queue.erase(handle); }

  Motorcycle_priority_queue() : queue(), handles() { }

  void initialize(Motorcycle_container& motorcycles);

private:
  MPQ queue;
  std::vector<handle_type> handles; // maps motorcycle_ids to their handle in the PQ
};

template<typename K>
void
Motorcycle_priority_queue<K>::
initialize(Motorcycle_container& motorcycles)
{
  typename Motorcycle_container::iterator m_it = motorcycles.begin();
  typename Motorcycle_container::iterator last = motorcycles.end();
  handles.resize(std::distance(m_it, last));

  for(; m_it != last; ++m_it)
  {
    const int motorcycle_id = m_it->id();
    CGAL_precondition(motorcycle_id >= 0 && motorcycle_id < handles.size());
    handles[motorcycle_id] = queue.push(*m_it);
  }
}

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_MOTORCYCLE_PRIORITY_QUEUE_H
