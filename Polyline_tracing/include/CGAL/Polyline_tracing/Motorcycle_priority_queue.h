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
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_POLYLINE_TRACING_MOTORCYCLE_PRIORITY_QUEUE_H
#define CGAL_POLYLINE_TRACING_MOTORCYCLE_PRIORITY_QUEUE_H

#include <CGAL/Polyline_tracing/Motorcycle.h>
#include <CGAL/Polyline_tracing/Motorcycle_graph.h>
#include <CGAL/Polyline_tracing/Motorcycle_graph_node_dictionary.h>
#include <CGAL/Polyline_tracing/Motorcycle_priority_queue_entry.h>

#include <boost/heap/fibonacci_heap.hpp>
#include <boost/shared_ptr.hpp>

#include <utility>

namespace CGAL {

namespace Polyline_tracing {

template<typename MotorcycleGraph>
class Motorcycle_priority_queue
{
  typedef Motorcycle_priority_queue<MotorcycleGraph>             Self;

public:
  typedef typename MotorcycleGraph::Geom_traits                  Geom_traits;

  typedef typename MotorcycleGraph::Motorcycle                   Motorcycle;
  typedef typename MotorcycleGraph::Motorcycle_container         Motorcycle_container;

  // Picked a fibonacci_heap for now. Would it be better to simply use
  // std::priority_queue and ignore+pop values that are meaningless ?
  typedef Motorcycle_priority_queue_entry<MotorcycleGraph>       MPQ_entry;
  typedef boost::heap::fibonacci_heap<MPQ_entry>                 MPQ;
  typedef typename MPQ::handle_type                              handle_type;

  bool empty() const { return queue.empty(); }
  std::size_t size() const { return queue.size(); }
  const MPQ_entry& top() const { return queue.top(); }
  handle_type handle(const Motorcycle& m) const { return handles[m.id()]; }
  handle_type push(const Motorcycle& m) { queue.push(&m); }
  void update(const Motorcycle& m) { return queue.update(handles[m.id()]); }
  void erase(const Motorcycle& m) { return queue.erase(handles[m.id()]); }

  Motorcycle_priority_queue() : queue(), handles() { }

  void initialize(Motorcycle_container& motorcycles);

  // output
  friend std::ostream& operator<<(std::ostream& out, const Self& mpq)
  {
    out << "###################################################" << std::endl;
    out << "---------------------------------------------------" << std::endl;
    if(mpq.queue.empty())
    {
      out << "Empty !" << std::endl;
      return out;
    }

    typename MPQ::ordered_iterator pq_it = mpq.queue.ordered_begin();
    typename MPQ::ordered_iterator end = mpq.queue.ordered_end();
    for(; pq_it!=end; ++pq_it)
    {
      const Motorcycle& mc = pq_it->motorcycle();

      out << "  MC#" << std::setw(4) << mc.id();

      if(mc.is_initialized())
      {
        out << " with position: " << &*(mc.current_position())
            << " (" << mc.current_position()->point() << ")"
            << " at: " << mc.current_time()
            << " and target: " << &*(mc.closest_target())
            << " (" << mc.closest_target()->point() << ")";
      }
      else
      {
        out << " (uninitialized)";
      }

      out << " at: " << pq_it->time_at_closest_target() << std::endl;
    }

    out << "---------------------------------------------------" << std::endl;
    out << "###################################################" << std::endl;
    return out;
  }

private:
  MPQ queue;
  std::vector<handle_type> handles; // maps motorcycle_ids to their handle in the PQ
};

template<typename MotorcycleGraphTraits>
void
Motorcycle_priority_queue<MotorcycleGraphTraits>::
initialize(Motorcycle_container& motorcycles)
{
  typename Motorcycle_container::iterator m_it = motorcycles.begin();
  typename Motorcycle_container::iterator last = motorcycles.end();
  handles.resize(std::distance(m_it, last));

  for(; m_it!=last; ++m_it)
  {
    const std::size_t motorcycle_id = m_it->id();
    CGAL_precondition(motorcycle_id >= 0 && motorcycle_id < handles.size());

    if(m_it->status() == Motorcycle::IN_MOTION)
      handles[motorcycle_id] = queue.push(&*m_it);
  }

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "Initial queue:" << std::endl << *this << std::endl;
#endif
}

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_MOTORCYCLE_PRIORITY_QUEUE_H
