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

#ifndef CGAL_POLYLINE_TRACING_MOTORCYCLE_PRIORITY_QUEUE_ENTRY_H
#define CGAL_POLYLINE_TRACING_MOTORCYCLE_PRIORITY_QUEUE_ENTRY_H

#include <CGAL/Polyline_tracing/Motorcycle_graph_node_dictionary.h>
#include <CGAL/Polyline_tracing/Motorcycle.h>

#include <utility>

namespace CGAL {

namespace Polyline_tracing {

template<typename MotorcycleGraph>
class Motorcycle_priority_queue_entry
{
  typedef Motorcycle_priority_queue_entry<MotorcycleGraph>        Self;

public:
  typedef typename MotorcycleGraph::Geom_traits                   Geom_traits;

  typedef typename Geom_traits::FT                                FT;
  typedef typename MotorcycleGraph::Motorcycle                    Motorcycle;
  typedef Motorcycle*                                             Motorcycle_ptr;

  FT time_at_closest_target() const { return mc->time_at_closest_target(); }
  FT is_initialized() const { return mc->is_initialized(); }

  Motorcycle& motorcycle() { return *mc; }
  const Motorcycle& motorcycle() const { return *mc; }

  Motorcycle_priority_queue_entry(Motorcycle_ptr mc);

  friend bool operator<(const Self& lhs, const Self& rhs)
  {
    // If the times are equal, give priority to uninitialized motorcycles
    if(lhs.time_at_closest_target() == rhs.time_at_closest_target())
      return (!lhs.is_initialized());

    // '>' because we want the priority queue to output the element with smallest time
    return lhs.time_at_closest_target() > rhs.time_at_closest_target();
  }

private:
  Motorcycle_ptr mc;
};

template<typename MotorcycleGraph>
Motorcycle_priority_queue_entry<MotorcycleGraph>::
Motorcycle_priority_queue_entry(Motorcycle_ptr mc)
  : mc(mc)
{ }

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_MOTORCYCLE_PRIORITY_QUEUE_ENTRY_H
