// Copyright (c) 2017, 2018 GeometryFactory (France).
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

#ifndef CGAL_POLYLINE_TRACING_TRACK_H
#define CGAL_POLYLINE_TRACING_TRACK_H

#include <CGAL/assertions.h>

#include <boost/graph/graph_traits.hpp>

#include <cstddef>
#include <fstream>
#include <iostream>
#include <list>

namespace CGAL {

namespace Polyline_tracing {

template<typename MotorcycleGraphTraits>
class Motorcycle_graph_node_dictionary;

template<typename MotorcycleGraphTraits>
class Motorcycle_track_segment
{
  typedef Motorcycle_track_segment<MotorcycleGraphTraits>                 Self;

public:
  typedef MotorcycleGraphTraits                                           Geom_traits;

  typedef typename Geom_traits::FT                                        FT;

  typedef std::size_t                                                     Motorcycle_ID;
  typedef Motorcycle_graph_node_dictionary<Geom_traits>                   Nodes;
  typedef typename Nodes::Node_ptr                                        Node_ptr;

  typedef typename Geom_traits::Halfedge_graph                            Halfedge_graph;
  typedef typename Geom_traits::Triangle_mesh                             Triangle_mesh;

  typedef typename boost::graph_traits<Halfedge_graph>::edge_descriptor   hg_edge_descriptor;
  typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor    face_descriptor;

  // Constructor
  Motorcycle_track_segment() { }
  Motorcycle_track_segment(const Motorcycle_ID mc_id,
                           const Node_ptr source,
                           const FT source_time,
                           const Node_ptr target,
                           const FT target_time)
    : mc_id(mc_id),
      source_node(source), target_node(target),
      source_time(source_time), target_time(target_time)
  {
    // some sanity checks
    CGAL_precondition(source_time <= target_time);
    CGAL_precondition(source_time != target_time || source == target);
    CGAL_precondition(source != target || source_time == target_time);
    CGAL_precondition(source_node->face() == target_node->face());
  }

  // Access
  Motorcycle_ID& motorcycle_id() { return mc_id; }
  const Motorcycle_ID& motorcycle_id() const { return mc_id; }
  Node_ptr& source() { return source_node; }
  const Node_ptr& source() const { return source_node; }
  Node_ptr& target() { return target_node; }
  const Node_ptr& target() const { return target_node; }
  FT& time_at_source() { return source_time; }
  const FT& time_at_source() const { return source_time; }
  FT& time_at_target() { return target_time; }
  const FT& time_at_target() const { return target_time; }

  hg_edge_descriptor& graph_edge() { return ed; }
  const hg_edge_descriptor& graph_edge() const { return ed; }

  bool is_degenerate() const { return source_node == target_node; }
  face_descriptor face() const { return source_node->face(); }

  bool operator==(const Self& ts) const
  {
    return (mc_id == ts.motorcycle_id() &&
            *source_node == *(ts.source()) &&
            *target_node == *(ts.target()) &&
            source_time == ts.time_at_source() &&
            target_time == ts.time_at_target());
  }

private:
  Motorcycle_ID mc_id;
  Node_ptr source_node, target_node;
  FT source_time, target_time;

  hg_edge_descriptor ed;
};

template<typename MotorcycleGraphTraits>
class Motorcycle_track
{
  typedef Motorcycle_track<MotorcycleGraphTraits>           Self;

public:
  typedef MotorcycleGraphTraits                             Geom_traits;

  typedef Motorcycle_track_segment<Geom_traits>             Track_segment;

  // This container must be stable because we store iterators to its elements
  typedef std::list<Track_segment>                          Track_segment_container;

  typedef typename Geom_traits::FT                          FT;

  typedef std::size_t                                       Motorcycle_ID;
  typedef Motorcycle_graph_node_dictionary<Geom_traits>     Nodes;
  typedef typename Nodes::Node_ptr                          Node_ptr;

  typedef typename Track_segment_container::iterator        iterator;
  typedef typename Track_segment_container::const_iterator  const_iterator;
  typedef typename Track_segment_container::reference       reference;
  typedef typename Track_segment_container::const_reference const_reference;

  Track_segment_container& container() { return tsc_; }
  const Track_segment_container& container() const { return tsc_; }

  inline iterator begin() { return tsc_.begin(); }
  inline const_iterator begin() const { return tsc_.begin(); }
  inline iterator end() { return tsc_.end(); }
  inline const_iterator end() const { return tsc_.end(); }

  inline void push_back(const Track_segment& ts) { return tsc_.push_back(ts); }
  inline reference front() { return tsc_.front(); }
  inline const_reference front() const { return tsc_.front(); }
  inline reference back() { return tsc_.back(); }
  inline const_reference back() const { return tsc_.back(); }

  inline std::size_t size() const { return tsc_.size(); }

  iterator split_track_segment(iterator it, const Node_ptr new_point, const FT time_at_new_point)
  {
    CGAL_precondition(it->time_at_source() <= time_at_new_point);
    CGAL_precondition(time_at_new_point <= it->time_at_target());

    const Node_ptr source = it->source();
    const FT time_at_source = it->time_at_source();

    Track_segment new_seg(it->motorcycle_id(), source, time_at_source, new_point, time_at_new_point);
    iterator new_it = tsc_.insert(it, new_seg); // std::list's 'insert()' insert before the iterator

    it->source() = new_point;
    it->time_at_source() = time_at_new_point;

    return new_it;
  }

  bool write_off(std::ofstream& os) const
  {
    const std::size_t pn = tsc_.size();
    const std::size_t fn = (pn == 0) ? 0 : pn - 1;

    os << "OFF" << '\n';
    os << pn << " " << fn << " 0" << '\n';

    if(pn == 0)
      return os.good();

    const_iterator tscit = begin();
    const_iterator tsend = end();
    for(; tscit!=tsend; ++tscit)
    {
      os << tscit->source()->point();

      if(Geom_traits::dimension() == 2) // The '.off' format expects 3D points
        os << " 0";
      os << '\n';
    }

    tsc_.back().target_point()->point(); // back() exists because pn > 0

    for(std::size_t j=0; j<fn; ++j)
      os << "3 " << j << " " << j+1 << " " << j << '\n';

    return os.good();
  }

  friend std::ostream& operator<<(std::ostream& out, const Self& track)
  {
    out << "<~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>" << std::endl;

    const_iterator tscit = track.begin();
    const_iterator tsend = track.end();
    for(; tscit!=tsend; ++tscit)
      out << "[" << tscit->source()->point() << " ----- " << tscit->target()->point() << "]" << std::endl;

    return out;
  }

private:
  Track_segment_container tsc_;
};

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_TRACK_H
