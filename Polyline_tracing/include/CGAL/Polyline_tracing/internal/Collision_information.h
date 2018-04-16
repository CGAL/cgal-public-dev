// Copyright (c) 2018 GeometryFactory (France).
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

#ifndef CGAL_POLYLINE_TRACING_COLLISION_INFORMATION_H
#define CGAL_POLYLINE_TRACING_COLLISION_INFORMATION_H

#include <boost/container/slist.hpp>
#include <boost/optional.hpp>
#include <boost/variant.hpp>

#include <limits>

namespace CGAL {

namespace Polyline_tracing {

namespace internal {

// This struct regroups all useful information on a potential intersection
template<typename MotorcycleGraph>
struct Collision_information
{
  typedef typename MotorcycleGraph::Triangle_mesh             Triangle_mesh;
  typedef typename MotorcycleGraph::FT                        FT;
  typedef typename MotorcycleGraph::Node_ptr                  Node_ptr;
  typedef typename MotorcycleGraph::face_descriptor           face_descriptor;
  typedef typename MotorcycleGraph::Face_location             Face_location;
  typedef typename MotorcycleGraph::Barycentric_coordinates   Barycentric_coordinates;

  typedef typename MotorcycleGraph::Track_segment_ptr         Track_segment_ptr;

  typedef boost::variant<Node_ptr, Face_location>             Node_ptr_or_Face_location;
  typedef boost::optional<Node_ptr_or_Face_location>          Collision;

  enum Collision_time_comparison_result
  {
    LATER_THAN_CURRENT_CLOSEST_TIME = 0,
    EQUAL_TO_CURRENT_CLOSEST_TIME,
    NEW_CLOSEST_TIME
  };

  // 0, 1, and 3 because we will use this enum in bitwise operations (see operator below)
  enum Collision_return
  {
    NO_COLLISION = 0,
    COLLISION = 1,
    SNAPPED_COLLISION_TO_EXISTING_POINT = 3
  };

  friend inline Collision_return operator|(Collision_return cr1, Collision_return cr2) {
    return static_cast<Collision_return>(static_cast<int>(cr1) | static_cast<int>(cr2));
  }

  // Foreign collision seen from the foreign motorcycle.
  // There is usually only one foreign motorcycle intersection 'mc' at a given time,
  // but in the case of multiple motorcycles, it is useful to add them all for robustness
  struct Foreign_collision_information
  {
    std::size_t fmc_id;
    FT foreign_time_at_closest_collision;
    Track_segment_ptr foreign_track_ptr; // if it's not a track, then it'll be the tentative track

    Foreign_collision_information(const std::size_t fmc_id, const FT foreign_time_at_collision)
      :
        fmc_id(fmc_id),
        foreign_time_at_closest_collision(foreign_time_at_collision),
        foreign_track_ptr()
    { }
  };

  typedef boost::container::slist<Foreign_collision_information> Foreign_collisions_container;

  // Constructor
  Collision_information(const FT max_time_at_collision)
    :
      maximum_time_at_collision(max_time_at_collision),
      closest_collision(boost::none),
      time_at_closest_collision(std::numeric_limits<FT>::max()),
      foreign_collisions()
  { }

  void reset()
  {
    // information relative to the current face
    closest_collision = boost::none;
    time_at_closest_collision = std::numeric_limits<FT>::max();
    foreign_collisions.clear();
  }

  // Functions
  bool found_collision() const { return (closest_collision != boost::none); }

  // Check if the times provided passed in arguments correspond to a collision
  // earlier than the current best, bounded by the maximum time (time at closest target)
  Collision_time_comparison_result compare_collision_time_to_closest(const FT time_at_collision) const
  {
    if(time_at_collision > maximum_time_at_collision)
      return LATER_THAN_CURRENT_CLOSEST_TIME;
    else if(time_at_collision == maximum_time_at_collision && !found_collision())
      return NEW_CLOSEST_TIME;
    else // (time_at_collision < maximum_time_at_collision) or (a collision with closest_time = max_time)
    {
      if(time_at_collision > time_at_closest_collision)
        return LATER_THAN_CURRENT_CLOSEST_TIME;
      else if(time_at_collision == time_at_closest_collision)
        return EQUAL_TO_CURRENT_CLOSEST_TIME;
      else // time_at_collision < time_at_closest_collision
        return NEW_CLOSEST_TIME;
    }
  }

  void set_new_collision(const FT time_at_collision, const Node_ptr_or_Face_location& collision)
  {
    reset();

    time_at_closest_collision = time_at_collision;
    closest_collision = collision;
  }

  void add_foreign_collision(const std::size_t fmc_id, const FT foreign_time_at_collision)
  {
    foreign_collisions.push_front(Foreign_collision_information(fmc_id, foreign_time_at_collision));
  }

  Collision_return treat_potential_collision(const FT time_at_collision,
                                             const Node_ptr_or_Face_location& collision,
                                             const std::size_t fmc_id,
                                             const FT foreign_time_at_collision)
  {
    Collision_time_comparison_result ctcr = compare_collision_time_to_closest(time_at_collision);

    if(ctcr != LATER_THAN_CURRENT_CLOSEST_TIME)
    {
      if(ctcr == NEW_CLOSEST_TIME)
        set_new_collision(time_at_collision, collision);

      add_foreign_collision(fmc_id, foreign_time_at_collision);
    }
  }

public:
  const FT maximum_time_at_collision;

  Collision closest_collision;
  FT time_at_closest_collision;

  Foreign_collisions_container foreign_collisions;
};

} // namespace internal

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_COLLISION_INFORMATION_H
