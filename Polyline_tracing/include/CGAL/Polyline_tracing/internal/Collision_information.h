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

#include <CGAL/assertions.h>

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
  typedef MotorcycleGraph                                      Motorcycle_graph;
  typedef typename Motorcycle_graph::Triangle_mesh             Triangle_mesh;
  typedef typename Motorcycle_graph::FT                        FT;
  typedef typename Motorcycle_graph::Node_ptr                  Node_ptr;
  typedef typename Motorcycle_graph::face_descriptor           face_descriptor;
  typedef typename Motorcycle_graph::Face_location             Face_location;
  typedef typename Motorcycle_graph::Barycentric_coordinates   Barycentric_coordinates;

  typedef typename Motorcycle_graph::Track_segment_ptr         Track_segment_ptr;

  typedef boost::variant<Node_ptr, Face_location>              Node_ptr_or_Face_location;
  typedef boost::optional<Node_ptr_or_Face_location>           Collision;

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
  };

  friend inline Collision_return operator|(Collision_return cr1, Collision_return cr2) {
    return static_cast<Collision_return>(static_cast<int>(cr1) | static_cast<int>(cr2));
  }

  // Foreign collision seen from the foreign motorcycle.
  // There is usually only one foreign motorcycle intersection 'mc' at a given time,
  // but in the case of multiple motorcycles, it is useful to add them all for robustness
  struct Foreign_collision_information
  {
    const std::size_t fmc_id;
    const FT foreign_time_at_closest_collision;
    const bool must_crash;
    Track_segment_ptr foreign_track_ptr; // if it's not a track, then it'll be the tentative track

    Foreign_collision_information(const std::size_t fmc_id, const FT foreign_time_at_collision,
                                  const bool must_foreign_motorcycle_crash = false)
      :
        fmc_id(fmc_id),
        foreign_time_at_closest_collision(foreign_time_at_collision),
        must_crash(must_foreign_motorcycle_crash),
        foreign_track_ptr()
    { }
  };

  typedef boost::container::slist<Foreign_collision_information> Foreign_collisions_container;

  // Constructor
  Collision_information(const Motorcycle_graph& mg, const FT max_time_at_collision)
    :
      mg(mg),
      maximum_time_at_collision(max_time_at_collision),
      must_crash(false),
      closest_collision(boost::none),
      time_at_closest_collision(std::numeric_limits<FT>::max()),
      foreign_collisions()
  { }

  void reset()
  {
    // information relative to the current face
    must_crash = false;
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

  void set_new_collision(const FT time_at_collision, const Node_ptr_or_Face_location& collision,
                         const bool must_motorcycle_crash = false)
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "Best new collision at time : " << time_at_collision << std::endl;
    std::cout << "Previous best time was: " << time_at_closest_collision << std::endl;
#endif

    reset();

    must_crash = must_motorcycle_crash;
    time_at_closest_collision = time_at_collision;
    closest_collision = collision;
  }

  Collision_return add_foreign_collision(const std::size_t fmc_id, const FT foreign_time_at_collision,
                                         const bool must_motorcycle_crash,
                                         const bool must_foreign_motorcycle_crash)
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "trying to add foreign collision with " << fmc_id << " at " << foreign_time_at_collision << std::endl;
#endif

#if 1 // @tmp
    // If multiple foreign motorcycles intersect 'mc' at the same point,
    // keep a foreign motorcycle that does not know about this intersection

    CGAL_assertion(foreign_collisions.size() == 0 || foreign_collisions.size() == 1);
    if(foreign_collisions.size() == 1)
    {
      // Due to snapping and other shenanigans, it should be enough to test if there is
      // a target at that time (and not furthermore test that it is the same collision point)
      const bool current_best_foreign_mc_knows_about_intersection =
        (mg.motorcycle(foreign_collisions.front().fmc_id).has_target_at_time(foreign_time_at_collision)).second;

      if(current_best_foreign_mc_knows_about_intersection)
      {
        const bool new_tentative_foreign_mc_knows_about_intersection =
          (mg.motorcycle(fmc_id).has_target_at_time(foreign_time_at_collision)).second;

        if(!new_tentative_foreign_mc_knows_about_intersection)
        {
          std::cout << "replaced collision with " << foreign_collisions.front().fmc_id
                    << " by collision with " << fmc_id << std::endl;
          foreign_collisions.pop_front();
          foreign_collisions.push_front(Foreign_collision_information(fmc_id, foreign_time_at_collision,
                                                                      must_foreign_motorcycle_crash));
          must_crash = must_motorcycle_crash;
          return COLLISION;
        }
      }

      return NO_COLLISION;
    }
    else
#endif
    {
      must_crash = must_crash || must_motorcycle_crash;
      foreign_collisions.push_front(Foreign_collision_information(fmc_id, foreign_time_at_collision,
                                                                  must_foreign_motorcycle_crash));
      return COLLISION;
    }
  }

  Collision_return treat_potential_collision(const Node_ptr_or_Face_location& collision,
                                             const FT time_at_collision,
                                             const std::size_t fmc_id,
                                             const FT foreign_time_at_collision,
                                             const bool crash_motorcycle = false,
                                             const bool crash_foreign_motorcycle = false)
  {
    Collision_time_comparison_result ctcr = compare_collision_time_to_closest(time_at_collision);

    if(ctcr != LATER_THAN_CURRENT_CLOSEST_TIME)
    {
      if(ctcr == NEW_CLOSEST_TIME)
        set_new_collision(time_at_collision, collision, crash_motorcycle);

      return add_foreign_collision(fmc_id, foreign_time_at_collision, crash_motorcycle, crash_foreign_motorcycle);
    }

    return NO_COLLISION;
  }

public:
  const Motorcycle_graph& mg;
  const FT maximum_time_at_collision;

  bool must_crash;
  Collision closest_collision;
  FT time_at_closest_collision;

  Foreign_collisions_container foreign_collisions;
};

} // namespace internal

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_COLLISION_INFORMATION_H
