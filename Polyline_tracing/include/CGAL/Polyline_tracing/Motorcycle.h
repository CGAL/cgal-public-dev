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

#ifndef CGAL_POLYLINE_TRACING_MOTORCYCLE_H
#define CGAL_POLYLINE_TRACING_MOTORCYCLE_H

#include <CGAL/Polyline_tracing/Dictionary.h>

#include <list>
#include <set>

namespace CGAL {

namespace Polyline_tracing {

template<typename K>
struct Target_point_set_comparer
{
  typedef typename K::FT                           FT;
  typedef typename Dictionary<K>::DEC_it           DEC_it;

  bool operator()(const std::pair<DEC_it, FT>& lhs,
                  const std::pair<DEC_it, FT>& rhs) const {
    return lhs.second < rhs.second;
  }
};

template<typename K>
class Motorcycle
{
public:
  typedef typename K::FT                                        FT;
  typedef typename K::Point_2                                   Point;
  typedef typename K::Vector_2                                  Vector;

  typedef Dictionary<K>                                         Dictionary;
  typedef Dictionary_entry<K>                                   Dictionary_entry;
  typedef typename Dictionary::DEC_it                           DEC_it;

  typedef std::pair<DEC_it, FT>                                 Target_point;
  typedef std::set<Target_point, Target_point_set_comparer<K> > Target_point_container;

  // Access
  int id() const { return i; }
  bool is_crashed() const { return crashed; }
  void crash() { crashed = true; }

  FT& distance() { return dist; }
  const FT& distance() const { return dist; }
  DEC_it& position() { return conf; }
  const DEC_it position() const { return conf; }
  DEC_it& destination() { return dest; }
  const DEC_it destination() const { return dest; }

  Target_point_container& targets() { return target_points; }
  const Target_point_container& targets() const { return target_points; }
  std::list<DEC_it>& path() { return path_points; }
  const std::list<DEC_it>& path() const { return path_points; }

  // Constructor
  Motorcycle(const int id, const DEC_it p, const DEC_it v, const FT dist_at_d);

  // Functions
  void add_new_target(const DEC_it target_point, const FT distance_at_target);
  const DEC_it closest_target() const;
  FT distance_at_closest_target() const;
  void erase_closest_target();
  bool has_reached_blocked_point() const;
  bool has_reached_simultaneous_collision_point() const;
  bool is_motorcycle_destination_final() const;
  void set_new_destination(const DEC_it new_dest, const FT dist);

  void output_path() const;

private:
  const int i;
  bool crashed;

  const DEC_it sour; // source point
  DEC_it dest; // destination, which might change
  FT dist; // distance travelled
  DEC_it conf; // last confirmed position

  Target_point_container target_points; // tentative targets (ordered)

  std::list<DEC_it> path_points;
};

template<typename K>
Motorcycle<K>::
Motorcycle(const int id,
           const DEC_it s, const DEC_it d,
           const FT dist_at_d)
  : i(id), crashed(false),
    sour(s), dest(d), dist(0.), conf(s),
    target_points(Target_point_set_comparer<K>()), path_points()
{
  target_points.insert(std::make_pair(sour, 0.));
  target_points.insert(std::make_pair(dest, dist_at_d));
}

template<typename K>
void
Motorcycle<K>::
add_new_target(const DEC_it target_point, const FT distance_at_target)
{
  target_points.insert(std::make_pair(target_point, distance_at_target));
}

template<typename K>
const typename Motorcycle<K>::DEC_it
Motorcycle<K>::
closest_target() const
{
  CGAL_assertion(!target_points.empty());
  return target_points.begin()->first;
}

template<typename K>
typename Motorcycle<K>::FT
Motorcycle<K>::
distance_at_closest_target() const
{
  assert(!target_points.empty());
  return closest_target()->second;
}

template<typename K>
void
Motorcycle<K>::
erase_closest_target()
{
  CGAL_assertion(!target_points.empty());
  return target_points.erase(target_points.begin());
}

template<typename K>
bool
Motorcycle<K>::
has_reached_blocked_point() const
{
  return conf->second.is_blocked();
}

template<typename K>
bool
Motorcycle<K>::
has_reached_simultaneous_collision_point() const
{
  return conf->second.has_simultaneous_collision();
}

template<typename K>
bool
Motorcycle<K>::
is_motorcycle_destination_final() const
{
  return true; // @todo
}

template<typename K>
void
Motorcycle<K>::
set_new_destination(const DEC_it new_dest, const FT new_dist)
{
  assert(target_points.empty());

  dest = new_dest;

  // Putting 'conf' again in the queue means that this will be immediately treated
  // as next point in the queue, allowing to then insert points if necessary
  // between 'conf' and 'new_dest'
  target_points.insert(std::make_pair(conf, dist));
  target_points.insert(std::make_pair(dest, new_dist));
}


template<typename K>
void
Motorcycle<K>::
output_path() const
{

}

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_MOTORCYCLE_H
