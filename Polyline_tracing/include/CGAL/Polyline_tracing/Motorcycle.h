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

  int id() const { return i; }
  void set_id(int id) { i = id; }
  bool is_crashed() const { return crashed; }
  void crash() { crashed = true; }

  const DEC_it position() const { return conf; }
  DEC_it& position() { return conf; }
  const DEC_it destination() const { return dest; }
  DEC_it& destination() { return dest; }
  const DEC_it closest_target() const {
    CGAL_assertion(!target_points.empty());
    return target_points.begin()->first;
  }

  Target_point_container targets() { return target_points; }
  const Target_point_container& targets() const { return target_points; }

  Motorcycle(const int id, const DEC_it p, const DEC_it v, const FT dist_at_d);

  void add_new_target(const Target_point& target);
  void remove_closest_target();

  bool has_reached_end();
  void trace();

private:
  int i;
  bool crashed;

  const DEC_it sour; // source point
  DEC_it dest; // destination, might change
  FT distance; // distance travelled
  DEC_it conf; // confirmed position
  DEC_it tent; // tentative target
  Target_point_container target_points;

  std::list<DEC_it> path;
};

template<typename K>
Motorcycle<K>::
Motorcycle(const int id,
           const DEC_it s, const DEC_it d,
           const FT dist_at_d)
  : i(id), crashed(false),
    sour(s), dest(d), distance(0.),
    conf(s), tent(s),
    target_points(Target_point_set_comparer<K>()), path()
{
  target_points.insert(std::make_pair(sour, 0.));
  target_points.insert(std::make_pair(dest, dist_at_d));
}

template<typename K>
void
Motorcycle<K>::
add_new_target(const Target_point& target)
{
  target_points.insert(target);
}

template<typename K>
void
Motorcycle<K>::
remove_closest_target()
{
  CGAL_assertion(!target_points.empty());
  return target_points.erase(target_points.begin());
}

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_MOTORCYCLE_H
