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

#include <fstream>
#include <iostream>
#include <list>
#include <set>
#include <sstream>

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
  typedef Motorcycle<K>                                         Self;

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

  const DEC_it& source() const { return sour; }
  DEC_it& destination() { return dest; }
  const DEC_it destination() const { return dest; }
  const FT speed() const { return v; }
  FT& current_time() { return time; }
  const FT& current_time() const { return time; }
  DEC_it& position() { return conf; }
  const DEC_it position() const { return conf; }

  Target_point_container& targets() { return target_points; }
  const Target_point_container& targets() const { return target_points; }
  std::list<DEC_it>& track() { return track_points; }
  const std::list<DEC_it>& track() const { return track_points; }

  // Constructor
  Motorcycle(const int id, const DEC_it p, const DEC_it v, FT speed,
             const FT dist_at_s, const FT dist_at_d);

  // Functions
  void add_new_target(const DEC_it target_point, const FT time_at_target);
  const DEC_it closest_target() const;
  FT time_at_closest_target() const;
  void erase_closest_target();
  bool has_reached_blocked_point() const;
  bool has_reached_simultaneous_collision_point() const;
  bool is_motorcycle_destination_final() const;
  void set_new_destination(const DEC_it new_dest, const FT dist);

  // output
  friend std::ostream& operator<<(std::ostream& out, const Self& mc) {
    out << "Motorcycle: " << mc.id() << " (crashed? " << mc.is_crashed() << ") "
        << "going from source: (" << mc.source()->point() << ")"
        << " to destination: (" << mc.destination()->point() << ")" << std::endl
        << "  currently at position: (" << mc.position()->point() << ")" << std::endl
        << "  with targets: " << std::endl;
    typename Target_point_container::const_iterator tpc_it = mc.targets().begin();
    typename Target_point_container::const_iterator end = mc.targets().end();
    for(; tpc_it!=end; ++tpc_it)
      out << "\t Point: (" << tpc_it->first->point() << ") time: " << tpc_it->second << std::endl;

    return out;
  }

  void output_intended_track() const;
  void output_track() const;

private:
  const int i;
  bool crashed;

  const DEC_it sour; // source point
  DEC_it dest; // destination, which might change
  const FT v; // speed of the motorcycle, 'const' for now
  FT time; // current time at the motorcycle
  DEC_it conf; // last confirmed position

  Target_point_container target_points; // tentative targets (ordered)

  std::list<DEC_it> track_points;
};

template<typename K>
Motorcycle<K>::
Motorcycle(const int id,
           const DEC_it source, const DEC_it destination, const FT speed,
           const FT dist_at_s, const FT dist_at_d)
  : i(id), crashed(false),
    sour(source), dest(destination), v(speed), time(dist_at_s), conf(source),
    target_points(Target_point_set_comparer<K>()), track_points()
{
  CGAL_assertion(dist_at_s <= dist_at_d);

  if(sour->point() == dest->point()) {
    std::cerr << "Warning: source and destination (" << sour->point() << ") "
              << "are equal for motorcycle: " << id << std::endl;
  }

  target_points.insert(std::make_pair(sour, dist_at_s));
  target_points.insert(std::make_pair(dest, dist_at_d));
}

template<typename K>
void
Motorcycle<K>::
add_new_target(const DEC_it target_point, const FT time_at_target)
{
  target_points.insert(std::make_pair(target_point, time_at_target));
}

template<typename K>
const typename Motorcycle<K>::DEC_it
Motorcycle<K>::
closest_target() const
{
  CGAL_precondition(!target_points.empty());
  return target_points.begin()->first;
}

template<typename K>
typename Motorcycle<K>::FT
Motorcycle<K>::
time_at_closest_target() const
{
  CGAL_precondition(!target_points.empty());
  return target_points.begin()->second;
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
  return conf->is_blocked();
}

template<typename K>
bool
Motorcycle<K>::
has_reached_simultaneous_collision_point() const
{
  return conf->has_simultaneous_collision();
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
set_new_destination(const DEC_it new_dest, const FT new_time)
{
  assert(target_points.empty());

  dest = new_dest;

  // Putting 'conf' again in the queue means that this will be immediately treated
  // as next point in the queue, allowing to then insert points if necessary
  // between 'conf' and 'new_dest'
  target_points.insert(std::make_pair(conf, time));
  target_points.insert(std::make_pair(dest, new_time));
}

template<typename K>
void
Motorcycle<K>::
output_track() const
{
  std::ostringstream out_filename;
  out_filename << "out_motorcycle_track_" << i << ".off" << std::ends;
  std::ofstream os(out_filename.str().c_str());

  const std::size_t pn = track_points.size();
  const std::size_t fn = pn - 1;

  os << "OFF" << '\n';
  os << pn << " " << fn << " 0" << '\n';

  typename std::list<DEC_it>::const_iterator tit = track_points.begin();
  typename std::list<DEC_it>::const_iterator end = track_points.end();
  for(; tit!=end; ++tit)
    os << (*tit)->point() << " 0" << '\n'; // the '0' is because OFF is a 3D format

  for(std::size_t j=0; j<fn; ++j)
    os << "3 " << j << " " << j+1 << " " << j << '\n';
}

template<typename K>
void
Motorcycle<K>::
output_intended_track() const
{
  // must be adapted to surface @todo

  std::ostringstream out_filename;
  out_filename << "out_motorcycle_intended_track_" << i << ".off" << std::ends;
  std::ofstream os(out_filename.str().c_str());

  os << "OFF" << '\n';
  os << "2 1 0" << '\n';
  os << sour->point() << " 0" << '\n';
  os << dest->point() << " 0" << '\n';
  os << "3 0 1 0" << std::endl;
}

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_MOTORCYCLE_H
