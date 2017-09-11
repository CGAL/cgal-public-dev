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

#include <boost/optional.hpp>
#include <boost/parameter.hpp>

#include <fstream>
#include <iostream>
#include <list>
#include <set>
#include <sstream>

namespace CGAL {

namespace parameters {

  BOOST_PARAMETER_NAME( (source, tag) source_ )
  BOOST_PARAMETER_NAME( (destination, tag) destination_ )
  BOOST_PARAMETER_NAME( (speed, tag) speed_ )
  BOOST_PARAMETER_NAME( (direction, tag) direction_ )
  BOOST_PARAMETER_NAME( (initial_time, tag) initial_time_ )

} // end namespace parameters

namespace Polyline_tracing {

template<typename K, typename PolygonMesh>
struct Target_point_set_comparer
{
  typedef typename K::FT                                    FT;
  typedef typename Dictionary<K, PolygonMesh>::DEC_it       DEC_it;

  bool operator()(const std::pair<DEC_it, FT>& lhs,
                  const std::pair<DEC_it, FT>& rhs) const {
    return lhs.second < rhs.second;
  }
};

// -----------------------------------------------------------------------------

template<typename K, typename PolygonMesh>
class Motorcycle_impl
{
  typedef Motorcycle_impl<K, PolygonMesh>                       Self;

public:
  typedef typename K::FT                                        FT;
  typedef typename K::Point_2                                   Point;
  typedef typename K::Vector_2                                  Vector;

  typedef Dictionary<K, PolygonMesh>                            Dictionary;
  typedef Dictionary_entry<K, PolygonMesh>                      Dictionary_entry;
  typedef typename Dictionary::DEC_it                           DEC_it;
  typedef typename Dictionary_entry::Face_location              Face_location;

  typedef std::pair<DEC_it, FT>                                 Target_point;
  typedef std::set<Target_point,
                   Target_point_set_comparer<K, PolygonMesh> >  Target_point_container;

  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor    face_descriptor;

  // Access
  std::size_t id() const { return i; }
  void set_id(std::size_t id) { i = id; }
  bool is_crashed() const { return crashed; }
  void crash() { crashed = true; }

  const Point& initial_source_point() const { return ini_sour_pt; }
  boost::optional<Point>& initial_destination_point() { return ini_dest_pt; }
  const boost::optional<Point>& initial_destination_point() const { return ini_dest_pt; }

  DEC_it& source() { return sour; }
  const DEC_it& source() const { return sour; }
  DEC_it& destination() { return dest; }
  const DEC_it destination() const { return dest; }

  void set_location(const Face_location& loc) { conf->set_location(loc); }
  const Face_location& current_location() const { return conf->location(); }

  const FT speed() const { return spee; }
  boost::optional<Vector>& direction() { return dir; }
  const boost::optional<Vector>& direction() const { return dir; }
  FT& current_time() { return time; }
  const FT& current_time() const { return time; }
  DEC_it& position() { return conf; }
  const DEC_it position() const { return conf; }

  Target_point_container& targets() { return target_points; }
  const Target_point_container& targets() const { return target_points; }
  std::list<DEC_it>& track() { return track_points; }
  const std::list<DEC_it>& track() const { return track_points; }

  // Constructor
  template<typename ArgumentPack>
  Motorcycle_impl(const ArgumentPack& args);

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
  std::size_t i;
  bool crashed;

  // the very first source and destination points, before insertion in the dictionary
  const Point ini_sour_pt;
  boost::optional<Point> ini_dest_pt;

  // below might change when we change faces
  DEC_it sour; // source point
  DEC_it dest; // destination
  DEC_it conf; // last confirmed position

  const FT spee; // speed of the motorcycle, 'const' for now
  boost::optional<Vector> dir; // direction of the motorcycle
  FT time; // current time at the position of the motorcycle

  face_descriptor face; // current face

  // tentative targets (ordered by increasing distance to 'conf')
  Target_point_container target_points;

  std::list<DEC_it> track_points;
};

template<typename K, typename PolygonMesh>
class Motorcycle
  : public Motorcycle_impl<K, PolygonMesh>
{
  typedef Motorcycle_impl<K, PolygonMesh>               Base;

public:
  BOOST_PARAMETER_CONSTRUCTOR(Motorcycle, (Base), parameters::tag,
                              (required (source_, *))
                              (optional (destination_, *)
                                        (speed_, *)
                                        (direction_, *)
                                        (initial_time_, *))
                             )
};

// -----------------------------------------------------------------------------

template<typename K, typename PolygonMesh>
template <class ArgumentPack>
Motorcycle_impl<K, PolygonMesh>::
Motorcycle_impl(const ArgumentPack& args)
  :
    i(-1),
    crashed(false),
    ini_sour_pt(args[parameters::source]),
    ini_dest_pt(args[parameters::destination|boost::none]),
    sour(), dest(), conf(),
    spee(args[parameters::speed|1.]),
    dir(args[parameters::direction|boost::none]),
    time(args[parameters::initial_time|0.]),
    target_points(Target_point_set_comparer<K, PolygonMesh>()),
    track_points()
{
  CGAL_precondition(speed() > 0.);
  CGAL_precondition(ini_dest_pt || dir);

  if(ini_dest_pt != boost::none && ini_sour_pt == *ini_dest_pt)
  {
    std::cerr << "Warning: Creation of a motorcycle with identical source "
              << "and destination: (" << ini_sour_pt << ") " << std::endl;
  }

  if(dir != boost::none && *dir == CGAL::NULL_VECTOR)
  {
    std::cerr << "Warning: Creation of a motorcycle with null direction" << std::endl;
  }
}

template<typename K, typename PolygonMesh>
void
Motorcycle_impl<K, PolygonMesh>::
add_new_target(const DEC_it target_point, const FT time_at_target)
{
  target_points.insert(std::make_pair(target_point, time_at_target));
}

template<typename K, typename PolygonMesh>
const typename Motorcycle_impl<K, PolygonMesh>::DEC_it
Motorcycle_impl<K, PolygonMesh>::
closest_target() const
{
  CGAL_precondition(!target_points.empty());
  return target_points.begin()->first;
}

template<typename K, typename PolygonMesh>
typename Motorcycle_impl<K, PolygonMesh>::FT
Motorcycle_impl<K, PolygonMesh>::
time_at_closest_target() const
{
  CGAL_precondition(!target_points.empty());
  return target_points.begin()->second;
}

template<typename K, typename PolygonMesh>
void
Motorcycle_impl<K, PolygonMesh>::
erase_closest_target()
{
  CGAL_assertion(!target_points.empty());
  return target_points.erase(target_points.begin());
}

template<typename K, typename PolygonMesh>
bool
Motorcycle_impl<K, PolygonMesh>::
has_reached_blocked_point() const
{
  return conf->is_blocked();
}

template<typename K, typename PolygonMesh>
bool
Motorcycle_impl<K, PolygonMesh>::
has_reached_simultaneous_collision_point() const
{
  return conf->has_simultaneous_collision();
}

template<typename K, typename PolygonMesh>
bool
Motorcycle_impl<K, PolygonMesh>::
is_motorcycle_destination_final() const
{
  return true; // @todo
}

template<typename K, typename PolygonMesh>
void
Motorcycle_impl<K, PolygonMesh>::
set_new_destination(const DEC_it new_dest, const FT new_time)
{
  CGAL_assertion(target_points.empty());

  dest = new_dest;

  // Putting 'conf' again in the queue means that this will be immediately treated
  // as next point in the queue, allowing to then insert points if necessary
  // between 'conf' and 'new_dest'
  target_points.insert(std::make_pair(conf, time));
  target_points.insert(std::make_pair(dest, new_time));
}

template<typename K, typename PolygonMesh>
void
Motorcycle_impl<K, PolygonMesh>::
output_track() const
{
  std::ostringstream out_filename;
  out_filename << "out_motorcycle_track_" << i << ".off" << std::ends;
  std::ofstream os(out_filename.str().c_str());

  const std::size_t pn = track_points.size();
  CGAL_assertion(pn != 0);
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

template<typename K, typename PolygonMesh>
void
Motorcycle_impl<K, PolygonMesh>::
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
