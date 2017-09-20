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
#include <CGAL/Polyline_tracing/Uniform_direction_tracer_visitor.h>
#include <CGAL/Polyline_tracing/Tracer.h>

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

template<typename K, typename TriangleMesh>
struct Target_point_set_comparer
{
  typedef typename K::FT                                    FT;
  typedef typename Dictionary<K, TriangleMesh>::DEC_it      DEC_it;

  bool operator()(const std::pair<DEC_it, FT>& lhs,
                  const std::pair<DEC_it, FT>& rhs) const {
    // don't want to insert the same point multiple times
    CGAL_assertion(lhs.first != rhs.first);
    return lhs.second < rhs.second;
  }
};

// -----------------------------------------------------------------------------

// Having a "_impl" class is only done because it's needed for BOOST_PARAMETER_CONSTRUCTOR
template<typename K, typename TriangleMesh>
class Motorcycle_impl
{
  typedef Motorcycle_impl<K, TriangleMesh>                      Self;
  typedef Motorcycle<K, TriangleMesh>                           Derived;

public:
  typedef typename K::FT                                        FT;
  typedef typename K::Point_2                                   Point;
  typedef typename K::Vector_2                                  Vector;

  typedef Dictionary<K, TriangleMesh>                           Dictionary;
  typedef Dictionary_entry<K, TriangleMesh>                     Dictionary_entry;
  typedef typename Dictionary::DEC_it                           DEC_it;
  typedef typename Dictionary_entry::Face_location              Face_location;

  typedef std::pair<DEC_it, FT>                                 Target_point;
  typedef std::set<Target_point,
                   Target_point_set_comparer<K, TriangleMesh> > Target_point_container;

  typedef Uniform_direction_tracer_visitor<K, TriangleMesh>     Tracer_visitor;
  typedef Tracer<K, TriangleMesh, Tracer_visitor>               Tracer;

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

  const Face_location& current_location() const { return conf->location(); }

  void set_destination_finality(bool b) { is_dest_final = b; }
  const bool& is_destination_final() const { return is_dest_final; }

  const FT speed() const { return spee; }
  boost::optional<Vector>& direction() { return dir; }
  const boost::optional<Vector>& direction() const { return dir; }
  FT& current_time() { return time; }
  const FT& current_time() const { return time; }
  FT& time_at_source() { return time_at_sour; }
  const FT& time_at_source() const { return time_at_sour; }
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
  void add_target(const DEC_it target_point, const FT time_at_target);
  const DEC_it closest_target() const;

  boost::tuple<bool, DEC_it, DEC_it, FT, bool>
  compute_next_destination(Dictionary& points, const TriangleMesh& mesh);

  void erase_closest_target();
  bool has_reached_blocked_point() const;
  bool has_reached_simultaneous_collision_point() const;
  bool is_motorcycle_destination_final() const;
  FT time_at_closest_target() const;

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
  // ID and status
  std::size_t i;
  bool crashed;

  // The very first source and destination points, before insertion in the dictionary
  const Point ini_sour_pt;
  boost::optional<Point> ini_dest_pt;

  // Below might change when we move in the mesh
  DEC_it sour; // source
  DEC_it dest; // destination
  DEC_it conf; // current position (last confirmed position)

  // indicates whether we should stop at the destination or try to trace
  bool is_dest_final;

  const FT spee; // speed of the motorcycle, 'const' for now
  boost::optional<Vector> dir; // direction
  FT time; // time at the current position
  FT time_at_sour; // time at the current source

  // The tentative targets, ordered by increasing distance from 'conf'
  Target_point_container target_points;

  // Tracer (computes the next target when we reach the destination)
  Tracer tracer;

  std::list<DEC_it> track_points;
};

template<typename K, typename TriangleMesh>
class Motorcycle
  : public Motorcycle_impl<K, TriangleMesh>
{
  typedef Motorcycle_impl<K, TriangleMesh>               Base;

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

template<typename K, typename TriangleMesh>
template <class ArgumentPack>
Motorcycle_impl<K, TriangleMesh>::
Motorcycle_impl(const ArgumentPack& args)
  :
    i(-1),
    crashed(false),
    ini_sour_pt(args[parameters::source]),
    ini_dest_pt(args[parameters::destination|boost::none]),
    sour(), dest(), conf(),
    is_dest_final(false),
    spee(args[parameters::speed|1.]),
    dir(args[parameters::direction|boost::none]),
    time(args[parameters::initial_time|0.]),
    time_at_sour(time),
    target_points(Target_point_set_comparer<K, TriangleMesh>()),
    tracer(),
    track_points()
{
  // Reject null speed
  CGAL_precondition(spee > 0.);

  // Either the destination or the direction should be provided
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

template<typename K, typename TriangleMesh>
void
Motorcycle_impl<K, TriangleMesh>::
add_target(const DEC_it target_point, const FT time_at_target)
{
  target_points.insert(std::make_pair(target_point, time_at_target));
}

template<typename K, typename TriangleMesh>
const typename Motorcycle_impl<K, TriangleMesh>::DEC_it
Motorcycle_impl<K, TriangleMesh>::
closest_target() const
{
  CGAL_precondition(!target_points.empty());
  return target_points.begin()->first;
}

template<typename K, typename TriangleMesh>
boost::tuple<bool, // successfuly computed a next path or not
             typename Motorcycle_impl<K, TriangleMesh>::DEC_it, // next source
             typename Motorcycle_impl<K, TriangleMesh>::DEC_it, // next destination
             typename K::FT, // time at next destination
             bool> // whether the destination is final or not
Motorcycle_impl<K, TriangleMesh>::
compute_next_destination(Dictionary& points, const TriangleMesh& mesh)
{
  CGAL_precondition(target_points.empty());

  // that derived cast is so that tracer visitor can indeed take a Motorcycle
  // and not a motorcycle_impl. It's safe since we only deal manipulate "full"
  // motorcycles, but it's kinda ugly @fixme
  return tracer.trace(static_cast<Derived&>(*this), points, mesh);
}

template<typename K, typename TriangleMesh>
void
Motorcycle_impl<K, TriangleMesh>::
erase_closest_target()
{
  CGAL_assertion(!target_points.empty());
  return target_points.erase(target_points.begin());
}

template<typename K, typename TriangleMesh>
bool
Motorcycle_impl<K, TriangleMesh>::
has_reached_blocked_point() const
{
  return conf->is_blocked();
}

template<typename K, typename TriangleMesh>
bool
Motorcycle_impl<K, TriangleMesh>::
has_reached_simultaneous_collision_point() const
{
  return conf->has_simultaneous_collision();
}

template<typename K, typename TriangleMesh>
bool
Motorcycle_impl<K, TriangleMesh>::
is_motorcycle_destination_final() const
{
  return is_dest_final;
}

template<typename K, typename TriangleMesh>
typename Motorcycle_impl<K, TriangleMesh>::FT
Motorcycle_impl<K, TriangleMesh>::
time_at_closest_target() const
{
  CGAL_precondition(!target_points.empty());
  return target_points.begin()->second;
}

template<typename K, typename TriangleMesh>
void
Motorcycle_impl<K, TriangleMesh>::
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

template<typename K, typename TriangleMesh>
void
Motorcycle_impl<K, TriangleMesh>::
output_intended_track() const
{
  // must be adapted to multiple destinations @todo

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
