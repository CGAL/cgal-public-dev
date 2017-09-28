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

#include <CGAL/Polygon_mesh_processing/locate.h>

#include <boost/optional.hpp>
#include <boost/parameter.hpp>
#include <boost/variant.hpp>

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
  BOOST_PARAMETER_NAME( (tracer, tag) tracer_ )

} // end namespace parameters

namespace Polyline_tracing {

namespace internal {

template<typename MotorcycleGraphTraits>
struct Target_point_set_comparer
{
  typedef typename MotorcycleGraphTraits::FT                  FT;
  typedef typename Dictionary<MotorcycleGraphTraits>::DEC_it  DEC_it;

  bool operator()(const std::pair<DEC_it, FT>& lhs, const std::pair<DEC_it, FT>& rhs) const
  {
    // Don't want to insert the same point multiple times
    CGAL_assertion(lhs.first != rhs.first);
    return lhs.second < rhs.second;
  }
};

template<typename MotorcycleGraphTraits>
struct Track_comparer
{
  typedef typename MotorcycleGraphTraits::FT                  FT;
  typedef typename Dictionary<MotorcycleGraphTraits>::DEC_it  DEC_it;

  bool operator()(const std::pair<DEC_it, FT>& lhs, const std::pair<DEC_it, FT>& rhs) const
  {
    // Points at the same time from the source should be equal (or almost)
    if(lhs.second == rhs.second)
    {
      CGAL_assertion(
        CGAL::squared_distance(lhs.first->point(), rhs.first->point()) <
          std::numeric_limits<FT>::epsilon());
    }

    return lhs.second < rhs.second;
  }
};

} // namespace internal

// -----------------------------------------------------------------------------

template<typename MotorcycleGraphTraits>
class Motorcycle_impl_base
{
  typedef Motorcycle_impl_base<MotorcycleGraphTraits>         Self;

public:
  typedef MotorcycleGraphTraits                               Geom_traits;
  typedef typename Geom_traits::Triangle_mesh                 Triangle_mesh;

  typedef typename Geom_traits::FT                            FT;
  typedef typename Geom_traits::Point_d                       Point;
  typedef typename Geom_traits::Vector_d                      Vector;

  typedef Dictionary<Geom_traits>                             Dictionary;
  typedef Dictionary_entry<Geom_traits>                       Dictionary_entry;
  typedef typename Dictionary::DEC_it                         DEC_it;
  typedef typename Dictionary_entry::Face_location            Face_location;

  typedef std::pair<DEC_it, FT>                               Track_point;
  typedef std::set<Track_point, internal::Target_point_set_comparer<Geom_traits> >
                                                              Target_point_container;

  typedef std::multiset<Track_point, internal::Track_comparer<Geom_traits> >
                                                              Track;

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
  Track& track() { return track_points; }
  const Track& track() const { return track_points; }

  // Constructor
  virtual ~Motorcycle_impl_base() { }
  Motorcycle_impl_base(const Point& source, const boost::optional<Point>& destination,
                       const FT speed, const boost::optional<Vector>& direction,
                       const FT initial_time);

  // Functions
  void add_target(const DEC_it target_point, const FT time_at_target);
  const DEC_it closest_target() const;

  virtual boost::tuple<bool, DEC_it, DEC_it, FT, bool>
  compute_next_destination(Dictionary& points, const Triangle_mesh& mesh) = 0;

  void remove_closest_target_from_targets();
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
      out << "\t " << &*(tpc_it->first) << " Point: (" << tpc_it->first->point() << ") time: " << tpc_it->second << std::endl;

    return out;
  }

  void output_intended_track() const;
  void output_track() const;

protected:
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

  Track track_points;
};

template<typename MotorcycleGraphTraits>
Motorcycle_impl_base<MotorcycleGraphTraits>::
Motorcycle_impl_base(const Point& source, const boost::optional<Point>& destination,
                     const FT speed, const boost::optional<Vector>& direction,
                     const FT initial_time)
  :
    i(-1),
    crashed(false),
    ini_sour_pt(source),
    ini_dest_pt(destination),
    sour(), dest(), conf(),
    is_dest_final(false),
    spee(speed),
    dir(direction),
    time(initial_time),
    time_at_sour(time),
    target_points(internal::Target_point_set_comparer<MotorcycleGraphTraits>()),
    track_points(internal::Track_comparer<MotorcycleGraphTraits>())
{
  // Reject null speed
  CGAL_precondition(spee > 0.);

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

template<typename MotorcycleGraphTraits>
void
Motorcycle_impl_base<MotorcycleGraphTraits>::
add_target(const DEC_it target_point, const FT time_at_target)
{
  target_points.insert(std::make_pair(target_point, time_at_target));
}

template<typename MotorcycleGraphTraits>
const typename Motorcycle_impl_base<MotorcycleGraphTraits>::DEC_it
Motorcycle_impl_base<MotorcycleGraphTraits>::
closest_target() const
{
  CGAL_precondition(!target_points.empty());
  return target_points.begin()->first;
}

template<typename MotorcycleGraphTraits>
void
Motorcycle_impl_base<MotorcycleGraphTraits>::
remove_closest_target_from_targets()
{
  CGAL_assertion(!target_points.empty());
  return target_points.erase(target_points.begin());
}

template<typename MotorcycleGraphTraits>
bool
Motorcycle_impl_base<MotorcycleGraphTraits>::
has_reached_blocked_point() const
{
  return conf->is_blocked();
}

template<typename MotorcycleGraphTraits>
bool
Motorcycle_impl_base<MotorcycleGraphTraits>::
has_reached_simultaneous_collision_point() const
{
  return conf->has_simultaneous_collision();
}

template<typename MotorcycleGraphTraits>
bool
Motorcycle_impl_base<MotorcycleGraphTraits>::
is_motorcycle_destination_final() const
{
  return is_dest_final;
}

template<typename MotorcycleGraphTraits>
typename Motorcycle_impl_base<MotorcycleGraphTraits>::FT
Motorcycle_impl_base<MotorcycleGraphTraits>::
time_at_closest_target() const
{
  CGAL_precondition(!target_points.empty());
  return target_points.begin()->second;
}

template<typename MotorcycleGraphTraits>
void
Motorcycle_impl_base<MotorcycleGraphTraits>::
output_track() const
{
  std::ostringstream out_filename;
  out_filename  << "results_" << Geom_traits::dimension << "/track_" << i << ".off" << std::ends;
  std::ofstream os(out_filename.str().c_str());
  os.precision(17);

  const std::size_t pn = track_points.size();
  CGAL_assertion(pn != 0);
  const std::size_t fn = pn - 1;

  os << "OFF" << '\n';
  os << pn << " " << fn << " 0" << '\n';

  std::cout << "track of motorcycle " << i << std::endl;

  typename Track::const_iterator tit = track_points.begin();
  typename Track::const_iterator end = track_points.end();
  for(; tit!=end; ++tit)
  {
    std::cout << *(tit->first) << std::endl;

    os << tit->first->point();

    if(Geom_traits::dimension == 2) // The xyz format expects 3D points
      os << " 0";
    os << '\n';
  }

  for(std::size_t j=0; j<fn; ++j)
    os << "3 " << j << " " << j+1 << " " << j << '\n';
}

template<typename MotorcycleGraphTraits>
void
Motorcycle_impl_base<MotorcycleGraphTraits>::
output_intended_track() const
{
  // must be adapted to multiple destinations and 2D/surface @todo

  std::ostringstream out_filename;
  out_filename << "results_" << Geom_traits::dimension << "/intended_track_" << i << ".off" << std::ends;
  std::ofstream os(out_filename.str().c_str());

  os << "OFF" << '\n';
  os << "2 1 0" << '\n';
  os << sour->point() << " 0" << '\n';
  os << dest->point() << " 0" << '\n';
  os << "3 0 1 0" << std::endl;
}

// -----------------------------------------------------------------------------

// Having a "_impl" class is only done because it's needed for BOOST_PARAMETER_CONSTRUCTOR
template<typename MotorcycleGraphTraits, typename Tracer>
class Motorcycle;

template<typename MotorcycleGraphTraits, typename Tracer>
class Motorcycle_impl
  : public Motorcycle_impl_base<MotorcycleGraphTraits>
{
  typedef Motorcycle_impl<MotorcycleGraphTraits, Tracer>            Self;
  typedef Motorcycle_impl_base<MotorcycleGraphTraits>               Base;
  typedef Motorcycle<MotorcycleGraphTraits, Tracer>                 Derived;

public:
  typedef MotorcycleGraphTraits                                     Geom_traits;
  typedef typename MotorcycleGraphTraits::Triangle_mesh             Triangle_mesh;

  typedef typename Base::FT                                         FT;
  typedef typename Base::Dictionary                                 Dictionary;
  typedef typename Base::DEC_it                                     DEC_it;
  typedef typename Base::Face_location                              Face_location;

  typedef typename Geom_traits::vertex_descriptor                   vertex_descriptor;
  typedef typename Geom_traits::halfedge_descriptor                 halfedge_descriptor;
  typedef typename Geom_traits::face_descriptor                     face_descriptor;
  typedef boost::variant<vertex_descriptor,
                         halfedge_descriptor,
                         face_descriptor>                           descriptor_variant;

  // Constructor
  virtual ~Motorcycle_impl() { }

  template<typename ArgumentPack>
  Motorcycle_impl(const ArgumentPack& args);

  virtual boost::tuple<bool, DEC_it, DEC_it, FT, bool>
  compute_next_destination(Dictionary& points, const Triangle_mesh& mesh);

private:
  // Tracer (computes the next target when we reach the destination)
  Tracer tracer;
};

template<typename MotorcycleGraphTraits, typename Tracer>
template <class ArgumentPack>
Motorcycle_impl<MotorcycleGraphTraits, Tracer>::
Motorcycle_impl(const ArgumentPack& args)
  :
    Base(args[parameters::source],
         args[parameters::destination|boost::none],
         args[parameters::speed|1.],
         args[parameters::direction|boost::none],
         args[parameters::initial_time|0.]),
    tracer(args[parameters::tracer|Tracer()])
{ }

template<typename MotorcycleGraphTraits, typename Tracer>
boost::tuple<bool, // successfuly computed a next path or not
             typename Motorcycle_impl<MotorcycleGraphTraits, Tracer>::DEC_it, // next source
             typename Motorcycle_impl<MotorcycleGraphTraits, Tracer>::DEC_it, // next destination
             typename Motorcycle_impl<MotorcycleGraphTraits, Tracer>::FT, // time at next destination
             bool> // whether the destination is final or not
Motorcycle_impl<MotorcycleGraphTraits, Tracer>::
compute_next_destination(Dictionary& points, const Triangle_mesh& mesh)
{
  CGAL_precondition(this->target_points.empty());

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*" << std::endl;
  std::cout << "Computing the next path for motorcycle #" << this->i << std::endl;
  std::cout << "Current position: " << this->conf->point() << std::endl
            << "Location: " << this->current_location().first << " b: "
            << this->current_location().second[0] << " "
            << this->current_location().second[1] << " "
            << this->current_location().second[2] << std::endl;
#endif

  const Face_location& loc = this->current_location();
  descriptor_variant dv =
    CGAL::Polygon_mesh_processing::internal::get_descriptor_from_location(loc, mesh);

  // The derived cast is so that the tracer uses the Motorcycle type and not the
  // Motorcycle_impl_base type. It is safe since we only manipulate "full" motorcycles,
  // but it's a bit ugly @fixme

  if(const vertex_descriptor* v = boost::get<vertex_descriptor>(&dv))
  {
    return tracer(*v, static_cast<Derived&>(*this), points, mesh);
  }
  else if(const halfedge_descriptor* h = boost::get<halfedge_descriptor>(&dv))
  {
    return tracer(*h, static_cast<Derived&>(*this), points, mesh);
  }
  else
  {
    const face_descriptor* f = boost::get<face_descriptor>(&dv);
    CGAL_assertion(f);
    return tracer(*f, static_cast<Derived&>(*this), points, mesh);
  }
}

// -----------------------------------------------------------------------------

template<typename MotorcycleGraphTraits, typename Tracer>
class Motorcycle
  : public Motorcycle_impl<MotorcycleGraphTraits, Tracer>
{
  typedef Motorcycle_impl<MotorcycleGraphTraits, Tracer>           Base;

public:
  BOOST_PARAMETER_CONSTRUCTOR(Motorcycle, (Base), parameters::tag,
                              (required (source_, *))
                              (optional (destination_, *)
                                        (speed_, *)
                                        (direction_, *)
                                        (initial_time_, *)
                                        (tracer_, *))
                             )
};

// -----------------------------------------------------------------------------

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_MOTORCYCLE_H
