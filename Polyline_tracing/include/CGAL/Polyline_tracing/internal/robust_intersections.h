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

#ifndef CGAL_POLYLINE_TRACING_ROBUST_INTERSECTIONS_H
#define CGAL_POLYLINE_TRACING_ROBUST_INTERSECTIONS_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/intersection_2.h>
#include <CGAL/intersection_3.h>

#include <boost/optional.hpp>

namespace CGAL {

namespace Polyline_tracing {

namespace internal {

template<typename EK, int dimension>
class Exact_kernel_d;

template<typename EK>
class Exact_kernel_d<EK, 2>
  : public EK
{
  typedef EK                                                             Base;

public:
  typedef typename Base::Point_2                                         Exact_point_d;
  typedef typename Base::Ray_2                                           Exact_ray_d;
  typedef typename Base::Segment_2                                       Exact_segment_d;

  typedef typename Base::Do_intersect_2                                  Do_intersect_d;
  typedef typename Base::Intersect_2                                     Intersect_d;

  typedef typename CGAL::cpp11::result_of<
    typename Base::Intersect_2(Exact_ray_d, Exact_segment_d)>::type      Exact_RS_intersection_result;
  typedef typename CGAL::cpp11::result_of<
    typename Base::Intersect_2(Exact_segment_d, Exact_segment_d)>::type  Exact_SS_intersection_result;

  Exact_kernel_d(const EK& ek = EK()) : Base(ek) { }

  Do_intersect_d
  do_intersect_d_object() const
  { return Base::do_intersect_2_object(); }

  Intersect_d
  intersect_d_object() const
  { return Base::intersect_2_object(); }
};

template<typename EK>
class Exact_kernel_d<EK, 3>
  : public EK
{
  typedef EK                                                             Base;

public:
  typedef typename Base::Point_3                                         Exact_point_d;
  typedef typename Base::Ray_3                                           Exact_ray_d;
  typedef typename Base::Segment_3                                       Exact_segment_d;

  typedef typename Base::Do_intersect_3                                  Do_intersect_d;
  typedef typename Base::Intersect_3                                     Intersect_d;

  typedef typename CGAL::cpp11::result_of<
    typename Base::Intersect_3(Exact_ray_d, Exact_segment_d)>::type      Exact_RS_intersection_result;
  typedef typename CGAL::cpp11::result_of<
    typename Base::Intersect_3(Exact_segment_d, Exact_segment_d)>::type  Exact_SS_intersection_result;

  Exact_kernel_d(const EK& ek = EK()) : Base(ek) { }

  Do_intersect_d
  do_intersect_d_object() const
  { return Base::do_intersect_3_object(); }

  Intersect_d
  intersect_d_object() const
  { return Base::intersect_3_object(); }
};

template<typename Traits>
boost::optional<typename Traits::Point_d>
robust_intersection(const typename Traits::Ray_d& r,
                    const typename Traits::Segment_d& s,
                    const Traits& traits = Traits())
{
  // This function should only be called if there is an intersection between
  // the ray and the segment
  CGAL_precondition(traits.do_intersect_d_object()(r, s));

  typedef typename Traits::Point_d                                 Point;
  typedef typename Traits::Segment_d                               Segment;
  typedef typename Traits::Ray_d                                   Ray;

  const Point* pp;
  boost::optional<Point> intersection_point;

  typedef typename CGAL::cpp11::result_of<
      typename Traits::Intersect_d(Ray, Segment)>::type Intersection_result;
  Intersection_result res = traits.intersect_d_object()(r, s);

  bool need_to_use_exact = false;
  if(!res) // There should be an intersection
  {
    need_to_use_exact = true;
  }
  else
  {
    // @todo (?) instead of ignoring the degenerate 'segment' return type,
    // return the opposite extremity of the segment (the one that is not the ray's source)
    pp = boost::get<Point>(&*res);

    // another type of sanity check ? @todo
    if(!pp || !r.has_on(*pp) || !s.has_on(*pp))
      need_to_use_exact = true;
  }

  if(need_to_use_exact)
  {
    typedef CGAL::Exact_predicates_exact_constructions_kernel   EPECK;
    typedef Exact_kernel_d<EPECK, Traits::dimension>            EK_d;

    EK_d ek;

    typedef typename EK_d::Exact_point_d                        Exact_point;
    typedef typename EK_d::Exact_ray_d                          Exact_ray;
    typedef typename EK_d::Exact_segment_d                      Exact_segment;
    typedef typename EK_d::Exact_RS_intersection_result         Exact_intersection_result;

    CGAL::Cartesian_converter<Traits, EPECK> to_exact;
    CGAL::Cartesian_converter<EPECK, Traits> from_exact;
    Exact_ray exact_r = to_exact(r);
    Exact_segment exact_s = to_exact(s);

    Exact_intersection_result exact_res = ek.intersect_d_object()(exact_r, exact_s);
    CGAL_assertion(exact_res);

    const Exact_point* epp = boost::get<Exact_point>(&*exact_res);
    if(!epp)
    {
      std::cerr << "Warning: intersection is not a point" << std::endl;
      return intersection_point;
    }

    const Exact_point exact_intersection_point = *epp;
    CGAL_postcondition(exact_s.has_on(exact_intersection_point));
    intersection_point = from_exact(exact_intersection_point);
  }
  else
  {
    intersection_point = *pp;
  }

  return intersection_point;
}

template<typename Traits>
typename Traits::Point_d
robust_intersection(const typename Traits::Segment_d& s,
                    const typename Traits::Segment_d& t,
                    const Traits& traits = Traits())
{
  // this function should only be called in the case of non collinear segments
  // that are known to intersect
  CGAL_precondition(!traits.collinear_d_object()(s.source(), s.target(), t.source()) ||
                    !traits.collinear_d_object()(s.source(), s.target(), t.target()));
  CGAL_precondition(traits.do_intersect_d_object()(s, t));

  typedef typename Traits::Point_d                                 Point;
  typedef typename Traits::Segment_d                               Segment;

  const Point* pp;
  Point intersection_point;

  typedef typename CGAL::cpp11::result_of<
      typename Traits::Intersect_d(Segment, Segment)>::type Intersection_result;
  Intersection_result res = traits.intersect_d_object()(s, t);

  bool need_to_use_exact = false;
  if(!res) // There should be an intersection
  {
    need_to_use_exact = true;
  }
  else
  {
    // The intersection cannot be a segment since the tracks are not collinear
    pp = boost::get<Point>(&*res);
    if(!pp || !s.has_on(*pp) || !t.has_on(*pp)) // another type of sanity check ? @todo
      need_to_use_exact = true;
  }

  if(need_to_use_exact)
  {
    typedef CGAL::Exact_predicates_exact_constructions_kernel   EPECK;
    typedef Exact_kernel_d<EPECK, Traits::dimension>            EK_d;

    EK_d ek;

    typedef typename EK_d::Exact_point_d                        Exact_point;
    typedef typename EK_d::Exact_segment_d                      Exact_segment;
    typedef typename EK_d::Exact_SS_intersection_result         Exact_intersection_result;

    CGAL::Cartesian_converter<Traits, EPECK> to_exact;
    CGAL::Cartesian_converter<EPECK, Traits> from_exact;
    Exact_segment exact_s = to_exact(s);
    Exact_segment exact_t = to_exact(t);

    CGAL_precondition(!CGAL::collinear(exact_s.source(), exact_s.target(), exact_t.source()) ||
                      !CGAL::collinear(exact_s.source(), exact_s.target(), exact_t.target()));
    CGAL_precondition(ek.do_intersect_d_object()(exact_s, exact_t));

    Exact_intersection_result exact_res = ek.intersect_d_object()(exact_s, exact_t);
    CGAL_assertion(exact_res);

    const Exact_point* epp = boost::get<Exact_point>(&*exact_res);
    const Exact_point exact_intersection_point = *epp;
    CGAL_assertion(epp);
    CGAL_postcondition(exact_s.has_on(exact_intersection_point) &&
                       exact_t.has_on(exact_intersection_point));
    intersection_point = from_exact(exact_intersection_point);
  }
  else
  {
    intersection_point = *pp;
  }

  return intersection_point;
}

} // namespace internal

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_TRACER_H
