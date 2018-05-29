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
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_POLYLINE_ROBUST_COLLINEAR_H
#define CGAL_POLYLINE_ROBUST_COLLINEAR_H

#include <CGAL/assertions.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

namespace CGAL {

namespace Polyline_tracing {

namespace internal {

template <class Traits>
bool are_logically_collinear_on_border(const typename Traits::Face_location& first_loc,
                                       const typename Traits::Face_location& second_loc,
                                       const typename Traits::Face_location& third_loc)
{
//  std::cout << "are_logically_collinear_on_border: " << std::endl;
//  std::cout << first_loc.first << " bc: " << first_loc.second[0] << " " << first_loc.second[1] << " " << first_loc.second[2] << std::endl;
//  std::cout << second_loc.first << " bc: " << second_loc.second[0] << " " << second_loc.second[1] << " " << second_loc.second[2] << std::endl;
//  std::cout << third_loc.first << " bc: " << third_loc.second[0] << " " << third_loc.second[1] << " " << third_loc.second[2] << std::endl;

  if(first_loc.first != second_loc.first ||
     first_loc.first != third_loc.first ||
     second_loc.first != third_loc.first)
    return false;

  for(int i=0; i<3; ++i)
    if(first_loc.second[i] == 0. && second_loc.second[i] == 0. && third_loc.second[i] == 0.)
      return true;

  return false;
}

template <class Traits>
bool are_logically_collinear_on_border(const typename Traits::Face_location& first_loc,
                                       const typename Traits::Face_location& second_loc,
                                       const typename Traits::Face_location& third_loc,
                                       const typename Traits::Face_location& fourth_loc)
{
  //  std::cout << "are_logically_collinear_on_border: " << std::endl;
  //  std::cout << first_loc.first << " bc: " << first_loc.second[0] << " " << first_loc.second[1] << " " << first_loc.second[2] << std::endl;
  //  std::cout << second_loc.first << " bc: " << second_loc.second[0] << " " << second_loc.second[1] << " " << second_loc.second[2] << std::endl;
  //  std::cout << third_loc.first << " bc: " << third_loc.second[0] << " " << third_loc.second[1] << " " << third_loc.second[2] << std::endl;
  //  std::cout << fourth_loc.first << " bc: " << fourth_loc.second[0] << " " << fourth_loc.second[1] << " " << fourth_loc.second[2] << std::endl;

  if(first_loc.first != second_loc.first || first_loc.first != third_loc.first ||
     first_loc.first != fourth_loc.first || second_loc.first != third_loc.first ||
     second_loc.first != fourth_loc.first || third_loc.first != fourth_loc.first)
    return false;

    for(int i=0; i<3; ++i)
      if(first_loc.second[i] == 0. && second_loc.second[i] == 0. &&
         third_loc.second[i] == 0. && fourth_loc.second[i] == 0.)
        return true;

    return false;
}

template <class Traits>
bool are_segments_collinear(const typename Traits::Segment_2& s,
                            const typename Traits::Segment_2& t,
                            const typename Traits::FT tolerance =
    std::numeric_limits<typename Traits::FT>::epsilon())
{
  CGAL_precondition(!s.is_degenerate());

  typedef typename Traits::FT   FT;

  const FT px = s.source().x();
  const FT py = s.source().y();
  const FT qx = s.target().x();
  const FT qy = s.target().y();
  const FT rx = t.source().x();
  const FT ry = t.source().y();

  if(CGAL::abs(determinant(qx-px, qy-py, rx-px, ry-py)) > tolerance)
    return false;

  const FT sx = t.target().x();
  const FT sy = t.target().y();

  if(CGAL::abs(determinant(qx-px, qy-py, sx-px, sy-py)) > tolerance)
    return false;

  return true;
}


// The purpose of this predicate is to find out the order of three points on
// a halfedge when numerical errors sneak in and `Collinear_2` returns `false`
//
// \pre q is different from 'p' and 'r'
// \pre p, q, and r are known to be collinear but might not be numerically
//
template <class Traits>
class Robust_collinear_are_strictly_ordered_along_line_2
  : public Traits
{
  typedef Traits                                              Base;

public:
  typedef typename Base::Point_2                              Point_2;

  typedef typename Base::Angle_2                              Angle_2;
  typedef typename Base::Collinear_2                          Collinear_2;

  typedef typename Base::Boolean                              result_type;

  Robust_collinear_are_strictly_ordered_along_line_2(const Traits& tr)
    : Base(tr)
  { }

  result_type
  operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
  {
    CGAL_precondition(p != q);
    CGAL_precondition(p != r);
    CGAL_precondition(q != r);

    if(Base::collinear_2_object()(p, q, r))
      return Base::collinear_are_strictly_ordered_along_line_2_object()(p, q, r);

    // @todo don't like that very much, something better ?
    if(Base::angle_2_object()(q, p, r) != CGAL::ACUTE)
      return false;

    if(Base::angle_2_object()(q, r, p) != CGAL::ACUTE)
      return false;

    // some (also disliked) alternative
#if 0
    Vector_2 st(s, t);
    Vector_2 sct2(s, ct2);
    FT sp = gt.compute_scalar_product_2_object()(st, sct2);
    if(sp < 0)
      return;

    Vector_2 ts(t, s);
    Vector_2 tct2(s, ct2);
    sp = gt.compute_scalar_product_2_object()(ts, tct2);
    if(sp < 0)
      return;
#endif

    return true;
  }
};

} // namespace internal

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_ROBUST_COLLINEAR_H
