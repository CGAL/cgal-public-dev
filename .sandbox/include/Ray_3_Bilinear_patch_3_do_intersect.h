// Copyright (c) 2023  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jeffrey Cochran

#ifndef CGAL_INTERNAL_INTERSECTIONS_RAY_3_BILINEAR_PATCH_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_RAY_3_BILINEAR_PATCH_3_DO_INTERSECT_H

#include <iostream>

#include <Bilinear_patch_3.h>
#include <Collision_candidate_3.h>

#include <CGAL/enum.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/intersection_3.h>
#include <CGAL/Interval_nt.h>

namespace CGAL {
namespace Intersections {
namespace internal {

enum class RAY_3_BILINEAR_PATCH_3_INTERSECTION_TYPE{ON_SURFACE, PASSING_THROUGH, NONE};

template <class K>
bool
do_intersect_odd_parity(
  const typename CGAL::BilinearPatchC3<K> &bp,
  const typename K::Ray_3 &r
) {
  // NOTE: This assumes that the origin of the ray is not on the patch.

  using Triangle = K::Triangle_3;
  using Segment = K::Segment_3;
  using Ray = K::Ray_3;
  using Point = K::Point_3;
  using Interval = ::CGAL::Interval_nt_advanced;

  CGAL_kernel_precondition(!bp.is_degenerate());
  CGAL_kernel_precondition(!r.is_degenerate());

  Point ray_source = r.source();

  // Case 1
  // The bilinear patch degenerates to coplanar triangles
  if(bp.is_planar()) {
    bool does_intersect_odd_parity_{false};
    for(const auto& t : bp.triangles_)
    { 
      does_intersect_odd_parity_ = does_intersect_odd_parity_ || do_intersect(t, r);
    }
    return does_intersect_odd_parity_;
  };

  // Case 2
  // Origin lies inside bounding tetrahedron
  if (bp.tetrahedron().has_on_bounded_side(ray_source) || bp.tetrahedron().has_on_boundary(ray_source))
  {
    // Otherwise, check the sign of phi(origin). Two of the bounding
    // tetrahedron's four triangles lie on the positive side of phi()==0,
    // and two lie on the negative side. If the origin is on one side,
    // check the ray for intersection with the two triangles on the other side
    const Point & a = bp.vertex(0);
    const Point & c = bp.vertex(2);
    const Point mid_point = Point(
      (a.x() + c.x())/2.,
      (a.y() + c.y())/2.,
      (a.z() + c.z())/2.
    ); 

    //  This will determine which triangles are on the opposite side
    double phi_source   = bp.signed_scaled_patch_distance(ray_source);
    double phi_midpoint = bp.signed_scaled_patch_distance(mid_point);
    if ( (phi_midpoint > 0) == (phi_source > 0) ) {
      // The edge connecting 0--2 is on the same side as the ray's source
      return (
            do_intersect(Triangle(bp.vertex(1), bp.vertex(2), bp.vertex(3)), r)
        ||  do_intersect(Triangle(bp.vertex(0), bp.vertex(1), bp.vertex(3)), r)
      );
    } else {
      // The edge connecting 0--2 is on the opposite side as the ray's source
      return (
            do_intersect(Triangle(bp.vertex(0), bp.vertex(1), bp.vertex(2)), r) 
        ||  do_intersect(Triangle(bp.vertex(0), bp.vertex(2), bp.vertex(3)), r)
      );
    }
  } // End Case 2

  // Case 3
  // Origin lies outside the bounding tetrahedron
  if (
       !(do_intersect(Triangle(bp.vertex(0), bp.vertex(1), bp.vertex(2)), r)) 
    != !(do_intersect(Triangle(bp.vertex(0), bp.vertex(2), bp.vertex(3)), r))
  ) {
    // The ray intersects exactly one of the bounding tetrahedron's two
    // triangles on the positive/negative side of phi(*)==0 _if_and_only_if_ the ray
    // intersects the bilinear patch an odd number of times.
    return true;
  }

  return false;
}



} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_RAY_3_BILINEAR_PATCH_3_DO_INTERSECT_H
