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


template <class K>
bool
do_intersect_odd_parity(
  const typename CGAL::BilinearPatchC3<K> &bp,
  const typename K::Ray_3 &r
) {
  using Triangle = K::Triangle_3;
  using Ray = K::Ray_3;
  using Point = K::Point_3;
  using Interval = ::CGAL::Interval_nt_advanced;

  // Unfortunately, Bridson's methodology only works for determining the
  // parity of intersections. Returns true if parity is odd.

  if(bp.is_degenerate()) {
    return false;
  };

  CGAL_kernel_precondition(!r.is_degenerate());

  // Case 0
  // Bilinear patch is planar and can be treated as two triangles
  if (bp.is_planar()) {
    std::cout << "Case 0"<< std::endl;
    return (
          do_intersect(K::Triangle_3(bp.vertex(0), bp.vertex(1), bp.vertex(2)), r)
      ||  do_intersect(K::Triangle_3(bp.vertex(1), bp.vertex(2), bp.vertex(3)), r)
    );
  }

  // Case 1
  // Origin lies inside bounding tetrahedron
  Point ray_source = r.source();
  if (bp.tetrahedron().has_on_bounded_side(ray_source) || bp.tetrahedron().has_on_boundary(ray_source))
  {
    // TODO compute exact arithmetic if within the interval
    double phi_source = to_double(bp.aux_phi(ray_source));
    if (abs(phi_source) < 1e-14)
    {
        std::cout << "Case 1a\n";
      // If the ray's origin lies on the bilinear patch,
      // count that as an intersection. The function phi(x)
      // returns zero <==> x is on the patch.
      return true;
    } 
    else
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
      double phi_midpoint = to_double(bp.aux_phi(mid_point));
      if ( (phi_midpoint > 0) == (phi_source > 0) ) {
        std::cout << "Case 1b\n";
        // The edge connecting 0--2 is on the same side as the ray's source
        return (
             do_intersect(Triangle(bp.vertex(1), bp.vertex(2), bp.vertex(3)), r)
          || do_intersect(Triangle(bp.vertex(0), bp.vertex(1), bp.vertex(3)), r)
        );
      } else {
        std::cout << "Case 1c\n";
        // The edge connecting 0--2 is on the opposite side as the ray's source
        return (
             do_intersect(Triangle(bp.vertex(0), bp.vertex(1), bp.vertex(2)), r) 
          || do_intersect(Triangle(bp.vertex(0), bp.vertex(2), bp.vertex(3)), r)
        );
      }
    }
  } // End Case 1

  // Case 2
  // Origin lies outside the bounding tetrahedron
  if (
       !(do_intersect(Triangle(bp.vertex(0), bp.vertex(1), bp.vertex(2)), r)) 
    != !(do_intersect(Triangle(bp.vertex(0), bp.vertex(2), bp.vertex(3)), r))
  ) {
    std::cout << "Case 2\n";
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
