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

#include <CGAL/Bilinear_patch_3.h>
#include <CGAL/Collision_candidate_3.h>

#include <CGAL/enum.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/intersection_3.h>

namespace CGAL {



// TODO: filter this predicate. 
// TODO: make sure that bp is constructed exactly in 
//       the collision test boundary.
template <class K>
bool
do_intersect_odd_parity(
  const typename CGAL::BilinearPatchC3<K> &bp,
  const typename K::Ray_3 &r
) {

  using Triangle = typename K::Triangle_3;
  using Segment = typename K::Segment_3;
  using Ray = typename K::Ray_3;
  using Point = typename K::Point_3;
  
  // TODO: decide whether the to treat point-not-on-boundary
  //       as a precondition 
  CGAL_kernel_precondition(!bp.is_degenerate());
  CGAL_kernel_precondition(!r.is_degenerate());

  Point ray_source = r.source();

  // Case 1
  // The bilinear patch degenerates to coplanar triangles
  // TODO: make sure this works for the case where the 
  //       ray is coplanar with the bilinear patch
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
  // THIS IMPLEMENTS BROCHU & BRIDSON 2012
  //if (bp.tetrahedron().has_on_bounded_side(ray_source) || bp.tetrahedron().has_on_boundary(ray_source))
  if ( !bp.tetrahedron().has_on_unbounded_side(ray_source) )
  {
    // Otherwise, check the sign of phi(origin). Two of the bounding
    // tetrahedron's four triangles lie on the positive side of phi()==0,
    // and two lie on the negative side. If the origin is on one side,
    // check the ray for intersection with the two triangles on the other side
    Point mid_point = ::CGAL::midpoint(bp.vertex(0), bp.vertex(2));

    //  This will determine which triangles are on the opposite side
    auto orientation_ray_source = bp.orientation(ray_source);
    auto orientation_midpoint   = bp.orientation(mid_point);
    if ( orientation_ray_source == orientation_midpoint ) {
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
  // TODO: check to see if ray intersects vertex/edge
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



} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_RAY_3_BILINEAR_PATCH_3_DO_INTERSECT_H
