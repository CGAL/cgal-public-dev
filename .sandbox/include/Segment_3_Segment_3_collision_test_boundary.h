// Copyright (c) 2023
// INRIA Sophia-Antipolis (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jeffrey Cochran

#ifndef SEGMENT_3_SEGMENT_3_COLLISION_TEST_BOUNDARY_H
#define SEGMENT_3_SEGMENT_3_COLLISION_TEST_BOUNDARY_H

#include <type_traits>
#include <variant>
#include <vector>
#include <Segment_3_Segment_3_collision_function.h>
#include <Bilinear_patch_3.h>
#include <Ray_3_Bilinear_patch_3_do_intersect.h>

namespace CGAL {
namespace Collisions {
namespace internal {



// ====================================================
// TODO: document the concept for which this is a model
// ====================================================

template <class K>
class Segment_3_Segment_3_collision_test_boundary{

  public:
    typedef typename K::Point_3                       Point;
    typedef typename K::Segment_3                     Segment;
    typedef typename K::Vector_3                      Vector;
    typedef typename K::Ray_3                         Ray;
    typedef typename K::FT                            FT;
    typedef typename ::CGAL::BilinearPatchC3<K>       Bilinear_patch;
    // typedef typename ::CGAL::BilinearPatchC3<K>       Facet; // for consistency

  private:
    using Collision_function = Segment_3_Segment_3_collision_function<K>;

    Collision_function          collision_function_;
    std::vector<Bilinear_patch> bilinear_patches_;

    FT ONE{1.};
    FT ZERO{0.};

    Point A0 = Point(ZERO, ZERO, ZERO); // (t, u, v)
    Point B0 = Point(ZERO, ZERO, ONE );
    Point C0 = Point(ZERO, ONE , ONE );
    Point D0 = Point(ZERO, ONE , ZERO);
    Point A1 = Point(ONE,  ZERO, ZERO);
    Point B1 = Point(ONE,  ZERO, ONE );
    Point C1 = Point(ONE,  ONE , ONE );
    Point D1 = Point(ONE,  ONE , ZERO);

    Bilinear_patch make_bilinear_patch_facet(const Point& A, const Point& B, const Point& C, const Point& D) const;


  public:
    Segment_3_Segment_3_collision_test_boundary( 
      const Segment& s0_cur, 
      const Segment& s0_next, 
      const Segment& s1_cur, 
      const Segment& s1_next
    ) 
      : collision_function_(s0_cur, s0_next, s1_cur, s1_next)
    {
        bilinear_patches_.reserve(6);
        bilinear_patches_.push_back(make_bilinear_patch_facet(A0, B0, C0, D0));  // t = 0
        bilinear_patches_.push_back(make_bilinear_patch_facet(A1, B1, C1, D1));  // t = 1
        bilinear_patches_.push_back(make_bilinear_patch_facet(A0, B0, B1, A1));  // u = 0
        bilinear_patches_.push_back(make_bilinear_patch_facet(C0, D0, D1, C1));  // u = 1
        bilinear_patches_.push_back(make_bilinear_patch_facet(A0, D0, D1, A1));  // v = 0
        bilinear_patches_.push_back(make_bilinear_patch_facet(B0, C0, C1, B1));  // v = 1
    }

    size_t num_ray_intersections(Ray r) const;
    std::vector<Bilinear_patch> facets() const;
};

// ================
// Member functions
// ================
template <class K>
auto Segment_3_Segment_3_collision_test_boundary<K>::make_bilinear_patch_facet(
  const Point& A, 
  const Point& B, 
  const Point& C,  
  const Point& D
) const -> Bilinear_patch
{   
    return Bilinear_patch(
      collision_function_(A.x(), A.y(), A.z()), // t, u, v
      collision_function_(B.x(), B.y(), B.z()), // t, u, v
      collision_function_(C.x(), C.y(), C.z()), // t, u, v
      collision_function_(D.x(), D.y(), D.z())  // t, u, v
    );
}

template <class K>
size_t Segment_3_Segment_3_collision_test_boundary<K>::num_ray_intersections(Ray r) const {

  size_t num_intersections{0};

  for(const auto& bp : bilinear_patches_) {
    // If the ray's origin lies on any bilinear patch,
    // count that as a collision.
    // TODO: transition to exact arithmetic if 0 is in
    //       the interval.
    if (bp.has_on(r.source()))
    {
      std::cout << "Origin on patch...\n";
      std::cout << bp << "\n";
      return 1;
    } 

    // Otherwise, proceed with ray-intersection testing to
    // compute the parity of intersections with the boundary
    if( ::CGAL::Intersections::internal::do_intersect_odd_parity(bp, r) ) { ++num_intersections; }

  }

  return num_intersections;
}

template<class K>
auto Segment_3_Segment_3_collision_test_boundary<K>::facets() const -> std::vector<Bilinear_patch> {
  return this->bilinear_patches_;
}



}
}
}

#endif