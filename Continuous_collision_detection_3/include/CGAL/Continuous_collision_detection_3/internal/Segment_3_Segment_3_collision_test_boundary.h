// Copyright (c) 2023 GeometryFactory (France).
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jeffrey Cochran

#ifndef SEGMENT_3_SEGMENT_3_COLLISION_TEST_BOUNDARY_H
#define SEGMENT_3_SEGMENT_3_COLLISION_TEST_BOUNDARY_H

#include <type_traits>
#include <variant>
#include <vector>

#include <CGAL/Continuous_collision_detection_3/internal/Cull_test_boundary.h>
#include <CGAL/Bilinear_patch_3.h>
#include <CGAL/Ray_3_Bilinear_patch_3_do_intersect.h>

namespace CGAL {
namespace Collisions {
namespace internal {



template <class K>
class Segment_3_Segment_3_collision_test_boundary{

  public:
    class Collision_function;
    typedef typename K::Point_3                       Point;
    typedef typename K::Segment_3                     Segment;
    typedef typename K::Vector_3                      Vector;
    typedef typename K::Ray_3                         Ray;
    typedef typename K::FT                            FT;
    typedef typename ::CGAL::Bilinear_patch_3<K>       Bilinear_patch;

  private:
    Collision_function          collision_function_;
    std::vector<Bilinear_patch> bilinear_patches_;
    std::vector<Point*>         points_;

    FT ONE{1.};
    FT ZERO{0.};

    bool IS_CULLED_{false};

    Point A0;
    Point B0;
    Point C0;
    Point D0;
    Point A1;
    Point B1;
    Point C1;
    Point D1;

  public:
    Segment_3_Segment_3_collision_test_boundary(
      const Segment& s0_cur,
      const Segment& s0_next,
      const Segment& s1_cur,
      const Segment& s1_next
    )
      : collision_function_(s0_cur, s0_next, s1_cur, s1_next)
    {
        // TODO: clean this up? I wrote them out
        //       explicitly to make it clear which
        //       vertices are which
        A0 = collision_function_(ZERO, ZERO, ZERO); // (t, u, v)
        B0 = collision_function_(ZERO, ZERO, ONE );
        C0 = collision_function_(ZERO, ONE , ONE );
        D0 = collision_function_(ZERO, ONE , ZERO);
        A1 = collision_function_(ONE,  ZERO, ZERO);
        B1 = collision_function_(ONE,  ZERO, ONE );
        C1 = collision_function_(ONE,  ONE , ONE );
        D1 = collision_function_(ONE,  ONE , ZERO);

        points_ = {&A0, &B0, &C0, &D0, &A1, &B1, &C1, &D1};

        IS_CULLED_ = ::CGAL::Collisions::internal::cull_test_boundary<K>(points_);

        if( !IS_CULLED_ )
        {
            bilinear_patches_.reserve(6);
            bilinear_patches_.push_back(Bilinear_patch(A0, B0, C0, D0));  // t = 0
            bilinear_patches_.push_back(Bilinear_patch(A1, B1, C1, D1));  // t = 1
            bilinear_patches_.push_back(Bilinear_patch(A0, B0, B1, A1));  // u = 0
            bilinear_patches_.push_back(Bilinear_patch(C0, D0, D1, C1));  // u = 1
            bilinear_patches_.push_back(Bilinear_patch(A0, D0, D1, A1));  // v = 0
            bilinear_patches_.push_back(Bilinear_patch(B0, C0, C1, B1));  // v = 1
        }
        // else {
        //     std::cout << "...CULLED EDGE-EDGE\n";
        // }
    }

    size_t num_ray_intersections(Ray r) const;

    const std::vector<Bilinear_patch>& bilinear_patches() const;
};

template <class K>
size_t Segment_3_Segment_3_collision_test_boundary<K>::num_ray_intersections(Ray r) const {


  if( IS_CULLED_ ) { return 0; }

  size_t num_intersections{0};

  for(const auto& bp : bilinear_patches_) {
    // If the ray's origin lies on any bilinear patch,
    // count that as a collision.
    // TODO: transition to exact arithmetic if 0 is in
    //       the interval.
    if (bp.has_on(r.source()))
    {
      // std::cout << "Point on bilinear patch: \n" << bp << std::endl;
      return 1;
    }

    // If the patch is degenerate and does not contain
    // the origin, ignore it.
    if (bp.is_degenerate())
    {
      // std::cout << "Bilinear patch is degenerate: \n" << bp << std::endl;
      continue;
    }

    // Otherwise, proceed with ray-intersection testing to
    // compute the parity of intersections with the boundary
    if( ::CGAL::do_intersect_odd_parity(bp, r) )
    {
      // std::cout << "Ray intersects bp: \n" << bp << std::endl;
      ++num_intersections;
    }

  }

  return num_intersections;
}

template <class K>
auto Segment_3_Segment_3_collision_test_boundary<K>::bilinear_patches() const
    -> const std::vector<Bilinear_patch>&
{
    return bilinear_patches_;
}

template <class K>
class Segment_3_Segment_3_collision_test_boundary<K>::Collision_function {

  public:
    typedef typename K::Point_3   Point;
    typedef typename K::Vector_3  Vector;
    typedef typename K::Segment_3 Segment;
    typedef typename K::FT        FT;
    typedef typename CGAL::Origin Origin;

  private:
    Origin origin = ::CGAL::ORIGIN;
    FT     ONE{1.};
    Vector x0;
    Vector x1;
    Vector x2;
    Vector x3;
    Vector x0_next;
    Vector x1_next;
    Vector x2_next;
    Vector x3_next;

  public:
    Collision_function(
      const Segment& s0_cur,
      const Segment& s0_next,
      const Segment& s1_cur,
      const Segment& s1_next
    )
      : x0{s0_cur.source() - origin}
      , x1{s0_cur.target() - origin}
      , x2{s1_cur.source() - origin}
      , x3{s1_cur.target() - origin}
      , x0_next{s0_next.source() - origin}
      , x1_next{s0_next.target() - origin}
      , x2_next{s1_next.source() - origin}
      , x3_next{s1_next.target() - origin}
      {}

    Point operator() (const FT& t, const FT& u, const FT& v) const {
      Vector interpolant = (
          (ONE - t)*((x0      + u*(x1      - x0     )) - (x2      + v*(x3      - x2     )))
        +       (t)*((x0_next + u*(x1_next - x0_next)) - (x2_next + v*(x3_next - x2_next)))
      );
      return origin + interpolant;

    };
};



}
}
}

#endif