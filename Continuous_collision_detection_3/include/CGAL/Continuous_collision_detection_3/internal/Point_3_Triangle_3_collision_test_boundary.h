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

#ifndef POINT_3_TRIANGLE_3_COLLISION_TEST_BOUNDARY_H
#define POINT_3_TRIANGLE_3_COLLISION_TEST_BOUNDARY_H

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
class Point_3_Triangle_3_collision_test_boundary {

  public:
    class Collision_function;
    typedef typename K::Point_3                             Point;
    typedef typename K::Triangle_3                          Triangle;
    typedef typename K::Segment_3                           Segment;
    typedef typename K::Vector_3                            Vector;
    typedef typename K::Ray_3                               Ray;
    typedef typename K::FT                                  FT;
    typedef typename ::CGAL::Bilinear_patch_3<K>             Bilinear_patch;
    typedef          std::variant<Triangle, Bilinear_patch> Facet;

  private:
    Collision_function          collision_function_;
    std::vector<Point*>         points_;
    std::vector<Triangle>       triangles_;
    std::vector<Bilinear_patch> bilinear_patches_;

    FT ONE{1.};
    FT ZERO{0.};

    bool IS_CULLED_{false};

    Point A0; // (t, u, v)
    Point B0;
    Point D0;
    Point A1;
    Point B1;
    Point D1;

  public:
    Point_3_Triangle_3_collision_test_boundary(
      const Point& p_cur,
      const Point& p_next,
      const Triangle& t_cur,
      const Triangle& t_next
    )
      : collision_function_(p_cur, p_next, t_cur, t_next)
    {

        // TODO: clean this up? I wrote them out
        //       explicitly to make it clear which
        //       vertices are which
        A0 = collision_function_(ZERO, ZERO, ZERO); // (t, u, v)
        B0 = collision_function_(ZERO, ZERO, ONE );
        D0 = collision_function_(ZERO, ONE , ZERO);
        A1 = collision_function_(ONE,  ZERO, ZERO);
        B1 = collision_function_(ONE,  ZERO, ONE );
        D1 = collision_function_(ONE,  ONE , ZERO);
        points_ = {&A0, &B0, &D0, &A1, &B1, &D1};

        IS_CULLED_ = ::CGAL::Collisions::internal::cull_test_boundary<K>(points_);

        if( !IS_CULLED_ ) {
            triangles_.reserve(2);
            triangles_.push_back(Triangle(A0, B0, D0));// t = 0
            triangles_.push_back(Triangle(A1, B1, D1));// t = 1

            bilinear_patches_.reserve(3);
            bilinear_patches_.push_back(Bilinear_patch(A0, B0, B1, A1));  // u = 0
            bilinear_patches_.push_back(Bilinear_patch(A0, D0, D1, A1));  // v = 0
            bilinear_patches_.push_back(Bilinear_patch(B0, D0, D1, B1));  // extruded hypotenuse
        }
        // else {
        //     std::cout << "...CULLED POINT-TRIANGLE\n";
        // }
    }

    size_t num_ray_intersections(Ray r) const;

    const std::vector<Bilinear_patch>& bilinear_patches() const;
    const std::vector<Triangle>& triangles() const;
};

template <class K>
size_t Point_3_Triangle_3_collision_test_boundary<K>::num_ray_intersections(Ray r) const {

  // If all the points lie one side of a plane
  // that passes through the origin, then
  // it's impossible for the origin to be inside
  // of the boundary
  if( IS_CULLED_ ) { return 0; }

  size_t num_intersections{0};

  for(const auto& t : triangles_) {
    // If the ray's origin lies on any triangle,
    // count that as a collision.
    // TODO: transition to exact arithmetic if 0 is in
    //       the interval.
    if (t.has_on(r.source()))
    {
      return 1;
    }

    if( t.is_degenerate() )
    {
        // std::cout << "Triangle is degenerate...\n";
        continue;
    }

    if( do_intersect(r, t) ) {
      ++num_intersections;
    }
  }

  for(const auto& bp : bilinear_patches_) {
    // If the ray's origin lies on any bilinear patch,
    // count that as a collision.
    // TODO: transition to exact arithmetic if 0 is in
    //       the interval.
    if (bp.has_on(r.source()))
    {
      return 1;
    }

    // If the patch is degenerate and does not contain
    // the origin, ignore it.
    if (bp.is_degenerate())
    {
        // std::cout << "Bilinear patch is degenerate...\n";
        continue;
    }

    // Otherwise, proceed with ray-intersection testing to
    // compute the parity of intersections with the boundary
    if( ::CGAL::do_intersect_odd_parity(bp, r) ) {
      ++num_intersections;
    }

  }

  return num_intersections;
}

template <class K>
auto Point_3_Triangle_3_collision_test_boundary<K>::bilinear_patches() const
    -> const std::vector<Bilinear_patch>&
{
    return bilinear_patches_;
}

template <class K>
auto Point_3_Triangle_3_collision_test_boundary<K>::triangles() const
    -> const std::vector<Triangle>&
{
    return triangles_;
}

template <class K>
class Point_3_Triangle_3_collision_test_boundary<K>::Collision_function {

  public:
    typedef typename K::Point_3     Point;
    typedef typename K::Vector_3    Vector;
    typedef typename K::Triangle_3  Triangle;
    typedef typename K::FT          FT;
    typedef typename ::CGAL::Origin Origin;

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
      const Point& p_cur,
      const Point& p_next,
      const Triangle& t_cur,
      const Triangle& t_next
    )
      : x0{p_cur - origin}
      , x1{t_cur[0] - origin}
      , x2{t_cur[1] - origin}
      , x3{t_cur[2] - origin}
      , x0_next{p_next - origin}
      , x1_next{t_next[0] - origin}
      , x2_next{t_next[1] - origin}
      , x3_next{t_next[2] - origin}
      {}

    Point operator() (const FT& t, const FT& u, const FT& v) const{
      FT complement{ONE - t};

      Vector interpolant = (
                    complement*x0 + t*x0_next
        -(ONE-u-v)*(complement*x1 + t*x1_next)
        -      (u)*(complement*x2 + t*x2_next)
        -      (v)*(complement*x3 + t*x3_next)
      );

      return origin + interpolant;
    };
};



}
}
}

#endif