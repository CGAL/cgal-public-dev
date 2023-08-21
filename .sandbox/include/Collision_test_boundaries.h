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

#ifndef COLLISION_TEST_BOUNDARIES_H
#define COLLISION_TEST_BOUNDARIES_H

#include <type_traits>
#include <variant>
#include <vector>
#include <Collision_functions.h>
#include <Bilinear_patch_3.h>
#include <Ray_3_Bilinear_patch_3_do_intersect.h>

namespace CGAL {
namespace Collisions {
namespace internal {



// ===========================================
// Shared functions
// ===========================================
template <class K>
bool cull_test_boundary(const std::vector<typename K::Point_3*>& points_) 
{
  // TODO: instantiate these normal vectors once at 
  //       compile time
  for (int i = -1; i < 1; ++i) {
  for (int j = -1; j < 1; ++j) {
  for (int k = -1; k < 1; ++k) {
    if(
      all_on_one_side_of_plane<K>(
        points_, 
        K::Vector_3(K::FT(i), K::FT(j), K::FT(k))
      )
    ) { 
      return true;
    }
  }
  }
  }
  return false;
}

template <class K>
bool all_on_one_side_of_plane(
  const std::vector<typename K::Point_3*>& points_, 
  const typename K::Vector_3& plane_normal
) {

  auto desired_orientation = ::CGAL::sign(
    ::CGAL::scalar_product( plane_normal, *(points_.front()) - ::CGAL::ORIGIN )
  );

  // If some of the points are on the plane,
  // it's possible that a line or patch contains the origin
  if( desired_orientation == ::CGAL::ZERO ) { return false; }

  return std::all_of(
    (++points_.begin()), // Don't need to check the first one
    points_.end(),
    [&plane_normal, &desired_orientation](const auto& m_p) {
      return (
            ::CGAL::sign(::CGAL::scalar_product( plane_normal, *m_p - ::CGAL::ORIGIN ))
        ==  desired_orientation // check same side of plane
      );
    }
  );
}



// ===========================================
// Point_3_Triangle_3_collision_test_boundary
// ===========================================
template <class K>
class Point_3_Triangle_3_collision_test_boundary{

  public:
    typedef typename K::Point_3                             Point;
    typedef typename K::Triangle_3                          Triangle;
    typedef typename K::Segment_3                           Segment;
    typedef typename K::Vector_3                            Vector;
    typedef typename K::Ray_3                               Ray;
    typedef typename K::FT                                  FT;
    typedef typename ::CGAL::BilinearPatchC3<K>             Bilinear_patch;
    typedef          std::variant<Triangle, Bilinear_patch> Facet;

  private:
    using Collision_function = Point_3_Triangle_3_collision_function<K>;
  
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

        // IS_CULLED_ = ::CGAL::Collisions::internal::cull_test_boundary<K>(points_);

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
    if( ::CGAL::Intersections::internal::do_intersect_odd_parity(bp, r) ) { 
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



// ===========================================
// Segment_3_Segment_3_collision_test_boundary
// ===========================================
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

        // IS_CULLED_ = ::CGAL::Collisions::internal::cull_test_boundary<K>(points_);

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
    if( ::CGAL::Intersections::internal::do_intersect_odd_parity(bp, r) ) 
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


}
}
}

#endif