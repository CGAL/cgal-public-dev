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

#ifndef POINT_3_TRIANGLE_3_COLLISION_TEST_BOUNDARY_H
#define POINT_3_TRIANGLE_3_COLLISION_TEST_BOUNDARY_H

#include <type_traits>
#include <variant>
#include <vector>
#include <Point_3_Triangle_3_collision_function.h>
#include <Bilinear_patch_3.h>
#include <Ray_3_Bilinear_patch_3_do_intersect.h>

namespace CGAL {
namespace Collisions {
namespace internal {

// ====================================================
// TODO: document the concept for which this is a model
// ====================================================

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
    // namespace Intersections  = ;

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

    Triangle make_triangle_facet(const Point& A, const Point& B, const Point& D);
    Bilinear_patch make_bilinear_patch_facet(const Point& A, const Point& B, const Point& C, const Point& D);

  public:
    Point_3_Triangle_3_collision_test_boundary( 
      const Point& p_cur, 
      const Point& p_next,
      const Triangle& t_cur,
      const Triangle& t_next
    ) 
      : collision_function_(p_cur, p_next, t_cur, t_next)
    {

      A0 = collision_function_(ZERO, ZERO, ZERO); // (t, u, v)
      B0 = collision_function_(ZERO, ZERO, ONE );
      D0 = collision_function_(ZERO, ONE , ZERO);
      A1 = collision_function_(ONE,  ZERO, ZERO);
      B1 = collision_function_(ONE,  ZERO, ONE );
      D1 = collision_function_(ONE,  ONE , ZERO);
      points_ = {&A0, &B0, &D0, &A1, &B1, &D1};

      IS_CULLED_ = false;//are_cullable(points_);

      if( !IS_CULLED_ ) {
        triangles_.reserve(2);
        triangles_.push_back(Triangle(A0, B0, D0));// t = 0
        triangles_.push_back(Triangle(A1, B1, D1));// t = 1
        
        bilinear_patches_.reserve(3);
        bilinear_patches_.push_back(Bilinear_patch(A0, B0, B1, A1));  // u = 0
        bilinear_patches_.push_back(Bilinear_patch(A0, D0, D1, A1));  // v = 0
        bilinear_patches_.push_back(Bilinear_patch(B0, D0, D1, B1));  // extruded hypotenuse
      }
    }

    size_t num_ray_intersections(Ray r) const;
    std::vector<Facet> facets() const;
};


// ================
// Member functions
// ================
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
auto Point_3_Triangle_3_collision_test_boundary<K>::facets() const -> std::vector<Facet> {
  std::vector<Facet> tmp;
  tmp.reserve(5);

  tmp.push_back(triangles_[0]);
  tmp.push_back(triangles_[1]);
  tmp.push_back(bilinear_patches_[0]);
  tmp.push_back(bilinear_patches_[1]);
  tmp.push_back(bilinear_patches_[2]);

  return tmp;
}



}
}
}

#endif