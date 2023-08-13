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
    std::vector<Triangle>       triangles_;
    std::vector<Bilinear_patch> bilinear_patches_;

    FT ONE{1.};
    FT ZERO{0.};

    Point A0{ZERO, ZERO, ZERO}; // (t, u, v)
    Point B0{ZERO, ZERO, ONE };
    Point C0{ZERO, ONE , ONE };
    Point D0{ZERO, ONE , ZERO};
    Point A1{ONE,  ZERO, ZERO};
    Point B1{ONE,  ZERO, ONE };
    Point C1{ONE,  ONE , ONE };
    Point D1{ONE,  ONE , ZERO};

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
      , triangles_{
          make_triangle_facet(A0, B0, D0),            // t = 0
          make_triangle_facet(A1, B1, D1)             // t = 1
        }
      , bilinear_patches_{
          make_bilinear_patch_facet(A0, B0, B1, A1),  // u = 0
          make_bilinear_patch_facet(A0, D0, D1, A1),  // v = 0
          make_bilinear_patch_facet(B0, D0, D1, B1)   // extruded hypotenuse
        }
    {}

    size_t num_ray_intersections(Ray r) const;
    std::vector<Facet> facets() const;
};


// ================
// Member functions
// ================
template <class K>
auto Point_3_Triangle_3_collision_test_boundary<K>::make_triangle_facet(const Point& A, const Point& B, const Point& D) -> Triangle
{              
    return Triangle( 
        collision_function_(A.x(), A.y(), A.z()), // collision_function(t, u, v)
        collision_function_(B.x(), B.y(), B.z()),
        collision_function_(D.x(), D.y(), D.z())  
    );
}

template <class K>
auto Point_3_Triangle_3_collision_test_boundary<K>::make_bilinear_patch_facet(const Point& A, const Point& B, const Point& C, const Point& D) -> Bilinear_patch
{
    return Bilinear_patch(
        collision_function_(A.x(), A.y(), A.z()), // t, u, v
        collision_function_(B.x(), B.y(), B.z()), // t, u, v
        collision_function_(C.x(), C.y(), C.z()), // t, u, v
        collision_function_(D.x(), D.y(), D.z())  // t, u, v
    );
}

template <class K>
size_t Point_3_Triangle_3_collision_test_boundary<K>::num_ray_intersections(Ray r) const {

  size_t num_intersections{0};

  for(const auto& t : triangles_) { 
    if( do_intersect(r, t) ) { ++num_intersections; }
  }

  for(const auto& bp : bilinear_patches_) {
    if( ::CGAL::Intersections::internal::do_intersect_odd_parity(bp, r) ) { ++num_intersections; }
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