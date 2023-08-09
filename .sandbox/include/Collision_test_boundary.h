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

#ifndef COLLISION_TEST_BOUNDARY_H
#define COLLISION_TEST_BOUNDARY_H

#include <type_traits>
#include <variant>
#include <vector>
#include <Collision_type.h>
#include <Collision_function.h>
#include <Bilinear_patch_3.h>

namespace CGAL {
namespace Collisions {
namespace internal {

template <class K, COLLISION_TYPE C>
class Collision_test_boundary{

  public:
    typedef typename K::Point_3                             Point;
    typedef typename K::Vector_3                            Vector;
    typedef typename K::Triangle_3                          Triangle;
    typedef typename K::FT                                  FT;
    typedef typename ::CGAL::BilinearPatchC3<K>             Bilinear_patch;
    typedef          std::variant<Triangle, Bilinear_patch> Facet;
    typedef          std::vector<Facet>                     Facets;

  private:
    Collision_function<K,C> collision_function;
    FT ONE{1.};
    FT ZERO{0.};
    //       t     u     v
    Point A0{ZERO, ZERO, ZERO};
    Point B0{ZERO, ZERO, ONE };
    Point C0{ZERO, ONE , ONE };
    Point D0{ZERO, ONE , ZERO};
    Point A1{ONE,  ZERO, ZERO};
    Point B1{ONE,  ZERO, ONE };
    Point C1{ONE,  ONE , ONE };
    Point D1{ONE,  ONE , ZERO};

  public:
    Facets facets;

    Collision_test_boundary( 
      const Point& x0, 
      const Point& x1,
      const Point& x2,
      const Point& x3,
      const Point& x0_next, 
      const Point& x1_next,
      const Point& x2_next,
      const Point& x3_next
    ) 
      : collision_function(x0, x1, x2, x3, x0_next, x1_next, x2_next, x3_next)
      , facets{compute_boundary_facets()}
    {}

    template <COLLISION_TYPE C_>
    using EnableIfEdgeEdge = std::enable_if_t<(C_ == COLLISION_TYPE::EDGE_EDGE), bool>;
    
    template <COLLISION_TYPE C_>
    using EnableIfPointTriangle = std::enable_if_t<(C_ == COLLISION_TYPE::POINT_TRIANGLE), bool>;
    
    template <COLLISION_TYPE C_=C, EnableIfEdgeEdge<C_> = true> 
    Facets compute_boundary_facets() {
      return compute_edge_edge_boundary_facets();
    };

    template <COLLISION_TYPE C_=C, EnableIfPointTriangle<C_> = true>
    Facets compute_boundary_facets() {
      return compute_point_triangle_boundary_facets();
    };

    Facets compute_point_triangle_boundary_facets();
    Facets compute_edge_edge_boundary_facets();

    Facet make_facet(const Point& A, const Point& B, const Point& D);
    Facet make_facet(const Point& A, const Point& B, const Point& C, const Point& D);
};

template <class K, COLLISION_TYPE C>
auto Collision_test_boundary<K, C>::compute_point_triangle_boundary_facets() -> Facets
{
    Facets tmp;
    tmp.reserve(5);

    tmp.push_back(make_facet(A0, B0, D0));     // t = 0
    tmp.push_back(make_facet(A1, B1, D1));     // t = 1
    tmp.push_back(make_facet(A0, B0, B1, A1)); // u = 0
    tmp.push_back(make_facet(A0, D0, D1, A1)); // v = 0
    tmp.push_back(make_facet(B0, D0, D1, B1)); // extruded hypotenuse

    return tmp;
}

template <class K, COLLISION_TYPE C>
auto Collision_test_boundary<K, C>::compute_edge_edge_boundary_facets() -> Facets
{
    Facets tmp;
    tmp.reserve(6);

    tmp.push_back(make_facet(A0, B0, C0, D0)); // t = 0
    tmp.push_back(make_facet(A1, B1, C1, D1)); // t = 1
    tmp.push_back(make_facet(A0, B0, B1, A1)); // u = 0
    tmp.push_back(make_facet(C0, D0, D1, C1)); // u = 1
    tmp.push_back(make_facet(A0, D0, D1, A1)); // v = 0
    tmp.push_back(make_facet(B0, C0, C1, B1)); // v = 1

    return tmp;
}

template <class K, COLLISION_TYPE C>
auto Collision_test_boundary<K,C>::make_facet(const Point& A, const Point& B, const Point& D) -> Facet
{
    Facet f = Triangle(
        collision_function(A.x(), A.y(), A.z()), // t, u, v
        collision_function(B.x(), B.y(), B.z()), // t, u, v
        collision_function(D.x(), D.y(), D.z())  // t, u, v
    );
    return f;
}

template <class K, COLLISION_TYPE C>
auto Collision_test_boundary<K,C>::make_facet(const Point& A, const Point& B, const Point& C,  const Point& D) -> Facet
{
    Facet f = Bilinear_patch(
        collision_function(A.x(), A.y(), A.z()), // t, u, v
        collision_function(B.x(), B.y(), B.z()), // t, u, v
        collision_function(C.x(), C.y(), C.z()), // t, u, v
        collision_function(D.x(), D.y(), D.z())  // t, u, v
    );
    return f;
}


}
}
}

#endif