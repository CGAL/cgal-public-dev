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

#ifndef COLLISION_FUNCTION_H
#define COLLISION_FUNCTION_H

#include <type_traits>
#include <Collision_type.h>
#include <CGAL\Origin.h>

namespace CGAL {
namespace Collisions {
namespace internal {

template <class K, COLLISION_TYPE C>
class Collision_function{

  public:
    typedef typename K::Point_3   Point;
    typedef typename K::Vector_3  Vector;
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

    Point point_triangle_collision_function(const FT& t, const FT& u, const FT& v);
    Point edge_edge_collision_function(const FT& t, const FT& u, const FT& v);

  public:
    Collision_function( 
      const Point& x0, 
      const Point& x1,
      const Point& x2,
      const Point& x3,
      const Point& x0_next, 
      const Point& x1_next,
      const Point& x2_next,
      const Point& x3_next
    ) 
      : x0{x0 - origin}
      , x1{x1 - origin}
      , x2{x2 - origin}
      , x3{x3 - origin}
      , x0_next{x0_next - origin}
      , x1_next{x1_next - origin}
      , x2_next{x2_next - origin}
      , x3_next{x3_next - origin}
      {}

    template <COLLISION_TYPE C_>
    using EnableIfEdgeEdge = std::enable_if_t<(C_ == COLLISION_TYPE::EDGE_EDGE), bool>;
    
    template <COLLISION_TYPE C_>
    using EnableIfPointTriangle = std::enable_if_t<(C_ == COLLISION_TYPE::POINT_TRIANGLE), bool>;
    
    template <COLLISION_TYPE C_=C, EnableIfEdgeEdge<C_> = true> 
    Point operator() (const FT& t, const FT& u, const FT& v){
      return edge_edge_collision_function(t, u, v);
    };

    template <COLLISION_TYPE C_=C, EnableIfPointTriangle<C_> = true>
    Point operator() (const FT& t, const FT& u, const FT& v){
      return point_triangle_collision_function(t, u, v);
    };
};

template <class K, COLLISION_TYPE C>
auto Collision_function<K, C>::point_triangle_collision_function(const FT& t, const FT& u, const FT& v) -> Point
{
    FT complement{ONE - t};
    
    Vector interpolant = (
                  complement*x0 + t*x0_next
      -(ONE-u-v)*(complement*x1 + t*x1_next)
      -        u*(complement*x2 + t*x2_next)
      -        v*(complement*x3 + t*x3_next)
    );

    return origin + interpolant;
}

template <class K, COLLISION_TYPE C>
auto Collision_function<K, C>::edge_edge_collision_function(const FT& t, const FT& u, const FT& v) -> Point
{
    FT complement_u{ONE - u};
    FT complement_v{ONE - v};
    
    Vector interpolant = (
       (ONE - t)*(complement_u*x0      + u*x1      - complement_v*x2      - v*x3     )
      -        t*(complement_u*x0_next + u*x1_next - complement_v*x2_next - v*x3_next)
    );
    
    return origin + interpolant;
}

}
}
}

#endif