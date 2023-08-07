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

#include <Collision_type.h>
#include <CGAL\Origin.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
class Collision_function{
  public:
    typedef typename K::Point_3   Point;
    typedef typename K::Vector_3  Vector;
    typedef typename K::FT        FT;

  private:
    Point  o(::CGAL::ORIGIN);
    FT     ONE{1.};
    Vector x0; 
    Vector x1;
    Vector x2;
    Vector x3;
    Vector x0_next; 
    Vector x1_next;
    Vector x2_next;
    Vector x3_next;
    COLLISION_TYPE ct;

    Point point_triangle_collision_function(const FT& t, const FT& u, const FT& v);
    Point edge_edge_collision_function(const FT& t, const FT& u, const FT& v);
    Point interpolate

  public:
    Collision_function( 
      const Point& x0, 
      const Point& x1,
      const Point& x2,
      const Point& x3,
      const Point& x0_next, 
      const Point& x1_next,
      const Point& x2_next,
      const Point& x3_next,
      Collision_type ct
    ) 
      : x0(o, x0)
      , x1(o, x1)
      , x2(o, x2)
      , x3(o, x3)
      , x0_next(o, x0)
      , x1_next(o, x1)
      , x2_next(o, x2)
      , x3_next(o, x3)
      , ct{ct}
      {}

    Point operator() (const FT& t, const FT& u, const FT& v)
    {
      switch(this->ct) 
      {
        case COLLISION_TYPE::edge_edge:
          return edge_edge_collision_function(t, u, v);
        case COLLISION_TYPE::point_triangle:
          return point_triangle_collision_function(t, u, v);
      }
    }
};


template <class K>
auto Collision_function<K>::point_triangle_collision_function(const FT& t, const FT& u, const FT& v) -> Point
{
    FT complement{ONE - t};
    
    Vector interpolant = (
                  complement*x0 + t*x0_next
      -(ONE-u-v)*(complement*x1 + t*x1_next)
      -        u*(complement*x2 + t*x2_next)
      -        v*(complement*x3 + t*x3_next)
    );

    return Point(interpolant);
}

template <class K>
auto Collision_function<K>::edge_edge_collision_function(const FT& t, const FT& u, const FT& v) -> Point
{
    FT complement_u{ONE - u};
    FT complement_v{ONE - v};
    
    Vector interpolant = (
       (ONE - t)*(complement_u*x0      + u*x1      - complement_v*x2      - v*x3     )
      -        t*(complement_u*x0_next + u*x1_next - complement_v*x2_next - v*x3_next)
    );
    
    return Point(interpolant);
}

}
}
}

#endif