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

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
class Collision_function{
  public:
    typedef typename K::Point_3 Point;
    typedef typename K::FT      FT;

  private:
    Point* x0; 
    Point* x1;
    Point* x2;
    Point* x3;
    Point* x0_next; 
    Point* x1_next;
    Point* x2_next;
    Point* x3_next;
    COLLISION_TYPE ct;

    Point point_triangle_collision_function(const FT& t, const FT& u, const FT& v);
    Point edge_edge_collision_function(const FT& t, const FT& u, const FT& v);
    Point interpolate

  public:
    Collision_function( 
      Point* x0, 
      Point* x1,
      Point* x2,
      Point* x3,
      Point* x0_next, 
      Point* x1_next,
      Point* x2_next,
      Point* x3_next,
      Collision_type ct
    ) 
      : x0{x0}
      , x1{x1}
      , x2{x2}
      , x3{x3}
      , x0_next{x0_next}
      , x1_next{x1_next}
      , x2_next{x2_next}
      , x3_next{x3_next}
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

    return 
}

template <class K>
auto Collision_function<K>::edge_edge_collision_function(const FT& t, const FT& u, const FT& v) -> Point
{

}

}
}
}

#endif