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

#ifndef POINT_3_TRIANGLE_3_COLLISION_FUNCTION_H
#define POINT_3_TRIANGLE_3_COLLISION_FUNCTION_H

#include <CGAL/Origin.h>



namespace CGAL {
namespace Collisions {
namespace internal {
    
template <class K>
class Point_3_Triangle_3_collision_function {

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
    Point_3_Triangle_3_collision_function( 
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