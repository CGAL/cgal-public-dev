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

#ifndef SEGMENT_3_SEGMENT_3_COLLISION_FUNCTION_H
#define SEGMENT_3_SEGMENT_3_COLLISION_FUNCTION_H

namespace CGAL {
namespace Collisions {
namespace internal {
    


template <class K>
class Segment_3_Segment_3_collision_function {

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
    Segment_3_Segment_3_collision_function( 
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