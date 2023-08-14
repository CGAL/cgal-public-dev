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

#ifndef SEGMENT_3_SEGMENT_3_DO_COLLIDE_H
#define SEGMENT_3_SEGMENT_3_DO_COLLIDE_H

#include <Segment_3_Segment_3_collision_test_boundary.h>
#include <CGAL\Origin.h>
#include <Trajectories.h>
#include <cstdlib>
#include <ctime>
#include <iostream>

namespace CGAL{

  template <class K>
  bool do_collide(
      Segment_3_trajectory<K> s0,
      Segment_3_trajectory<K> s1
  ){
    using Test_boundary = ::CGAL::Collisions::internal::Segment_3_Segment_3_collision_test_boundary<K>;
    using Ray = K::Ray_3;
    using Point = K::Point_3;
    using FT = K::FT;

    Test_boundary test_boundary(
      s0.current(),
      s0.next(),
      s1.current(),
      s1.next()
    );

    FT random_number = FT(std::rand());
    Ray r(
        Point(::CGAL::ORIGIN), 
        Point(random_number, random_number, FT(0.))
    );
    size_t num_intersections = test_boundary.num_ray_intersections(r);

    return (num_intersections % 2) == 1; //
  }

} // end CGAL

#endif

