// Copyright (c) 2023 GeometryFactory (France).
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jeffrey Cochran

#ifndef SEGMENT_3_SEGMENT_3_DO_COLLIDE_H
#define SEGMENT_3_SEGMENT_3_DO_COLLIDE_H

#include <CGAL/Continuous_collision_detection_3/internal/Segment_3_Segment_3_collision_test_boundary.h>

#include <CGAL/Origin.h>
#include <CGAL/Trajectories.h>

#include <cstdlib>
#include <ctime>
#include <iostream>

namespace CGAL{



  // TODO: speed this up.
  template <class K>
  bool do_collide(
      const Segment_3_trajectory<K>& s0,
      const Segment_3_trajectory<K>& s1,
      typename K::Ray_3       test_ray
  ){
    using Test_boundary = ::CGAL::Collisions::internal::Segment_3_Segment_3_collision_test_boundary<K>;
    using Ray = typename K::Ray_3;
    using Point = typename K::Point_3;
    using FT = typename K::FT;

    Test_boundary test_boundary(
      s0.current(),
      s0.next(),
      s1.current(),
      s1.next()
    );

    std::size_t num_intersections = test_boundary.num_ray_intersections( test_ray );

    return (num_intersections % 2) == 1; //
  }

  template <class K>
  bool do_collide(
      const Segment_3_trajectory<K>& s0,
      const Segment_3_trajectory<K>& s1
  ){
    using FT = typename K::FT;
    using Point = typename K::Point_3;

    //TODO: use CGAL::get_default_random() or ask CGAL::Random as a parameter
    FT random_max = FT(RAND_MAX);
    FT x = FT(std::rand()) / random_max;
    FT y = FT(std::rand()) / random_max;
    FT z = FT(std::rand()) / random_max;

    typename K::Ray_3 test_ray(
        Point(::CGAL::ORIGIN),
        Point(x, y, z)
    );

    return do_collide(s0, s1, test_ray); //
  }



} // end CGAL

#endif
