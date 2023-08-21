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

#ifndef POINT_3_TRIANGLE_3_DO_COLLIDE_H
#define POINT_3_TRIANGLE_3_DO_COLLIDE_H

#include <Collision_test_boundaries.h>
#include <CGAL\Origin.h>
#include <Trajectories.h>
#include <cstdlib>
#include <ctime>
#include <iostream>

namespace CGAL{

  // TODO: speed this up.
  template <class K>
  bool do_collide(
      Point_3_trajectory<K> p,
      Triangle_3_trajectory<K> t,
      typename K::Ray_3       test_ray
  ){
    using Test_boundary = ::CGAL::Collisions::internal::Point_3_Triangle_3_collision_test_boundary<K>;
    using Ray = K::Ray_3;
    using Point = K::Point_3;
    using FT = K::FT;

    Test_boundary test_boundary(
      p.current(),
      p.next(),
      t.current(),
      t.next()
    );
    
    size_t num_intersections = test_boundary.num_ray_intersections( test_ray );
    
    return (num_intersections % 2) == 1; //
  }

  template <class K>
  bool do_collide(
      Point_3_trajectory<K> p,
      Triangle_3_trajectory<K> t
  ){
    using FT = K::FT;

    FT random_max = FT(RAND_MAX);
    FT x = FT(std::rand()) / random_max;
    FT y = FT(std::rand()) / random_max;
    FT z = FT(std::rand()) / random_max;

    typename K::Ray_3 test_ray = K::Ray_3(
        Point(::CGAL::ORIGIN), 
        Point(x, y, z)
    );

    return do_collide(p, t, test_ray); //
  }

} // end CGAL

#endif

