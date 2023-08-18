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

#ifndef TRIANGLE_3_TRIANGLE_3_DO_COLLIDE_H
#define TRIANGLE_3_TRIANGLE_3_DO_COLLIDE_H

#include <Segment_3_Segment_3_do_collide.h>
#include <Trajectories.h>
#include <vector>

namespace CGAL{

  template <class K>
  bool do_collide(
      Triangle_3_trajectory<K> t0,
      Triangle_3_trajectory<K> t1
  ){
    using Point_trajectory = Point_3_trajectory<K>;
    using Segment_trajectory = Segment_3_trajectory<K>;

    std::vector<Segment_trajectory> t0_segments{
      Segment_3_trajectory(std::get<0>(t0), std::get<1>(t0)),
      Segment_3_trajectory(std::get<1>(t0), std::get<2>(t0)),
      Segment_3_trajectory(std::get<2>(t0), std::get<0>(t0)),
    };

    std::vector<Segment_trajectory> t1_segments{
      Segment_3_trajectory(std::get<0>(t1), std::get<1>(t1)),
      Segment_3_trajectory(std::get<1>(t1), std::get<2>(t1)),
      Segment_3_trajectory(std::get<2>(t1), std::get<0>(t1)),
    };

    // Check edge-edge collision
    for( const Segment_trajectory& s0 : t0_segments ) {
      for ( const Segment_trajectory& s1 : t1_segments ) {
        if( do_collide(s0, s1) ) {
          std::cout << "edge-edge collision\n: " << s0 << "\n" << s1 << "\n";
          return true;
        }
      }
    }

    // Check point-triangle collision
    // for( const auto& p0 : t0 ) {
    //   if( do_collide(p0, t1)) {
    //     std::cout << "point-triangle collision\n";
    //     return true;
    //   }
    // }
    // for( const auto& p1 : t1 ) {
    //   if( do_coolide(p1, t0)) {
    //     std::cout << "point-triangle collision\n";
    //     return true;
    //   }
    // }

    return false; //
  }

} // end CGAL

#endif

