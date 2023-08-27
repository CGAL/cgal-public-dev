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

#include <CGAL/Point_3_Triangle_3_do_collide.h>
#include <CGAL/Segment_3_Segment_3_do_collide.h>
#include <CGAL/Trajectories.h>
#include <vector>

namespace CGAL{



  // TODO: minimize the number of copies occuring
  template <class K>
  bool do_collide(
      Triangle_3_trajectory<K> t0,
      Triangle_3_trajectory<K> t1
  ){
    using Point_trajectory      = Point_3_trajectory<K>;
    using Segment_trajectory    = Segment_3_trajectory<K>;
    using Triangle_trajectory   = Triangle_3_trajectory<K>;
    using Point_triangle_pair   = std::pair<Point_trajectory*, Triangle_trajectory*>;

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
          // std::cout << "...segment-segment collision: \n..." << s0 << "\n..." << s1 << "\n";
          return true;
        }
      }
    }

    std::vector<Point_triangle_pair> point_triangle_pairs{
      Point_triangle_pair(&std::get<0>(t0), &t1),
      Point_triangle_pair(&std::get<1>(t0), &t1),
      Point_triangle_pair(&std::get<2>(t0), &t1),
      Point_triangle_pair(&std::get<0>(t1), &t0),
      Point_triangle_pair(&std::get<1>(t1), &t0),
      Point_triangle_pair(&std::get<2>(t1), &t0),
    };

    // Check point-triangle collision
    for( const auto& pt : point_triangle_pairs ) {
      if( do_collide( *pt.first, *pt.second ) ) 
      {
        // std::cout << "...point-triangle collision: \n..." << *pt.first << "\n..." << *pt.second << "\n"; 
        return true;
      }
    }

    return false; //
  }



} // end CGAL

#endif

