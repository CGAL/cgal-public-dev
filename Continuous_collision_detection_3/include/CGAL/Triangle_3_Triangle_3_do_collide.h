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

#ifndef TRIANGLE_3_TRIANGLE_3_DO_COLLIDE_H
#define TRIANGLE_3_TRIANGLE_3_DO_COLLIDE_H

#include <CGAL/Point_3_Triangle_3_do_collide.h>
#include <CGAL/Segment_3_Segment_3_do_collide.h>
#include <CGAL/Trajectories.h>
#include <vector>

namespace CGAL{

/// \defgroup do_collide_grp CGAL::do_collide()
/// \ingroup PkgCollisions3Predicates
/// @{

/// \brief Returns true if a collision occurs between the two triangle trajectories provided
/// \details This function determines whether a collision occurs by checking all possible segment-segment and point-triangle collisions that can be inferred from the provided triangle trajectories.
template <class K>
bool do_collide(
    const Triangle_3_trajectory<K>& t0,
    const Triangle_3_trajectory<K>& t1
){
  // TODO: minimize the number of copies occuring
  using Point_trajectory      = Point_3_trajectory<K>;
  using Segment_trajectory    = Segment_3_trajectory<K>;
  using Triangle_trajectory   = Triangle_3_trajectory<K>;
  using Point_triangle_pair   = decltype(std::make_pair(&std::get<0>(t0), &t1));

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

  std::vector<Point_triangle_pair> point_triangle_pairs;
  point_triangle_pairs.reserve(6);
  point_triangle_pairs.push_back(std::make_pair(&std::get<0>(t0), &t1));
  point_triangle_pairs.push_back(std::make_pair(&std::get<1>(t0), &t1));
  point_triangle_pairs.push_back(std::make_pair(&std::get<2>(t0), &t1));
  point_triangle_pairs.push_back(std::make_pair(&std::get<0>(t1), &t0));
  point_triangle_pairs.push_back(std::make_pair(&std::get<1>(t1), &t0));
  point_triangle_pairs.push_back(std::make_pair(&std::get<2>(t1), &t0));

  // Check point-triangle collision
  for( const auto& pt : point_triangle_pairs ) {
    if( do_collide( *pt.first, *pt.second ) )
    {
      return true;
    }
  }

  return false; //
}

/// @}


} // end CGAL

#endif

