// Copyright (c) 2011 Andreas von Dziegielewski and Michael Hemmer (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id:$
//
// Author(s)     : Michael Hemmer (mhsaar@googlemail.com)
//                 Andreas von Dziegielewski (dziegiel@uni-mainz.de)
//
// ================================================================================


#ifndef SV_COMPUTE_SCOUT_H
#define SV_COMPUTE_SCOUT_H

#include <SV/Octree_3.h>

namespace SV {

template <class INT>
typename Octree_3<INT>::Voxel 
compute_scout(
    const Octree_3<INT>& oct_hull0, 
    const typename Octree_3<INT>::Voxel& seed ){

  CGAL_precondition(oct_hull0.contains(seed));
  
  using CGAL::cpp0x::get;
  typedef typename Octree_3<INT>::Voxel Voxel;

  INT scout_x,scout_y,scout_z; 
  INT x,y,z;
  x = y = z = 0;
  
  CGAL_precondition(!oct_hull0.empty());
  CGAL_precondition(oct_hull0.contains(get<0>(seed),get<1>(seed),get<2>(seed)));

  CGAL_precondition(!oct_hull0.contains(x,y,z));
  while(x < get<0>(seed)){
    if(oct_hull0.contains(x+1,y,z)){
      scout_x = x; scout_y=y; scout_z=z;
      break;
    }
    x++;
  }

  CGAL_precondition(!oct_hull0.contains(x,y,z));
  while(y < get<1>(seed)){
    if(oct_hull0.contains(x,y+1,z)){
      scout_x = x; scout_y=y; scout_z=z;
      break;
    }
    y++;
  }

  CGAL_precondition(!oct_hull0.contains(x,y,z));
  while(z < get<2>(seed)){
    if(oct_hull0.contains(x,y,z+1)){
      scout_x = x; scout_y=y; scout_z=z;
      break;
    }
    z++;
  }

  return Voxel(scout_x,scout_y,scout_z);
}

} // namespace SV 

#endif // SV_COMPUTE_SCOUT_H
