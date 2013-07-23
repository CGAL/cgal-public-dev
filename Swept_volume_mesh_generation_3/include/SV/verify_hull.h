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
//
// ================================================================================


#ifndef SV_VERIFY_HULL_H
#define SV_VERIFY_HULL_H

namespace SV { 

template < class INT, class InputIterator, class OutputIterator>
OutputIterator verify_hull(
    InputIterator begin, InputIterator end,
    const Octree_3<INT>& octree, 
    OutputIterator oit){
  typedef typename Octree_3<INT>::Voxel Voxel; 
  
  using CGAL::cpp0x::get; 
  
  for(InputIterator it = begin; it != end; it++){  
    Voxel v = *it; 
#define SV_TEST_VOXEL(i,j,k)                                            \
    Voxel vv(get<0>(v)+i,get<1>(v)+j,get<2>(v)+k);                      \
    if(octree.contains(vv)){                                            \
      *oit++ = vv;                                                      \
    }                                                                   \
    
    {SV_TEST_VOXEL(-1,0,0);}   
    {SV_TEST_VOXEL(+1,0,0);}   
    {SV_TEST_VOXEL(0,-1,0);}   
    {SV_TEST_VOXEL(0,+1,0);}   
    {SV_TEST_VOXEL(0,0,+1);}   
    {SV_TEST_VOXEL(0,0,-1);}   
#undef SV_TEST_VOXEL 
  }
  return oit; 
}

} // namespace SV

#endif // SV_VERIFY_HULL_H
