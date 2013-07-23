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


#ifndef SV_COMPUTE_HULL_H
#define SV_COMPUTE_HULL_H

#include <SV/Octree_3.h>

namespace SV{


template < class INT, class OutputIterator>
OutputIterator compute_hull(
    const Octree_3<INT>& oct_hull0,
    const typename Octree_3<INT>::Voxel& scout, 
    OutputIterator oit){
  HashSet_Voxel<INT> set = compute_hull(oct_hull0,scout);
  return std::copy(set.begin(),set.end(),oit);
}

template < class INT >
inline bool is_in_hull_of(
    const typename Octree_3<INT>::Voxel& cv,
    const Octree_3<INT>& oct_hull0){
  using CGAL::cpp0x::get; 
  return ( 
      oct_hull0.contains(get<0>(cv)+1,get<1>(cv)+1,get<2>(cv)+1) ||
      oct_hull0.contains(get<0>(cv)-1,get<1>(cv)-1,get<2>(cv)-1) ||
      oct_hull0.contains(get<0>(cv)+1,get<1>(cv)+1,get<2>(cv)-1) ||
      oct_hull0.contains(get<0>(cv)+1,get<1>(cv)-1,get<2>(cv)+1) ||
      oct_hull0.contains(get<0>(cv)-1,get<1>(cv)+1,get<2>(cv)+1) ||
      oct_hull0.contains(get<0>(cv)-1,get<1>(cv)-1,get<2>(cv)+1) ||
      oct_hull0.contains(get<0>(cv)-1,get<1>(cv)+1,get<2>(cv)-1) ||
      oct_hull0.contains(get<0>(cv)+1,get<1>(cv)-1,get<2>(cv)-1) ||
      oct_hull0.contains(get<0>(cv)+1,get<1>(cv)+1,get<2>(cv)  ) ||
      oct_hull0.contains(get<0>(cv)+1,get<1>(cv)  ,get<2>(cv)+1) ||
      oct_hull0.contains(get<0>(cv)+1,get<1>(cv)  ,get<2>(cv)  ) ||
      oct_hull0.contains(get<0>(cv)+1,get<1>(cv)  ,get<2>(cv)-1) ||
      oct_hull0.contains(get<0>(cv)+1,get<1>(cv)-1,get<2>(cv)  ) ||
      oct_hull0.contains(get<0>(cv)  ,get<1>(cv)+1,get<2>(cv)+1) ||
      oct_hull0.contains(get<0>(cv)  ,get<1>(cv)+1,get<2>(cv)  ) ||
      oct_hull0.contains(get<0>(cv)  ,get<1>(cv)+1,get<2>(cv)-1) ||
      oct_hull0.contains(get<0>(cv)  ,get<1>(cv)  ,get<2>(cv)+1) ||
      //oct_hull0.contains(get<0>(cv)  ,get<1>(cv)  ,get<2>(cv)  ) ||
      oct_hull0.contains(get<0>(cv)  ,get<1>(cv)  ,get<2>(cv)-1) ||
      oct_hull0.contains(get<0>(cv)  ,get<1>(cv)-1,get<2>(cv)+1) ||
      oct_hull0.contains(get<0>(cv)  ,get<1>(cv)-1,get<2>(cv)  ) ||
      oct_hull0.contains(get<0>(cv)  ,get<1>(cv)-1,get<2>(cv)-1) ||
      oct_hull0.contains(get<0>(cv)-1,get<1>(cv)+1,get<2>(cv)  ) ||
      oct_hull0.contains(get<0>(cv)-1,get<1>(cv)  ,get<2>(cv)+1) ||
      oct_hull0.contains(get<0>(cv)-1,get<1>(cv)  ,get<2>(cv)  ) ||
      oct_hull0.contains(get<0>(cv)-1,get<1>(cv)  ,get<2>(cv)-1) ||
      oct_hull0.contains(get<0>(cv)-1,get<1>(cv)-1,get<2>(cv)));
}

template <class INT>  
HashSet_Voxel<INT>  compute_hull(
    const Octree_3<INT>& oct_hull0,
    const typename HashSet_Voxel<INT>::Voxel& scout,
    HashSet_Voxel<INT>& hull){
  std::cerr << "hull computation start (HULL A)" << std::endl; 
   
  CGAL::Real_timer timer;
  timer.start(); 
  using CGAL::cpp0x::get; 
  
  typedef HashSet_Voxel<INT> VoxelHashSet;
  typedef typename VoxelHashSet::Voxel Voxel; 
    
  CGAL_precondition(!oct_hull0.contains(scout));
  
  std::deque<Voxel> queue; 
//  HashSet_Voxel<INT> hull; 
  hull.clear(); 
  
  queue.push_back(scout);
  
  //int NN = CGAL::ipower(2,14);
   
  while(!queue.empty()){
    Voxel v = queue.front(); queue.pop_front();

    if(hull.contains(v)) continue;
    if(oct_hull0.contains(v)) continue;

    if(is_in_hull_of(v,oct_hull0)){
      assert(!oct_hull0.contains(v));
      hull.insert(v);
      queue.push_back(Voxel(get<0>(v)+1,get<1>(v)  ,get<2>(v)  ));
      queue.push_back(Voxel(get<0>(v)-1,get<1>(v)  ,get<2>(v)  ));
      queue.push_back(Voxel(get<0>(v)  ,get<1>(v)+1,get<2>(v)  ));
      queue.push_back(Voxel(get<0>(v)  ,get<1>(v)-1,get<2>(v)  ));
      queue.push_back(Voxel(get<0>(v)  ,get<1>(v)  ,get<2>(v)+1));
      queue.push_back(Voxel(get<0>(v)  ,get<1>(v)  ,get<2>(v)-1));
    }  

    if(hull.size()%10000==0)
      std::cerr << " queue size " << queue.size() 
                << " hull size " << hull.size() 
                << " Memory usage " << memory_usage()*100 << " %\r"; 
  }
  std::cerr << std::endl;

  timer.stop();
  std::cerr << "hull computation time:" << timer.time() << std::endl; 
  return hull; 
}

















#if 0
#if 1
template <class INT>  
HashSet_Voxel<INT>  compute_hull(
    const Octree_3<INT>& oct_hull0,
    const typename HashSet_Voxel<INT>::Voxel& scout
){
  std::cerr << "hull computation start (HULL B)" << std::endl; 
   
  CGAL::Real_timer timer;
  timer.start(); 
  using CGAL::cpp0x::get; 
  
  //std::cerr << get<0>(scout) << " " << get<1>(scout) << " " << get<2>(scout) << " " << std::endl;
   
  typedef HashSet_Voxel<INT> VoxelHashSet;
  typedef typename VoxelHashSet::Voxel Voxel; 
 
  Octree_3<INT>  oct_hull1 = oct_hull0.offset(); 
 
  CGAL_precondition(!oct_hull0.contains(scout));
  
  std::queue<Voxel> voxel_queue; 
  VoxelHashSet visited_voxels; 
  HashSet_Voxel<INT> hull_voxels; 
  
  voxel_queue.push(scout);
  hull_voxels.insert(scout);
  visited_voxels.insert(scout);
  
  //int NN = CGAL::ipower(2,14);
   
  while(!voxel_queue.empty()){
    Voxel v = voxel_queue.front(); 
    voxel_queue.pop();

    for ( short i = -1; i<=1;i++){
      for ( short j = -1; j<=1;j++){
        for ( short k = -1; k<=1;k++){
          if (i==0 && j == 0&& k==0) continue;

          Voxel cv(get<0>(v)+i,get<1>(v)+j,get<2>(v)+k);
          if(get<0>(cv) < 0 ) continue;
          if(get<1>(cv) < 0 ) continue;
          if(get<2>(cv) < 0 ) continue;
          
          if(oct_hull0.contains(cv)) continue;
          if(!oct_hull1.contains(cv)) continue;
          if(!visited_voxels.insert(cv)) continue; 
          
//          visited_voxels.insert(cv);
          voxel_queue.push(cv);
          hull_voxels.insert(cv);
        }
      }
    }
    if(hull_voxels.size()%10000==1)
      std::cerr << " queue size " << voxel_queue.size() 
                << " hull size " << hull_voxels.size() 
                << " Memory usage " << memory_usage()*100 << " %\r"; 
  }
  std::cerr << std::endl;
  timer.stop();
  std::cerr << "hull computation time:" << timer.time() << std::endl; 
  return hull_voxels; 
}


#else
template <class INT>  
HashSet_Voxel<INT>  compute_hull(
    const Octree_3<INT>& oct_hull0,
    const typename HashSet_Voxel<INT>::Voxel& scout
){

  std::cerr << "hull computation start (HULL  )" << std::endl; 
  CGAL::Real_timer timer;
  timer.start(); 
  using CGAL::cpp0x::get; 
  
  //std::cerr << get<0>(scout) << " " << get<1>(scout) << " " << get<2>(scout) << " " << std::endl;
   
  typedef HashSet_Voxel<INT> VoxelHashSet;
  typedef typename VoxelHashSet::Voxel Voxel; 
  
  
  CGAL_precondition(!oct_hull0.contains(scout));
  
  std::queue<Voxel> voxel_queue; 
  VoxelHashSet visited_voxels; 
  HashSet_Voxel<INT> hull_voxels; 
  
  voxel_queue.push(scout);
  hull_voxels.insert(scout);
  visited_voxels.insert(scout);
  
  //int NN = CGAL::ipower(2,14);
   
  while(!voxel_queue.empty()){
    Voxel v = voxel_queue.front(); 
    voxel_queue.pop();

    for ( short i = -1; i<=1;i++){
      for ( short j = -1; j<=1;j++){
        for ( short k = -1; k<=1;k++){
          if (i==0 && j == 0&& k==0) continue;

          Voxel cv(get<0>(v)+i,get<1>(v)+j,get<2>(v)+k);
          
          if(oct_hull0.contains(cv)) continue;
          if(visited_voxels.contains(cv)) continue; 
          visited_voxels.insert(cv);

          //if(!(get<0>(cv) >= 0   && get<1>(cv) < NN ) continue;
          //if(!(get<1>(cv) >= 0  && get<1>(cv) < NN) ) continue;
          //if(!(get<2>(cv) >= 0  && get<2>(cv) < NN) ) continue;
                   
          // put it into hull if one of its neighbours is in oct_hull0
          if( is_in_hull_of(cv,oct_hull0){
             voxel_queue.push(cv);
             hull_voxels.insert(cv);
          }
        }
      }
    }
    if(hull_voxels.size()%10000==1)
      std::cerr << " queue size " << voxel_queue.size() 
                << " hull size " << hull_voxels.size() 
                << " Memory usage " << memory_usage()*100 << " %\r"; 
  }
  std::cerr << std::endl;

  timer.stop();
  std::cerr << "hull computation time:" << timer.time() << std::endl; 
  return hull_voxels; 
}
#endif
#endif 


}// namespace SV 

#endif // SV_COMPUTE_HULL_H
