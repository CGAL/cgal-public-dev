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

#ifndef SV_FILL_HULL_H
#define SV_FILL_HULL_H

#include <SV/Octree_3.h>
#include <SV/HashSet_Voxel.h>
#include <SV/Filtered_iterator.h>
#include <boost/optional.hpp>
#include <vector>
#include <CGAL/Real_timer.h> 

namespace SV{

//#if 1template <class INT, class InputIterator>  
template <class INT>
Octree_3<INT> fill_hull(
    int n,
    const HashSet_Voxel<INT>& hull,
    const typename HashSet_Voxel<INT>::Voxel& seed){
  return fill_hull(n,hull,&seed, (&seed)+1); 
}

template <class INT, class InputIterator>  
Octree_3<INT> fill_hull(
    int n,
    const HashSet_Voxel<INT>& hull,
    InputIterator begin_seeds, InputIterator end_seeds){
  
  //std::cerr << "fill_hull with hierarchy start" << std::endl; 
  CGAL::Real_timer timer; 
  timer.start(); 



  typedef HashSet_Voxel<INT> VoxelHashSet; 
  typedef typename VoxelHashSet::Voxel Voxel; 
  typedef Octree_3<INT> Octree; 
  
  //std::cerr << "Hull hierarchy .. "; 
  // ======== Compute Voxel Hull Hierarchy======== 
  std::vector<VoxelHashSet> voxel_hull_hierarchy(n,VoxelHashSet());
  std::vector<VoxelHashSet> seed_hierarchy(n,VoxelHashSet());
  
  //std::cerr << n-1 << " " << std::flush; 
  
  voxel_hull_hierarchy[n-1]=hull;
  seed_hierarchy[n-1]=VoxelHashSet(begin_seeds,end_seeds);
  for(int L = n-2; L >=0 ; L--){
    
    //std::cerr << L << " " << std::flush; 
    
    BOOST_FOREACH(Voxel v, voxel_hull_hierarchy[L+1]){
      voxel_hull_hierarchy[L].insert(Voxel( 
                                         CGAL::cpp0x::get<0>(v)>>1,
                                         CGAL::cpp0x::get<1>(v)>>1,
                                         CGAL::cpp0x::get<2>(v)>>1));       
    }
      
    // assert(!seed_hierarchy[L]);

    BOOST_FOREACH(const Voxel& s, seed_hierarchy[L+1]){
      INT x =  CGAL::cpp0x::get<0>(s);
      INT y =  CGAL::cpp0x::get<1>(s);
      INT z =  CGAL::cpp0x::get<2>(s);
      for(int i = x-1; i<= x+1; i++){
        for(int j = y-1; j<= y+1; j++){
          for(int k = z-1; k<= z+1; k++){
            int X = i>>1;
            int Y = j>>1;
            int Z = k>>1;
            Voxel candidate(X,Y,Z);
            if(!voxel_hull_hierarchy[L].contains(candidate)){ 
              seed_hierarchy[L].insert(candidate);
            }
          }
        }
      }
    } 
  }
  //std::cerr <<" done! " << std::endl; 
  
  // ==== FILL HULL AND STORE IT AS A OCTREE 
    
  //std::cerr << "Fill hierarchy .. "; 
  Octree offset0(0); 
  int ubound = 1;
  std::queue<Voxel> queue;
    
  for(int L = 0; L < n; L++){  
  
    
    offset0.set_resolution(offset0.resolution()+1);
    ubound <<=1;

    {
      std::list<Voxel> ohull;
      // get all voxel on the hull of the current offset0 
      // insert those into the offset and queue that are not in the hull
      offset0.hull(std::back_inserter(ohull));
      // also put the seed in case offset is empty 
      BOOST_FOREACH(const Voxel& seed,seed_hierarchy[L]){
        CGAL_precondition(!voxel_hull_hierarchy[L].contains(seed));
        ohull.push_back(seed);  
      }
      BOOST_FOREACH(Voxel v, ohull){
        if(!voxel_hull_hierarchy[L].contains(v)){
          queue.push(v);
        }
      }
    };
    
//    if(!queue.empty()) std::cerr << L << " " << std::flush; 
    
    while(!queue.empty()){
      INT x =  CGAL::cpp0x::get<0>(queue.front());
      INT y =  CGAL::cpp0x::get<1>(queue.front());
      INT z =  CGAL::cpp0x::get<2>(queue.front());
               
#if 1
      assert( x >= 0      && y >= 0      && z >= 0);      // is not water tight? 
      assert( x <= ubound && y <= ubound && z <= ubound); // is not water tight? 
#else
      if(!( x >= 3      && y >= 3      && z >= 3))  { queue.pop(); continue;}     // is not water tight? 
      if(!( x <= ubound && y <= ubound && z <= ubound))  {queue.pop(); continue;} // is not water tight? 
#endif
   
      // we only need to look at the cross and not at the cube 
      // since the volume is conservative  
      std::vector<Voxel> voxels; 
#if 0 
      voxels.push_back(Voxel(x+1,y,z));
      voxels.push_back(Voxel(x-1,y,z));
      voxels.push_back(Voxel(x,y+1,z));
      voxels.push_back(Voxel(x,y-1,z));
      voxels.push_back(Voxel(x,y,z+1));
      voxels.push_back(Voxel(x,y,z-1));
#else
      for(int i = -1; i <= 1; i++){
        for(int j = -1; j <= 1; j++){
          for(int k = -1; k <= 1; k++){
            voxels.push_back(Voxel(x+i,y+j,z+k));
          }
        }
      }
      
#endif
        
      BOOST_FOREACH(Voxel v, voxels){
        if(!voxel_hull_hierarchy[L].contains(v)){
          if(offset0.insert(Voxel(v))){
            queue.push(v);
          }
        }
      }        
      queue.pop();
    } 
    // std::cerr <<  offset0.size() << " " <<  queue.size()<< std::endl; 
  }
  //std::cerr <<" done! " << std::endl; 
 

  BOOST_FOREACH(Voxel v, hull){
    assert(!offset0.contains(v));
  }
  //std::cerr << " offset0.size()" << offset0.size() << std::endl; 

  timer.stop(); 
  std::cerr << " fill_hull time: " << timer.time() << std::endl;
  return offset0;

}
//#else

template <class INT>  
void fill_hull(
    int n,
    const HashSet_Voxel<INT>& hull,
    Octree_3<INT>& offset0
    //const typename HashSet_Voxel<INT>::Voxel seed
){
  //std::cerr << "fill_hull using octree as seed:  start" << std::endl; 
  CGAL::Real_timer timer; 
  timer.start(); 

  using CGAL::cpp0x::get;
  typedef HashSet_Voxel<INT> VoxelHashSet;
  typedef typename Octree_3<INT>::Voxel Voxel; 
  
  //std::cerr << " offset0.resolution() " <<  offset0.resolution() << std::endl; 
  int ubound = 1 << offset0.resolution();
  std::deque<Voxel> queue;
 
  internal::Contains<VoxelHashSet> contains(&hull);
  offset0.hull(internal::make_filter_iterator(std::back_inserter(queue),contains));

  while(!queue.empty()){
    Voxel v = queue.front(); queue.pop_front();
    if(hull.contains(v))   continue;
    if(!offset0.insert(v)) continue;
    
    INT x =  CGAL::cpp0x::get<0>(v);
    INT y =  CGAL::cpp0x::get<1>(v);
    INT z =  CGAL::cpp0x::get<2>(v);
    
#if 1
    assert( x >= 0      && y >= 0      && z >= 0);      // is not water tight? 
    assert( x <= ubound && y <= ubound && z <= ubound); // is not water tight? 
#else
    if(!( x >= 3      && y >= 3      && z >= 3))  { continue;}     // is not water tight? 
    if(!( x <= ubound && y <= ubound && z <= ubound))  {continue;} // is not water tight? 
#endif
  
#if 0
    queue.push_back(Voxel(x+1,y,z));
    queue.push_back(Voxel(x-1,y,z));
    queue.push_back(Voxel(x,y+1,z));
    queue.push_back(Voxel(x,y-1,z));
    queue.push_back(Voxel(x,y,z+1));
    queue.push_back(Voxel(x,y,z-1));
#else
    for(int i = -1; i <= 1; i++){
      for(int j = -1; j <= 1; j++){
        for(int k = -1; k <= 1; k++){
          queue.push_back(Voxel(x+i,y+j,z+k));
        }
      }
    }
#endif
   
    if(offset0.size() % 10000 == 0){
      //std::cerr << " queue size   " << queue.size() 
      //<< " octree size  " << offset0.size() 
      //        << " Memory usage " << memory_usage()*100 << " %\r"; 
    }
  } 
  
  timer.stop(); 
  std::cerr << " fill_hull time: " << timer.time() << std::endl;
}

//#endif 



}// namespace SV 

#endif // SV_FILL_HULL_H
