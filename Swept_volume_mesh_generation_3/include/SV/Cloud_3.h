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

#ifndef CLOUD_3_H_
#define CLOUD_3_H_

#include <SV/Octree_3.h>
#include <boost/foreach.hpp>
#include <SV/memory.h>


namespace SV {

template < class Kernel , typename INT >
class Cloud_3 {
public:
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point_3;
  
  typedef SV::Octree_3<INT> Octree_3; 
  

public:
  Octree_3 offset0, offset1, offset2;
  int resolution() const {return offset0.resolution();}
  double scale() const {return offset0.scale();}
//  int resolution;
//  double scale;

public:
  // constructors
  Cloud_3(int resolution_ = 0) :
    offset0(resolution_), offset1(resolution_),offset2(resolution_)
  {
  }
  
  Cloud_3(const Octree_3& offset0_):offset0(offset0_){
    init_offsets();
  }

  void
  clear()
  {
    offset0.clear();
    offset1.clear();
    offset2.clear();
  }

  void init_offsets(){
//     std::cerr << "offset0 insertion timer:   " <<  offset0.insertion_timer.time() << std::endl;
//     std::cerr << "offset0 compression timer: " <<  offset0.compression_timer.time() << std::endl;
//     std::cerr << "offset0 erasing timer:     " <<  offset0.erase_timer.time() << std::endl;
//     std::cerr << "contains timer:           " <<  offset0.contains_timer.time() << std::endl;
  
    {
      std::cerr << "\n Memory before offset1 " << SV::memory_usage()*100 << " %\n" ;
      
      CGAL::Real_timer timer; 
      timer.start();
      offset1.clear();
      offset1 = offset0.offset();
      timer.stop();
      std::cerr << "Offset1 computation time: " <<  timer.time() << std::endl;
    }
    {
      std::cerr << "\n Memory before offset2 " << SV::memory_usage()*100 << " %\n" ;
         
      CGAL::Real_timer timer; 
      timer.start();
      offset2.clear();
      offset2 = offset1.offset();
      timer.stop();
      std::cerr << "Offset2 computation time: " <<  timer.time() << std::endl;
    }
    std::cerr << "\n Memory after offset2 " << SV::memory_usage()*100 << " %\n" ;

    std::cerr << "offset1.size() " <<  offset1.size() << std::endl; 
    std::cerr << "offset2.size() " <<  offset2.size() << std::endl; 
  }

  void
  compute_key(const Point_3& p, INT&x, INT&y, INT&z) const
  {
    double dx = p.x();
    double dy = p.y();
    double dz = p.z();
    
    //std::cerr << x << " " << y << " " << z << std::endl;
    //std::cerr << "resolution: " << resolution() << std::endl;
    //std::cerr << "scale:      " << scale() << std::endl;

    x = INT(dx * scale());
    y = INT(dy * scale());
    z = INT(dz * scale());
    // std::cerr << x << " " << y << " " << z << std::endl;
  }

  bool insert(const Point_3& p)
  {
    if(!offset1.empty()){
      offset1.clear();
      offset2.clear();
    }
    INT x, y, z;
    compute_key(p, x, y, z);
    return offset0.insert(x,y,z);
  }

  unsigned int
  size() const
  {
//     std::cerr << "max load_factor offset0:  " << offset0.set().max_load_factor() << std::endl; 
//     std::cerr << "    load factor offset0:  " << offset0.set().load_factor() << std::endl; 
//     std::cerr << "max load_factor offset1: " << offset1.set().max_load_factor() << std::endl; 
//     std::cerr << "    load factor offset1: " << offset1.set().load_factor() << std::endl; 
//     std::cerr << "max load_factor offset2: " << offset2.set().max_load_factor() << std::endl; 
//     std::cerr << "    load factor offset2: " << offset2.set().load_factor() << std::endl; 
//     std::cerr <<" ------------  "<< std::endl;
       
   
//    int max_size=0; 
//     for(int i=0; i <offset0.set().bucket_count();i++){
//       max_size=std::max(int(offset0.set().bucket_size(i)),max_size);
//     }
//     std::cerr<< "Max size bucket offset0: " << max_size << std::endl;   


    //std::cerr << (offset0.offset().size() == offset1.size()) << std::endl; 
    //std::cerr << (offset0.offset().offset().size() == offset2.size()) << std::endl; 
    return offset0.size();
  } 
  
  
  int
  classify(const Point_3& p) const
  {
    if(p.x()<0 || p.y()<0 || p.z()<0) return 3; 
    if(p.x()>1 || p.y()>1 || p.z()>1) return 3; 
    
    CGAL_precondition(!offset0.set().empty());
    CGAL_precondition(!offset1.set().empty());
    CGAL_precondition(!offset2.set().empty());
    
    INT x, y, z;
    compute_key(p, x, y, z);
    // std::cerr << "key " << x << " " << y << " " << z << std::endl;  
    return classify(x,y,z);
  }

  int
  classify(INT x, INT y, INT z) const
  {

    if (offset0.contains(x, y, z)) return 0;
    if (offset1.contains(x, y, z)) return 1;
#if SV_USE_TWO_OFFSETS
    if (offset2.contains(x, y, z)) return 2;
#endif
    return 3;
  }

};

} // namespace SV
#endif /* CLOUD_3_H_ */
