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

#ifndef SV_IO_H
#define SV_IO_H

#include <iostream> 
#include <fstream> 
#include <boost/foreach.hpp> 
//#include <vecmath.h>


#include <CGAL/Exact_predicates_exact_constructions_kernel.h>


namespace SV{

    typedef CGAL::Exact_predicates_inexact_constructions_kernel    Kernel; 
	  typedef typename Kernel::Point_3 Point_3;
    typedef typename Kernel::Aff_transformation_3 AT_3;

template < class C3t3, class Stream> 
Stream& save_as_off(const C3t3& c3t3, Stream& os)
{ 
  AT_3 I(1,0,0,0, 0,1,0,0, 0,0,1,0, 1);
  //I.scale(1,1,1);
  return save_as_off(c3t3, os, I); 
}

template < class C3t3, class Stream>  
Stream& save_as_off(const C3t3& c3t3, Stream& os, const AT_3& MT){
 


  //	std::cerr<<MT[0][0]<<" "<<MT[1][0]<<" "<<MT[2][0]<<std::endl; 
  //	std::cerr<<MT[1][0]<<" "<<MT[1][1]<<" "<<MT[1][2]<<std::endl; 
  //	std::cerr<<MT[2][0]<<" "<<MT[2][1]<<" "<<MT[2][2]<<std::endl; 
  Point_3 e1 = MT.transform(Point_3(1,0,0));
  Point_3 e2 = MT.transform(Point_3(0,1,0));
  Point_3 e3 = MT.transform(Point_3(0,0,1));
  std::cerr << "Witing to off with back_trafo:"  << std::endl;
  //std::cerr << e1 << " " << e2 << " " << e3 << std::endl;


  typedef typename C3t3::Facet_iterator Facet_iterator ;
  typedef typename C3t3::Facet Facet;
  typedef typename C3t3::Cell_iterator Cell_iterator ;
  typedef typename C3t3::Cell_handle Cell_handle;
  typedef typename C3t3::Vertex_handle Vertex_handle;
  typedef typename C3t3::Subdomain_index Subdomain_index;

  // collect all vertices 
  std::set<Vertex_handle> vertex_set; 
  for(Facet_iterator fit = c3t3.facets_begin(); 
      fit != c3t3.facets_end(); fit++){
    vertex_set.insert(fit->first->vertex((1+fit->second)%4));
    vertex_set.insert(fit->first->vertex((2+fit->second)%4));
    vertex_set.insert(fit->first->vertex((3+fit->second)%4));
  }

  // write header for OFF file 
  os << "OFF" << std::endl;
  os << vertex_set.size() << " "
     << c3t3.number_of_facets() << " "
     << "0" << std::endl; 


  // write all vertecies and assign index to each vertex 
  std::map<Vertex_handle,int> vertex_index_map; 
  int vcount = 0;
  for(typename std::set<Vertex_handle>::iterator vit = vertex_set.begin(); 
      vit != vertex_set.end(); vit++){
    vertex_index_map[*vit]=vcount;

   // Vector3d v = Vector3d( (*vit)->point().x(),  
   //     (*vit)->point().y(),
   //     (*vit)->point().z());
    Point_3 v = Point_3( (*vit)->point().x(),  
        (*vit)->point().y(),
        (*vit)->point().z());
    //std::cerr << v.x() << " " <<v.y()<<" " <<v.z()<<"  --> " ;
    v = MT.transform(v);
    //v.x() *= 0.001;
    //v.y() *= 0.001;
    //v.z() *= 0.001;
    //std::cerr << v.x() << " " <<v.y()<<" " <<v.z()<<std::endl;
    os  << v.x() << " "
        << v.y() << " " 
        << v.z() << std::endl;
    //	os  << (*vit)->point().x() 
    //			<< (*vit)->point().y() << " " 
    //			<< (*vit)->point().z() << std::endl;
    vcount++;
  }

  typedef std::pair<Cell_iterator,int> Bfacet; 
  std::vector<Bfacet> bfacets; //  boundary facets 

  for(Cell_iterator cit = c3t3.cells_begin(); cit != c3t3.cells_end(); cit++){
    for(int i = 0; i<4;i++){
      if(c3t3.is_in_complex(cit,i)){
        bfacets.push_back(std::make_pair(cit,i));
      }
    }
  }

  BOOST_FOREACH(Bfacet f, bfacets){
    Vertex_handle v0 = f.first->vertex((0+f.second)%4);
    Vertex_handle v1 = f.first->vertex((1+f.second)%4);
    Vertex_handle v2 = f.first->vertex((2+f.second)%4);
    Vertex_handle v3 = f.first->vertex((3+f.second)%4);

    if(CGAL::orientation(v1->point(),v2->point(),v3->point(),v0->point()) == CGAL::NEGATIVE){
      os << 3 << " "
         <<  vertex_index_map[f.first->vertex((1+f.second)%4)] << " "
         <<  vertex_index_map[f.first->vertex((2+f.second)%4)] << " "
         <<  vertex_index_map[f.first->vertex((3+f.second)%4)] << std::endl;
    }else{
      os << 3 << " "
         <<  vertex_index_map[f.first->vertex((3+f.second)%4)] << " "
         <<  vertex_index_map[f.first->vertex((2+f.second)%4)] << " "
         <<  vertex_index_map[f.first->vertex((1+f.second)%4)] << std::endl;
    }

  } 

  return os; 
}

}//namespace SV 

#endif // SV_IO_H
