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


#include <SV/Triangle_voxelizer.h> 
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>


namespace SV { 

  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel; 
  typedef Kernel::Point_3 Point_3;
  typedef Kernel::Vector_3 Vector_3;



const int Triangle_voxelizer::Bits[16]={1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768};


//static inline void Triangle_voxelizer::normal(const Point_3& a,const Point_3& b, const Point_3& c, Point_3& n){
//
//		Vector_3 ab = b - a;			
//		Vector_3 ac = c - a;			
//		Vector_3 r = CGAL::cross_product(ab,ac);
//		Point_3 n(r.x(),r.y(),r.z());
//}


Triangle_voxelizer::Hint::Hint(const Point_3& a_, const Point_3& b_, const Point_3& c_)
  :a(a_),b(b_),c(c_){

  // bbox of triangle 
  tMin = Point_3(
      std::min(std::min(a.x(),b.x()),c.x()),
      std::min(std::min(a.y(),b.y()),c.y()),
      std::min(std::min(a.z(),b.z()),c.z()));
  tMax = Point_3(
      std::max(std::max(a.x(),b.x()),c.x()),
      std::max(std::max(a.y(),b.y()),c.y()),
      std::max(std::max(a.z(),b.z()),c.z()));
      
  { // generate all relevant directions 
    
    Point_3 d; 
    double tA,tB,tC;
    
#define SV_INIT_DIRECTION                                       \
    directions.push_back(d);                                    \
    tA = ( d.x()*(a.x()) + d.y()*(a.y()) + d.z()*(a.z()) );     \
    tB = ( d.x()*(b.x()) + d.y()*(b.y()) + d.z()*(b.z()) );     \
    tC = ( d.x()*(c.x()) + d.y()*(c.y()) + d.z()*(c.z()) );     \
    tmins.push_back(std::min(std::min(tA,tB),tC));              \
    tmaxs.push_back(std::max(std::max(tA,tB),tC));              \
                                                                \
    // the normal of the triangle
    //d = normal(a,b,c) ;
		Vector_3 ppab = b - a;			
		Vector_3 ppac = c - a;			
		Vector_3 ppr = CGAL::cross_product(ppab,ppac);
		Point_3 ppn(ppr.x(),ppr.y(),ppr.z());
  	d = Point_3(ppn.x(), ppn.y(), ppn.z()); 
	 
		SV_INIT_DIRECTION
          
      //d = ( b-a ) x (1 0 0)	
      //d.x() = 0.0;  d.y() = b.z()-a.z();  d.z() = -(b.y()-a.y()); 
      d = Point_3( 0.0,  b.z()-a.z(),  - (b.y()-a.y()) ) ; 
    SV_INIT_DIRECTION
      
      //d = ( c-b ) x (1 0 0)	
      //d.x() = 0.0;  d.y() = c.z()-b.z();  d.z() = -(c.y()-b.y()); 
      d = Point_3( 0.0,  c.z()-b.z(),   -(c.y()-b.y()) ); 
    SV_INIT_DIRECTION

      //d = ( a-c ) x (1 0 0)	
      //d.x() = 0.0;  d.y() = a.z()-c.z();  d.z() = -(a.y()-c.y()); 
      d = Point_3( 0.0,  a.z()-c.z(),  -(a.y()-c.y()) ); 
    SV_INIT_DIRECTION
      
      //d = ( b-a ) x (0 1 0)	
      //d.x() = -(b.z()-a.z());  d.y() = 0.0;  d.z() = b.x()-a.x(); 
      d = Point_3( -(b.z()-a.z()),  0.0,  b.x()-a.x() ); 
    SV_INIT_DIRECTION
      
      //d = ( c-b ) x (0 1 0)	
      //d.x() = -(c.z()-b.z());  d.y() = 0.0;  d.z() = c.x()-b.x();
      d = Point_3 (-(c.z()-b.z()),  0.0,  c.x()-b.x() );
    SV_INIT_DIRECTION
      
      //d = ( a-c ) x (0 1 0)	
      //d.x() = -(a.z()-c.z());  d.y() = 0.0;  d.z() = a.x()-c.x(); 
      d = Point_3( -(a.z()-c.z()), 0.0, a.x()-c.x() ); 
    SV_INIT_DIRECTION
      
      //d = ( b-a ) x (0 0 1)	
      //d.x() = b.y()-a.y(); d.y() = -(b.x()-a.x()); d.z() = 0.0;  
      d = Point_3( b.y()-a.y(),  -(b.x()-a.x()),  0.0 );  

    SV_INIT_DIRECTION
      
      //d = ( c-b ) x (0 0 1)	
      //d.x() = c.y()-b.y(); d.y() = -(c.x()-b.x()); d.z() = 0.0;  
      d = Point_3( c.y()-b.y(),  -(c.x()-b.x()), 0.0 );  
    SV_INIT_DIRECTION
      
      //d = ( a-c ) x (0 0 1)	
      //d.x() = a.y()-c.y(); d.y() = -(a.x()-c.x()); d.z() = 0.0; 
      d = Point_3( a.y()-c.y(),  -(a.x()-c.x()), 0.0 ) ;
    SV_INIT_DIRECTION
      }
}

#if 1
int Triangle_voxelizer::box_triangle_intersect(
    //const Vector3d& bMin, const Vector3d& bMax, const Hint& hint){
    const Point_3& bMin, const Point_3& bMax, const Hint& hint){
  
  //const Vector3d& a = hint.a; 
  //const Vector3d& b = hint.b; 
  //const Vector3d& c = hint.c; 
  
  // Dreiecks-AABB mit Box schneiden:
  //const Vector3d& tMin = hint.tMin; 
  //const Vector3d& tMax = hint.tMax; 
  const Point_3& tMin = hint.tMin; 
  const Point_3& tMax = hint.tMax; 
  
  if( tMin.x() > bMax.x() ) return 0;
  if( tMin.y() > bMax.y() ) return 0;
  if( tMin.z() > bMax.z() ) return 0;
  if( tMax.x() < bMin.x() ) return 0;
  if( tMax.y() < bMin.y() ) return 0;
  if( tMax.z() < bMin.z() ) return 0;
  
  Point_3 o = Point_3((bMax.x()+bMin.x())*0.5, 
											(bMax.y()+bMin.y())*0.5,
											(bMax.z()+bMin.z())*0.5);
  Vector_3 lv = (bMax-bMin)*0.5;
  Point_3 l(lv.x(), lv.y(), lv.z());
  Point_3 d; 
  double cbox, lbox;
  
#define SV_TEST_AXIS                                                    \
  cbox =  d.x()*o.x() + d.y()*o.y() + d.z()*o.z();                      \
  lbox =  fabs(d.x()*l.x()) + fabs(d.y()*l.y()) + fabs(d.z()*l.z());    \
  if( cbox+lbox < tmin ) return 0;                                      \
  if( tmax < cbox-lbox ) return 0;                                      \
                                                                        \
  // for each direction d do SV_TEST_AXIS
  for(int i = 0; i < hint.directions.size();i++){
    const Point_3& d= hint.directions[i];
    const double& tmin = hint.tmins[i];
    const double& tmax = hint.tmaxs[i];
    SV_TEST_AXIS; 
  }
  return 1;  
}
#endif

} // namespace SV 
