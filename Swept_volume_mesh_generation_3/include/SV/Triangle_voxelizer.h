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

#ifndef SV_TRIANGLE_VOXELIZER_H
#define SV_TRIANGLE_VOXELIZER_H

#include <CGAL/basic.h>
#include <CGAL/tuple.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <vecmath.h>

namespace SV { 

struct Triangle_voxelizer{




  static const int Bits[16];
  
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel; 
  typedef Kernel::Point_3 Point_3;
  typedef Kernel::Vector_3 Vector_3;

  typedef CGAL::cpp0x::tuple<int,int,int>                   Voxel;
  typedef CGAL::cpp0x::tuple<int,int,int,Point_3,Point_3> VB; 
 
 //typedef CGAL::cpp0x::tuple<int,int,int,Vector3d,Vector3d> VB; 



  static inline double fabs(double a){
    return a > 0 ? a : -a;
  }

//static inline void normal(const Point_3& a,const Point_3& b, const Point_3& c, Point_3& n);
    
  struct Hint{
    //const Vector3d& a,b,c;            // points of triangle 
    //Vector3d tMin,tMax;               // bbox of triangle   (3 directions)
    //std::vector<Vector3d> directions; // all other relevant directions 
		const Point_3 a,b,c;
		Point_3 tMin, tMax;
		std::vector<Point_3> directions;		

    std::vector<double> tmins;        // min value of triangle interval in direction
    std::vector<double> tmaxs;        // max value of triangle interval in direction
  
    //Hint(const Vector3d& a_, const Vector3d& b_, const Vector3d& c_);
    Hint(const Point_3& a_, const Point_3& b_, const Point_3& c_);
  };
    
  static int box_triangle_intersect( 
      const Point_3& bMin, const Point_3& bMax, const Hint& hint);
      //const Vector3d& bMin, const Vector3d& bMax, const Hint& hint);

 template <class OutputIterator> inline 
 OutputIterator recurse(
     int res, 
     int I ,int J ,int K,
     //const Vector3d& bmin,
     //const Vector3d& bmax,
     const Point_3& bmin,
     const Point_3& bmax,
     const Hint& hint, 
     OutputIterator oit) const {
   
   using CGAL::cpp0x::get; 

   if(res<=0) {
     //std::cerr << I<<" , " << J<<" , " << K << std::endl; 
     return *oit++=Voxel(I,J,K);
   }
   
   float B[][3]={{bmin.x(),(bmin.x()+bmax.x())*0.5,bmax.x()},
                 {bmin.y(),(bmin.y()+bmax.y())*0.5,bmax.y()},
                 {bmin.z(),(bmin.z()+bmax.z())*0.5,bmax.z()}};
#if 1
   res -=1; 
   int II = I|Bits[res];
   int JJ = J|Bits[res];
   int KK = K|Bits[res];
   
   //std::cerr << res << " , "<< I << " , " << J << " , " << K << " , "  << II << " , " << std::endl; 
   
#define SV_TEST_BOX(i,j,k)                                              \
   Point_3 bmin(B[0][i]  ,B[1][j]  ,B[2][k]);                          \
   Point_3 bmax(B[0][i+1],B[1][j+1],B[2][k+1]);                        \
   if(box_triangle_intersect(bmin,bmax,hint)){                          \
     if(!res){                                                           \
       *oit++=Voxel(i?II:I,j?JJ:J,k?KK:K);                              \
     }else{                                                             \
       oit = recurse(res, i?II:I,j?JJ:J,k?KK:K, bmin, bmax ,hint,oit);  \
     }                                                                  \
   }                                                                    \
   
   {SV_TEST_BOX(0,0,0);}                                                   
   {SV_TEST_BOX(0,0,1);}                                                   
   {SV_TEST_BOX(0,1,0);}                                                   
   {SV_TEST_BOX(0,1,1);}   
   {SV_TEST_BOX(1,0,0);}                                                   
   {SV_TEST_BOX(1,0,1);}                                                   
   {SV_TEST_BOX(1,1,0);}                                                   
   {SV_TEST_BOX(1,1,1);}

#undef SV_TEST_BOX

#else
   for(int i = 0; i <=1 ; i++){
     for(int j = 0; j <=1 ; j++){          
       for(int k = 0; k <=1 ; k++){
         Vector3d bmin(B[0][i]  ,B[1][j]  ,B[2][k]);
         Vector3d bmax(B[0][i+1],B[1][j+1],B[2][k+1]);
         if(box_triangle_intersect(bmin,bmax,hint)){ 
           oit = recurse(res-1, (I<<1)+i, (J<<1)+j, (K<<1)+k, bmin, bmax ,hint,oit); 
         }
       }
     }
   } 
#endif
   return oit; 
 }

  template <class OutputIterator>
  OutputIterator operator()(
      int res,
//      const Vector3d& a, const Vector3d& b, const Vector3d& c,
//      const Vector3d& bmin, const Vector3d& bmax,
      const Point_3& A, const Point_3& B, const Point_3& C,
      const Point_3& bmin, const Point_3& bmax,
      OutputIterator oit) const {
	
	//	Vector3d va(A.x(),A.y(),A.z());	
	//	Vector3d vb(B.x(),B.y(),B.z());	
	//	Vector3d vc(C.x(),C.y(),C.z());	

	//	Vector3d vbmin(bmin.x(),bmin.y(),bmin.z());	
	//	Vector3d vbmax(bmax.x(),bmax.y(),bmax.z());	
		
    Hint hint(A,B,C);
    return recurse(res,0,0,0,bmin,bmax,hint,oit);
  }
  
  
};

} // namespace SV 
#endif // SV_TRIANGLE_VOXELIZER_H
