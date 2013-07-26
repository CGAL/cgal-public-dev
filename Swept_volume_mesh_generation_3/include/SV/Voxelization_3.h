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

#ifndef SV_VOXELIZATION_3_H
#define SV_VOXELIZATION_3_H


#ifdef _OPENMP
#ifndef USE_OMP
#define USE_OMP 0
#include <omp.h>
#endif
#endif 

#include <vector>
#include <set>
#include <list>

#include <tr1/tuple>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
//#include </usr/include/GL/glut.h>
//#include "vecmath.h"
#include <stdlib.h>
#include <algorithm>


#include <CGAL/tuple.h>
#include <CGAL/Real_timer.h> 
#include <CGAL/Random.h> 

#include <SV/HashSet_Voxel.h>
#include <SV/Octree_3.h>
#include <SV/memory.h>
#include <SV/Triangle_voxelizer.h>

#include <boost/foreach.hpp>
#include <boost/optional.hpp>

#ifndef USE_EXACT_CULLING
#define USE_EXACT_CULLING 1
#endif 

#ifndef DONT_USE_CULLING_PATCH 
#define DONT_USE_CULLING_PATCH 0
#endif
#ifndef DONT_USE_CULLING_FACET 
#define DONT_USE_CULLING_FACET 0
#endif
#ifndef DONT_USE_CULLING_FACET_NORMAL 
#define DONT_USE_CULLING_FACET_NORMAL 0
#endif

#if USE_EXACT_CULLING
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#endif

namespace SV { 

class Voxelization_3
{
public:
  typedef unsigned int uint;
  typedef std::pair<uint, uint> edge;
  typedef std::pair<edge, uint> wing;
  typedef std::vector<uint> wingedEdge;

public: 
  typedef short INT; 
  typedef CGAL::cpp0x::tuple<INT,INT,INT> Voxel;
  typedef std::set<Voxel> VoxelSet;
  typedef std::vector<Voxel> VoxelVec;
  typedef SV::HashSet_Voxel<INT> VoxelHashSet;
  typedef SV::Octree_3<INT> Octree_3;
  typedef std::vector<Voxel> Voxels;

  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel; 
  typedef Kernel::Aff_transformation_3 AT_3;
  typedef Kernel::Point_3 Point_3;

public:
  typedef int UINT; 		

  //std::vector<Vector3d> P;     // Points 
 	//std::vector<Matrix4d> Track; 

  std::vector<Point_3> P;     // Points 
  std::vector<int> ind;        // index chain 
  std::vector<AT_3> Track;

  std::vector<wing> Wings; 
  std::vector<wingedEdge> wingedEdges;
  std::vector<edge> non_manifold_edges; 
  
  int m_resolution;
  int m_downstep;
  int m_num_threads; 
  int m_epsilon;
  
  mutable int cull_candidate;
  mutable int cull_success;

  //Matrix4d m_back_transformation;
  AT_3 m_back_transformation;

private:
  mutable boost::optional<Octree_3> opt_volume;  
  
public:
  mutable CGAL::Real_timer timer_total; 
  mutable CGAL::Real_timer timer_hull; 
  mutable CGAL::Real_timer timer_voxelize; 
  mutable CGAL::Real_timer timer_voxelize_lookup; 
  mutable CGAL::Real_timer timer_voxelize_hash_lookup; 
  mutable CGAL::Real_timer timer_hull_lookup; 
  mutable CGAL::Real_timer timer_bt_intersect; 
  mutable CGAL::Real_timer timer_box_gen; 
  
  mutable CGAL::Real_timer timer_emit_boundary;
  mutable CGAL::Real_timer timer_compute_hull;
  mutable CGAL::Real_timer timer_fill_hull;


private:
  
  int& resolution(){return m_resolution;}
  const int& resolution() const {return m_resolution;}

public:
  Octree_3& volume(){
    generate_voxelization();
    return *opt_volume;
  }
  const Octree_3& volume() const {
    generate_voxelization();
    return *opt_volume;
  }
  
  //Matrix4d& back_transformation(){return m_back_transformation; }
  //const Matrix4d& back_transformation() const {return m_back_transformation; }
  AT_3& back_transformation(){return m_back_transformation; }
  const AT_3& back_transformation() const {return m_back_transformation; }

public:
 
  Voxelization_3(int n = 1):m_resolution(n){};
  
  Voxelization_3(
      const char* filename_track, const char* filename_off, 
      int n = 10, int downstep = 1 , int num_threads = 1 );
  
  Voxelization_3(
      const char* filename_track, const char* filename_off, 
      double eps, int downstep = 1 , int num_threads = 1 );
  
  void clear(){   
    //P = std::vector<Vector3d>(); 
    P = std::vector<Point_3>(); 

    ind = std::vector<int>();   
    Track = std::vector< AT_3 >();
    //Track = std::vector< Matrix4d >();
    Wings = std::vector<wing>();
    wingedEdges = std::vector<wingedEdge>(); 
    opt_volume = boost::optional<Octree_3>();
  }
  
  void printTimings() const {
    std::cerr << "# Voxelization_3 total  : " << timer_total.time() << std::endl;
    std::cerr << "-> time voxelization    : " << timer_voxelize.time() << std::endl
              << "---> gen boxes            : " << timer_box_gen.time() << std::endl
              << "---> box/tri intersect    : " << timer_bt_intersect.time() << std::endl
              << "-> ohull0 insertion     : " << timer_voxelize_hash_lookup.time() << std::endl;
    std::cerr << "-> time hull            : " << timer_hull.time() << std::endl;
  }

private:
  void generate_voxelization() const;
  void initialize_winged_edge();
  void normalize(std::vector<Point_3>&, std::vector< AT_3 >&);
  void computeResAndScale(std::vector<Point_3>& P, std::vector< AT_3 >& Track,double eps );
  void genObjectList(); 	
  void compress(int current_res, bool compute_inner_offset) const;
  
  void loadFile(const char* filename, std::vector< AT_3 >& T) ;
  void loadPathFile(const char* filename, std::vector< AT_3 >& T) ;
	void LoadOffFile(const char* as_FileName, std::vector<Point_3>& ak_Vertices, 
	std::vector<int>& ak_Indices);
 
  double pointsOnSameSide(Point_3 a, Point_3 b, Point_3 n, Point_3 p) const;
  //bool pointsOnSameSide(Vector3d a, Vector3d b, Vector3d c, Vector3d n, Vector3d p) const ;
  void emit_boundary(int current_res, Octree_3& octree) const;  
  void emit_boundary_serial(int current_res, Octree_3& octree) const;    
#if USE_OMP
  void emit_boundary_omp(int current_res, Octree_3& octree) const;  
#endif
  
public:

  template  <class OutputIterator> 
  OutputIterator get_voxelization( 
      int res,
      const float ax, const float ay, const float az, 
      const float bx, const float by, const float bz,
      const float cx, const float cy, const float cz,
      // const int n,
      OutputIterator oit) const {    using CGAL::cpp0x::get;
    
    typedef CGAL::cpp0x::tuple<int,int,int,Point_3,Point_3> VB; 
    
    Point_3 a(ax,ay,az);
    Point_3 b(bx,by,bz);
    Point_3 c(cx,cy,cz); 
    Point_3 bmin(0,0,0); 
    Point_3 bmax(1,1,1);    
    
    return Triangle_voxelizer()(res,a,b,c,bmin,bmax,oit);
  }
  
  template <class OuputIterator > 
  OuputIterator emitQuad( 
      int res, const Point_3& a0, const Point_3& a1, const Point_3& b0, const Point_3& b1,
      OuputIterator oit) const {
    oit = get_voxelization(
        res,
	a0.x(),a0.y(),a0.z(),
        b0.x(),b0.y(),b0.z(),
        b1.x(),b1.y(),b1.z(),
        oit);
    
    oit = get_voxelization(
        res, 
	a0.x(),a0.y(),a0.z(),
        b1.x(),b1.y(),b1.z(),
        a1.x(),a1.y(),a1.z(),
        oit);
    return oit; 
  }
  
  template <class OuputIterator > 
  OuputIterator emitFacet(
      int res, const Point_3& a, const Point_3& b, const Point_3& c, OuputIterator oit) const {
    return get_voxelization(
        res, 
	a.x(),a.y(),a.z(),
        b.x(),b.y(),b.z(),
        c.x(),c.y(),c.z(),
        oit);
  }
  
 
#if USE_EXACT_CULLING
  inline bool needPatch(
       AT_3 D0, AT_3 D1, wingedEdge we,
      Point_3 &a, Point_3 &b, Point_3 &c, Point_3 &d) const {    


    cull_candidate += 2;
    //Vector3d a0 = D0*P[we.at(0)]; 	
    //Vector3d a1 = D1*P[we.at(0)]; 	
    //Vector3d b0 = D0*P[we.at(1)]; 	
    //Vector3d b1 = D1*P[we.at(1)]; 	

    //Vector3d c0 = D0*P[we.at(2)]; 
    //Vector3d c1 = D1*P[we.at(2)]; 	
    //Vector3d c2 = D0*P[we.at(3)]; 	
    //Vector3d c3 = D1*P[we.at(3)]; 	

	Point_3 A0 = D0.transform(P[we.at(0)]);
	Point_3 A1 = D1.transform(P[we.at(0)]);
	Point_3 B0 = D0.transform(P[we.at(1)]);
	Point_3 B1 = D1.transform(P[we.at(1)]);
     
	Point_3 C0 = D0.transform(P[we.at(2)]);
	Point_3 C1 = D1.transform(P[we.at(2)]);
	Point_3 C2 = D0.transform(P[we.at(3)]);
	Point_3 C3 = D1.transform(P[we.at(3)]);

    
    //Point_3 A0(a0.x(),a0.y(),a0.z()); 	
    //Point_3 A1(a1.x(),a1.y(),a1.z()); 
    //Point_3 B0(b0.x(),b0.y(),b0.z()); 	
    //Point_3 B1(b1.x(),b1.y(),b1.z());  
    //Point_3 C0(c0.x(),c0.y(),c0.z()); 	
    //Point_3 C1(c1.x(),c1.y(),c1.z());  
    //Point_3 C2(c2.x(),c2.y(),c2.z()); 	
    //Point_3 C3(c3.x(),c3.y(),c3.z());  
    
#define SV_ON_SAME_SIDE(P1,P2,P3,Q1,Q2)         \
    (CGAL::orientation(P1,P2,P3,Q1)==            \
      CGAL::orientation(P1,P2,P3,Q2))

//#define SV_ON_SAME_SIDE(P1,P2,P3,Q1,Q2)         \
//    (CGAL::orientation(P1,P2,P3,Q1)==            \
//      CGAL::orientation(P1,P2,P3,Q2)&& \
//      CGAL::orientation(P1,P2,P3,Q2)!=CGAL::COPLANAR) 
    
    // we are supposed to use the triangle pair that formes the largest angle 
    // this it done by chosing the 
    
    if(CGAL::compare_distance(A1,B0,A0,B1)==CGAL::SMALLER){
      //a = b0; b = a0; c = b1; d = a1;
      a = B0; b = A0; c = B1; d = A1;
      if(!SV_ON_SAME_SIDE(A0,A1,B0,C0,C1)) { return true; } // emits T(a0,b0,b1) and T(a0,b1,a1)	
      if(!SV_ON_SAME_SIDE(A1,B0,B1,C0,C1)) { return true; }	
      if(!SV_ON_SAME_SIDE(A0,A1,B0,C2,C3)) { return true; }	
      if(!SV_ON_SAME_SIDE(A1,B0,B1,C2,C3)) { return true; }	

      if( SV_ON_SAME_SIDE(A0,A1,B0,C0,C2)) { return true; }	
      if( SV_ON_SAME_SIDE(A1,B0,B1,C0,C2)) { return true; }	
      if( SV_ON_SAME_SIDE(A0,A1,B0,C1,C3)) { return true; }	
      if( SV_ON_SAME_SIDE(A1,B0,B1,C1,C3)) { return true; }
    }else{
      //a = a0; b = b0; c = a1; d =b1;  
      a = A0; b = B0; c = A1; d =B1;  
      if(!SV_ON_SAME_SIDE(A0,A1,B1,C0,C1)) { return true; } // emits T(1,3,4) and T(1,4,2)	
      if(!SV_ON_SAME_SIDE(A0,B0,B1,C0,C1)) { return true; } //   T (a0,b0,a1) and T(a0,b1,a1))
      if(!SV_ON_SAME_SIDE(A0,A1,B1,C2,C3)) { return true; }	
      if(!SV_ON_SAME_SIDE(A0,B0,B1,C2,C3)) { return true; }	

      if( SV_ON_SAME_SIDE(A0,A1,B1,C0,C2)) { return true; }	
      if( SV_ON_SAME_SIDE(A0,B0,B1,C0,C2)) { return true; }	
      if( SV_ON_SAME_SIDE(A0,A1,B1,C1,C3)) { return true; }	
      if( SV_ON_SAME_SIDE(A0,B0,B1,C1,C3)) { return true; }	
    }
#undef SV_ON_SAME_SIDE


#if DONT_USE_CULLING_PATCH
	return true;	
#endif
    // code reaches here: Quad discarded!
    cull_success += 2;
    return false; 
  }
  
  inline bool needFacet(uint i, uint j, Point_3 &A, Point_3 &B, Point_3 &C) const {

		
    cull_candidate++;
    //const Matrix4d& D   = Track[i];
    //const Matrix4d& D_1 = Track[i-1];
    //const Matrix4d& D1  = Track[i+1];
	AT_3 D   = Track[i];
	AT_3 D_1 = Track[i-1];
	AT_3 D1  = Track[i+1];

    A   = D.transform(P[ind[j]]  ) ;
    B   = D.transform(P[ind[j+1]]) ;
    C   = D.transform(P[ind[j+2]]) ;
    //Vector3d a_1 = D_1 *P[ind[j]]   ;
    //Vector3d b_1 = D_1 *P[ind[j+1]] ;
    //Vector3d c_1 = D_1 *P[ind[j+2]] ;
    //Vector3d a1  = D1  *P[ind[j]]   ;
    //Vector3d b1  = D1  *P[ind[j+1]] ;
    //Vector3d c1  = D1  *P[ind[j+2]] ;

    Point_3 A_1 = D_1.transform( P[ind[j]]  ) ;
    Point_3 B_1 = D_1.transform( P[ind[j+1]]) ;
    Point_3 C_1 = D_1.transform( P[ind[j+2]]) ;
    Point_3 A1  =  D1.transform( P[ind[j]]  ) ;
    Point_3 B1  =  D1.transform( P[ind[j+1]]) ;
    Point_3 C1  =  D1.transform( P[ind[j+2]]) ;


    //typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel; 
    //typedef Kernel::Point_3 Point_3;
     
    //Point_3 A(a.x(),a.y(),a.z()); 
    //Point_3 B(b.x(),b.y(),b.z());  
    //Point_3 C(c.x(),c.y(),c.z());  

    //Point_3 A1(a1.x(),a1.y(),a1.z()); 
    //Point_3 B1(b1.x(),b1.y(),b1.z());  
    //Point_3 C1(c1.x(),c1.y(),c1.z());  

    //Point_3 A_1(a_1.x(),a_1.y(),a_1.z()); 
    //Point_3 B_1(b_1.x(),b_1.y(),b_1.z());  
    //Point_3 C_1(c_1.x(),c_1.y(),c_1.z());  


#if DONT_USE_CULLING_FACET
	return true;	
#endif


    // for orientable meshes (correct normals)
    // if triangle normal points inwards for at least one prism (of timestep +/-)
    // and prism is well formed - we can cull
    if( CGAL::orientation(A,B,C,A1) == CGAL::POSITIVE &&
        CGAL::orientation(A,B,C,B1) == CGAL::POSITIVE &&
        CGAL::orientation(A,B,C,C1) == CGAL::POSITIVE ){
      cull_success++; return false;
    }
    if( CGAL::orientation(A,B,C,A_1) == CGAL::POSITIVE &&
        CGAL::orientation(A,B,C,B_1) == CGAL::POSITIVE &&
        CGAL::orientation(A,B,C,C_1) == CGAL::POSITIVE ){
      cull_success++; return false;
    } 
     
    // otherwise triangle is local boundary or something very strange happens (i.e. both
    // prisms are malformed)
    
    return true; 
  }
 
  
  //When mesh is oriented for begin and end config
  
  inline bool needFacetForwardNormal( uint i0, uint i1, uint j, Point_3 &A, Point_3 &B, Point_3 &C) const {  
    cull_candidate++;
   // const Matrix4d& D   = Track[i0];
   // const Matrix4d& D1  = Track[i1];
    AT_3 D   = Track[i0];
    AT_3 D1  = Track[i1];
    A   = D.transform(P[ind[j]]  ) ;
    B   = D.transform(P[ind[j+1]]) ;
    C   = D.transform(P[ind[j+2]]) ;
    Point_3 A1  = D1.transform(P[ind[j]]  ) ;
    Point_3 B1  = D1.transform(P[ind[j+1]]) ;
    Point_3 C1  = D1.transform(P[ind[j+2]]) ;


    //typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel; 
    //typedef Kernel::Point_3 Point_3;
    
    //Point_3 A(a.x(),a.y(),a.z()); 
    //Point_3 B(b.x(),b.y(),b.z());  
    //Point_3 C(c.x(),c.y(),c.z());  
    //
    //Point_3 A1(a1.x(),a1.y(),a1.z()); 
    //Point_3 B1(b1.x(),b1.y(),b1.z());  
    //Point_3 C1(c1.x(),c1.y(),c1.z());  


#if DONT_USE_CULLING_FACET_NORMAL
	return true;	
#endif
        
    if( CGAL::orientation(A,B,C,A1) == CGAL::POSITIVE &&
        CGAL::orientation(A,B,C,B1) == CGAL::POSITIVE &&
        CGAL::orientation(A,B,C,C1) == CGAL::POSITIVE ){
      cull_success++; return false;
    }
    
    return true;
  }
#endif
#if 0
  inline bool needPatch(
      const Matrix4d& D0, const Matrix4d& D1, wingedEdge we, 
      Vector3d& a, Vector3d& b, Vector3d& c,Vector3d& d ) const {    
    cull_candidate += 2;
    Vector3d a0 = D0*P[we.at(0)]; 	
    Vector3d a1 = D1*P[we.at(0)]; 	
    Vector3d b0 = D0*P[we.at(1)]; 	
    Vector3d b1 = D1*P[we.at(1)]; 	

    Vector3d c0 = D0*P[we.at(2)]; 
    Vector3d c1 = D1*P[we.at(2)]; 	
    Vector3d c2 = D0*P[we.at(3)]; 	
    Vector3d c3 = D1*P[we.at(3)]; 	
    
    Vector3d n0, n1, p0, p1;
    Vector3d n0b, n1b, p0b, p1b;
    
    n0 =  (a0-a1)^(b1-a1);  n0.normalize(); p0 = a1;
    n1 = -(a0-b0)^(b1-b0);  n1.normalize(); p1 = b0;
    n0b = -(a1-a0)^(b0-a0); n0b.normalize(); p0b = a0;
    n1b =  (a1-b1)^(b0-b1); n1b.normalize(); p1b = b1;
    
    // epsilon for save culling:
    double eps_cull = 0;//1E-5;
    // choose diagonal with lower dihedral angle 
    // = greater dot-product of normals 
    if (n0*n1 < n0b*n1b){ // diagonal is a0-b1
      n0=n0b;
      n1=n1b;
      p0=p0b;
      p1=p1b;
      
      a = b0; b = a0; c = b1; d = a1;
#if DONT_USE_CULLING_PATCH
    return true;
#endif
      if( pointsOnSameSide(c0,c1,n0,p0) < eps_cull ) { return true; } //  emit T(a0,b0,b1) and T(a0,b1,a1)	
      if( pointsOnSameSide(c0,c1,n1,p1) < eps_cull ) { return true; }	
      if( pointsOnSameSide(c2,c3,n0,p0) < eps_cull ) { return true; }	
      if( pointsOnSameSide(c2,c3,n1,p1) < eps_cull ) { return true; }	

      if( pointsOnSameSide(c0,c2,n0,p0) > -eps_cull ) { return true; }	
      if( pointsOnSameSide(c0,c2,n1,p1) > -eps_cull ) { return true; }	
      if( pointsOnSameSide(c1,c3,n0,p0) > -eps_cull ) { return true; }	
      if( pointsOnSameSide(c1,c3,n1,p1) > -eps_cull ) { return true; }	
      // degenerate situation?:
      if(n0*n1 < 0) return true;

    }else{
      a = a0; b = b0; c = a1; d =b1;  
#if DONT_USE_CULLING_PATCH
    return true;
#endif
      if( pointsOnSameSide(c0,c1,n0,p0) < eps_cull) { return true; } // emit T(1,3,4) and T(1,4,2)	
      if( pointsOnSameSide(c0,c1,n1,p1) < eps_cull) { return true; } //   T (a0,b0,a1) and T(b0,b1,a1))
      if( pointsOnSameSide(c2,c3,n0,p0) < eps_cull) { return true; }	
      if( pointsOnSameSide(c2,c3,n1,p1) < eps_cull) { return true; }	

      if( pointsOnSameSide(c0,c2,n0,p0) > -eps_cull) { return true; }	
      if( pointsOnSameSide(c0,c2,n1,p1) > -eps_cull) { return true; }	
      if( pointsOnSameSide(c1,c3,n0,p0) > -eps_cull) { return true; }	
      if( pointsOnSameSide(c1,c3,n1,p1) > -eps_cull) { return true; }	
      // degenerate situation?:
      if(n0*n1 < 0.0) return true;
    }

    // code reaches here: Quad discarded!
    cull_success += 2;
    return false;                                                                           // return 0 => discard
  }
  
  inline bool needFacet(uint i, uint j, Vector3d& a, Vector3d& b, Vector3d& c ) const {
		
    cull_candidate++;
    const Matrix4d& D   = Track[i];
    const Matrix4d& D_1 = Track[i-1];
    const Matrix4d& D1  = Track[i+1];
    a   = D   *P[ind[j]]   ;
    b   = D   *P[ind[j+1]] ;
    c   = D   *P[ind[j+2]] ;
#if DONT_USE_CULLING_FACET
return true;
#endif
    Vector3d a_1 = D_1 *P[ind[j]]   ;
    Vector3d b_1 = D_1 *P[ind[j+1]] ;
    Vector3d c_1 = D_1 *P[ind[j+2]] ;
    Vector3d a1  = D1  *P[ind[j]]   ;
    Vector3d b1  = D1  *P[ind[j+1]] ;
    Vector3d c1  = D1  *P[ind[j+2]] ;
    Vector3d n = (b-a)^(c-a); n.normalize();
    
    // for safe face culling
    
    double eps_cf = 0;//1E-5;
    
 
    // for orientable meshes (correct normals)
    // if triangle normal points inwards for at least one prism (of timestep +/-)
    // and prism is well formed - we can cull
    if ((a1 -a)*n > eps_cf && (b1 -b)*n > eps_cf && (c1 -c)*n > eps_cf  ) { return false; /*cull_success++; return oit;*/} // return 0 => discard tri
    if ((a_1-a)*n > eps_cf && (b_1-b)*n > eps_cf && (c_1-c)*n > eps_cf  ) { return false; /*cull_success++; return oit;*/}
    // otherwise triangle is local boundary or something very strange happens (i.e. both
    // prisms are malformed)
    return true; /*emitFacet(res,a,b,c, oit); */                                                            // return 1 => emit facet      
  }
 
  
  //When mesh is oriented for begin and end config
  
  inline bool needFacetForwardNormal(uint i0, uint i1, uint j, Vector3d &a, Vector3d &b, Vector3d &c) const {  
    cull_candidate++;
    const Matrix4d& D   = Track[i0];
    const Matrix4d& D1  = Track[i1];
    a   = D   *P[ind[j]]   ;
    b   = D   *P[ind[j+1]] ;
    c   = D   *P[ind[j+2]] ;
#if DONT_USE_CULLING_FACET_NORMAL
    return true;
#endif
    Vector3d a1  = D1  *P[ind[j]]   ;
    Vector3d b1  = D1  *P[ind[j+1]] ;
    Vector3d c1  = D1  *P[ind[j+2]] ;
    Vector3d n = (b-a)^(c-a); n.normalize();
    
    if (!((a1-a)*n < 0 || (b1-b)*n<0 || (c1-c)*n<0)){ cull_success++; return false;}
    return true;
  }
#endif

  inline Voxel get_seed(int current_res) const {
    std::vector<Voxel> voxels;
    for (uint j = 0; j< ind.size(); j+=3){
      if(voxels.size()>0) break; 
      genFacetForwardNormal(current_res,0,1,j,std::back_inserter(voxels));
    }
    return *voxels.begin();
  }
  
  template<class OutputIterator>
  int get_seed(int current_res, OutputIterator oit) const {
    // get at least the old seed 
    *oit++=get_seed(current_res);
    // thats why counter already = 1 
    int count = 1;
    
    // get more voxels 
    std::vector<Voxel> voxels;
    CGAL::Random random(244);
    while(count<20){
      voxels.clear();
      int index = random.get_int(0,ind.size()/3)*3;
      genFacetForwardNormal(current_res,0,1,index,std::back_inserter(voxels));
      if(voxels.size()>0){
        count++;
        *oit++=*voxels.begin();
      }
    }
    return count;     
  }


  template<class OutputIterator>
  inline OutputIterator genPatch(int res, const AT_3 A, const AT_3 B, const wingedEdge &edge, OutputIterator oit) const{
    Point_3 a,b,c,d;  
    if(needPatch(A,B,edge,a,b,c,d)){
      oit = emitFacet(res, a,b,c, oit);
      oit = emitFacet(res, b,c,d, oit);
    }
    return oit;  
  }

  template<class OutputIterator>
  inline OutputIterator genFacet(int res, int i, int j, OutputIterator oit) const{
    Point_3 a,b,c;  
    if(needFacet(i,j,a,b,c)){
      oit = emitFacet(res, a,b,c, oit);
    }
    return oit;  
  }

  template<class OutputIterator>
  inline OutputIterator genFacetForwardNormal(int res, int i0, int i1, int j, OutputIterator oit) const{
    Point_3 a,b,c;  
    if(needFacetForwardNormal(i0,i1,j,a,b,c)){
      oit = emitFacet(res, a,b,c, oit);
    }
    return oit;  
  }
  
  template<class OutputIterator> 
  inline OutputIterator genNonManifoldEdge(int res, int i, int j, OutputIterator oit) const{

    //const Matrix4d& D0   = Track[i];
    //const Matrix4d& D1   = Track[i+1];
    //Vector3d a0  = D0   *P[non_manifold_edges[j].first ] ;
    //Vector3d b0  = D0   *P[non_manifold_edges[j].second] ;
    //Vector3d a1  = D1   *P[non_manifold_edges[j].first ] ;
    //Vector3d b1  = D1   *P[non_manifold_edges[j].second] ;
    AT_3 D0   = Track[i];
    AT_3 D1   = Track[i+1];
    Point_3 a0  = D0.transform(P[non_manifold_edges[j].first ]) ;
    Point_3 b0  = D0.transform(P[non_manifold_edges[j].second]) ;
    Point_3 a1  = D1.transform(P[non_manifold_edges[j].first ]) ;
    Point_3 b1  = D1.transform(P[non_manifold_edges[j].second]) ;

    // just draw all possible triangles .-) 
    oit = emitFacet(res, a0,a1,b0,oit);
    oit = emitFacet(res, a0,a1,b1,oit);
    oit = emitFacet(res, b0,b1,a0,oit);
    oit = emitFacet(res, b0,b1,a1,oit); 
    return oit;  
  }

};

} // namespace SV 

#endif // SV_VOXELIZATION_3_H 
