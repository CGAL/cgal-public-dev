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
//                 Andreas von Dziegielewski (anvdzi@gmail.com)
//
// ================================================================================

#ifndef SV_SWEEPT_VOLUME_WITH_VHULL_3_H
#define SV_SWEEPT_VOLUME_WITH_VHULL_3_H


#ifndef SV_HAS_INITIAL_POINTS
#define SV_HAS_INITIAL_POINTS 1
#endif

#ifndef SV_HAS_INTERSECT
#define SV_HAS_INTERSECT 1
#endif

#ifndef SV_USE_TWO_OFFSETS
#define SV_USE_TWO_OFFSETS 1
#endif 

#ifdef _OPENMP
#ifndef USE_OMP
#define USE_OMP 0
#endif
#endif 

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Real_timer.h>
#include <CGAL/IO/Scanner_OFF.h> 

#include <boost/optional.hpp>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

//#include <defines.h>
#include <stack>
#include <queue>
#include <algorithm>

//#include <timing.h>

//#include <SV/geometry.h>
#include <SV/Cloud_3.h>

#include <SV/Voxelization_3.h>



namespace SV{ 

// Swept_volume models interface of a implicit function
template < class Kernel_ = CGAL::Exact_predicates_inexact_constructions_kernel >
class Swept_volume_with_vhull_3 {
public:
  // TYPEDEFS
  typedef Kernel_ Kernel;
  typedef Swept_volume_with_vhull_3 <Kernel> Self;
   
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Sphere_3 Sphere_3;
  typedef typename Kernel::Aff_transformation_3 AT_3;

private:
  // TYPEDEFS
  typedef typename Voxelization_3::Voxel Voxel; 
  typedef typename Voxelization_3::VoxelHashSet VoxelHashSet; 
  typedef typename Voxelization_3::INT INT; 
  typedef SV::Cloud_3 <Kernel, INT> Cloud_3;
  typedef typename Cloud_3::Octree_3 Octree_3; 
  
 // typedef Kernel::Point_3 Point_3;
    
private:
  // MEMBER
  mutable std::map <std::vector <Point_3>, bool> is_accepted_facet;
  // std::vector <Matrix4 <double> > trajectory;
  Cloud_3 cloud;

public:
  Voxelization_3 vhull; 

public:
  mutable CGAL::Real_timer compute_facet_badness_timer;
  mutable CGAL::Real_timer intersect_timer;
  mutable CGAL::Real_timer operator_timer;
  mutable CGAL::Real_timer initialization_timer;
  mutable CGAL::Real_timer offset_timer; 
  mutable unsigned int compute_facet_badness_counter;

private:
  double m_delta,m_squared_delta,m_epsilon, m_squared_epsilon;
	AT_3 backtrafo;

public:
  void printHullTimings(){
    vhull.printTimings();
  }

public:
  // CONSTRUCTORS
  Swept_volume_with_vhull_3(
      const std::string& filename_track, const std::string& filename_off, 
      const int n, const int downstep = 0 , const int num_threads = 1) 
    :vhull(n), compute_facet_badness_counter(1){
    this->initialize(filename_track,filename_off, n, downstep, num_threads);
  }
  
  Swept_volume_with_vhull_3(): compute_facet_badness_counter(1){}
private:
//  Swept_volume_with_vhull_3(): compute_facet_badness_counter(1){} // not sure this is nice 
  //Swept_volume_with_vhull_3(const Self&) {}
  //Self operator=(Self&) {}

public:
// MAIN FUNCTION !!! 
  void initialize(const std::string& filename_track, const std::string& filename_off, 
      const int n, const int downstep, const int num_threads){
    std::cerr << "\n =================== \n " << n << std::endl;
    std::cerr << "initialize with resolution:   " << n << std::endl;

    m_delta    = (1.0/(2<<(n)));
    m_squared_delta   = (m_delta*m_delta);
    m_epsilon  = (m_delta/2.0f);
    m_squared_epsilon = (m_epsilon*m_epsilon);
    
    
    initialization_timer.start();

#if 0
    vhull = 
      Voxelization_3(
          filename_track.c_str(),filename_off.c_str(),
          0.01, downstep, num_threads);
#else
    vhull = 
      Voxelization_3(
          filename_track.c_str(),filename_off.c_str(),
          n, downstep, num_threads);
#endif
    backtrafo = vhull.back_transformation();
    Octree_3 volume = vhull.volume(); 
    // vhull.clear();
    
    std::cerr << "\n Memory before vhull clear " << SV::memory_usage()*100 << " %\n" ;
    // vhull = Voxelization_3(n);
    // std::cerr << "\n Memory after  vhull clear " << SV::memory_usage()*100 << " %\n" ;
    
    std::cerr << "Offset0 size: " <<  vhull.volume().size() << std::endl; 
    // ====== COMPUTE TWO OFFSETS FROM THIS OCTREE ============
    offset_timer.start();
    cloud = Cloud_3(volume);
    offset_timer.stop(); 
    

    std::cerr << "\n Memory after cloud " << SV::memory_usage()*100 << " %\n" ;
    
    std::cerr << "offset_time:   " << offset_timer.time() << std::endl ; 
    initialization_timer.stop();
  }
    

  ~Swept_volume_with_vhull_3() {
    std::cerr << std::endl;
    std::cerr << "initialization_timer.time()       : "
              << initialization_timer.time() << std::endl;
    std::cerr << "\t timer_compute_hull.time()          : "
              << vhull.timer_compute_hull.time() << std::endl;
    std::cerr << "\t timer_fill_hull.time()          : "
              << vhull.timer_fill_hull.time() << std::endl;
    std::cerr << "\t offset_timer.time()        : "
              << offset_timer.time() << std::endl;
    
    
    std::cerr << "compute_facet_badness_timer.time(): "
              << compute_facet_badness_timer.time() << " ... "
              << compute_facet_badness_timer.time() / compute_facet_badness_counter
              << " per call" << std::endl;
    std::cerr << "intersect_timer.time()            : "
              << intersect_timer.time() << std::endl;
    std::cerr << "operator_timer.time()             : "
              << operator_timer.time() << std::endl;
  }
    
   
public:

//#if SV_HAS_INTERSECT
  Point_3
  intersect(Point_3 source,Point_3 target) const
  {
    //intersect_timer.start();
    //std::cerr << " s " << source << std::endl;
    //std::cerr << " t " << target << std::endl;
    int source_is_contained = this->operator()(source);
    int target_is_contained = this->operator()(target);
      
    //std::cerr << source << " | " << target << "   " 
    //          << source_is_contained << " | " << target_is_contained << std::endl;
 
    if(source_is_contained == target_is_contained) {
      std::cerr << "ERROR: Segment does not intersect Swept Volume "
                << std::endl;
      exit(0);
    }
   
    //while (CGAL::squared_distance(source,target) >= m_squared_delta / 2000) {
    int count = 0;
    Point_3 m;
    while(true){
      m = CGAL::midpoint(source , target);
      //std::cerr << "m: " << m << " " << this->classify(m) << std::endl; 
      if((this->classify(m)==1 || this->classify(m)==2)) count++;
      if(count > 15) break; 
      // std::cerr << m << std::endl;
      CGAL_precondition(this->operator()(source) != this->operator()(target));
      if (this->operator()(m) == source_is_contained) {
        source = m;
      }
      else {
        target = m;
      }
    }
    
    // std::cerr << "intersect result: " << m << " " << this->classify(m) << std::endl; 
    verify_point(m);
    //intersect_timer.stop();
    return m;
  }  
//#endif

  template  <class OutputIterator> 
  void construct_initial_points(OutputIterator oit, int n) const {
    std::cerr <<" construct_initial_points begin " << std::endl; 
    CGAL_precondition(!cloud.offset0.empty());
    CGAL_precondition(!cloud.offset1.empty());
    CGAL_precondition(!cloud.offset2.empty());

    typename Octree_3::USet::iterator off0_it = cloud.offset0.set().begin();
    while(CGAL::cpp0x::get<3>(*off0_it) != cloud.offset0.resolution()-1){off0_it++;}
        
    Point_3 p( 
        (CGAL::cpp0x::get<0>(*off0_it)+0.5)/cloud.scale ,  
        (CGAL::cpp0x::get<1>(*off0_it)+0.5)/cloud.scale ,
        (CGAL::cpp0x::get<2>(*off0_it)+0.5)/cloud.scale );

    assert(cloud.classify(p) == 0); // point must be inside 
      

    typename Octree_3::USet::iterator off2_it = cloud.offset2.set().begin();

    int count = 0; 
    Point_3 last(0,0,0); 
    while(count < n ){
      while(CGAL::cpp0x::get<3>(*off2_it) != cloud.offset2.resolution()-1){off2_it++;}
      Point_3 q( 
          (CGAL::cpp0x::get<0>(*off2_it)+0.5)/cloud.scale ,  
          (CGAL::cpp0x::get<1>(*off2_it)+0.5)/cloud.scale ,
          (CGAL::cpp0x::get<2>(*off2_it)+0.5)/cloud.scale );
      //std::cerr << q << " " << cloud.classify(q) << std::endl; 
      assert(cloud.classify(q) == 2); // point must be outside 
      Point_3 r = intersect(p,q);
      // std::cerr << r << " " << cloud.classify(r) ;
      if(CGAL::squared_distance(last,r)>m_squared_delta){
        *oit++ = last = r ;
        count++;   
        //std::cerr << " accepted " << std::endl;
      }else{
        //std::cerr << " rejected " << std::endl;
      }
      off2_it++;
   }
    std::cerr <<" construct_initial_points end " << std::endl; 
  }
    
  typedef int return_type;
  return_type
  operator()(const Point_3& p, const bool = true) const{
    
    operator_timer.start();
    return_type r = (this->classify(p)<=1)?-1:1;
    operator_timer.stop();
    // std::cerr << " " << r << std::flush; 
    return r; 
  }

  // 0,1 is in, 2 or 3 is out 
  return_type
  classify(const Point_3& p) const
  {
#if 0
    return cloud.classify(p);   
#endif 
    
#if SV_USE_TWO_OFFSETS
       
    if(p.x()<0 || p.y()<0 || p.z()<0) return 3; 
    if(p.x()>1 || p.y()>1 || p.z()>1) return 3; 
    
    INT i,j,k; 
    cloud.compute_key(p,i,j,k);
    int c = cloud.classify(i,j,k);
    if (c == 0) return c;
    if (c == 3) return c;    

    double dx = p.x()*cloud.scale()-i-0.5; 
    double dy = p.y()*cloud.scale()-j-0.5; 
    double dz = p.z()*cloud.scale()-k-0.5; 

    int sx = (dx < 0)?-1:1; 
    int sy = (dy < 0)?-1:1; 
    int sz = (dz < 0)?-1:1; 
        
    double ax = CGAL::abs(dx); 
    double ay = CGAL::abs(dy); 
    double az = CGAL::abs(dz); 


    int ncx   = cloud.classify(i+sx,j,k);
    int ncy   = cloud.classify(i,j+sy,k);    
    int ncz   = cloud.classify(i,j,k+sz);
    int ncxy  = cloud.classify(i+sx,j+sy,k);
    int ncxz  = cloud.classify(i+sx,j,k+sz);    
    int ncyz  = cloud.classify(i,j+sy,k+sz);  
    int ncxyz = cloud.classify(i+sx,j+sy,k+sz);

//#define SVS(X,Y,Z) (1-(std::max)((std::max)(X,Y),Z))
//#define SVS(W,X,Y,Z) W = std::max(0.0,1-CGAL::sqrt(X*X+Y*Y+Z*Z))

#define SVS(W,X,Y,Z){                                      \
        int k = 10;                                         \
        double t(CGAL::kth_root(k,                         \
                CGAL::ipower(X,k)+                         \
                CGAL::ipower(Y,k)+                         \
                CGAL::ipower(Z,k)));                       \
            t *= 0.5;                                      \
            t += 0.5;                                      \
            t = std::min(1.0,t);                           \
            W= t*t*(1-t)*(1-t);                            \
            }                                              \
      
    double wo,wx,wy,wz,wxy,wxz,wyz,wxyz;
    SVS(wo,ax,ay,az);
    SVS(wx,(1-ax),ay,az);
    SVS(wy,ax,(1-ay),az);
    SVS(wz,ax,ay,(1-az));
    SVS(wxy,(1-ax),(1-ay),az);
    SVS(wxz,(1-ax),ay,(1-az));
    SVS(wyz,ax,(1-ay),(1-az));
    SVS(wxyz,(1-ax),(1-ay),(1-az));             
#undef SVS

    double C = c*wo
      +ncx*wx
      +ncy*wy
      +ncz*wz
      +ncxy*wxy
      +ncxz*wxz
      +ncyz*wyz
      +ncxyz*wxyz;
    C /= (wo+wx+wy+wz+wxy+wxz+wyz+wxyz);
    
    // std::cerr << C << " " << std::endl; 

    return (C <= 1.2)?1:2; 
    

    

//     if( ax+ay+az > 1){
//       if (c == 1){
//         if (ncx == 2 && ncy == 2 && ncz == 2){
//           return 2 ;
//         }
//       }else{
//         if (ncx == 1 && ncy == 1 && ncz == 1){
//           return 1 ;
//         }
//       }
//     }
//     return c;
#else
    int c = cloud.classify(p);
    if (c == 0) return c;
    if (c == 3) return c;    
    
    CGAL_precondition(c == 1);
    
    INT i,j,k; 
    cloud.compute_key(p,i,j,k);

    double dx = p.x()*cloud.scale()-i-0.5; 
    double dy = p.y()*cloud.scale()-j-0.5; 
    double dz = p.z()*cloud.scale()-k-0.5; 

    int sx = (dx < 0)?-1:1; 
    int sy = (dy < 0)?-1:1; 
    int sz = (dz < 0)?-1:1; 
    
    double ax = CGAL::abs(dx); 
    double ay = CGAL::abs(dy); 
    double az = CGAL::abs(dz); 
    // std::cerr << sx << " " << sy << " " << sz << std::endl ;
   
   
    // consider all neighbours of point 
    c = cloud.classify(i+sx,j,k);
    if (c == 0) return c;
    c = cloud.classify(i,j+sy,k);
    if (c == 0) return c;
    c = cloud.classify(i,j,k+sz);
    if (c == 0) return c;
    
    if (ay+az > 0.5){
      c = cloud.classify(i,j+sy,k+sz);
      if (c == 0) return c;
    }else{
      c = cloud.classify(i,j-sy,k)
        + cloud.classify(i,j,k-sz);
      if (c == 0) return c;
    }

    if (ax+az > 0.5){
      c = cloud.classify(i+sx,j,k+sz);
      if (c == 0) return c;
    }else{
      c = cloud.classify(i-sx,j,k)
        + cloud.classify(i,j,k-sz);
      if (c == 0) return c;
    }

    if (ax+ay > 0.5){
      c = cloud.classify(i+sx,j+sy,k);
      if (c == 0) return c;
    }else{
      c = cloud.classify(i-sx,j,k)
        + cloud.classify(i,j-sy,k);
      if (c == 0) return c;
    }

    if (ax+ay+az > 1){
      c = cloud.classify(i+sx,j+sy,k+sz);
      if (c == 0) return c;
    }else{
      c = cloud.classify(i-sx,j,k)
        + cloud.classify(i,j-sy,k)
        + cloud.classify(i,j,k-sz);
      if (c == 0) return c;
    }


    c = cloud.classify(i+sx,j+sy,k)
      + cloud.classify(i,j+sy,k+sz)
      + cloud.classify(i+sx,j,k+sz);
    if (c == 0) return c;
   
    if (ax+ay-az < 0){
      c = cloud.classify(i-sx,j-sy,k)
        + cloud.classify(i,j-sy,k+sz)
        + cloud.classify(i-sx,j,k+sz);
      if (c == 0) return c;
    }

    if (ax-ay+az < 0){
      c = cloud.classify(i-sx,j+sy,k)
        + cloud.classify(i,j+sy,k-sz)
        + cloud.classify(i-sx,j,k-sz);
      if (c == 0) return c;
    } 
    
    if (-ax+ay+az < 0){
      c = cloud.classify(i+sx,j-sy,k)
        + cloud.classify(i,j-sy,k-sz)
        + cloud.classify(i+sx,j,k-sz);
      if (c == 0) return c;
    }

    if (ax-ay-az < 0){
      c = cloud.classify(i-sx,j+sy,k)
        + cloud.classify(i,j+sy,k+sz)
        + cloud.classify(i-sx,j,k+sz);
      if (c == 0) return c;
    }
        
    if (-ax+ay-az < 0){
      c = cloud.classify(i+sx,j-sy,k)
        + cloud.classify(i,j-sy,k+sz)
        + cloud.classify(i+sx,j,k+sz);
      if (c == 0) return c;
    }
    if (-ax-ay+az < 0){
      c = cloud.classify(i+sx,j+sy,k)
        + cloud.classify(i,j+sy,k-sz)
        + cloud.classify(i+sx,j,k-sz);
      if (c == 0) return c;
    }
    // with non of them gives a decision we count it as out side 
    return 3; 
#endif
  }

  Sphere_3
  bounding_sphere() const
  {
    return Sphere_3(Point_3(1.0, 1.0, 1.0), 4.0);
  }

  void
  verify_point(const Point_3& p) const
  {
    int classification = cloud.classify(p);
    if (classification == 0)
      std::cerr << "point inside  " << std::endl; 
    //if (classification == 1)
    //  std::cerr << "point in off 1 " << std::endl; 
    //if (classification == 2)
    //  std::cerr << "point in off 2 " << std::endl; 
    if (classification == 3)
      std::cerr << "point too far  " << std::endl; 
  }

  template < class OutputIterator >
  OutputIterator
  refine_triangle_sampling(int n, const Point_3& p1, const Point_3& p2,  const Point_3& p3, OutputIterator oit) const
  {
    if (n <= 0) {
      // std::cerr <<"  reached recursion depth !! " << std::endl;
      return oit;
    }

    // std::cerr << "refine_triangle_sampling depth: "<< n << std::endl;
    double d12 = CGAL::squared_distance(p1,p2);
    double d13 = CGAL::squared_distance(p1,p3);
    double d23 = CGAL::squared_distance(p2,p3);
    double d = std::max(std::max(d12, d13), d23);
    if (d <= m_squared_epsilon) {
      // std::cerr << "edge length2: " << d << " is  <" << m_squared_epsilon << std::endl;
      return oit;
    }
    //std::cerr << "edge length2: " << d << " is! <" << m_squared_epsilon << std::endl;
      
    if (d12 >= d13 && d12 >= d23) {
      //std:: cout << "D12" << std::endl;
      Point_3 c = CGAL::midpoint(p1 , p2);
      *oit++ = c;
      oit = refine_triangle_sampling(n - 1, p1, c, p3, oit);
      oit = refine_triangle_sampling(n - 1, p2, c, p3, oit);
      return oit;
    }
    if (d13 >= d12 && d13 >= d23) {
      //std:: cout << "D13" << std::endl;
      Point_3 c = CGAL::midpoint(p1 , p3);
      *oit++ = c;
      oit = refine_triangle_sampling(n - 1, p1, c, p2, oit);
      oit = refine_triangle_sampling(n - 1, p3, c, p2, oit);
      return oit;
    }
    if (d23 >= d12 && d23 >= d13) {
      //std:: cout << "D14" << std::endl;
      Point_3 c = CGAL::midpoint(p2 , p3);
      *oit++ = c;
      oit = refine_triangle_sampling(n - 1, p2, c, p1, oit);
      oit = refine_triangle_sampling(n - 1, p3, c, p1, oit);
      return oit;
    }
    exit(0);
  }
    
  typedef boost::optional <double> Facet_badness;
  Facet_badness
  compute_facet_badness(const Point_3& p1, const Point_3& p2, const Point_3& p3) const
  {
    compute_facet_badness_timer.start();
    Facet_badness result = compute_facet_badness_(p1, p2, p3);
    compute_facet_badness_timer.stop();
    return result;
  }

  

  Facet_badness   compute_facet_badness_(const Point_3& a, const Point_3& b,   const Point_3& c) const
  {
    ++compute_facet_badness_counter;
    // std::cerr << is_accepted_facet.size() << std::endl;

    // returning the default constructed boost optional means that the facet is ok.
    // return a positive value if the facet is not ok. e.g. Facet_badness(1.0)
    // where 0 is the worst quality.

    std::vector <Point_3> key;
    key.push_back(a);
    key.push_back(b);
    key.push_back(c);
    std::sort(key.begin(), key.end());

    if (is_accepted_facet[key]) {
      //std::cerr << " old facet that was already accepted " << std::endl;
      return Facet_badness();
    }

    //verify_point(p1);
    //verify_point(p2);
    //verify_point(p3);

    double maxedge = CGAL::squared_distance(a,b);
    maxedge = std::max( maxedge, CGAL::squared_distance(a,c));
    maxedge = std::max( maxedge, CGAL::squared_distance(b,c));
    
  

    // I keep this part for fast discartds
#if 1
    {
      std::vector <Point_3> samples;
      refine_triangle_sampling(4, a, b, c, std::back_inserter(samples));
      
      for (uint i = 0; i < samples.size(); i++) {
        int classification = cloud.classify(samples[i]);
        if(classification == 0 || classification == 3){
          return Facet_badness(1.0f / maxedge);
        }
      }
      samples.clear();
    }
#endif 
    std::vector<Voxel> voxelization;
    vhull.get_voxelization( 
        cloud.resolution(),
        float(a.x()), float(a.y()), float(a.z()), 
        float(b.x()), float(b.y()), float(b.z()), 
        float(c.x()), float(c.y()), float(c.z()), 
        std::back_inserter(voxelization));
    
    std::random_shuffle(voxelization.begin(), voxelization.end());
    
//    {   // this is all debug stuff 
//      std::cerr << "voxelize:\n"  
//                << "a: "<< a  << std::endl
//                << "b: "<< b  << std::endl
//                << "c: "<< c  << std::endl
//                << std::endl;
//
//      std::cerr << "resolution: "<< cloud.resolution() << std::endl;
//      std::cerr << "scale:      "<< cloud.scale()      << std::endl;
//      std::cerr << "voxel values should be between:"   << std::endl;
//      INT x,y,z; 
//      cloud.compute_key(a,x,y,z);
//      std::cerr << "key for a " << x << " " << y << " " << z << std::endl;  
//      cloud.compute_key(b,x,y,z);
//      std::cerr << "key for b " << x << " " << y << " " << z << std::endl;  
//      cloud.compute_key(c,x,y,z);
//      std::cerr << "key for c " << x << " " << y << " " << z << std::endl;  
//      
//      std::cerr << "this is what vhull produces:" << std::endl;
//    
//      BOOST_FOREACH(Voxel v, voxelization){
//        std::cerr <<  CGAL::cpp0x::get<0>(v)<<" "  << CGAL::cpp0x::get<1>(v) <<" " << CGAL::cpp0x::get<2>(v) << " is " 
//                  <<  cloud.classify(CGAL::cpp0x::get<0>(v),CGAL::cpp0x::get<1>(v),CGAL::cpp0x::get<2>(v)) <<  std::endl;
//      }
//      std::cerr <<std::endl; 
//    }
   
    int classification;
    BOOST_FOREACH(Voxel v, voxelization){     
      classification = cloud.classify(
          CGAL::cpp0x::get<0>(v),CGAL::cpp0x::get<1>(v),CGAL::cpp0x::get<2>(v));      
      if(classification == 0 || classification == 3){
        return Facet_badness(1.0f / maxedge);
      }
    }
        
    is_accepted_facet[key] = true;
    return Facet_badness();
  }


AT_3& getBackTrafo(){
			return backtrafo;
}

};
  

  template < class Kernel >
  std::ostream&
  operator<<(std::ostream &os, Swept_volume_with_vhull_3 <Kernel> const &n)
  {
    // TODO
    return os << "write Swept_volume_with_vhull_3 into a .SV3 file here"
              << std::endl;
  }

  template < class Kernel >
  std::istream&
  operator>>(std::istream &is, Swept_volume_with_vhull_3 <Kernel> &n)
  {
    std::cerr << "In instream operator  ... " << std::endl;
    n.initialize("current.tarr","current.off",10);
    std::cerr << "...done." << std::endl;
    return is;
  }

} // namespace SV 

#endif // SV_SWEEPT_VOLUME_WITH_VHULL_3_H
