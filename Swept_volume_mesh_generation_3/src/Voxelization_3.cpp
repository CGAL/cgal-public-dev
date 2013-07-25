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


#include <SV/Voxelization_3.h>

#if USE_OMP
#define USE_NO_PARALLEL 0
#else
#define USE_NO_PARALLEL 1
#endif

#include <vector>
#include <set>
#include <stack>
#include <queue>
#include <CGAL/tuple.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
//#include </usr/include/GL/glut.h>
//#include "vecmath.h"
//#include "OffReader.h"
#include <stdlib.h>
#include <algorithm>

#include <SV/compute_hull.h>
#include <SV/compute_scout.h>
#include <SV/fill_hull.h>
#include <SV/verify_hull.h>
#include <boost/foreach.hpp>
#include <CGAL/iterator.h>


//#define NAME_TARR_FILE "../../data/Swept_volume/track.tarr"
//#define NAME_OFF_FILE "../../data/Swept_volume/newDragon.off"
#define NAME_TARR_FILE "../../data/Swept_volume/Ausbau_M274_Langner.tarr"
#define NAME_OFF_FILE "../../data/Swept_volume/M274.off"
//#define NAME_TARR_FILE "../../data/Swept_volume/motor_ausbau.tarr"
//#define NAME_OFF_FILE "../../../../andreasd/OffFiles/Daimler/M272_NAG2_Silhouette_1mm.off"

/* 
 * Implementation 
 * 
 */
namespace SV{  



  Voxelization_3::Voxelization_3(
      const char* filename_track, const char* filename_off, 
      int n, int downstep, int num_threads) 
    :
      m_resolution(n),
      m_downstep(downstep), 
      m_num_threads(num_threads),
      cull_success(0),
      cull_candidate(0)
  {
    // timer_total.start();

    char* name_off = (char*) filename_off;	
    char* name_track = (char*) filename_track;

    std::cout <<"Loading off-File: " << filename_off << std::endl;
    //Off-File laden:
    LoadOffFile(filename_off,P,ind);
    std::cout <<"Size of mesh: " << ind.size()/3 << std::endl;


    loadFile(filename_track, Track);     // use tarr file
#if 0
    loadPathFile(filename_track, Track); // use teraport path file
#endif

    normalize(P, Track);

    //Koordinaten normalisieren
    std::cerr << "Anzahl Dreiecke: " << ind.size() << std::endl;
    std::cerr << "Anzahl Vertices: " << P.size() << std::endl;

    initialize_winged_edge();
    generate_voxelization();

    // timer_total.stop();
  }  

  void Voxelization_3::initialize_winged_edge() {
    // Generate Winged edge Structure
    for (uint i = 0; i< ind.size(); i+=3){
      for (uint j = 0; j<3; j++){
        edge ea;
        uint e1    = ind[i+((j+0)%3)];
        uint e2    = ind[i+((j+1)%3)];
        ea.first   = std::min(e1,e2); 
        ea.second  = std::max(e1,e2); 
        wing wa;
        wa.first  = ea;
        wa.second = ind[i+((j+2)%3)]; 
        Wings.push_back(wa);
      }
    }

    std::cerr << " wings " << Wings.size() << std::endl; 
    std::sort(Wings.begin(),Wings.end());
    for (uint i = 0; i< Wings.size(); i+=2){

      // manifold edge check
      if (Wings[i  ].first.first  != Wings[i+1].first.first  ||
          Wings[i  ].first.second != Wings[i+1].first.second){ 
        non_manifold_edges.push_back(Wings[i].first);
        i-=1;
        continue;
      }
      wingedEdge we;
      we.push_back(Wings[i  ].first.first); 
      we.push_back(Wings[i  ].first.second);  
      we.push_back(Wings[i  ].second);  
      we.push_back(Wings[i+1].second);  
      wingedEdges.push_back(we);
    }
    std::cerr  << " # expected edges----: " << Wings.size()/2 << std::endl; 
    std::cerr  << " # non-manifold edges: " << non_manifold_edges.size()/2 << std::endl; 
    std::cerr  << " # winged edges------: " << wingedEdges.size() << std::endl;
  }

  void Voxelization_3::generate_voxelization() const {
    if(opt_volume) return; 

    timer_total.start();
    opt_volume = Octree_3(m_resolution);

    using CGAL::cpp0x::get;  

    for(int current_res = resolution() - m_downstep; current_res < resolution(); current_res++){
      opt_volume->set_resolution(current_res);
      std::cerr << " // ------------------------------ // " << std::endl;
      std::cerr << " opt_volume->resolution() " <<  opt_volume->resolution() << std::endl; 

      timer_emit_boundary.start();
      emit_boundary(current_res,*opt_volume);  
      timer_emit_boundary.stop();
      compress(current_res,true); 
      std::cerr << " reduced volume size (octree): " << opt_volume->size() << std::endl; 
    }


    int current_res = resolution();
    opt_volume->set_resolution(current_res);
    Octree_3 old_octree = *opt_volume; 

    while(true){  
      std::cerr << " // ------------------------------ // " << std::endl;
      std::cerr << " opt_volume->resolution() " <<  opt_volume->resolution() << std::endl; 

      timer_emit_boundary.start(); 
      emit_boundary(current_res,*opt_volume);  
      timer_emit_boundary.stop(); 

#if 1 
      compress(current_res,false);
      break; // ignore conflicts 
#else 
      VoxelHashSet hull_voxels;  
      std::vector<Voxel> seeds; 
      get_seed(current_res,std:back_inserter(seeds));
      Voxel scout = SV::compute_scout(*opt_volume,seeds[0]); 
      timer_compute_hull.start(); 
      SV::compute_hull(*opt_volume, scout, hull_voxels);
      timer_compute_hull.stop();
      std::cerr << "Hull size: " << hull_voxels.size() << std::endl; 

      // fill hull, use hierrachcal fill in the first run 
      timer_fill_hull.start(); 
      if(old_octree.size() == 0){
        *opt_volume = SV::fill_hull(current_res, hull_voxels, seeds.begin(),seeds.end());
      }else{
        SV::fill_hull(current_res, hull_voxels, *opt_volume);
      }
      std::cerr << "Volume size (octree): " << opt_volume->size() << std::endl; 
      timer_fill_hull.stop(); 

      if(old_octree.size() == 0) break; // no conflicts possible 


      // check if the old octree did not occupy too much 
      // that is, if the hull can see a voxel of the octree of the last run
      // if so, we remove the relevant voxels from the old octree 
      // and restart the last run, that is, current_res--, opt_volume=old_octree-conflicts  
      HashSet_Voxel<INT> conflicts;
      verify_hull(hull_voxels.begin(), hull_voxels.end(),old_octree,CGAL::inserter(conflicts));
      if(conflicts.size() == 0) break; 

      std::cerr << " !!!!!! THERE WHERE " << conflicts.size() << " CONFLICTS " << std::flush;

      BOOST_FOREACH(const Voxel& v, conflicts){
        // std::cerr << get<0>(v)<<" , "<<get<1>(v)<<" , "<<get<2>(v) << " " 
        // 		<< opt_volume->contains(v) << " "  
        // 		<< old_octree.contains(v) << " "
        //	<< std::flush;
        for(int i = -1; i <= 1; i++){
          for(int j = -1; j <= 1; j++){
            for(int k = -1; k <= 1; k++){
              opt_volume->erase((get<0>(v))+i,(get<1>(v))+j,(get<2>(v))+k);
              old_octree.erase((get<0>(v))+i,(get<1>(v))+j,(get<2>(v))+k);
            }
          }
        }
        //opt_volume->erase(v);
        //old_octree.erase(v);
        // std::cerr << opt_volume->contains(v) << " "  
        // 		<< old_octree.contains(v) << " " 
        // 		<< std::endl;
      }


      // if(conflicts.size() < 30){
      //   BOOST_FOREACH(const Voxel& v, conflicts){
      // 	std::cerr << get<0>(v)<<" , "<<get<1>(v)<<" , "<<get<2>(v)<< std::endl;
      //   }
      // }

      hull_voxels.clear();
#endif
      std::cerr << " Resolution level: " << current_res << " opt_volume->size() " << opt_volume->size() << std::endl;  
      std::cerr << " Memory usage    : " << SV::memory_usage()*100 << " %\n" ; 
      std::cerr << std::endl;
    }


    std::cerr  << "\n  # culling ratio -----: " 
      << int(cull_success*100.0/cull_candidate)<<" %" << std::endl ; 

    std::cerr << " Memory usage " << SV::memory_usage()*100 << " %\n" ;

    // CLEAN UP
    timer_total.stop();
    this->printTimings();  
  }



  void Voxelization_3::compress(int current_res, bool compute_inverse_offset) const {

    std::cerr << "\n Octree size before compression:      " << opt_volume->size() << "\n" ;
    using CGAL::cpp0x::get;

    VoxelHashSet hull_voxels;  
    std::vector<Voxel> seeds;
    get_seed(current_res,std::back_inserter(seeds));
    Voxel scout = SV::compute_scout(*opt_volume,seeds[0]); 

    timer_compute_hull.start(); 
    SV::compute_hull(*opt_volume, scout, hull_voxels);
    timer_compute_hull.stop();
    std::cerr << "Hull size: " << hull_voxels.size() << std::endl; 


    // use hierrachcal fill in the first run 
    timer_fill_hull.start(); 
    if(current_res == resolution() - m_downstep){
      *opt_volume = SV::fill_hull(current_res, hull_voxels, seeds.begin(),seeds.end());
    }else{
      SV::fill_hull(current_res, hull_voxels, *opt_volume);
    }        
    timer_fill_hull.stop(); 


    if(compute_inverse_offset){
      std::cerr << "compute inverse offset" << std::endl; 
      BOOST_FOREACH(Voxel v, hull_voxels){
        for(int i = -1; i <= 1; i++){
          for(int j = -1; j <= 1; j++){
            for(int k = -1; k <= 1; k++){
              opt_volume->erase(get<0>(v)+i,get<1>(v)+j,get<2>(v)+k);
            }
          }
        }
      }
    }
    hull_voxels.clear();
    std::cerr << "Volume size (octree) after compression: " << opt_volume->size() << std::endl; 
    // std::cerr << opt_volume->size() << std::endl;
  }

  // returns (p-a).normalize.dot(n) * (p-b).normalize.dot(n)
  // that is, return value is positive, iff a and b are on 
  // the same side of the plane defined by p and n 
  // dziegiel - 17.6.13
  double Voxelization_3::pointsOnSameSide(Point_3 a, Point_3 b, Point_3 n, Point_3 p) const {
#if 0
// #if TERAPORT_BUILT
    Vector3d pa = a-p;
    pa.normalize();
    double aSide = pa*n;
    Vector3d pb = b-p;
    pb.normalize();
    double bSide = pb*n;
    return aSide*bSide;  
//   if (aSide*bSide > 0.001) return true;
 //   return false;
//#else 
//  // with no eps
//    double aSide = (a-p)*n;
//    double bSide = (b-p)*n;
//    if (aSide*bSide > 0) return true;
//    return false;
//#endif

#endif
	return 0;  // TODO _ANDI: is this needed for EXACT?
  }
  /*bool Voxelization_3::pointsOnSameSide(Vector3d a, Vector3d b, Vector3d c, Vector3d n, Vector3d p) const {
    double aSide = (a-p)*n;
    double bSide = (b-p)*n;
    if (aSide*bSide < 0) return false;
    double cSide = (c-p)*n;
    if (aSide*cSide < 0) return false;
    return true;
  }*/

  void Voxelization_3::loadPathFile(const char* filename, std::vector< AT_3 >& T) 
  {

#if 0
    std::cout << " Loading path-File: "<<filename  << std::endl;
    FILE* in;
    double v;
    uint k,zeilen=0;
    in = fopen(filename, "r");
    while( fscanf( in, "%lf", &v ) != EOF )
      ++zeilen;
    uint steps = zeilen/16;
    fseek(in, 0L, SEEK_SET);

    //printf("%d", steps);
    for(k=0;k<steps;k++)
    { 
      double  m00=0,m01=0,m02=0,m03=0,
              m10=0,m11=0,m12=0,m13=0,
              m20=0,m21=0,m22=0,m23=0,
              m30=0,m31=0,m32=0,m33=0;

      int dummy;
      dummy=  fscanf(in, "%lf", &m00);
      dummy=  fscanf(in, "%lf", &m01);
      dummy=  fscanf(in, "%lf", &m02);
      dummy=  fscanf(in, "%lf", &m03);
      dummy=  fscanf(in, "%lf", &m10);
      dummy=  fscanf(in, "%lf", &m11);
      dummy=  fscanf(in, "%lf", &m12);
      dummy=  fscanf(in, "%lf", &m13);
      dummy=  fscanf(in, "%lf", &m20);
      dummy=  fscanf(in, "%lf", &m21);
      dummy=  fscanf(in, "%lf", &m22);
      dummy=  fscanf(in, "%lf", &m23);
      dummy=  fscanf(in, "%lf", &m30);
      dummy=  fscanf(in, "%lf", &m31);
      dummy=  fscanf(in, "%lf", &m32);
      dummy=  fscanf(in, "%lf", &m33);

      dummy+=2.0; //prevent warning
      T.push_back(Matrix4d( 
            m00, m01, m02, m03,
            m10, m11, m12, m13,
            m20, m21, m22, m23,
            m30, m31, m32, m33 ))  ;
    }	
    fclose(in);
    std::cout << "Size of track: " << steps << std::endl;
  
	#endif // TODO _ANDI this is not needed
}




void Voxelization_3::LoadOffFile(
	const char* as_FileName, 
	std::vector<Point_3>& ak_Vertices, 
	std::vector<int>& ak_Indices) 
{
  std::cerr << "LoadOffFile(\"" << as_FileName << "\");" << std::endl;
  std::vector<float> lk_Faces;
  int li_TriangleCounter = 0;
  int li_VerticeLength = 0;
  int li_FaceCounter = 0;
	
  std::ifstream lk_InStream(as_FileName);
  if(!lk_InStream)
    {
      std::cerr << "Off-Datei nicht gefunden!" << std::endl;
      return;
    }
  std::string ls_FileHeader;
  std::getline(lk_InStream, ls_FileHeader);
  if (ls_FileHeader.length() > 0 && ls_FileHeader[ls_FileHeader.length()-1] == 0xd) 
    ls_FileHeader.erase (ls_FileHeader.length()-1,1);

  if(ls_FileHeader.compare ("OFF") == 0)
    {
      std::cerr << "Lese OFF-Datei..." << std::endl;
    }
  else 
    {
      std::cerr << "Keine gueltige OFF-Datei!" << std::endl;
      return;
    }

  lk_InStream >> li_VerticeLength;
  std::cerr << "Knoten: " << li_VerticeLength << std::endl;
  lk_InStream >> li_FaceCounter;
  std::cerr << "Flaechen: " << li_FaceCounter << std::endl;
  int li_Edges;
  lk_InStream >> li_Edges;
  std::cerr << "Kanten: " << li_Edges << std::endl;

  for(int i = 0; i < li_VerticeLength; i++) 
    {
      float x, y, z;
      lk_InStream >> x;
      lk_InStream >> y;
      lk_InStream >> z;
      ak_Vertices.push_back(Point_3(x, y, z));
    }

  int li_Temp1, li_Temp2, li_Temp3;

  for (int i = 0; i < li_FaceCounter; i++) 
    {
      int li_Kanten;
      lk_InStream >> li_Kanten;
      for(int j = 0; j < li_Kanten - 2; j++)
        {
          if (j == 0)
            {
              lk_InStream >> li_Temp1;
              lk_InStream >> li_Temp2;
              lk_InStream >> li_Temp3;
            }
          if (j > 0)
            {
              li_Temp2 = li_Temp3;
              lk_InStream >> li_Temp3;
            }
          li_TriangleCounter++;
          lk_Faces.push_back(float(li_Temp1));
          lk_Faces.push_back(float(li_Temp2));
          lk_Faces.push_back(float(li_Temp3));
        }
    }

  for (int i = 0; i < li_TriangleCounter * 3; i++)
    {
      ak_Indices.push_back(int(lk_Faces[i]));
    }

  //OUT: lk_Vertices, lk_Indices
}








  void Voxelization_3::loadFile(const char* filename, std::vector< AT_3 >& T) 
  {
    FILE* in;
    double v;
    uint k,zeilen=0;
    in = fopen(filename, "r");
    while( fscanf( in, "%lf", &v ) != EOF )
      ++zeilen;
    uint steps = zeilen / 12;
    fseek(in, 0L, SEEK_SET);
    for(k=0;k<steps;k++) 
    {
      double t0=0,t1=0,t2=0;
      double m00=0,m01=0,m02=0,m10=0,m11=0,m12=0,m20=0,m21=0,m22=0;

      int dummy;
      dummy= 	fscanf(in, "%lf", &m00);
      dummy= 	fscanf(in, "%lf", &m10);
      dummy= 	fscanf(in, "%lf", &m20);
      dummy= 	fscanf(in, "%lf", &m01);
      dummy= 	fscanf(in, "%lf", &m11);
      dummy= 	fscanf(in, "%lf", &m21);
      dummy= 	fscanf(in, "%lf", &m02);
      dummy= 	fscanf(in, "%lf", &m12);
      dummy= 	fscanf(in, "%lf", &m22);
      dummy= 	fscanf(in, "%lf", &t0);	
      dummy= 	fscanf(in, "%lf", &t1);	
      dummy= 	fscanf(in, "%lf", &t2);	

      dummy+=2.0; //prevent warning
      T.push_back(AT_3( 
            m00, m01, m02, t0,
            m10, m11, m12, t1,
            m20, m21, m22, t2,
                          1.0 ))  ;
    }	
    fclose(in);
    std::cerr << "Tarr mit " << steps << " Punkten geladen ..." << std::endl;
  }


  //Koordinaten auf [-1,1]^3 normalisieren:
  void Voxelization_3::normalize(std::vector<Point_3>& P, std::vector< AT_3 >& Track ) {

    double x = P[0].x(); double X = P[0].x();
    double y = P[0].y(); double Y = P[0].y();
    double z = P[0].z(); double Z = P[0].z();

    BOOST_FOREACH(Point_3 v, P){
      BOOST_FOREACH(AT_3 M, Track){
        Point_3 Mv = M.transform(v);
        x= (std::min)(x,Mv.x());  
        y= (std::min)(y,Mv.y());  
        z= (std::min)(z,Mv.z()); 
        X= (std::max)(X,Mv.x());  
        Y= (std::max)(Y,Mv.y());  
        Z= (std::max)(Z,Mv.z());                              
      }
    }


		// TODO _ANDI 0.7 seems a bit overkill to scale!
    double s =(0.7)/(std::max)((std::max)(X-x,Y-y),Z-z);
    std::cerr << "s: " << s << std::endl;
    // move lower left corner to origin 
    // invert scale to get (0,0,0)-(1,1,1)cube 
    // In order to write the original back, take this matrix 
    // and multiply each point with its inverse 
    //Matrix4d T1,T2,T3;
    //T3.makeTranslate(0.125,0.125,0.125);
    //T2.makeScale(s,s,s);
    //T1.makeTranslate(-x,-y,-z);
    //Matrix4d T = T3*T2*T1;
		AT_3 T1(1,0,0,0.125,
						0,1,0,0.125,
						0,0,1,0.125, 1);
		AT_3 T2(s,0,0,0,
						0,s,0,0,
						0,0,s,0, 1);
		AT_3 T3(1,0,0,-x,
						0,1,0,-y,
						0,0,1,-z, 1);
		AT_3 T = T3*T2*T1;

    BOOST_FOREACH(AT_3& M, Track){
      M=T*M;
    }

    // global MT for backtrafo
    //Matrix4d TT1,TT2,TT3;
    //TT3.makeTranslate(-0.125,-0.125,-0.125);
    //TT2.makeScale(1/s,1/s,1/s);
    //TT1.makeTranslate(x,y,z);
    //m_back_transformation = TT1*TT2*TT3;
		AT_3 TT1(1,0,0,-0.125,
					   0,1,0,-0.125,
						 0,0,1,-0.125, 1);
		AT_3 TT2(1/s,0,0,0,
					 	 0,1/s,0,0,
						 0,0,1/s,0, 1);
		AT_3 TT3(1,0,0,x,
						 0,1,0,y,
						 0,0,1,z, 1);
    m_back_transformation = TT1*TT2*TT3;
  }


#if USE_NO_PARALLEL || USE_OMP 
  void Voxelization_3::emit_boundary(int current_res, Octree_3& octree) const 
  {
#if USE_OMP
    if(m_num_threads == 1)
      emit_boundary_serial(current_res, octree);
    else
      emit_boundary_omp(current_res, octree);
#else
    emit_boundary_serial(current_res, octree);
#endif
  }
#endif

  void Voxelization_3::emit_boundary_serial(int current_res, Octree_3& octree) const {
    CGAL::Real_timer timer; 
    timer.start();

    std::cerr << "Start emit boundary: no multi threading" << std::endl;
    int i = 0; 

    for (uint j = 0; j< ind.size(); j+=3)
      genFacetForwardNormal(current_res,0,1,j,octree.inserter());

    while(i < Track.size()){
      // Edge Edge Patches
      // puts voxelization of edge-edge-patches-trianlges in voxSet
      if(0 <= i && i < Track.size()-1){
        for (uint j = 0; j< non_manifold_edges.size();j++){
          genNonManifoldEdge(current_res,i, j, octree.inserter());
        }
        for (uint j = 0; j< wingedEdges.size(); j++){
          genPatch(current_res,Track.at(i),Track.at(i+1),wingedEdges.at(j),octree.inserter());         
        }
      }       
      // Face-Vertex-Facets
      // puts voxelization of Facet-triangles in voxSet    
      if(0 < i && i < Track.size()-1){
        for (uint j = 0; j< ind.size(); j+=3){
          genFacet(current_res,i,j,octree.inserter());		
        }
      }
      i++;
      std::cerr << " TRACK " << i << " Memory usage " << SV::memory_usage()*100 << " %\r" << std::flush;
    }
    // put the full solid into octhull to close swept volume hull for compression 
    for (uint j = 0; j< ind.size(); j+=3)
      genFacetForwardNormal(current_res,i-1,i-2,j,octree.inserter());

    //m_seed = get_seed(current_res);

    timer.stop();
    std::cerr << "\ntotal time emit boundary:" << timer.time() << std::endl;
    std::cerr<< std::endl; 

  }





#if USE_OMP
  void  Voxelization_3::emit_boundary_omp(int current_res, Octree_3& octree) const {
    CGAL::Real_timer timer; 
    timer.start();

    std::cerr << "Start emit boundary with " << m_num_threads -1<<"+1 threads" << std::endl;

    for (uint j = 0; j< ind.size(); j+=3)
      genFacetForwardNormal(current_res,0,1,j,octree.inserter());


    std::vector<VoxelHashSet > buffered_voxels(m_num_threads-1,VoxelHashSet());
    std::vector<VoxelHashSet > buffered_voxels_(m_num_threads-1,VoxelHashSet());

    VoxelHashSet buffer; 
    omp_lock_t lock; 
    omp_init_lock(&lock);
    int tid, i;
    int I = 0; 
    int bid = 0; 
    int finished_threads = 0; 
    bool need_compression = false; 
    bool use_filter = m_downstep; 
    bool update_octree = false; 
    double max_memory = 0.35;
    while(true){    
      update_octree = false;
      Octree_3 old_octree_for_threads(*opt_volume); // for thread safety  
      internal::Contains<Octree_3> contains(&old_octree_for_threads);

#pragma omp parallel private(tid) private(i) num_threads(m_num_threads)
      { // begin pragma 
        tid = omp_get_thread_num(); 
        while(true){ 

          if(need_compression || update_octree ) break; 

          omp_set_lock(&lock); // LOCK START =============


          if(SV::memory_usage()> max_memory) {
            need_compression = true; 
            max_memory+=0.05;
          }
          if(tid!=m_num_threads-1){
            if(I%1000==999) {
              std::cerr << "update octree !" << std::endl;
              update_octree = true;
            } 
            //if(buffered_voxels_[tid].size()==0){         
            // while(I<active_track.size()&&!active_track[I]){I++;}
            i = I;
            I++;
            if(buffered_voxels_[tid].size()==0){
              std::swap(buffered_voxels[tid],buffered_voxels_[tid]);
            }else{
              if(buffered_voxels[tid].size()*100>opt_volume->size()){
                std::cerr << "update octree !" << std::endl;
                update_octree=true; 
              }
            }

            if(i < Track.size())
              std::cerr 
                << " Thread " <<  tid 
                << " taking track " << I 
                << " buffer sizes [" << buffered_voxels[tid].size() 
                << " | " << buffered_voxels_[tid].size() << "]"
                << " Memory usage "  << SV::memory_usage()*100 
                << "%\r" << std::flush; 
            else{
              if(!finished_threads){
                std::cerr << "\n finishing thread: " << tid << std::flush;
              }else{
                std::cerr << " " << tid << std::flush; 
              }
              finished_threads++; 
              omp_unset_lock(&lock);  // LOCK STOP  =============
              break; 
            }
            //}
            //<< " %\r" << std::flush;
          }else{
            bid = (bid+1)%(m_num_threads-1);
            std::swap(buffered_voxels_[bid],buffer); 
          }

          omp_unset_lock(&lock); // LOCK STOP  =============


          if(tid!=m_num_threads-1){
            {
              // Edge Edge Patches
              // puts voxelization of edge-edge-patches-trianlges in voxSet
              if(0 <= i && i < Track.size()-1){

                for (uint j = 0; j< non_manifold_edges.size();j++){
                  if( use_filter ){
                    genNonManifoldEdge(current_res,i, j,  
                        internal::make_filter_iterator(CGAL::inserter(buffered_voxels[tid]),contains));
                  }else{
                    genNonManifoldEdge(current_res,i, j, 
                        CGAL::inserter(buffered_voxels[tid]));
                  }
                }

                for (uint j = 0; j< wingedEdges.size(); j++){

                  if( use_filter ){
                    genPatch(current_res,Track.at(i),Track.at(i+1),wingedEdges.at(j), 
                        internal::make_filter_iterator(CGAL::inserter(buffered_voxels[tid]),contains));
                  }else{
                    genPatch(current_res,Track.at(i),Track.at(i+1),wingedEdges.at(j), 
                        CGAL::inserter(buffered_voxels[tid]));              
                  }
                  // std::back_inserter(buffered_voxels[tid]));
                  // internal::make_filter_iterator(std::back_inserter(buffered_voxels[tid]),contains));  	
                }
              }       
              // Face-Vertex-Facets
              // puts voxelization of Facet-triangles in voxSet    
              if(0 < i && i < Track.size()-1){
                for (uint j = 0; j< ind.size(); j+=3){

                  if( use_filter){
                    genFacet(current_res,i,j,
                        internal::make_filter_iterator(CGAL::inserter(buffered_voxels[tid]),contains));
                  }else{
                    genFacet(current_res,i,j,
                        CGAL::inserter(buffered_voxels[tid]));
                  }
                  //std::back_inserter(buffered_voxels[tid]));
                  //internal::make_filter_iterator(std::back_inserter(buffered_voxels[tid]),contains));		
                }
              }
            }
            // if(buffered_voxels[tid].size()==0) active_track[i] = false;
          }else{
            if(buffer.size()){
              int count = 0; 
              BOOST_FOREACH(const Voxel& v, buffer){
                if(opt_volume->insert(v)) count++;
              }

              if((timer_fill_hull.time() + timer_compute_hull.time()) *20 < timer_emit_boundary.time()){
                //if(double(count)/buffer.size()< 0.01 && double(buffer.size())/opt_volume->size() > 0.1){ 
                if(double(count)/buffer.size()< 0.01){
                  std::cerr << "\n" << count 
                    << "  " << buffer.size() 
                    << "  " << double(count)/buffer.size()
                    << "  " << double(buffer.size())/opt_volume->size()
                    << std::endl; 
                  need_compression = true; 
                }
              }
              buffer.clear();

              }
              if(finished_threads  >= m_num_threads-1) {
                std::cerr << " finishing master thread: " << tid << std::endl;
                break;    
              }
            }
          }
        } // end pragma 


        // collect all remaining voxels 
        BOOST_FOREACH(VoxelHashSet& voxels, buffered_voxels){
          BOOST_FOREACH(const Voxel& v, voxels){
            opt_volume->insert(v); 
          }
          voxels.clear();
        }    
        BOOST_FOREACH(VoxelHashSet& voxels, buffered_voxels_){
          BOOST_FOREACH(const Voxel& v, voxels){ 
            opt_volume->insert(v); 
          }
          voxels.clear();
        }    

        if(I >= Track.size()){
          for (uint j = 0; j< ind.size(); j+=3)
            genFacetForwardNormal(current_res,Track.size()-1,Track.size()-2,j,octree.inserter());
          break; 
        }else{
          if(need_compression){
            for (uint j = 0; j< ind.size(); j+=3)
              genFacetForwardNormal(current_res,I,I-1,j,octree.inserter()); 
            timer_emit_boundary.stop();
            compress(current_res, false);
            timer_emit_boundary.start();
            need_compression = false; // continue with all threads
          }
          use_filter = true; 
        }
      } // outer while loop that may launch compression 


      // m_seed = get_seed(current_res);
      timer.stop();
      std::cerr << "\n opt_volume->size(): " <<  opt_volume->size() << std::endl;
      std::cerr << "\n total time emit boundary:" << timer.time() << std::endl;
      std::cerr<< std::endl; 
    }
#endif

  } // namespace 
