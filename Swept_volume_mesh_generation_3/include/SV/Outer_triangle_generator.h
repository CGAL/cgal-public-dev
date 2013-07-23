#ifndef SV_OUTER_TRIANGLE_GENERATOR_3_H
#define SV_OUTER_TRIANGLE_GENERATOR_3_H

#include <boost/foreach.hpp>
#include "vecmath.h"


namespace SV{


template< class Voxelization_> 
class Outer_triangle_generator{
public:
  typedef Voxelization_ Voxelization; 
  typedef typename Voxelization::wingedEdge wingedEdge; 
  typedef typename Voxelization::edge Edge; 
  
  Voxelization m_voxelization;
  std::map<Vector3d,int> index_map; 
  std::vector<int> triangles; // three consecutive ints make up a triangle 
  
  Outer_triangle_generator(Voxelization voxelizaion)
    :m_voxelization(voxelizaion){
    initialize();
  }
  
  int get_index(const Vector3d& p){
    int old_size = index_map.size();
    int& index = index_map[p]; 
    if(old_size != index_map.size()){
      index = index_map.size()-1;
    }
    return index; 
  }
  
  
  void genPatch(const Matrix4d& A, const Matrix4d B, const wingedEdge &edge){
    Vector3d a,b,c,d;  
    if(m_voxelization.needPatch(A,B,edge,a,b,c,d)){
      triangles.push_back(get_index(a));
      triangles.push_back(get_index(b));
      triangles.push_back(get_index(c));

      triangles.push_back(get_index(b));
      triangles.push_back(get_index(c));
      triangles.push_back(get_index(d));
    }
  }

  void genFacet(int i, int j){
    Vector3d a,b,c;  
    if(m_voxelization.needFacet(i,j,a,b,c)){
      triangles.push_back(get_index(a));
      triangles.push_back(get_index(b));
      triangles.push_back(get_index(c));
    }
  }
  
 
  void genFacetForwardNormal(int i0, int i1, int j){
    Vector3d a,b,c;  
    if(m_voxelization.needFacetForwardNormal(i0,i1,j,a,b,c)){
      triangles.push_back(get_index(a));
      triangles.push_back(get_index(b));
      triangles.push_back(get_index(c));
    }
  }
  
  void genNonManifoldEdge(int i, int j){

    const Matrix4d& D0 = m_voxelization.Track[i];
    const Matrix4d& D1 = m_voxelization.Track[i+1];
    const Edge& edge   = m_voxelization.non_manifold_edges[j]; 
    const Vector3d& a  = m_voxelization.P[edge.first ]; 
    const Vector3d& b  = m_voxelization.P[edge.second]; 
    
    Vector3d a0  = D0*a ;
    Vector3d b0  = D0*b ;
    Vector3d a1  = D1*a ;
    Vector3d b1  = D1*b ;

    // just draw all possible triangles .-) 
     
    //oit = emitFacet(res, a0,a1,b0,oit);
    triangles.push_back(get_index(a0));
    triangles.push_back(get_index(a1));
    triangles.push_back(get_index(b0));
    
    // oit = emitFacet(res, a0,a1,b1,oit);
    triangles.push_back(get_index(a0));
    triangles.push_back(get_index(a1));
    triangles.push_back(get_index(b1));

    // oit = emitFacet(res, b0,b1,a0,oit);
    triangles.push_back(get_index(b0));
    triangles.push_back(get_index(b1));
    triangles.push_back(get_index(a0));

    // oit = emitFacet(res, b0,b1,a1,oit); 
    triangles.push_back(get_index(b0));
    triangles.push_back(get_index(b1));
    triangles.push_back(get_index(a1));
  }
  

  void initialize(){
    std::cerr << "Start emit boundary: no multi threading" << std::endl;
    int i = 0; 
      
    for (uint j = 0; j< m_voxelization.ind.size(); j+=3)
      genFacetForwardNormal(0,1,j);

    while(i < m_voxelization.Track.size()){
      // Edge Edge Patches
      // puts voxelization of edge-edge-patches-trianlges in voxSet
      if(0 <= i && i < m_voxelization.Track.size()-1){
        for (uint j = 0; j< m_voxelization.non_manifold_edges.size();j++){
          genNonManifoldEdge(i, j);
        }
        for (uint j = 0; j< m_voxelization.wingedEdges.size(); j++){
          genPatch(
              m_voxelization.Track.at(i),
              m_voxelization.Track.at(i+1),
              m_voxelization.wingedEdges.at(j));         
        }
      }       
      // Face-Vertex-Facets
      // puts voxelization of Facet-triangles in voxSet    
      if(0 < i && i < m_voxelization.Track.size()-1){
        for (uint j = 0; j< m_voxelization.ind.size(); j+=3){
          genFacet(i,j);		
        }
      }
      i++;
      std::cerr << " TRACK " << i << " Memory usage " << SV::memory_usage()*100 << " %\r" << std::flush;
    }
    // put the full solid into octhull to close swept volume hull for compression 
    for (uint j = 0; j< m_voxelization.ind.size(); j+=3)
      genFacetForwardNormal(i-1,i-2,j);
  
  }

  
  void write_filtered_triangle_soup(const char* filename) const {
    std::ofstream os(filename);
    os << "OFF" << std::endl;
    os << index_map.size() << " " << triangles.size()/3 << " 0 " <<std::endl;

    std::vector<Vector3d> points(index_map.size()); 
    typedef std::pair<Vector3d,int> VI_PAIR; 
    BOOST_FOREACH(VI_PAIR vi, index_map){
      points[vi.second]=vi.first;
    }
    BOOST_FOREACH(Vector3d v, points){
      os << v[0] << " " << v[1] << " " << v[2] << std::endl; 
    }
    int counter = 0; 
    BOOST_FOREACH(int t, triangles){ 
      if (counter == 0) os << 3 << " "; 
      os << t << " ";
      counter=(counter +1)%3; 
      if (counter == 0) os << std::endl; 
    }
  }
};

} // namespace SV 

#endif // SV_OUTER_TRIANGLE_GENERATOR_3_H 
