#ifndef SV_SURFACE_MESH_TRAITS_3_H
#define SV_SURFACE_MESH_TRAITS_3_H


#include <SV/Mesh_domain_3.h> 

namespace SV{

template<class Kernel> class Swept_volume_with_vhull_3; 
template<class SweptVolume> class Surface_mesh_traits_3; 

template<class Kernel_>
class Surface_mesh_traits_3<Swept_volume_with_vhull_3<Kernel_> >{

public:
  typedef Kernel_ Kernel; 
  typedef Swept_volume_with_vhull_3<Kernel> Swept_volume; 
  typedef Mesh_domain_3<Swept_volume> Mesh_domain; 
  typedef Surface_mesh_traits_3<Swept_volume> Self; 

private:  
  const Swept_volume  m_mesh_domain 
public:
  const Swept_volume& mesh_domain() const {return m_mesh_domain;}
  
  Surface_mesh_traits_3(Swept_volume& swept_volume):
    m_mesh_domain(swept_volume){
  }  
  
public:

  typedef Mesh_domain::Intersect_3 Intersect_3;
  typedef Mesh_domain::Construct_initial_points Construct_initial_points; 
  
  Intersect_3  intersect_3_object(){
    me
  }
    Construct_initial_points 	mesh_domain.construct_initial_points_object () 
    
    

};



} // namespace SV


#endif // SV_SURFACE_MESH_TRAITS_3_H
