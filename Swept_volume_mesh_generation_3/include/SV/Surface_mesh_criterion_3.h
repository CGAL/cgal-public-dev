#ifndef SV_SURFACE_MESH_CRITERION_3_H
#define SV_SURFACE_MESH_CRITERION_3_H


#include <CGAL/Surface_mesh_default_criteria_3.h>

namespace SV{

template <class Tr, class Volume_> 
class Surface_mesh_criterion_3{

public:  
  typedef Volume_ Volume;
  typedef CGAL::Surface_mesh_default_criteria_3<Tr> Default_criteria_3;

private:
  Volume m_volume; 
  Default_criteria_3 m_default_criteria_3; 
  
public:
  Surface_mesh_criterion_3(
      const Volume& volume, 
      const Default_criteria_3& default_criteria)
    :m_volume(volume), m_default_criteria_3(default_criteria)
  {}
  
  typedef typename Default_criteria_3::Facet Facet; 
  typedef typename Default_criteria_3::Quality Quality; 
  
  bool is_bad ( Facet f, Quality& q) const {
    typedef typename Tr::Point Point_3;
    const Point_3& p1 = f.first->vertex((f.second+1)&3)->point();
    const Point_3& p2 = f.first->vertex((f.second+2)&3)->point();
    const Point_3& p3 = f.first->vertex((f.second+3)&3)->point();
    
    // if it is bad (return in a boost optional )
    if(m_volume.compute_facet_badness(p1,p2,p3)){
      return true;
    }
    return m_default_criteria_3.is_bad(f,q);
  }
};



} // namespace SV
#endif // SV_SURFACE_MESH_CRITERION_3_H
