#ifndef SV_SURFACE_MESH_CRITERION_H
#define SV_SURFACE_MESH_CRITERION_H


#include <CGAL/Surface_mesh_default_criteria_3.h>

namespace SV{

template <class TR, class Volume> 
class Surface_mesh_criterion{
  
  typedef CGAL::Surface_mesh_default_criteria_3<Tr> Default_criteria_3;
  typedef Facet_criterion; 
  Volume m_volume; 
  Default_criteria_3 m_default_criteria_3; 
  
  Surface_mesh_criterion(const Volume& volume, const Default_criteria_3& default_criteria)
    :m_volume(volume), m_default_criteria_3(default_criteria)
  {}

  typedef Default_criteria_3::Facet Facet; 
  typedef Default_criteria_3::Quality Quality; 
  
  bool 	criteria.is_bad ( Facet f, Quality& q) const {
    typedef typename Tr::Point Point_3;
    const Point_3& p1 = f.first->vertex((f.second+1)&3)->point();
    const Point_3& p2 = f.first->vertex((f.second+2)&3)->point();
    const Point_3& p3 = f.first->vertex((f.second+3)&3)->point();
    
    // if it is bad (return in a boost optional )
    if(m_volume.compute_facet_badness(p1,p2,p3)){
      return true;
    }
#if 1
    q = Quality(std::max(
        std::max(CGAL::squared_distance(p1,p2),CGAL::squared_distance(p2,p3)),
        CGAL::squared_distance(p1,p2)));
    return false; 
#else
    return m_default_criteria_3(f,q);
#endif

  }
}



} // namespace SV
#endif // SV_SURFACE_MESH_CRITERION_H
