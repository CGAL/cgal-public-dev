#ifndef CGAL_VDOL_APPROXIMATION_INFO_H
#define CGAL_VDOL_APPROXIMATION_INFO_H

#include <CGAL/basic.h>

namespace CGAL {
namespace VDOL_3 { 

struct Approximation_info{
  double sradius_bounding_sphere;  
  double sradius_clipping_sphere;
  double sdistance_generation; 
  double sdistance_final;
  double point_weight;    
  
  Approximation_info(){
    sradius_bounding_sphere = 110; 
    sradius_clipping_sphere = 100; 
    sdistance_generation = 0.5;
    sdistance_final = 0.61;
    point_weight   = 0.6; 
  }

  
};

} // namespace VDOL_3
} // namespace CGAL 


#endif // CGAL_VDOL_APPROXIMATION_INFO_H
