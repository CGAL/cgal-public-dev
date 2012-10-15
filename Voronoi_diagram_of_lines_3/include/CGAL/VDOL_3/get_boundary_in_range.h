// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Michael Hemmer <hemmer@mpi-inf.mpg.de>
//
// =============================================================================

#ifndef CGAL_GET_BOUNDARY_IN_X_RANGE_H
#define CGAL_GET_BOUNDARY_IN_X_RANGE_H

#include <CGAL/basic.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>


namespace CGAL {


namespace VDOL_3 {

// This function is a template in CurveKernel_2 (e.g. Algebraic_kernel_d_2)
// see documetation of Curved_kernel_via_analysis_2 


template <class CurveKernel_2> 
typename CurveKernel_2::Boundary 
get_boundary_in_x_range_interior(
    const Curved_kernel_via_analysis_2<CurveKernel_2>* ckva, 
    const typename Curved_kernel_via_analysis_2<CurveKernel_2>::X_monotone_curve_2& arc){
  return arc.boundary_in_x_range_interior();
}

template <class CurveKernel_2> 
typename CurveKernel_2::Boundary 
get_boundary_in_y_range_interior(
    const Curved_kernel_via_analysis_2<CurveKernel_2>* ckva, 
    const typename Curved_kernel_via_analysis_2<CurveKernel_2>::X_monotone_curve_2& arc){
  
  typedef typename CurveKernel_2::Boundary Boundary;
  
  // may also be implemented for other arcs but not needed within this context 
  CGAL_precondition(arc.is_vertical()); 
  
  CurveKernel_2 ck = ckva->kernel(); 
  Boundary y_coord(0);

  if (arc.is_finite(ARR_MIN_END)){
      if (arc.is_finite(ARR_MAX_END)){
        y_coord = ck.boundary_between_y_2_object() (
            arc.curve_end(ARR_MIN_END).xy(),
            arc.curve_end(ARR_MAX_END).xy());
      }else{
        y_coord = 
          ck.upper_boundary_y_2_object()(arc.curve_end(ARR_MIN_END).xy())+Boundary(1);
      }
  } else {
    if (arc.is_finite(ARR_MAX_END)){
      y_coord = 
        ck.lower_boundary_y_2_object()(arc.curve_end(ARR_MAX_END).xy())-Boundary(1);
    }
  }
  return y_coord;
}

} // namespace VDOL_3
} //namespace CGAL

#endif // CGAL_GET_BOUNDARY_IN_X_RANGE_H
