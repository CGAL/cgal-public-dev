// Copyright (c) 2005-2009 Tel-Aviv University (Israel).
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
// $URL: $
// $Id:
// 
//
// Author(s): Ophir Setter          <ophirset@post.tau.ac.il>
//            Michael Hemmer        <hemmer@mpi-inf.mpg.de>
//


// SVCET_3 = Single_voronoi_cell_envelope_traits_3

#ifndef CGAL_VDL3_SVCET_33_INTERSECT_2
#define CGAL_VDL3_SVCET_33_INTERSECT_2

#include <CGAL/VDOL_3/SVCET_3_functors/SVCET_3_functor_base.h>
#include <CGAL/Arr_enums.h>

namespace CGAL{ 
namespace VDOL_3 {
namespace SVCET_3_functors {

template < class SVCET_3 >
class Intersect_2 : public 
Single_voronoi_cell_envelope_traits_3_functor_base< SVCET_3 > {

public:
  typedef SVCET_3 Single_voronoi_cell_envelope_traits_3;
  typedef Single_voronoi_cell_envelope_traits_3_functor_base< SVCET_3 > Base;
  
  CGAL_SNAP_SVCET_3_TYPEDEFS;
  
  Intersect_2(const SVCET_3 *traits): Base(traits) {};

  //! Overloading Intersect_2 to throw exceptions on discontinuity events.
  /*! We overload the Intersect_2 of the arrangement traits to throw an
    exception when there are intersection points at the discontinuity line.
    We know that there is an event at the discontinuity line if there are
    two vertical asymtotes with the same x-value.
  */
  template<typename OutputIterator>
  OutputIterator operator()(const X_monotone_curve_2& xcv1,
      const X_monotone_curve_2& xcv2,
      OutputIterator oi)
  {
    
    // The boundary curves of the $xy$-monotone surfaces are also considered
    // here. In that case, it is clear that there are no intersections on the
    // boundary. We know that a curve is a boundary curve if it is a linear 
    // curve.
    // \todo Check with Michael.
    
    if(xcv1.id() == xcv2.id() || xcv1.do_overlap(xcv2)){
        return this->svcet_3()->Curved_kernel_2::intersect_2_object() 
        (xcv1, xcv2, oi);
    }
    
    if (CGAL::total_degree(xcv1.curve().polynomial_2()) > 1 &&
        CGAL::total_degree(xcv2.curve().polynomial_2()) > 1)
    {
      bool xcv1_min_asymptot = false, xcv1_max_asymptot = false;
      bool xcv2_min_asymptot = false, xcv2_max_asymptot = false;
      
      typename Curved_kernel_2::Parameter_space_in_x_2 parameter_space_in_x_2 =
        this->svcet_3()->parameter_space_in_x_2_object();
      typename Curved_kernel_2::Parameter_space_in_y_2 parameter_space_in_y_2 =
        this->svcet_3()->parameter_space_in_y_2_object();
      
      
      if ((parameter_space_in_x_2(xcv1, ARR_MIN_END) == ARR_INTERIOR) &&
          (parameter_space_in_y_2(xcv1, ARR_MIN_END) != ARR_INTERIOR))
        xcv1_min_asymptot = true;
      
      if ((parameter_space_in_x_2(xcv1, ARR_MAX_END) == ARR_INTERIOR) &&
          (parameter_space_in_y_2(xcv1, ARR_MAX_END) != ARR_INTERIOR))
        xcv1_max_asymptot = true;
      
      if ((parameter_space_in_x_2(xcv2, ARR_MIN_END) == ARR_INTERIOR) &&
          (parameter_space_in_y_2(xcv2, ARR_MIN_END) != ARR_INTERIOR))
        xcv2_min_asymptot = true;
      
      if ((parameter_space_in_x_2(xcv2, ARR_MAX_END) == ARR_INTERIOR) &&
          (parameter_space_in_y_2(xcv2, ARR_MAX_END) != ARR_INTERIOR))
        xcv2_max_asymptot = true;
      
      if (xcv1_min_asymptot && xcv2_min_asymptot)
      {
        if (this->svcet_3()->kernel().compare_x_2_object()(xcv1.curve_end_x(ARR_MIN_END),
                                     xcv2.curve_end_x(ARR_MIN_END)) == EQUAL)
          throw Event_at_discontinuity_exception("two curves intersect on a "
                                                 "discontinuity line "
                                                 "(min, min)");
      }
      
      if (xcv1_min_asymptot && xcv2_max_asymptot)
      {
        if (this->svcet_3()->kernel().compare_x_2_object()(xcv1.curve_end_x(ARR_MIN_END),
                                     xcv2.curve_end_x(ARR_MAX_END)) == EQUAL)
          throw Event_at_discontinuity_exception("two curves intersect on a "
                                                 "discontinuity line "
                                                 "(min, max)");
      }
      
      if (xcv1_max_asymptot && xcv2_min_asymptot)
      {
        if (this->svcet_3()->kernel().compare_x_2_object()(xcv1.curve_end_x(ARR_MAX_END),
                                     xcv2.curve_end_x(ARR_MIN_END)) == EQUAL)
          throw Event_at_discontinuity_exception("two curves intersect on a "
                                                 "discontinuity line "
                                                 "(max, min)");
      }
      
      if (xcv1_max_asymptot && xcv2_max_asymptot)
      {
        if (this->svcet_3()->kernel().compare_x_2_object()(xcv1.curve_end_x(ARR_MAX_END),
                                     xcv2.curve_end_x(ARR_MAX_END)) == EQUAL)
          throw Event_at_discontinuity_exception("two curves intersect on a "
                                                 "discontinuity line "
                                                 "(max, max)");
      }
    }
      
    return this->svcet_3()->Curved_kernel_2::intersect_2_object() 
      (xcv1, xcv2, oi);
  }
};

} // namespace SVCET_3_functors
} // namespace VDOL_3 
} // namespace CGAL

#endif // CGAL_VDL3_SVCET_33_INTERSECT_2
