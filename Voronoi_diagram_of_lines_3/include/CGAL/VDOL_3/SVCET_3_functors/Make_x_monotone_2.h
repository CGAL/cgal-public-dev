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

#ifndef CGAL_VDL3_SVCET_33_MAKE_X_MONOTONE_2
#define CGAL_VDL3_SVCET_33_MAKE_X_MONOTONE_2

#include <CGAL/VDOL_3/SVCET_3_functors/SVCET_3_functor_base.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Object.h>
#include <boost/foreach.hpp>

namespace CGAL{ 
namespace VDOL_3 {
namespace SVCET_3_functors {

template < class SVCET_3 >
class Make_x_monotone_2 : public 
Single_voronoi_cell_envelope_traits_3_functor_base< SVCET_3 > {

public:
  typedef SVCET_3 Single_voronoi_cell_envelope_traits_3;
  typedef Single_voronoi_cell_envelope_traits_3_functor_base< SVCET_3 > Base;
  
  CGAL_SNAP_SVCET_3_TYPEDEFS;
  
  Make_x_monotone_2(const SVCET_3 *traits): Base(traits) {};

  template<typename OutputIterator>
  OutputIterator operator()(const Curve_2& cv, OutputIterator oi)
  {
#if 0
#warning Make_x_monotone_2 splits for insertion bug 
    // std::cout << " called Make_x_monotone_2 "<< std::endl; 
    
    std::vector<CGAL::Object> objects; 
    this->svcet_3()->Curved_kernel_2::make_x_monotone_2_object()
      (cv,std::back_inserter(objects));
    
    X_monotone_curve_2 xcv; 
    X_monotone_curve_2 xcv_a,xcv_b; 
    
    std::vector<X_monotone_curve_2> xcurves;  
    BOOST_FOREACH(cosnt CGAL::Object& obj, objects){
      if( CGAL::assign(xcv,obj) && 
          xcv.location(ARR_MIN_END) == xcv.location(ARR_MAX_END) && 
          xcv.location(ARR_MAX_END) != ARR_INTERIOR){
          
        Point_2 p = this->svcet_3()->construct_interior_vertex_2_object()(xcv);
        this->svcet_3()->split_2_object()(xcv,p,xcv_a, xcv_b);
        *oi++ = CGAL::make_object(xcv_a);
        *oi++ = CGAL::make_object(xcv_b);
        CGAL_postcondition(
            !(xcv_a.location(ARR_MIN_END) == xcv_a.location(ARR_MAX_END) && 
                xcv_a.location(ARR_MAX_END) != ARR_INTERIOR));
        CGAL_postcondition(
            !(xcv_b.location(ARR_MIN_END) == xcv_b.location(ARR_MAX_END) && 
                xcv_b.location(ARR_MAX_END) != ARR_INTERIOR));
      }else{
        *oi++=obj; 
      }
    }
    
    return oi; 
#else
    return this->svcet_3()->Curved_kernel_2::make_x_monotone_2_object()(cv, oi);
#endif
  }
};

} // namespace SVCET_3_functors
} // namespace VDOL_3 
} // namespace CGAL

#endif // CGAL_VDL3_SVCET_33_MAKE_X_MONOTONE_2
