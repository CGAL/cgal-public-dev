// Copyright (c) 2005-2007 Tel-Aviv University (Israel).
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
// Author(s): Michael Hemmer        <hemmer@mpi-inf.mpg.de>
//
//

// This files illustrates the TOP/BOTTOM BOUNDARY INSERTION BUG 
// The bug arises when inserting an x_monotone_curve that has both endpoints 
// at the TOP BOUNDARY (or alternatively both on the BOTTOM BOUNDARY)

// First analysis: 
// The problem seems to be due to the fact that the insertion of ARR_MIN_END 
// already splits an fictitious halfedge, but does not insert the new edge since 
// the ARR_MAX_END is not yet located. During the location of the ARR_MAX_END
// it happens that it is compared to the new vertex that corresponds to 
// ARR_MIN_END, which is not possible, since the halfedge that could provide 
// the required supporting curve is not available yet. 

#define CGAL_AK_ENABLE_DEPRECATED_INTERFACE 1 


#include <CGAL/basic.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_curve_kernel_2.h> // template arg to CKvA_2
#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_impl.h> // traits for Arr_2
#include <CGAL/Arrangement_2.h>
#include <CGAL/Algebraic_kernel_d_1.h>


template<typename Poly_> Poly_ from_string(const char* s) {
  std::stringstream ss(s);
  Poly_ f;
  ss >> f;
  return f;
}

int main(int argc, char** argv)
{
 
  // get default Integer / AK_1 / AK_2/ CK_2
  typedef CGAL::Arithmetic_kernel::Integer Integer;
  typedef CGAL::Algebraic_kernel_d_1<Integer>                  AK_1;
  typedef CGAL::Algebraic_curve_kernel_2<AK_1>                 AK_2;
  typedef CGAL::Curved_kernel_via_analysis_2< AK_2 >           CK_2;

  // define curved kernel 
  CK_2 ck; 
  
  typedef AK_2::Polynomial_2 Polynomial_2;   
  Polynomial_2 poly = from_string<Polynomial_2>(
      "P[4(0,P[4(0,-17161)(1,22532)(2,-8554)(3,1548)(4,-81)])(1,P[4(0,-68644)(1,79124)(2,-30460)(3,5124)(4,-216)])(2,P[4(0,-102966)(1,102180)(2,-27096)(3,1044)(4,558)])(3,P[4(0,-68644)(1,57116)(2,2972)(3,-8604)(4,1656)])(4,P[4(0,-17161)(1,11528)(2,8162)(3,-6072)(4,639)])]");
  
  
 
  // construct curve object
  typedef CK_2::Curve_2 Curve_2;
  Curve_2 curve = ck.kernel().construct_curve_2_object()(poly);

  // construct x monoton curves (in CGAL Object)
  std::list<CGAL::Object> x_mono_objects;
  ck.make_x_monotone_2_object()
    (curve, std::back_inserter(x_mono_objects));
  
  // take out x monotone curves 
  typedef CK_2::X_monotone_curve_2 X_monotone_curve_2;
  std::vector<X_monotone_curve_2> xcurves; 
 
  X_monotone_curve_2 x_monotone_curve_2;
  for (std::list<CGAL::Object>::iterator it = x_mono_objects.begin();
       it != x_mono_objects.end(); ++it){
    assign(x_monotone_curve_2, *it);
    xcurves.push_back(x_monotone_curve_2);
  }

  // define Arrangement 
  CGAL::Arrangement_2<CK_2> arr; 

  // insert exactly one xcurve into arrangement
  // namely an xcurve with both endpoints at ARR_TOP_BOUNDARY
  // that is, there are two vertical aymptotes. 
  std::vector<X_monotone_curve_2>::iterator it; 
  for(it = xcurves.begin(); it != xcurves.end(); it++){
    if( it->location(CGAL::ARR_MIN_END) == CGAL::ARR_TOP_BOUNDARY && 
        it->location(CGAL::ARR_MAX_END) == CGAL::ARR_TOP_BOUNDARY ){
      std::cout<< "insert at top top " << std::endl;
      CGAL::insert(arr,*it);
    }
     
  }


  
  return 0;
}
