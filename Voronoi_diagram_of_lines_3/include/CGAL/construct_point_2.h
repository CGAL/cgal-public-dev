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
// $Id: $
// 
//
// Author(s): Ophir Setter          <ophirset@post.tau.ac.il>
//            Michael Hemmer        <hemmer@mpi-inf.mpg.de>
//
//

// Helper function to construct a (ratinal) CK_2::Point_2 

#ifndef CGAL_CONSTRUCT_CK2_POINT_2_H
#define CGAL_CONSTRUCT_CK2_POINT_2_H

namespace CGAL {

template <class SVCET_3>
typename SVCET_3::Point_2 construct_point_2(
    const SVCET_3* svcet_3,
    const typename SVCET_3::FT& x, 
    const typename SVCET_3::FT& y){
  typedef typename SVCET_3::Point_2::Coordinate_1 Coordinate_1; 
  return construct_point_2(svcet_3,Coordinate_1(x),y);
}
    

template <class SVCET_3>
typename SVCET_3::Point_2 construct_point_2(
    const SVCET_3* svcet_3,
    const typename SVCET_3::Point_2::Coordinate_1& x, 
    const typename SVCET_3::FT& y){
    
  typedef typename SVCET_3::FT Rational;
  typedef typename SVCET_3::Curved_kernel_2 Curved_kernel_2;
  typedef typename SVCET_3::Algebraic_kernel_d_2 Algebraic_kernel_d_2;
  typedef typename SVCET_3::Algebraic_kernel_d_1 Algebraic_kernel_d_1;
  typedef typename SVCET_3::Point_2 Point_2;
  
  typedef typename Algebraic_kernel_d_2::Polynomial_2 Polynomial_2;
  typedef typename Algebraic_kernel_d_1::Polynomial_1 Polynomial_1;
  
  // TODO: use CGAL::Polynomial_traits_d
  //typedef CGAL::Polynomial_traits_d<Polynomial_2> PT_2;
  //typedef CGAL::Polynomial_traits_d<Polynomial_1> PT_1;
 
  typedef typename CGAL::Fraction_traits< Rational >::Numerator_type Integer; 
  typedef typename SVCET_3::Point_2::Coordinate_1 Coordinate_1;

  typename CGAL::Fraction_traits< Rational >::Decompose decompose;
    
  // Construct horizontal line that supports the point. 
  Integer yn,yd;
  decompose(y, yn, yd);
  Polynomial_2 poly_y(Polynomial_1(-yn), Polynomial_1(yd));
  typename Algebraic_kernel_d_2::Curve_analysis_2 horizontal_line =
    svcet_3->kernel().construct_curve_2_object()(poly_y);

  Point_2 point = svcet_3->construct_point_2_object()(x,horizontal_line,0);
  return point; 
}


template <class SVCET_3>
typename SVCET_3::Point_2 construct_point_2(
    const SVCET_3* svcet_3, 
    const Arr_parameter_space& psx, 
    const typename SVCET_3::FT& y){
  CGAL_precondition(psx == ARR_LEFT_BOUNDARY || psx == ARR_RIGHT_BOUNDARY);
  
  typedef typename SVCET_3::FT Rational;
  typedef typename SVCET_3::Algebraic_kernel_d_2 Algebraic_kernel_d_2;
  typedef typename SVCET_3::Algebraic_kernel_d_1 Algebraic_kernel_d_1;
  typedef typename SVCET_3::Point_2 Point_2;
  
  typedef typename Algebraic_kernel_d_2::Polynomial_2 Polynomial_2;
  typedef typename Algebraic_kernel_d_1::Polynomial_1 Polynomial_1;

  typename CGAL::Fraction_traits< Rational >::Decompose decompose;
  // Construct horizontal line that supports the point. 
  typename CGAL::Fraction_traits< Rational >::Numerator_type yn,yd;
  decompose(y, yn, yd);
  Polynomial_2 poly_y(Polynomial_1(-yn), Polynomial_1(yd));
  typename Algebraic_kernel_d_2::Curve_analysis_2 horizontal_line =
    svcet_3->kernel().construct_curve_2_object()(poly_y);
  Point_2 point(
      (psx == ARR_LEFT_BOUNDARY)?ARR_MIN_END:ARR_MAX_END,
      horizontal_line,0);
  CGAL_assertion(point.location() == psx);
  return point; 
}

template <class SVCET_3>
typename SVCET_3::Point_2 construct_point_2(
    const SVCET_3* svcet_3,  
    const typename SVCET_3::FT& x,
    const Arr_parameter_space& psy){
  CGAL_precondition(psy == ARR_TOP_BOUNDARY || psy == ARR_BOTTOM_BOUNDARY);
  
  typedef typename SVCET_3::FT Rational;
  typedef typename SVCET_3::Algebraic_kernel_d_2 Algebraic_kernel_d_2;
  typedef typename SVCET_3::Algebraic_kernel_d_1 Algebraic_kernel_d_1;
  typedef typename SVCET_3::Point_2 Point_2;
  typedef typename Point_2::Coordinate_1 Coordinate_1;
  
  typedef typename Algebraic_kernel_d_2::Polynomial_2 Polynomial_2;
  typedef typename Algebraic_kernel_d_1::Polynomial_1 Polynomial_1;

  typename CGAL::Fraction_traits< Rational >::Decompose decompose; 
  typedef typename CGAL::Fraction_traits< Rational >::Numerator_type Integer; 

  Integer xn,xd;
  decompose(x, xn, xd);
  Polynomial_1 poly_x(-xn, xd);
  
  std::vector<Coordinate_1> roots;  
  svcet_3->kernel().solve_1_object()(poly_x,std::back_inserter(roots),true);
  CGAL_postcondition(roots.size()==1); 
  typename Algebraic_kernel_d_2::Curve_analysis_2 vertical_line =
    svcet_3->kernel().construct_curve_2_object()(poly_x);
  Point_2 point( roots[0],vertical_line, 
      (psy == ARR_BOTTOM_BOUNDARY)?ARR_MIN_END:ARR_MAX_END);
  CGAL_assertion(point.location() == psy);
  return point; 
}

} //namespace CGAL

#endif // CGAL_CONSTRUCT_CK2_POINT_2_H
