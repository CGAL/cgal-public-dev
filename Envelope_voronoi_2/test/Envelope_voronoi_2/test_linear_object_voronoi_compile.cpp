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
// $Id:  $
// 
//
// Author(s): Ophir Setter          <ophir.setter@post.tau.ac.il>
//

/*!
  \file   test_linear_object_voronoi_compile.cpp
  \brief  The test file only checks that the segment Voronoi diagram properly
  compiles.
*/

#include <list>

#include <CGAL/Arithmetic_kernel.h>

#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_impl.h>
#include <CGAL/Algebraic_curve_kernel_2_generator.h>

#include <CGAL/envelope_3.h>
#include <CGAL/Envelope_voronoi_traits_2/Linear_objects_Voronoi_traits_2.h>


typedef CGAL::Arithmetic_kernel                       AK;
typedef AK::Integer                                   Integer;
typedef AK::Rational                                  Rational;

typedef Integer                                       Coefficient;
typedef CGAL::Algebraic_curve_kernel_2_generator<Coefficient>::
Algebraic_curve_kernel_with_qir_and_bitstream_2       Algebraic_curve_kernel_2;
typedef CGAL::Curved_kernel_via_analysis_2< Algebraic_curve_kernel_2 > 
                                                      Curved_kernel_2;
typedef CGAL::Linear_objects_Voronoi_traits_2< 
  Curved_kernel_2 >                                   Segment_traits_2;
typedef Segment_traits_2::Surface_3                   Surface_3;
typedef CGAL::Envelope_diagram_2<Segment_traits_2>    Segment_diagram_2;


int main( int argc, char **argv )
{
  Segment_diagram_2 segment_diagram;
  std::list<Surface_3> segments;

  CGAL::lower_envelope_3 (segments.begin(), segments.end(),
                          segment_diagram);
  return 0;

}


