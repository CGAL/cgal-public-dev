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
// This file is provided AS IS with NO
// WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id:
// 
//
// Author(s): Ophir Setter          <ophirset@post.tau.ac.il>
//            Michael Hemmer        <hemmer@mpi-inf.mpg.de>
//
//

#define CGAL_AK_ENABLE_DEPRECATED_INTERFACE 1 

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/macros.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d_1.h>

#include <CGAL/Algebraic_kernel_d/Algebraic_curve_kernel_2.h> // template arg to CKvA_2
#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_impl.h> // traits for Arr_2


#include <CGAL/VDOL_3/Single_voronoi_cell_envelope_traits_3.h>
#include <CGAL/envelope_3.h>
#include <CGAL/Timer.h>
#include <CGAL/Random.h>
 
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/test_envelop_for_lines.h>

int main(int argc, char** argv)
{
  CGAL::set_pretty_mode(std::cout);
  CGAL::set_pretty_mode(std::cerr);

  // get default arithmetic 
  typedef CGAL::Arithmetic_kernel AK;
  typedef AK::Integer Integer;
  typedef AK::Rational Rational;

  // define linear kernel 
  typedef CGAL::Cartesian<Rational>                       Linear_kernel_3;
  typedef Linear_kernel_3::Line_3                         Line_3;
  typedef Linear_kernel_3::Point_3                        Point_3;
  
  // define algebraic kernel / curved kernel 
  typedef CGAL::Algebraic_kernel_d_1<Integer>             AK_1;
  typedef CGAL::Algebraic_curve_kernel_2<AK_1>            Algebraic_curve_kernel_2;
  
  // define curved kernel 
  typedef CGAL::Curved_kernel_via_analysis_2< Algebraic_curve_kernel_2 > 
    Curved_kernel_2;

  // define traits for lower envelop diagram 
  typedef CGAL::VDOL_3::
    Single_voronoi_cell_envelope_traits_3
    <Linear_kernel_3, Curved_kernel_2>                    SVCET_3;
  typedef CGAL::Envelope_diagram_2<SVCET_3>          Envelope_diagram_2;
  



  CGAL::test_envelop_for_lines<SVCET_3>(
      Line_3 (Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
      Line_3 (Point_3( 0, 0, 0), Point_3( 2, -1, 1)),
      Line_3 (Point_3( 1, 0, 0), Point_3( 3, 2, -2)));

  std::cout<< " TEST: Generic tests" << std::endl; 
  CGAL::test_envelop_for_lines<SVCET_3>(
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
      Line_3(Point_3( 2,-8,-9), Point_3(-2,-1, 1)));  

  CGAL::test_envelop_for_lines<SVCET_3>(
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
      Line_3(Point_3( 1, 1, 1), Point_3( 2, 3, 4)),
      Line_3(Point_3( 2, 2, 2), Point_3( 2, 4, 2)));

  CGAL::test_envelop_for_lines<SVCET_3>( 
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
      Line_3(Point_3(-5, 7, 2), Point_3( 5,-2, 1)),
      Line_3(Point_3( 2,-8,-9), Point_3(-2,-1, 1)),
      Line_3(Point_3( 4,-3, 7), Point_3( 8, 4,-8)));
  


  // special configurations for one line 
  std::cout<< "TEST: line in direction of v2 " << std::endl; 
  CGAL::test_envelop_for_lines<SVCET_3>(
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
      Line_3(Point_3( 1, 1, 1), Point_3( 1, 2, 1)));
  
  std::cout<< "TEST: line in direction of v3 "  << std::endl; 
  CGAL::test_envelop_for_lines<SVCET_3>(
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
      Line_3(Point_3( 1, 1, 1), Point_3( 1, 1, 2)));

  // special configurations for one line 
  std::cout<< "TEST: intersecting line in direction of v2 " << std::endl; 
  CGAL::test_envelop_for_lines<SVCET_3>(
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
      Line_3(Point_3( 0, 0, 0), Point_3( 0, 1, 0)));

  std::cout<< "TEST: intersecting line in direction of v3 "  << std::endl; 
  CGAL::test_envelop_for_lines<SVCET_3>(
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
      Line_3(Point_3( 0, 0, 0), Point_3( 0, 0, 1)));
  
 
  // Special case ! 
  // base line and the following two lines are parallel to a common plane. 
  // according to Everett et al. this forms a nodal quartic with singular point at infinity.
  //  
  std::cout<< "TEST: Three lines two parallel planes: "  << std::endl; 
  std::cout<< "TEST: -- Boundaries on negative plane " << std::endl; 
  CGAL::test_envelop_for_lines<SVCET_3>(
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
      Line_3(Point_3( 1, 1, 1), Point_3( 2, 3, 1)),
      Line_3(Point_3( 2, 2, 2), Point_3( 2, 4, 2)));

  std::cout<< "TEST: -- Boundaries on different plane " << std::endl; 
  CGAL::test_envelop_for_lines<SVCET_3>(
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
      Line_3(Point_3( 1, 1, -1), Point_3( 2, 3, -1)),
      Line_3(Point_3( 2, 2, 2), Point_3( 2, 4, 2)));

  std::cout<< "TEST: -- Boundaries on positive plane " << std::endl; 
  CGAL::test_envelop_for_lines<SVCET_3>(
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
      Line_3(Point_3( 2, 2, -2), Point_3( 2, 4, -2)));
  CGAL::test_envelop_for_lines<SVCET_3>(
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
      Line_3(Point_3( 1, 1, -1), Point_3( 2, 3, -1)));
  CGAL::test_envelop_for_lines<SVCET_3>(
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
      Line_3(Point_3( 1, 1, -1), Point_3( 2, 3, -1)),
      Line_3(Point_3( 2, 2, -2), Point_3( 2, 4, -2)));


  std::cout<< " TEST: Parallel lines  " << std::endl; 
  // Special case parallel lines 
  std::cout<< " TEST: -- two parallel but not to base line " << std::endl; 
  // Bisector surface qith baseline are still hyperboloids of one sheet 
  CGAL::test_envelop_for_lines<SVCET_3>(
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
      Line_3(Point_3( 1, 2, 1), Point_3( 2, 1, 1)),
      Line_3(Point_3( 2, 3, 2), Point_3( 3, 2, 2)));
  
  std::cout<< " TEST: -- only one line that is parallel to base line " << std::endl; 
   CGAL::test_envelop_for_lines<SVCET_3>(
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
      Line_3(Point_3( 1, 1, 1), Point_3( 2, 1, 1)));

   std::cout<< " TEST: -- two lines one parallel line to base line " << std::endl; 
   CGAL::test_envelop_for_lines<SVCET_3>(
       Line_3(Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
       Line_3(Point_3( 1, 1, 1), Point_3( 2, 1, 1)),
       Line_3(Point_3( 2, 5, -1), Point_3( 2, 4, -2)));
   
  std::cout<< " TEST: -- two parallel lines to base line " << std::endl; 
   CGAL::test_envelop_for_lines<SVCET_3>(
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
      Line_3(Point_3( 1, 1, 1), Point_3( 2, 1, 1)),
      Line_3(Point_3( 2, 3, 2), Point_3( 3, 3, 2)));
   
  std::cout<< " TEST: -- two parallel lines to base line with same boundary " << std::endl; 
  CGAL::test_envelop_for_lines<SVCET_3>(
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
      Line_3(Point_3( 1, 1, 1), Point_3( 2, 1, 1)),
      Line_3(Point_3( 2, 2, 2), Point_3( 3, 2, 2)));

  std::cout<< " TEST: Intersecting line " << std::endl; 
  std::cout<< " TEST: one line intersecting the base line " << std::endl; 
  CGAL::test_envelop_for_lines<SVCET_3>(
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 1, 1)));
  

  std::cout<< " TEST: two lines both intersecting the base line " << std::endl; 
  CGAL::test_envelop_for_lines<SVCET_3>(
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
      Line_3(Point_3( 0, 0, 0), Point_3( 2, -1, 1)),
      Line_3(Point_3( 1, 0, 0), Point_3( 3, 2, -2)));

  std::cout<< " TEST: two lines both intersecting .. but strange case  " << std::endl; 
  CGAL::test_envelop_for_lines<SVCET_3>(
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
      Line_3(Point_3( 0, 0, 0), Point_3( 2, -1, 1)),
      Line_3(Point_3( 1, 0, 0), Point_3( 3, 2, -2)));

  std::cout<< " TEST: two parallel lines intersecting the base line" << std::endl; 
  CGAL::test_envelop_for_lines<SVCET_3>(
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
      Line_3(Point_3( 0, 0, 0), Point_3( 2, -1, 5)),
      Line_3(Point_3( 1, 0, 0), Point_3( 3, -1, 5)));

  std::cout<< " TEST: three lines intersecting in three points " << std::endl; 
  CGAL::test_envelop_for_lines<SVCET_3>(
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 1, 0)),
      Line_3(Point_3( 1, 1, 0), Point_3( 1, 0, 0)));

  std::cout<< " TEST: two lines both intersecting the base line in same point " << std::endl; 
  CGAL::test_envelop_for_lines<SVCET_3>(
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
      Line_3(Point_3( 0, 0, 0), Point_3( 2, -1, 1)),
      Line_3(Point_3( 0, 0, 0), Point_3( 3, 5, -2))); 
  
  std::cout<< " TEST: three lines in same plane intersect in one point" << std::endl; 
  CGAL::test_envelop_for_lines<SVCET_3>(
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
      Line_3(Point_3( 0, 0, 0), Point_3( 2,-1, 0)),
      Line_3(Point_3( 0, 0, 0), Point_3( 3, 5, 0)));

  std::cout<< " TEST: two lines one intersecting one parallel to base line" << std::endl; 
  CGAL::test_envelop_for_lines<SVCET_3>(
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
      Line_3(Point_3( 0, 0, 0), Point_3( 2,-1, 0)),
      Line_3(Point_3( 0, 1, 1), Point_3( 1, 1, 1)));
  
  std::cout<< " TEST: two lines one intersecting the base line - I" << std::endl;  
  CGAL::test_envelop_for_lines<SVCET_3>(
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
      Line_3(Point_3( 2, 5, -1), Point_3( 3, 2, -1)),
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 3, 2)));  
  std::cout<< " TEST: two lines one intersecting the base line - II" << std::endl; 
  CGAL::test_envelop_for_lines<SVCET_3>(
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 0, 0)),
      Line_3(Point_3( 2, 5, -1), Point_3( 3, 2, -1)),
      Line_3(Point_3( 0, 0, 0), Point_3( 1, 1, 1)));  


  std::cout<< "TEST: some generic instances: "  << std::endl; 
  // Generic position: 
  Line_3      line_0 (Point_3( 0, 0, 0), Point_3( 1, 0, 0)); // base line 
  std::vector<SVCET_3::Surface_3> lines;
  // two randome lines 
  lines.push_back(Line_3(Point_3(-5, 7, 2), Point_3( 5,-2, 1)));
  lines.push_back(Line_3(Point_3( 2,-8,-9), Point_3(-2,-1, 1)));
  lines.push_back(Line_3(Point_3( 4,-3, 7), Point_3( 8, 4,-8)));
  lines.push_back(Line_3(Point_3( 6, 8, 0), Point_3( 8,-8, 9)));
  lines.push_back(Line_3(Point_3( 3,-7,-9), Point_3( 6,-4,-2)));
  lines.push_back(Line_3(Point_3( 3, 7, 9), Point_3( 1, 1, 0)));


  std::cout << "################################  LINE 0 1 2" << std::endl;
  CGAL::test_envelop_for_lines<SVCET_3>(line_0,lines[0],lines[1],lines[2]);

  
  for (int i = 0; i<lines.size();i++){
    std::cout << "################################  LINE:  "<< i << std::endl;
    CGAL::test_envelop_for_lines<SVCET_3>(line_0,lines[i]);
  }
  for (int i = 0; i<lines.size();i++){
    for (int j = i+1; j<lines.size();j++){
      std::cout << "################################  LINES: "<< i << " " << j << std::endl;
    CGAL::test_envelop_for_lines<SVCET_3>(line_0,lines[i],lines[j]);
    }
  }
  for (int i = 0; i<lines.size();i++){
    for (int j = i+1; j<lines.size();j++){
      for (int k = j+1; k<lines.size();k++){
        std::cout << "################################  LINES: "
                  << i << " " << j  << " " << k << std::endl;
        CGAL::test_envelop_for_lines<SVCET_3>(line_0,lines[i],lines[j],lines[k]);
      }
    }   
  }
  
  // works but takes 343.796 in debug mode 
  CGAL::test_envelop_for_lines<SVCET_3>(line_0,lines);

  return 0; 
}


