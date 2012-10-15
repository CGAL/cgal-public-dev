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
// Author(s): Ophir Setter          <ophirset@post.tau.ac.il>
//
//

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_quadratic_refinement_rep_bfi.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>
#include <CGAL/Algebraic_kernel_1.h>

#include <CGAL/Algebraic_kernel_d/Algebraic_curve_kernel_2.h> // template arg to CKvA_2
#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_impl.h> // traits for Arr_2

#include <CGAL/Envelope_voronoi_traits_2/Algebraic_curve_envelope_traits_2.h>
#include <CGAL/envelope_3.h>


int main(int argc, char** argv)
{

  typedef CGAL::Arithmetic_kernel AK;

  typedef AK::Integer Integer;
  typedef AK::Rational Rational;


// Definition of Algebraic_kernel_1
  typedef Integer Coefficient;
  typedef Rational Boundary; 
  typedef CGAL::CGALi::Algebraic_real_quadratic_refinement_rep_bfi
    < Coefficient , Boundary >                            Rep_class;
  typedef CGAL::CGALi::Bitstream_descartes
    < CGAL::Polynomial< Coefficient >, Rational >         Isolator;
  typedef CGAL::Algebraic_kernel_1
    <Coefficient, Rational, Rep_class, Isolator>          Algebraic_kernel_1;

// Definition of Algebraic_kernel_2 (Algebraic_curve_kernel_2)
  typedef CGAL::Algebraic_curve_kernel_2<Algebraic_kernel_1>
    Algebraic_curve_kernel_2;

  typedef CGAL::Curved_kernel_via_analysis_2< Algebraic_curve_kernel_2 > 
    Curved_kernel_2;


  typedef CGAL::Algebraic_curve_envelope_traits_2<
      CGAL::Adaptable_points_traits< Curved_kernel_2 > >  Env_traits_3;
  typedef CGAL::Envelope_diagram_2<Env_traits_3>          Envelope_diagram_2;

  
  Env_traits_3 tr;
  Envelope_diagram_2 diag(&tr); 
  
  std::list<Env_traits_3::Surface_3> points;
  // add the infinite plane
  points.push_back(std::make_pair(Rational(0), Rational(0)));
  points.push_back(std::make_pair(Rational(1), Rational(1)));

  CGAL::lower_envelope_3(points.begin(), points.end(), diag);


}


