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
  \file   anisotropic_diagram.cpp
  \brief  anisotropic diagram example.
*/

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Cartesian.h>

#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_impl.h> // traits for Arr_2
#include <CGAL/Algebraic_curve_kernel_2_generator.h>

#include <CGAL/envelope_3.h>
#include <CGAL/Envelope_voronoi_traits_2/Anistropic_voronoi_traits_2.h>


using namespace std;

typedef CGAL::Arithmetic_kernel                       AK;
typedef AK::Integer                                   Integer;
typedef AK::Rational                                  Rational;
  


// Definition of Algebraic_kernel_2 (Algebraic_curve_kernel_2)
typedef Integer                                       Coefficient;
typedef CGAL::Algebraic_curve_kernel_2_generator<Coefficient>::
Algebraic_curve_kernel_with_qir_and_bitstream_2       Algebraic_curve_kernel_2;
typedef CGAL::Curved_kernel_via_analysis_2< Algebraic_curve_kernel_2 > 
                                                      Curved_kernel_2;

typedef Curved_kernel_2::Point_2                      Point_2;
typedef Curved_kernel_2::Curve_2                      Curve_2;
typedef Curved_kernel_2::X_monotone_curve_2           X_monotone_curve_2;

typedef CGAL::Anisotropic_voronoi_traits_2< 
  Curved_kernel_2 >                                   VD_traits_2;

typedef CGAL::Envelope_diagram_2<VD_traits_2>         Anisotropic_diagram_2;

typedef VD_traits::Site_2                             Site_2;

int main( int argc, char **argv )
{

  std::list<Site_2>               sites; // anisotropic sites

  int n;
  std::cin >> n;
  for (int i = 0; i < n; ++i)
  {
    Raional x, y, a, b, c, d;
    std::cin >> x >> y >> a >> b >> c >> d;
    std::cout << x << " " << y << " " << a << " " << b << " " << c << " " << d
              << endl;

    // making the matrix possitive.
    CGAL_assertion_msg(a*c - b*b > 0, "input invalid");

    sites.push_back(Site_2(make_pair(x, y), a, b, c, d));
  }  

  Anisotropic_diagram_2 VD;
  CGAL::lower_envelope_3 (sites.begin(), sites.end(), VD);

  std::cout << "Number of sites: " << sites.size() << std::endl;
  
  std::cout << "Anisotropic VD:" << std::endl <<
    "V = " << VD.number_of_vertices() << ", E = " << 
    VD.number_of_edges() << ", F = " << VD.number_of_faces() << 
    std::endl;
  
  
  return 0;
}
