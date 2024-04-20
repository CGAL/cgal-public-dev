
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
  \file   power_diagram.cpp
  \brief  Construct the power diagram of circles in the plane. The traits class
  also support the standard Voronoi diagram by adding a circle with radius of
  zero. The circles squared radii should be given instead of the radii 
  themselves. Only the squared radii should be rational.  
*/


//! \file examples/Envelope_3/ex_envelope_planes.cpp
// Constructing the lower and the upper envelope of a set of planes.

#include <iostream>
#include <vector>

#include <boost/program_options.hpp>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/envelope_voronoi_2.h>
#include <CGAL/Envelope_voronoi_traits_2/Spherical_power_diagram_traits_2.h>

namespace po = boost::program_options;

typedef CGAL::Exact_predicates_exact_constructions_kernel    Kernel;

typedef Kernel::Point_3                                      Point_3;
typedef Kernel::Plane_3                                      Plane_3;


typedef CGAL::Spherical_power_diagram_traits_2<Kernel>       Traits_2;
typedef CGAL::Envelope_voronoi_2::
          Spherical_voronoi_diagram_2<Traits_2>              Voronoi_diagram_2;

int main (int argc, char **argv)
{ 
  std::vector<Traits_2::Site_2>   sites;

  sites.push_back(Point_3());
  sites.push_back(Plane_3());

  Voronoi_diagram_2       VD;
  voronoi_2 (sites.begin(), sites.end(), VD);

  return (0);
}



