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
#include <CGAL/Envelope_voronoi_traits_2/Power_diagram_traits_2.h>

namespace po = boost::program_options;

typedef CGAL::Exact_predicates_exact_constructions_kernel    Kernel;

typedef Kernel::Circle_2                                     Circle_2;


typedef CGAL::Power_diagram_traits_2<Kernel>                 Traits_2;
typedef CGAL::Envelope_voronoi_2::
          Voronoi_diagram_2<Traits_2>                        Voronoi_diagram_2;

int main (int argc, char **argv)
{
  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("farthest", "Farthest-site voronoi diagram")
    ;
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);    
  
  if (vm.count("help")) 
  {
    std::cout << desc << "\n";
    return 1;
  }
  
  std::vector<Circle_2>     sites;

  int n_circles;
  std::cin >> n_circles;
  for (int i = 0; i < n_circles; ++i)
  {
    Kernel::Point_2 c;
    Kernel::FT r;
    std::cin >> c >> r;
    sites.push_back(Circle_2(c, r));
//    std::cout << c << std::endl;
  }
  
  // Compute and print the minimization diagram.
  Voronoi_diagram_2       VD;
  bool furthest = vm.count("furthest");

  if (furthest)
    farthest_voronoi_2 (sites.begin(), sites.end(), VD);
  else
    voronoi_2 (sites.begin(), sites.end(), VD);

  std::cout << "Number of sites " << n_circles << std::endl;
  std::cout << "V = " << VD.number_of_vertices()
            << ", E = " << VD.number_of_edges()
            << ", F = " << VD.number_of_faces() << std::endl;

  return (0);
}
