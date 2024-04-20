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
  \file   triangle_area_distance_voronoi.cpp
  \brief  Example of the 2-point triangle area distance function.
*/


#include <iostream>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/envelope_voronoi_2.h>
#include <CGAL/Envelope_voronoi_2/Voronoi_diagram_2.h>
#include <CGAL/Envelope_voronoi_traits_2/Triangle_area_distance_traits_2.h>

#include "cartesian_product.h"

using namespace std;

typedef CGAL::Exact_predicates_exact_constructions_kernel        Kernel;
typedef CGAL::Triangle_area_distance_traits_2<Kernel>            VD_traits;

typedef CGAL::Envelope_voronoi_2::Voronoi_diagram_2<VD_traits>   Triangle_VD;


typedef Kernel::Point_2                                          Point_2;
typedef VD_traits::Site_2                                        Site_2;

int main( int argc, char **argv )
{

  int n;
  std::cin >> n;
  
  std::list<Point_2> sites;
  for (int i = 0; i < n; ++i)
  {
    Point_2 p;
    cin >> p;
    sites.push_back(p);
  }

  std::list<Site_2> site_pairs;
  cartesian_square_without_order (sites.begin(), sites.end(),
                                  std::back_inserter(site_pairs));

  std::cout << "Number of sites: " << site_pairs.size() << std::endl;
  
  Triangle_VD diagram;
  CGAL::voronoi_2 (site_pairs.begin(), site_pairs.end(),
                   diagram);

  std::cout << "Triangle area VD:" << std::endl <<
    "V = " << diagram.number_of_vertices() << ", E = " << 
    diagram.number_of_edges() << ", F = " << diagram.number_of_faces() << 
    std::endl;
  
  
  return 0;
}
