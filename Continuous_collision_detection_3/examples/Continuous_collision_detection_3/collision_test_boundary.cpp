// Copyright (c) 2023
// INRIA Sophia-Antipolis (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jeffrey Cochran

#include <variant>
#include <iostream>
#include <Bilinear_patch_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Continuous_collision_detection_3/internal/Segment_3_Segment_3_collision_test_boundary.h>
#include <CGAL/Continuous_collision_detection_3/internal/Point_3_Triangle_3_collision_test_boundary.h>

typedef CGAL::Simple_cartesian<double>  Kernel;
typedef Kernel::Point_3                 Point;
typedef Kernel::Vector_3                Vector;
typedef Kernel::Triangle_3              Triangle;
typedef Kernel::Segment_3               Segment;
typedef CGAL::BilinearPatchC3<Kernel>   Bilinear_patch;

using Collision_test_boundary_PT = CGAL::Collisions::internal::Point_3_Triangle_3_collision_test_boundary<Kernel>;
using Collision_test_boundary_SS = CGAL::Collisions::internal::Segment_3_Segment_3_collision_test_boundary<Kernel>;

int main(int argc, char* argv[])
{

  Point p1(3.0, 2.0, 2.0);
  Point q1(2.0, 3.0, 2.0);
  Point r1(2.0, 2.0, 3.0);
  Point s1(2.0, 2.0, 2.0);  

  Point p2(1.0, 0.0, 0.0);
  Point q2(0.0, 1.0, 0.0);
  Point r2(0.0, 0.0, 1.0);
  Point s2(0.0, 0.0, 0.0);

  //
  // Check that the 
  Collision_test_boundary_SS CTB_SS(
    Segment(p1, q1), 
    Segment(p2, q2),
    Segment(r1, s1), 
    Segment(r2, s2)
  );

  size_t num_bp{std::size(CTB_SS.bilinear_patches())};

  std::cout << "Check that edge-edge collision boundary is a cuboid..." << std::endl;
  std::cout << "...Number of Facets = " << num_bp << std::endl;
  std::cout << "...Number of Bilinear Patches = " << num_bp << std::endl;



  Collision_test_boundary_PT CTB_PT(
    p1, 
    p2,
    Triangle(q1, r1, s1), 
    Triangle(q2, r2, s2)
  );

  size_t num_tr{std::size(CTB_PT.triangles())};
  num_bp = std::size(CTB_SS.bilinear_patches());

  std::cout << "\nCheck that point-triangle collision boundary is a triangular prism..." << std::endl;
  std::cout << "...Number of Facets = " << num_tr + num_bp<< std::endl;
  std::cout << "...Number of Triangles = " << num_tr << std::endl;
  std::cout << "...Number of Bilinear Patches = " << num_bp << std::endl;

  return EXIT_SUCCESS;
}
  