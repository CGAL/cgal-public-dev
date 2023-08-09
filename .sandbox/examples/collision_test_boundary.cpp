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
#include <Collision_test_boundary.h>

typedef CGAL::Simple_cartesian<double>  Kernel;
typedef Kernel::Point_3                 Point;
typedef Kernel::Vector_3                Vector;
typedef Kernel::Triangle_3              Triangle;
typedef ::CGAL::BilinearPatchC3<Kernel> Bilinear_patch;

typedef CGAL::Collisions::internal::COLLISION_TYPE COLLISION_TYPE;

template<COLLISION_TYPE C>
using Collision_test_boundary = CGAL::Collisions::internal::Collision_test_boundary<Kernel, C>;

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
  Collision_test_boundary<COLLISION_TYPE::EDGE_EDGE> CTB_EE(
    p1, q1, r1, s1, 
    p2, q2, r2, s2
  );

  std::cout << "Check that edge-edge collision boundary is a cuboid..." << std::endl;
  std::cout << "...Number of Facets = " << std::size(CTB_EE.facets) << std::endl;

  size_t num_bp{0};
  size_t num_tr{0};
  for(const auto& facet_ : CTB_EE.facets) {
    if( std::holds_alternative<Triangle>(facet_)) { ++num_tr; }
    if( std::holds_alternative<Bilinear_patch>(facet_)) { ++num_bp; }
  }

  std::cout << "...Number of Triangles = " << num_tr << std::endl;
  std::cout << "...Number of Bilinear Patches = " << num_bp << std::endl;



  Collision_test_boundary<COLLISION_TYPE::POINT_TRIANGLE> CTB_PT(
    p1, q1, r1, s1, 
    p2, q2, r2, s2
  );

  std::cout << "\nCheck that point-triangle collision boundary is a triangular prism..." << std::endl;
  std::cout << "...Number of Facets = " << std::size(CTB_PT.facets) << std::endl;

  num_bp = 0;
  num_tr = 0;
  for(const auto& facet_ : CTB_PT.facets) {
    if( std::holds_alternative<Triangle>(facet_)) { ++num_tr; }
    if( std::holds_alternative<Bilinear_patch>(facet_)) { ++num_bp; }
  }

  std::cout << "...Number of Triangles = " << num_tr << std::endl;
  std::cout << "...Number of Bilinear Patches = " << num_bp << std::endl;

  // Collision_function<COLLISION_TYPE::POINT_TRIANGLE> CF_PT(
  //   p1, q1, r1, s1, 
  //   p2, q2, r2, s2
  // );

  // double t{0.};
  // double u{0.};
  // double v{0.};

  // Point P_EE = CF_EE(t, u, v); 
  // Point P_PT = CF_PT(t, u, v); 

  // std::cout << "P1: " << p1 << std::endl;
  // std::cout << "Q1: " << q1 << std::endl;
  // std::cout << "R1: " << r1 << std::endl;
  
  // std::cout << "Should be P1 - R1: " << P_EE << std::endl;
  // std::cout << "Should be P1 - Q1: " << P_PT << std::endl;

  return EXIT_SUCCESS;
}
  