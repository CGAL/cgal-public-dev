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

#include <CGAL/Simple_cartesian.h>

#include <Collision_function.h>
#include <Collision_type.h>

typedef CGAL::Simple_cartesian<double>  Kernel;
typedef Kernel::Point_3                 Point;
typedef Kernel::Vector_3                Vector;

typedef CGAL::Collisions::internal::COLLISION_TYPE COLLISION_TYPE;

template<COLLISION_TYPE C>
using Collision_function = CGAL::Collisions::internal::Collision_function<Kernel, C>;

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

  Collision_function<COLLISION_TYPE::EDGE_EDGE> CF_EE(
    p1, q1, r1, s1, 
    p2, q2, r2, s2
  );

  Collision_function<COLLISION_TYPE::POINT_TRIANGLE> CF_PT(
    p1, q1, r1, s1, 
    p2, q2, r2, s2
  );

  double t{0.};
  double u{0.};
  double v{0.};

  Point P_EE = CF_EE(t, u, v); 
  Point P_PT = CF_PT(t, u, v); 

  std::cout << "P1: " << p1 << std::endl;
  std::cout << "Q1: " << q1 << std::endl;
  std::cout << "R1: " << r1 << std::endl;
  
  std::cout << "Should be P1 - R1: " << P_EE << std::endl;
  std::cout << "Should be P1 - Q1: " << P_PT << std::endl;

  return EXIT_SUCCESS;
}
  