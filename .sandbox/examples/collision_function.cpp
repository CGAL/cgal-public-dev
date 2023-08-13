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

#include <Point_3_Triangle_3_collision_function.h>
#include <Segment_3_Segment_3_collision_function.h>

typedef CGAL::Simple_cartesian<double>  Kernel;
typedef Kernel::Point_3                 Point;
typedef Kernel::Vector_3                Vector;
typedef Kernel::Triangle_3              Triangle;
typedef Kernel::Segment_3               Segment;

using Collision_function_PT = typename CGAL::Collisions::internal::Point_3_Triangle_3_collision_function<Kernel>;
using Collision_function_SS = typename CGAL::Collisions::internal::Segment_3_Segment_3_collision_function<Kernel>;

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

  Collision_function_SS CF_EE(
    Segment(p1, q1), 
    Segment(p2, q2),
    Segment(r1, s1), 
    Segment(r2, s2)
  );

  Collision_function_PT CF_PT(
    p1, 
    p2, 
    Triangle(q1, r1, s1), 
    Triangle(q2, r2, s2)
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
  