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

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <Bilinear_patch_3.h>
#include <Trajectories.h>
#include <Segment_3_Segment_3_do_collide.h>
#include <Ray_3_Bilinear_patch_3_do_intersect.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Origin.h>
#include <cstdlib>
#include <ctime>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef Kernel::Point_3                   Point;
typedef Kernel::Vector_3                  Vector;
typedef Kernel::Ray_3                     Ray;
typedef ::CGAL::Bilinear_patch_3<Kernel>   BilinearPatch;
typedef ::CGAL::Point_3_trajectory<Kernel> P_trajectory;
typedef ::CGAL::Segment_3_trajectory<Kernel> S_trajectory;
typedef Kernel::FT                        FT;

int main(int argc, char* argv[])
{

Point s10_past(-1,  1, 0);
Point s10_next(-1, -1, 0);

Point s11_past( 1,  1, 0);
Point s11_next( 1, -1, 0);

Point s00_past(0, -1, -1);
Point s00_next(0,  1, -1);

Point s01_past(0, -1,  1);
Point s01_next(0,  1,  1);

P_trajectory s00(s00_past, s00_next);
P_trajectory s01(s01_past, s01_next);
P_trajectory s10(s10_past, s10_next);
P_trajectory s11(s11_past, s11_next);

S_trajectory s0(s00, s01);
S_trajectory s1(s10, s11);

std::cout << "Confirm these edges collide: " << CGAL::do_collide(s0, s1) << std::endl;

s10_past = Point(1,  1, -1);
s10_next = Point(1, -1, -1);

s11_past = Point(1,  1, 1);
s11_next = Point(1, -1, 1);

s00_past = Point(0, -1, -1);
s00_next = Point(0,  1, -1);

s01_past = Point(0, -1,  1);
s01_next = Point(0,  1,  1);
 
s00 = P_trajectory(s00_past, s00_next);
s01 = P_trajectory(s01_past, s01_next);
s10 = P_trajectory(s10_past, s10_next);
s11 = P_trajectory(s11_past, s11_next);

s0 = S_trajectory(s00, s01);
s1 = S_trajectory(s10, s11);

std::cout << "Confirm these edges do not collide: " << CGAL::do_collide(s0, s1) << std::endl;

s10_past = Point(1,  0,  1);
s10_next = Point(1,  0, -1);

s11_past = Point(-1, 1,  1);
s11_next = Point(-1, 1, -1);

s00_past = Point(0, 0, 0.01);
s00_next = Point(0, 0, 0.00);

s01_past = Point(0.01, 0, 1.01);
s01_next = Point(0.01, 0, 1.00);
 
s00 = P_trajectory(s00_past, s00_next);
s01 = P_trajectory(s01_past, s01_next);
s10 = P_trajectory(s10_past, s10_next);
s11 = P_trajectory(s11_past, s11_next);

s0 = S_trajectory(s00, s01);
s1 = S_trajectory(s10, s11);

std::cout << "Confirm these edges do not collide: " << CGAL::do_collide(s0, s1) << std::endl;

return EXIT_SUCCESS;
}
  





