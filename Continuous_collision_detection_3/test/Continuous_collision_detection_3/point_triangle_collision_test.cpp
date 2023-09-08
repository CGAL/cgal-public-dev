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
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cassert>

#include <CGAL/CCD_3.h>

typedef ::CGAL::Simple_cartesian<double>  Kernel;
typedef Kernel::Point_3                   Point;
typedef Kernel::Vector_3                  Vector;
typedef Kernel::Ray_3                     Ray;
typedef ::CGAL::Bilinear_patch_3<Kernel>   BilinearPatch;
typedef ::CGAL::Point_3_trajectory<Kernel> P_trajectory;
typedef ::CGAL::Triangle_3_trajectory<Kernel> T_trajectory;
typedef Kernel::FT                        FT;

int main(int argc, char* argv[])
{

Point t0_past(-1, -1,  1);
Point t0_next(-1, -1, -1);
Point t1_past( 0,  1,  1);
Point t1_next( 0,  1, -1);
Point t2_past( 1, -1,  1);
Point t2_next( 1, -1, -1);

Point p_past(0, 0, -1);
Point p_next(0, 0,  1);

P_trajectory t0(t0_past, t0_next);
P_trajectory t1(t1_past, t1_next);
P_trajectory t2(t2_past, t2_next);
P_trajectory p(p_past,  p_next);

T_trajectory t(t0, t1, t2);

assert(CGAL::do_collide(p, t));


t0_past = Point(-1, -1,  -2);
t0_next = Point(-1, -1, -1);
t1_past = Point( 0,  1, -2);
t1_next = Point( 0,  1, -1);
t2_past = Point( 1, -1, -2);
t2_next = Point( 1, -1, -1);

t0 = P_trajectory(t0_past, t0_next);
t1 = P_trajectory(t1_past, t1_next);
t2 = P_trajectory(t2_past, t2_next);

t = T_trajectory(t0, t1, t2);

assert(CGAL::do_collide(p, t));


return EXIT_SUCCESS;
}
  





