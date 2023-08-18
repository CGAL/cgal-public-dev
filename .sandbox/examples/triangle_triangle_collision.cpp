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
#include <Bilinear_patch_3.h>
#include <Trajectories.h>
#include <Triangle_3_Triangle_3_do_collide.h>
#include <Ray_3_Bilinear_patch_3_do_intersect.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Origin.h>
#include <cstdlib>
#include <ctime>
#include <iostream>

typedef ::CGAL::Simple_cartesian<double>  Kernel;
typedef Kernel::Point_3                   Point;
typedef Kernel::Vector_3                  Vector;
typedef Kernel::Ray_3                     Ray;
typedef ::CGAL::BilinearPatchC3<Kernel>   BilinearPatch;
typedef ::CGAL::Point_3_trajectory<Kernel> P_trajectory;
typedef ::CGAL::Triangle_3_trajectory<Kernel> T_trajectory;
typedef Kernel::FT                        FT;

int main(int argc, char* argv[])
{

Point t00_past(-1, -1,  1);
Point t00_next(-1, -1, -1);
Point t01_past( 0,  1,  1);
Point t01_next( 0,  1, -1);
Point t02_past( 1, -1,  1);
Point t02_next( 1, -1, -1);

Point t10_past(0,  0,  0.01);
Point t10_next(0,  0,  0);
Point t11_past(0,  0.01,  1.01);
Point t11_next(0,  0.01,  1);
Point t12_past(0, -0.01,  1.01);
Point t12_next(0, -0.01,  1);

P_trajectory t00(t00_past, t00_next);
P_trajectory t01(t01_past, t01_next);
P_trajectory t02(t02_past, t02_next);

P_trajectory t10(t10_past, t10_next);
P_trajectory t11(t11_past, t11_next);
P_trajectory t12(t12_past, t12_next);

T_trajectory t0(t00, t01, t02);
T_trajectory t1(t10, t11, t12);

std::cout << "Confirm triangle and triangle do collide: " << CGAL::do_collide(t0, t1) << std::endl;


return EXIT_SUCCESS;
}
  





