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
typedef Kernel::FT                        FT;
typedef Kernel::Triangle_3                Triangle;

int main(int argc, char* argv[])
{

std::srand(std::time(nullptr));
FT rand_max{RAND_MAX};
FT u = FT(std::rand())/rand_max;
FT v = FT(std::rand())/rand_max;

// =============================
// Stuff for me to revisit later
// =============================

Point a(0, 0, 0);
Point b(1, 0, 1);
Point c(1, 1, 0);
Point d(0, 1, 1);

BilinearPatch bp = BilinearPatch(a, b, c, d);

Vector bump(0.,0.,.01);
Vector huge_bump(0.,0.1,1000.);

Point surface_point   = bp(0.5, 0.5);
Point surface_point2  = bp(u, v);
Point point_below     = ::CGAL::ORIGIN + ((surface_point-::CGAL::ORIGIN) + bump); // meaning phi > 0, not above the surface
Point point_above     = ::CGAL::ORIGIN + ((surface_point-::CGAL::ORIGIN) - bump);
Point point_way_above = ::CGAL::ORIGIN + ((surface_point-::CGAL::ORIGIN) - huge_bump);

std::cout << "Value of phi(surface_point) [should be ~zero]: "    << bp.signed_scaled_patch_distance(surface_point) << "\n";
std::cout << "Value of phi(surface_point2) [should be ~zero]: "   << bp.signed_scaled_patch_distance(surface_point2) << "\n";
std::cout << "Value of phi(point_below) [should be < zero]: "     << bp.signed_scaled_patch_distance(point_below) << "\n";
std::cout << "Value of phi(point_above) [should be > zero]: "     << bp.signed_scaled_patch_distance(point_above) << "\n";
std::cout << "Value of phi(point_way_above) [should be >> zero]: "<< bp.signed_scaled_patch_distance(point_way_above) << "\n\n";


// Case 1a:
// Ray r(surface_point, point_above);
// bool does_intersect = ::CGAL::Intersections::internal::do_intersect_odd_parity(bp, r);
// std::cout << "Ray originating on the surface does intersect: " << does_intersect << "\n";

// Case 1b:
Ray r = Ray(point_above, point_below);
bool does_intersect = ::CGAL::Intersections::internal::do_intersect_odd_parity(bp, r);
std::cout << "Ray originating in the bounding tetrahedron (above surface) and passing through does intersect: " << does_intersect << "\n";

// Case 1c:
r = Ray(point_below, point_above);
does_intersect = ::CGAL::Intersections::internal::do_intersect_odd_parity(bp, r);
std::cout << "Ray originating in the bounding tetrahedron (below surface) and passing through does intersect: " << does_intersect << "\n";

// Case 2:
r = Ray(point_way_above, point_below);
does_intersect = ::CGAL::Intersections::internal::do_intersect_odd_parity(bp, r);
std::cout << "Ray originating outside the bounding tetrahedron and passing through does intersect: " << does_intersect << "\n";

//
// Now we look at some edge cases...
//
// Planar patch, but not a parallelogram
a = Point(-1, -1, 0);
b = Point(-1,  1, 0);
c = Point( 1, -1, 0);
d = Point( 1,  1, 0);
bp = BilinearPatch(a, b, c, d);
point_above = Point(0, 1, 0);
surface_point = Point(-.5, 0, 0);

std::cout << "Testing planar case I, decomposable to triangles...\n";
std::cout << "...confirm planar: " << bp.is_planar() << "\n";
std::cout << "...confirm non-degenerate: " << !bp.is_degenerate() << "\n";
std::cout << "...point_above is not on the patch: " << bp.has_on(point_above) << "\n";
std::cout << "...surface_point is on the patch: " << bp.has_on(surface_point) << "\n";
//
// Planar patch, but parallelogram
a = Point(-1, -1, 0);
b = Point(-1,  1, 0);
c = Point( 1, -1, 0);
d = Point( 1,  1, 0);
bp = BilinearPatch(a, b, d, c);
point_above = Point(0, 1, 0);
surface_point = Point(-.5, 0, 0);

std::cout << "Testing planar case II, decomposable to triangles...\n";
std::cout << "...confirm planar: " << bp.is_planar() << "\n";
std::cout << "...confirm non-degenerate: " << !bp.is_degenerate() << "\n";
std::cout << "...point_above is now on the patch: " << bp.has_on(point_above) << "\n";
std::cout << "...surface_point is still on the patch: " << bp.has_on(surface_point) << "\n";
//
// Planar patch, but degenerate
a = Point(-2, 0, 0);
b = Point(-1, 0, 0);
c = Point( 1, 0, 0);
d = Point( 2, 0, 0);
bp = BilinearPatch(a, b, d, c);
point_above = Point(3, 0, 0);
surface_point = Point(-.5, 0, 0);

std::cout << "Testing planar case III, degenerate and not decomposable to triangles...\n";
std::cout << "...confirm planar: " << bp.is_planar() << "\n";
std::cout << "...confirm degenerate: " << bp.is_degenerate() << "\n";
std::cout << "...point_above is not on the patch: " << bp.has_on(point_above) << "\n";
std::cout << "...surface_point is still  on the patch: " << bp.has_on(surface_point) << "\n";
  
return EXIT_SUCCESS;
}
  





