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

std::cout << "Value of phi(surface_point) [should be ~zero]: "    << CGAL::to_double(bp.aux_phi(surface_point)) << "\n";
std::cout << "Value of phi(surface_point2) [should be ~zero]: "   << CGAL::to_double(bp.aux_phi(surface_point2)) << "\n";
std::cout << "Value of phi(point_below) [should be < zero]: "     << CGAL::to_double(bp.aux_phi(point_below)) << "\n";
std::cout << "Value of phi(point_above) [should be > zero]: "     << CGAL::to_double(bp.aux_phi(point_above)) << "\n";
std::cout << "Value of phi(point_way_above) [should be >> zero]: "<< CGAL::to_double(bp.aux_phi(point_way_above)) << "\n\n";


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

// Tetrahedron t = bp.tetrahedron();

// Point outside_point = Point(1, 1, 1);
// Point inside_point = Point(0.25, 0.25, 0.25);
// std::cout << "Orientation is positive: " << t.orientation() << std::endl;
// std::cout << "Inside point is inside: " << t.has_on_bounded_side(inside_point) << std::endl;
// std::cout << "Outside point is outside: " << t.has_on_unbounded_side(outside_point) << std::endl;
// std::cout << "Should be zero: " << bp.phi(d) << std::endl;
// std::cout << "Should be false: " << bp.is_planar() << std::endl;
// Point a = Point(0, 0, 0);
// Point b = Point(1, 0, 0);
// Point c = Point(0, 1, 0);
// Point d = Point(0, 0, 1);
// Point outside_point = Point(1, 1, 1);
// Point inside_point = Point(0.25, 0.25, 0.25);

// BilinearPatch bp = BilinearPatch(a, b, c, d);
// Tetrahedron t = bp.tetrahedron();
// Ray r = Ray(a, outside_point);

// std::cout << "Orientation is positive: " << t.orientation() << std::endl;
// std::cout << "Inside point is inside: " << t.has_on_bounded_side(inside_point) << std::endl;
// std::cout << "Outside point is outside: " << t.has_on_unbounded_side(outside_point) << std::endl;
// std::cout << "Should be true: " << bp.aux_phi(d).do_overlap(0) << std::endl;
// std::cout << "Should be true: " << bp.has_on(d) << std::endl;
// std::cout << "Should be false: " << bp.has_on(outside_point) << std::endl;
// std::cout << "Should be false: " << bp.is_planar() << std::endl;

// auto vcm = a.add_property_map<Mesh::Vertex_index, CGAL::IO::Color>("v:color").first;
// auto ecm = a.add_property_map<Mesh::Edge_index, CGAL::IO::Color>("e:color").first;
// auto fcm = a.add_property_map<Mesh::Face_index>("f:color", CGAL::IO::white() /*default*/).first; 

// for(auto v : vertices(a))
// {
//   if(v.idx()%2) 
//   {
//     Point& p = a.point(v);
//     Vector& q = a.velocity(v);

//     Transform transformation(CGAL::TRANSLATION, q);

//     // swap(q,p.transform(transformation));

//     std::cout << "P.x current: " << a.point(v).x() << std::endl;
//     std::cout << "V.x current: " << a.velocity(v).x() << std::endl;
//   }
// }

// for(auto f : a.faces())
// {
//   // or the same again, but directly with a range based loop
//   for(SMesh::Vertex_index vi : a.vertices_around_face(a.halfedge(f)))
//   {
//     std::cout << vi << std::endl;
//   }
// }

// CGAL_USE(fcm);



// CGAL::draw_color(a);

// Case 0: planar patch
// Point a = Point(0, 0, 0);
// Point b = Point(1, 0, 0);
// Point c = Point(0, 1, 0);
// Point d = Point(1, 1, 0);
// Point outside_point = Point(0.25, 0.25, 1);
// Point outside_point_2 = Point(-0.25, -0.25, 1);
// Point inside_point = Point(0.25, 0.25, 0);

// BilinearPatch bp = BilinearPatch(a, b, c, d); 

// Ray ray_true = Ray(outside_point, inside_point);
// Ray ray_false = Ray(outside_point, outside_point_2);

// bool intersects = CGAL::Intersections::internal::do_intersect_odd_parity(bp, ray_true);
// std::cout << "Case 0, True: " << intersects << std::endl;

// intersects = CGAL::Intersections::internal::do_intersect_odd_parity(bp, ray_false);
// std::cout << "Case 0, False: " << intersects << std::endl;

// // Case 1a: non-planar, ray origin on the bilinear patch
// a = Point(0, 0, 0);
// b = Point(1, 0, 0);
// c = Point(0, 1, 0);
// d = Point(0, 0, 1);
// outside_point = Point(-0.1, -0.1, -0.1);
// outside_point_2 = Point(1, 1, 1);
// inside_point = Point(0.333, 0.333, 0.333);

// bp = BilinearPatch(a, b, c, d); 

// ray_true = Ray(a, outside_point);

// intersects = CGAL::Intersections::internal::do_intersect_odd_parity(bp, ray_true);
// std::cout << "Case 1a, True: " << intersects << std::endl;

// // Case 1b: non-planar, ray origin inside bounding tetrahedron, but not on patch
// ray_true = Ray(inside_point, outside_point);
// ray_false = Ray(inside_point, outside_point_2);

// intersects = CGAL::Intersections::internal::do_intersect_odd_parity(bp, ray_true);
// std::cout << "Case 1b, True: " << intersects << std::endl;

// intersects = CGAL::Intersections::internal::do_intersect_odd_parity(bp, ray_false);
// std::cout << "Case 1b, False: " << intersects << std::endl;

  return EXIT_SUCCESS;
}
  





