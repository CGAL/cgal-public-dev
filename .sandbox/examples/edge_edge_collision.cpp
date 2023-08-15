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
#include <Segment_3_Segment_3_do_collide.h>
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
  





