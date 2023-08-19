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

#include <iostream>
#include <fstream>
#include <utility>
#include <iterator>
#include <vector>
#include <conio.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/intersections.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/draw_surface_mesh.h>

#include <Ray_3_Bilinear_patch_3_do_intersect.h>
#include <AABB_triangle_trajectory_primitive.h>
#include <Collision_mesh_3.h>
#include <Collision_scene_3.h>
#include <Bilinear_patch_3.h>
#include <Collision_candidate_3.h>
#include <Collisions_3.h>

typedef CGAL::Simple_cartesian<double>  Kernel;
typedef Kernel::Point_3                 Point;
typedef Kernel::Vector_3                Vector;
typedef Kernel::Tetrahedron_3           Tetrahedron;
typedef Kernel::Ray_3                   Ray;
typedef Kernel::Aff_transformation_3    Transform;                        
typedef CGAL::Surface_mesh<Point>       Mesh;                    
typedef CGAL::Collision_mesh<Kernel>    SMesh;
typedef CGAL::Collision_scene<Kernel>   Scene;
typedef Scene::Mesh_index               Mesh_index;
typedef Scene::Face_index               Face_index;
typedef CGAL::BilinearPatchC3<Kernel>   BilinearPatch;
typedef Scene::Primitive_id             Primitive_id;
typedef Scene::Scene_face_index         Scene_face_index;
typedef Scene::Trajectory               Trajectory;

typedef CGAL::Collision_candidate<Trajectory> Collision_candidate;
typedef std::vector<Collision_candidate>      OutputIterator;

int main(int argc, char* argv[])
{

  // const std::string inner_sphere_filename = CGAL::data_file_path("meshes/cactus.off");
  // const std::string outer_sphere_filename = CGAL::data_file_path("meshes/camel.off");

  // Mesh inner_sphere_mesh;
  // Mesh outer_sphere_mesh;

  // if(!CGAL::IO::read_polygon_mesh(inner_sphere_filename, inner_sphere_mesh))
  // {
  //   std::cerr << "Invalid input file for internal sphere." << std::endl;
  //   return EXIT_FAILURE;
  // }

  // if(!CGAL::IO::read_polygon_mesh(outer_sphere_filename, outer_sphere_mesh))
  // {
  //   std::cerr << "Invalid input file." << std::endl;
  //   return EXIT_FAILURE;
  // }

  
    Point p1(1.0, 0.0, 0.0);
    Point q1(0.0, 1.0, 0.0);
    Point r1(0.0, 0.0, 1.0);
    Point s1(0.0, 0.0, 0.0);


    double t{.3};
    Point p2(1.0 + t, t,       t      );
    Point q2(t,       1.0 + t, t      );
    Point r2(t,       t,       1.0 + t);
    Point s2(t,       t,       t      );

    Mesh m1, m2;
    CGAL::make_tetrahedron(p1, q1, r1, s1, m1);
    CGAL::make_tetrahedron(p2, q2, r2, s2, m2);

  std::vector<SMesh> meshes;
  meshes.reserve(2);
  meshes.push_back(m1);
  meshes.push_back(m2);

  Scene scene(meshes);

  OutputIterator collision_candidates = CGAL::get_collision_candidates<Kernel>(scene);
  std::set<Scene_face_index> filtered_indices;

  Scene_face_index ti(Mesh_index(1), Face_index(0));


  std::cout << "\nCollision candidates containing Face 0 of Mesh 0:\n===================" << std::endl;
  std::for_each(
    collision_candidates.begin(), 
    collision_candidates.end(), 
    [&ti, &filtered_indices](const auto& candidate){
      if(candidate.first->index == ti || candidate.second->index == ti)
      {
        std::cout << candidate << std::endl;
        filtered_indices.insert(candidate.first->index);
        filtered_indices.insert(candidate.second->index);
      };
    }
  );
  std::cout << "===================" << std::endl;
  
  std::cout << "\nPress <m> to color the faces:" << std::endl;
  std::cout << "  (blue)  : Face 0 of Mesh 0" << std::endl;
  std::cout << "  (red)   : Candidate for collision with Face 0 of Mesh 0" << std::endl;
  std::cout << "  (white) : Does not collide with Face 0 of Mesh 0\n" << std::endl;


  int k{0};
  for( const auto& fi : filtered_indices )
  {
    scene.color(fi, CGAL::IO::red());
  }

  scene.color(ti, CGAL::IO::blue());

  CGAL::draw_color(scene.joined_meshes());

  return EXIT_SUCCESS;
}
  





// =============================
// Stuff for me to revisit later
// =============================

// Point a = Point(0, 0, 0);
// Point b = Point(1, 0, 0);
// Point c = Point(0, 1, 0);
// Point d = Point(0, 0, 1);

// BilinearPatch bp = BilinearPatch(a, b, c, d);
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

