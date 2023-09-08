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

#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/CCD_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef Kernel::Point_3                 Point;
typedef Kernel::Vector_3                Vector;
typedef Kernel::Tetrahedron_3           Tetrahedron;
typedef Kernel::Ray_3                   Ray;
typedef Kernel::Aff_transformation_3    Transform;
typedef CGAL::Surface_mesh<Point>       Mesh;
typedef CGAL::Collision_mesh<Kernel>    SMesh;
typedef CGAL::Collision_scene<Kernel>   Scene;
typedef CGAL::Point_3_trajectory<Kernel>Point_trajectory;
typedef Scene::Vertex_index             Vertex_index;
typedef Scene::Halfedge_index           Halfedge_index;
typedef Scene::Mesh_index               Mesh_index;
typedef Scene::Face_index               Face_index;
typedef CGAL::Bilinear_patch_3<Kernel>   BilinearPatch;
typedef Scene::Primitive_id             Primitive_id;
typedef Scene::Scene_face_index         Scene_face_index;

typedef CGAL::Collisions::internal::Swap_current_next_functor<Scene>  Swap_functor;
typedef CGAL::Collisions::internal::Translate_functor<Scene>          Translate_functor;
typedef CGAL::Collisions::internal::Contraction_functor<Scene>        Contraction_functor;


typedef CGAL::Collision_candidate<Scene::Trajectory> Collision_candidate;
typedef std::vector<Collision_candidate>            OutputIterator;

int main(int argc, char* argv[])
{

  const std::string outer_sphere_filename = CGAL::data_file_path("meshes/Sphere_External.off");

  Point p2(1, 0, 0);
  Point q2(0, 1, 0);
  Point r2(0, 0, 1);
  Point s2(0, 0, 0);

  Mesh m1;
  CGAL::make_tetrahedron(p2, q2, r2, s2, m1);

  Mesh outer_sphere_mesh;

  if(!CGAL::IO::read_polygon_mesh(outer_sphere_filename, outer_sphere_mesh))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  std::vector<SMesh> meshes;
  meshes.reserve(2);
  meshes.push_back(SMesh(m1));
  meshes.push_back(SMesh(outer_sphere_mesh));

  Scene scene(meshes);

  Swap_functor swap_functor = Swap_functor();

  // Translate next_point, mesh
  Vector v(  -3,      0, 0);
  Point  p(   3.333,  0, 0);

  Translate_functor set_mesh  = Translate_functor(Mesh_index(0), v);
  Contraction_functor collide_meshes = Contraction_functor(Mesh_index(0), p, 0.9);

  scene.update_state(set_mesh, true);
  scene.update_state(swap_functor);
  scene.update_state(collide_meshes, true);

  std::cout << "Scene before collision..." << std::endl;
  ::CGAL::draw_collision_scene(scene);

  std::cout << "Scene after collision\n(press 'm' to color the colliding faces white)..." << std::endl;
  ::CGAL::draw_collision_scene(scene, true);

  collide_meshes = Contraction_functor(Mesh_index(0), p, 0.2);
  scene.update_state(collide_meshes, true);

  std::cout << "Scene during collision, t=0.2\n(press 'm' to color the colliding faces white)..." << std::endl;
  ::CGAL::draw_collision_scene(scene, true);

  collide_meshes = Contraction_functor(Mesh_index(0), p, 0.25);
  scene.update_state(collide_meshes, true);

  std::cout << "Scene during collision, t=0.25\n(press 'm' to color the colliding faces white)..." << std::endl;
  ::CGAL::draw_collision_scene(scene, true);

  collide_meshes = Contraction_functor(Mesh_index(0), p, 0.3);
  scene.update_state(collide_meshes, true);

  std::cout << "Scene during collision, t=0.3\n(press 'm' to color the colliding faces white)..." << std::endl;
  ::CGAL::draw_collision_scene(scene, true);

  collide_meshes = Contraction_functor(Mesh_index(0), p, 0.33);
  scene.update_state(collide_meshes, true);

  std::cout << "Scene during collision, t=0.33\n(press 'm' to color the colliding faces white)..." << std::endl;
  ::CGAL::draw_collision_scene(scene, true);

  collide_meshes = Contraction_functor(Mesh_index(0), p, 0.6);
  scene.update_state(collide_meshes, true);

  std::cout << "Scene during collision, t=0.6\n(press 'm' to color the colliding faces white)..." << std::endl;
  ::CGAL::draw_collision_scene(scene, true);

  collide_meshes = Contraction_functor(Mesh_index(0), p, 0.64);
  scene.update_state(collide_meshes, true);

  std::cout << "Scene during collision, t=0.64\n(press 'm' to color the colliding faces white)..." << std::endl;
  ::CGAL::draw_collision_scene(scene, true);

  collide_meshes = Contraction_functor(Mesh_index(0), p, 0.66);
  scene.update_state(collide_meshes, true);

  std::cout << "Scene during collision, t=0.66\n(press 'm' to color the colliding faces white)..." << std::endl;
  ::CGAL::draw_collision_scene(scene, true);

  return EXIT_SUCCESS;
}
