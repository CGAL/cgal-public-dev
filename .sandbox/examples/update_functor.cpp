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

typedef CGAL::Collision_candidate<Scene_face_index> Collision_candidate;
typedef std::vector<Collision_candidate>            OutputIterator;

struct Translate_functor {

  Transform translation;
  Mesh_index mi{0};

  Translate_functor(Vector translation_vector) : translation(CGAL::TRANSLATION, translation_vector) {}

  void operator() (SMesh* mesh, typename Scene::Scene_vertex_index svi) {
    if( mi == svi.mesh_index()) 
    {
      swap(
        mesh->point(svi.local_index()).transform(this->translation),
        mesh->next_point(svi.local_index())
      );
    }
  }
};

struct Swap_current_next_functor {

    Transform translation;
    Mesh_index mi{0};

    Swap_current_next_functor(Mesh_index mi) : mi{mi} {}

    void operator() (SMesh* mesh, typename Scene::Scene_vertex_index svi) {
      if( mi == svi.mesh_index()) 
      {
        swap(
          mesh->point(svi.local_index()),
          mesh->next_point(svi.local_index())
        );
      }
    }
  };

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

  Mesh m1;
  Mesh m2;
  CGAL::make_tetrahedron(p1, q1, r1, s1, m1);
  CGAL::make_tetrahedron(p2, q2, r2, s2, m2);

  std::vector<SMesh> meshes;
  meshes.reserve(2);
  meshes.push_back(m1);
  meshes.push_back(m2);

  Scene scene(meshes);

  // Translate next_point
  Vector v(-1.9, -1.9, -1.9);
  Translate_functor t(v);
  scene.update_state(t);

  // Draw before
  std::cout << "Visualize current state..." std::endl;
  CGAL::draw_color(scene.joined_meshes());

  // Swap current and next points for visual
  Swap_current_next_functor scn(Mesh_index(0));
  scene.update_state(scn);

  // Draw after
  std::cout << "Visualize next state..." std::endl;
  CGAL::draw_color(scene.joined_meshes());

  // Return to current state.
  scene.update_state(scn);

  return EXIT_SUCCESS;
}
  