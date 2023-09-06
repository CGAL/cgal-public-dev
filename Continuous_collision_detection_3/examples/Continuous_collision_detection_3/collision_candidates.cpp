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
// #include <conio.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/intersections.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/draw_surface_mesh.h>

#include <CGAL/CCD_3.h>

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

  CGAL::draw(scene.joined_meshes());

  return EXIT_SUCCESS;
}

