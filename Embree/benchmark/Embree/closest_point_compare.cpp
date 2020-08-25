#include <iostream>
#include <fstream>

#include <CGAL/Embree/AABB_tree.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/Real_timer.h>

#include "../../../AABB_tree/benchmark/AABB_tree/include/RaysGenerate.h"

typedef CGAL::Simple_cartesian<float> K;
// typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_3 Point;
typedef K::Ray_3 Ray;
typedef CGAL::Surface_mesh<Point> Mesh;

typedef CGAL::Embree::Triangle_mesh_geometry<Mesh, K> Geometry;
typedef CGAL::Embree::AABB_tree<Geometry, K> TreeEmbree;

typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> TreeCgal;

int main(int argc, char const *argv[])
{
  const char* filename =  (argc > 1) ? argv[1] :  "../../test/Embree/data/bunny00.off" ;
  std::ifstream input(filename);

  Mesh triangle_mesh;
  input >> triangle_mesh;

  int num_points = 1000; 
  RaysGenerate rg(num_points);

  CGAL::Real_timer time;

  time.start();
  TreeEmbree tree1;
  tree1.insert(triangle_mesh);
  time.stop();
  std::cout<<"Embree AABB construction time : "<< time.time()<<std::endl;

  Point centre((tree1.bbox().min(0)+tree1.bbox().max(0))/2, 
               (tree1.bbox().min(1)+tree1.bbox().max(1))/2,
               (tree1.bbox().min(2)+tree1.bbox().max(2))/2); 

  time.reset();

  time.start();
  TreeCgal tree2(faces(triangle_mesh).first, faces(triangle_mesh).second, triangle_mesh);
  tree2.accelerate_distance_queries(); //accelerated structure.
  time.stop();
  std::cout<<"CGAL AABB construction time : "<< time.time()<<std::endl;
  time.reset();

  {
    time.start();
    for (int i=0; i<num_points;i++){
      Point p (centre.x()+500*rg.normalisedRayDirections[i]._x, 
               centre.y()+500*rg.normalisedRayDirections[i]._y,
               centre.z()+500*rg.normalisedRayDirections[i]._z);
      auto p1 = tree1.closest_point(p);
    }
    time.stop();
    std::cout<<"Embree Closest_point : "<<time.time()<<std::endl;
    time.reset();
  }
  
  {
    time.start();
        for (int i=0; i<num_points;i++){
      Point p (centre.x()+500*rg.normalisedRayDirections[i]._x, 
               centre.y()+500*rg.normalisedRayDirections[i]._y,
               centre.z()+500*rg.normalisedRayDirections[i]._z);
      auto p2 = tree2.closest_point(p); 
    } 
    time.stop();
    std::cout<<"CGAL Closest_point : "<<time.time()<<std::endl;
    time.reset();
  }

  return 0;
}
