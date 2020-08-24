#include <iostream>
#include <fstream>

#include <CGAL/Embree/AABB_tree.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/Real_timer.h>

typedef CGAL::Simple_cartesian<float> K;
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
  const char* filename =  (argc > 1) ? "../../test/Embree/data/bunny00.off" : argv[1];
  std::ifstream input(filename);

  Mesh triangle_mesh;
  input >> triangle_mesh;

  TreeEmbree tree1;
  tree1.insert(triangle_mesh);

  TreeCgal tree2(faces(triangle_mesh).first, faces(triangle_mesh).second, triangle_mesh);
  tree2.build();

  Point origin(0.1f,0.1f,-1.0f);
  {
    CGAL::Real_timer time;
    time.start();
    auto p1 = tree1.closest_point(origin);
    time.stop();
    std::cout<<"Embree Closest_point : "<<time.time()<<std::endl;
    time.reset();
  }
  
  {
    CGAL::Real_timer time;
    time.start();
    auto p2 = tree2.closest_point(origin);  
    time.stop();
    std::cout<<"CGAL Closest_point : "<<time.time()<<std::endl;
    time.reset();
  }

  return 0;
}
