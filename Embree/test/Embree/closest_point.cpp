#include <iostream>
#include <fstream>

#include <CGAL/Embree/AABB_tree.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Simple_cartesian<float> K;
typedef K::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef CGAL::Embree::Triangle_mesh_geometry<Mesh, K> Geometry;
typedef CGAL::Embree::AABB_tree<Geometry, K> Tree;

int main(int argc, char const *argv[])
{
  const char* filename =  (argc > 1)? "data/bunny00.off" : argv[1];
  std::ifstream input(filename);

  Mesh triangle_mesh;
  input >> triangle_mesh;

  Tree tree;
  tree.insert(triangle_mesh);

  Point origin(0,0,0);

  Tree::Point_and_primitive_id ppid = tree.closest_point_and_primitive(origin);

  std::cout << ppid.first << std::endl;

  return 0;
}
