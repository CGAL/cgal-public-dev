#include <fstream>
#include <iostream>

#include <CGAL/Embree/AABB_tree.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Simple_cartesian<float> K;
typedef K::Point_3 Point;
typedef K::Ray_3 Ray;
typedef K::Vector_3 Vector;
typedef CGAL::Surface_mesh<Point> Mesh;

typedef CGAL::Embree::Triangle_mesh_geometry<Mesh, K> TriangleMesh;
typedef CGAL::Embree::AABB_tree<TriangleMesh, K> Tree;

typedef boost::optional<Tree::Intersection_and_primitive_id> Ray_intersection;

int main(int argc, char const *argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "../../../AABB_tree/examples/AABB_tree/data/tetrahedron.off";
  std::ifstream input(filename);
  Mesh mesh;
  input >> mesh;

  Tree tree;
  tree.insert(mesh);

  Point rayOrigin(0.5f, 0.5f, -1.0f);
  Vector rayDirection(0.0f, 0.0f, 1.0f); /*Direction need not be normalized.*/
  Ray ray(rayOrigin, rayDirection);

  Ray_intersection intersection = tree.first_intersection(ray);
  if(intersection){
    Point p = intersection->first;
    std::cout<<"Point of intersection : "<<p<<std::endl;
  }

  return 0;
}
