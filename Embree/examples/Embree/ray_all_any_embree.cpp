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

  Point rayOrigin(0.1f, 0.1f, -1.0f);
  Vector rayDirection(0.0f, 0.0f, 1.0f);
  Ray ray(rayOrigin, rayDirection);

  std::cout<<"Any intersection."<<std::endl;
  Ray_intersection intersection = tree.any_intersection(ray);
  if(intersection){
    Point p = intersection->first;
    std::cout<<"Intersection point : "<<p<<std::endl;
  }
  std::cout<<"---------------------"<<std::endl;
  std::vector<Ray_intersection> intersections;

  std::cout<<"All intersections"<<std::endl;
  tree.all_intersections(ray, std::back_inserter(intersections));
    for (int i=0;i<intersections.size();i++){
      if(intersections[i]){
        Point p = intersections[i]->first;
        std::cout<<"Intersection point : "<<p<<std::endl;
      }
    }        
  
  return 0;
}
