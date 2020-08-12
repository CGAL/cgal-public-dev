#include <fstream>
#include <iostream>

#include <CGAL/Embree/AABB_tree.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Simple_cartesian<float> K;
typedef K::Point_3 Point;
typedef K::Segment_3 Segment;
typedef K::Vector_3 Vector;
typedef CGAL::Surface_mesh<Point> Mesh;

typedef CGAL::Embree::Triangle_mesh_geometry<Mesh, K> TriangleMesh;
typedef CGAL::Embree::AABB_tree<TriangleMesh, K> Tree;

typedef boost::optional<Tree::Intersection_and_primitive_id> Segment_intersection;

int main(int argc, char const *argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "../../../AABB_tree/examples/AABB_tree/data/tetrahedron.off";
  std::ifstream input(filename);
  Mesh mesh;
  input >> mesh;

  Tree tree;
  tree.insert(mesh);

  Point A(0.5f, 0.5f, -1.0f);
  Point B(0.5f, 0.5f, 0.5f);
  Point C(0.5f, 0.5f, -0.5f);

  Segment segment(A, B);
  Segment segment2(A, C);  
  
  Segment_intersection intersection = tree.first_intersection(segment);
  Segment_intersection intersection2 = tree.first_intersection(segment2);
    
  if(intersection){
    Point p = intersection->first;
    std::cout<<"Point of intersection : "<<p<<std::endl;
  }
  else
    std::cout<<"No intersection."<<std::endl;
  if(intersection2){
    Point p = intersection2->first;
    std::cout<<"Point of intersection : "<<p<<std::endl;
  }
  else
    std::cout<<"No intersection."<<std::endl;
    
  return 0;
}
