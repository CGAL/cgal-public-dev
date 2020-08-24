#include <fstream>
#include <iostream>

#include <CGAL/Embree/AABB_tree.h>
#include <CGAL/Real_timer.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

template <class K>
double test (int argc, char const *argv[]){
  // typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef typename K::Point_3 Point;
  typedef typename K::Ray_3 Ray;
  typedef typename K::Vector_3 Vector;
  typedef typename CGAL::Surface_mesh<Point> Mesh;

  typedef CGAL::Embree::Triangle_mesh_geometry<Mesh, K> TriangleMesh;
  typedef CGAL::Embree::AABB_tree<TriangleMesh, K> Tree;

  const char* filename = (argc > 1) ? argv[1] : "../../../AABB_tree/examples/AABB_tree/data/tetrahedron.off";

  std::ifstream input(filename);
  Mesh mesh;
  input >> mesh;

  Tree tree;
  tree.insert(mesh);

  Point rayOrigin(0.5f, 0.5f, -1.0f);
  Vector rayDirection(0.0f, 0.0f, 1.0f); /*Direction need not be normalized.*/
  Ray ray(rayOrigin, rayDirection);
  
  std::cout<<std::endl;
  std::cout<<"Ray information : "<<ray<<std::endl;
  std::cout<<"Bounding Box : "<<tree.bbox()<<std::endl;

  CGAL::Real_timer time;
  {
    time.start();
    auto intersection = tree.first_intersection(ray);
    time.stop();
  }
  std::cout<<"Time : ";
  return time.time();
}

int main(int argc, char const *argv[])
{
  std::cout<<"Exact_predicates_inexact_constructions_kernel"<<test<CGAL::Exact_predicates_inexact_constructions_kernel>(argc, argv)<<std::endl;
  std::cout<<std::endl;
  
  std::cout<<"CGAL::Simple_cartesian<float>"<<test<CGAL::Simple_cartesian<float>>(argc, argv)<<std::endl;
  
  return 0;
}

