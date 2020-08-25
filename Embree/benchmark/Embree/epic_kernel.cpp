#include <fstream>
#include <iostream>

#include <CGAL/Embree/AABB_tree.h>
#include <CGAL/Real_timer.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include "../../../AABB_tree/benchmark/AABB_tree/include/RaysGenerate.h"

template <class K>
double test (int argc, char const *argv[], const int numberOfRays){
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

  Point rayOrigin(-70.0f, -20.0f, 50.0f); //Point for gargoyle dataset. change it accordingly.
  RaysGenerate rg(numberOfRays);
  
  std::cout<<std::endl;
  std::cout<<"Bounding Box : "<<tree.bbox()<<std::endl;

  CGAL::Real_timer time;
  {
    time.start();
    for (int i=0; i<numberOfRays;i++){
        Vector v(rg.rayDirections[i]._x, rg.rayDirections[i]._y, rg.rayDirections[i]._z);
        Ray ray(rayOrigin, v);

        auto intersection = tree.first_intersection(ray);
    }
    time.stop();
  }
  std::cout<<"Time : ";
  return time.time();
}

int main(int argc, char const *argv[])
{
  int numRays = 1000000;
  std::cout<<"Exact_predicates_inexact_constructions_kernel"<<test<CGAL::Exact_predicates_inexact_constructions_kernel>(argc, argv, numRays)<<std::endl;
  std::cout<<std::endl;
  
  std::cout<<"CGAL::Simple_cartesian<float>"<<test<CGAL::Simple_cartesian<float>>(argc, argv, numRays)<<std::endl;
  
  return 0;
}

