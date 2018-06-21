#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Point_3.h>
#include "triangulation.h"
#include "grad_fit.h"
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Point_3<K> Point;
typedef CGAL::Vector_3<K> Vector;



int main(int argc, char** argv){
  typedef VB<K> Vb;
  typedef CB<K> Cb;
  typedef CGAL::Triangulation_data_structure_3<Vb, Cb> TDS;
  typedef Min_triangulation_3D<K, TDS> Triangulation;
  typedef Triangulation::Vertex_handle Vertex_handle;

  if(argc != 2)
  {
    std::cout << "Usage: ./func <input file name> " << std::endl;
    return 0;
  }

  Triangulation tr;
  std::cout << "reading file..." << std::endl;
  std::string str(argv[1]);
  std::size_t pos = str.find(".");
  if(str.substr(pos) == ".xyz")
  {
    tr.read_xyz(argv[1]);
    std::cout << "num vertices: " << tr.number_of_vertices() << std:: endl;
    std::cout << "done" << std::endl;
  }

  std::set<Vertex_handle> vertices;
  Vertex_handle v;
  for( auto it = tr.finite_vertices_begin(); it != tr.finite_vertices_end();
   it++ )
  {
    if(!(it->point() == Point(1, 0, 0)))
      vertices.insert(it);
    else v = it;
  }
  //std::cout << vertices.size() << std::endl;
  grad_fit(vertices, v);
  return 0;
}
