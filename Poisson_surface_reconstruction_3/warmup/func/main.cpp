#include <iostream>
#include <fstream>
#include <string>
#include "triangulation.h"
#include "function.h"
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_cell_base_3.h>

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>


// TODO: use argc and argv to call directly eg
// func test.xyz 0.01 0.001 (sizing, approximation)
int main(int argc, char** argv){

  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef K::FT FT;
  typedef K::Point_3 Point;

  typedef VB<K> Vb;
  typedef CB<K> Cb;
  typedef CGAL::Triangulation_data_structure_3<Vb, Cb> TDS;
  typedef Min_triangulation_3D<K, TDS> Triangulation;
  typedef Triangulation::Vertex_handle Vertex_handle;

  typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
  typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
  typedef Tr::Geom_traits GT;
  typedef GT::Sphere_3 Sphere_3;
  typedef GT::Point_3 Point_3;
  typedef GT::FT FT;

  typedef Func<K, Point, Triangulation> Function;

  typedef typename CGAL::Implicit_surface_3<GT, Function> Surface_3;

  if(argc != 4){
    std::cout << "Usage: ./func <filename> <sizing> <approximation>" << std::endl;
    return 0;
  }
  double sizing = std::stod(argv[2]); double approximation = std::stod(argv[3]);

  Triangulation tr;
  std::cout << "reading file...";
  tr.read_xyz(argv[1]);
  std::cout << "done" << std::endl;

  tr.compute_grad_per_cell();
  tr.compute_grad_per_vertex();

  Tr t;
  C2t3 c2t3(t);
  Sphere_3 bounding_sphere(CGAL::ORIGIN, 25.0);

  Function function(&tr);

  const FT dichotomy = 1e-10;
  Surface_3 surface(function, bounding_sphere, dichotomy);

  std::cout << "meshing...";
  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30, sizing, approximation);
  make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_with_boundary_tag());
  std::cout << "done (" << c2t3.number_of_facets() << " facets)" << std::endl;

  if (c2t3.number_of_facets() > 0)
  {
	  std::ofstream out("output.off");
	  CGAL::output_surface_facets_to_off(out, c2t3);
  }

  return 0;
}
