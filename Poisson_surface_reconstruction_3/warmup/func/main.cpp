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
#include"random.h"

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
  typedef GT::Sphere_3 Sphere;
  typedef GT::Point_3 Point_3;
  typedef GT::FT FT;

  typedef Func<K, Point, Triangulation> Function;
  typedef FuncSmooth<K, Point, Triangulation> SmoothFunction;
  typedef typename CGAL::Implicit_surface_3<GT, Function> Surface_3;
  typedef typename CGAL::Implicit_surface_3<GT, SmoothFunction> Smooth_Surface_3;

  
  
  if(argc != 6){
    std::cout << "Usage: ./func <input file name> <isovalue ><sizing> <approximation> <output file name (without extension)>" << std::endl;
    return 0;
  }
  
  /*
  //creating tests with normals as input
  std::ofstream ofile("input1.xyzn");
  for(int i = 0; i < 50; i++){
    Point_3 p1 = random_point_on_sphere<Point_3>(1.0);
    Point_3 p2 = random_point_on_sphere<Point_3>(2.0);
    ofile << p1[0] << " " << p1[1] << " " << p1[2] << " " << 2.0 * p1[0] << " " << 2.0 * p1[1] << " " << 2.0 * p1[2] << " " << 1.0 << std::endl;
    ofile << p2[0] << " " << p2[1] << " " << p2[2] << " " << 2.0 * p2[0] << " " << 2.0 * p2[1] << " " << 2.0 * p2[2] << " " << 4.0 << std::endl;
         
  }
  */
  double isovalue = std::stod(argv[2]);
  double sizing = std::stod(argv[3]);
  double approximation = std::stod(argv[4]);
  std::cout << "isovalue: " << isovalue << std::endl;
  std::cout << "sizing: " << sizing << std::endl;
  std::cout << "approximation: " << approximation << std::endl;

  std::string output_filename(argv[5]);

  Triangulation tr;
  std::cout << "reading file..." << std::endl;
  std::string str(argv[1]);
  std::size_t pos = str.find(".");
  if(str.substr(pos) == ".xyz"){
    tr.read_xyz(argv[1]);
    
    tr.compute_grad_per_cell();
    tr.compute_grad_per_vertex();

    std::cout << "num vertices: " << tr.number_of_vertices() << std:: endl;
    std::cout << "done" << std::endl;
  }
  
  else if(str.substr(pos) == ".xyzn"){
    tr.read_xyzn(argv[1]);
    
    std::cout << "num vertices: " << tr.number_of_vertices() << std:: endl;
    std::cout << "done" << std::endl;
  }
  
  
  tr.output_grads_to_off();

  Tr t1, t2;
  C2t3 c2t3(t1);
  C2t3 c2t3_smooth(t2);
  Sphere bounding_sphere(CGAL::ORIGIN, 25.0);

  Function function(&tr, isovalue);
  SmoothFunction smooth_function(&tr, isovalue);

  const FT dichotomy = 1e-10;
  Surface_3 surface(function, bounding_sphere, dichotomy);
  Smooth_Surface_3 smooth_surface(smooth_function, bounding_sphere, dichotomy);

  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30, sizing, approximation);

  std::cout << "meshing...";
  make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_with_boundary_tag());
  std::cout << "done (" << c2t3.number_of_facets() << " facets)" << std::endl;

  std::cout << "smooth meshing...";
  make_surface_mesh(c2t3_smooth, smooth_surface, criteria, CGAL::Manifold_with_boundary_tag());
  std::cout << "done (" << c2t3_smooth.number_of_facets() << " facets)" << std::endl;

  if (c2t3.number_of_facets() > 0)
  {
	  std::ofstream out(output_filename + ".off");
	  CGAL::output_surface_facets_to_off(out, c2t3);
  }

  if (c2t3_smooth.number_of_facets() > 0)
  {
	  std::ofstream out(output_filename + "_smooth.off");
	  CGAL::output_surface_facets_to_off(out, c2t3_smooth);
  }

  return 0;
}
