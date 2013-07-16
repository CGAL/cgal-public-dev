
#include "convertors.h"

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
namespace po = boost::program_options;

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/IO/Arr_iostream.h>

#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Linear_cell_complex_constructors.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel   Kernel;
typedef Kernel::FT                                          Number_type;
typedef CGAL::Arr_linear_traits_2<Kernel>                   Traits_2;
typedef Traits_2::Point_2                                   Point_2;
typedef Traits_2::Segment_2                                 Segment;
typedef Traits_2::Ray_2                                     Ray;
typedef Traits_2::Line_2                                    Line;
typedef Traits_2::X_monotone_curve_2                        X_monotone_curve;
typedef CGAL::Arrangement_2<Traits_2>                       Arrangement;
typedef Arrangement::Vertex_handle                          Vertex_handle;
typedef Arrangement::Halfedge_handle                        Halfedge_handle;
typedef Arrangement::Face_handle                            Face_handle;

typedef CGAL::Linear_cell_complex_traits<2, Kernel>         Traits;
typedef CGAL::Linear_cell_complex<2, 2, Traits>             LCC;
typedef LCC::Dart_handle                                    Dart_handle;
typedef LCC::Point                                          Point;
typedef LCC::FT                                             FT;

int main(int argc, char* argv[])
{
  LCC lcc;
  Arrangement arr;
  Dart_handle dh;
  Halfedge_handle he1;
   
  //Boost-options to override the input or output files
  boost::program_options::options_description desc;
  desc.add_options()
    ("help", "produce help")
    ("input,i", po::value< std::string>(), "input file")
    ("output,o", po::value<std::string>(),"output file")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  //Check whether there exists the input file
  std::string filename;
  if (vm.count("input")) {
    filename = vm["input"].as<std::string>();
  } else {
    filename = "lcc.dat";
  }
  //Read the input file
  std::ifstream lcc_file(filename.c_str());
  //Check whether the inpur file is valid
  if (lcc_file.is_open()) {
    dh = import_from_plane_graph(lcc, lcc_file);
  } else {
    std::cout << "Invalid input file" << std::endl;
    return EXIT_SUCCESS;
  }
  lcc_file.close();
  //Call the function to make the transfer
  he1 = CGAL::lcc2arr<LCC, Arrangement>(lcc, arr);
  //Write file
  std::ofstream myfile;
  //Check whether there exists the output file
  if (vm.count("output")) {
    myfile.open((vm["output"].as<std::string>()).c_str());
  }
  else myfile.open("arr.dat");
  //Write the associate coordinates corresponding to the segments
    std::cout<<"#vertices: " << arr.number_of_vertices()<<std::endl;
    std::cout<<"#edges: " << arr.number_of_edges()<<std::endl;
    std::cout<<"#faces: " << arr.number_of_faces()<<std::endl;

    /*
  int num_edge = arr.number_of_edges();
  int count = 0;
  Arrangement::Edge_const_iterator he;
  for (he = arr.edges_begin(); he != arr.edges_end(); ++he) {
    Point_2 p1 = he->source()->point();
    Point_2 p2 = he->target()->point();
    ++count;
    if (count == num_edge) myfile << p1 << "  " << p2;
    else myfile << p1 << "  " << p2 << std::endl;
  }
     */
    myfile << arr;
  myfile.close();
  return 0;
}
