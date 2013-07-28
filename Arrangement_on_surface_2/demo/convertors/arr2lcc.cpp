//This program is to read input file of Arrangement and transfer
//it into Linear_cell_complex. An output file of Linear_cell_complex
//will be created
#include "convertors.h"

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/timer.hpp>
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
#include <CGAL/Combinatorial_map_operations.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel   Kernel;
typedef Kernel::FT                                          Number_type;
typedef CGAL::Arr_linear_traits_2<Kernel>                   Traits_2;
typedef Traits_2::Point_2                                   Point_2;
typedef Traits_2::Segment_2                                 Segment;
typedef Traits_2::Ray_2                                     Ray;
typedef Traits_2::Line_2                                    Line;
typedef Traits_2::X_monotone_curve_2                        X_monotone_curve;
typedef CGAL::Arrangement_2<Traits_2>                       Arrangement;

typedef CGAL::Linear_cell_complex_traits<2, Kernel>         Traits;
typedef CGAL::Linear_cell_complex<2, 2, Traits>             LCC;
typedef LCC::Dart_handle                                    Dart_handle;
typedef LCC::Point                                          Point;
typedef LCC::FT                                             FT;

int main(int argc, char* argv[])
{
  //declare the objects
  Arrangement arr;
  LCC lcc;
  Dart_handle dart;
    
  //Boost-options to override the input or output files
  boost::program_options::options_description desc;
  desc.add_options()
    ("help", "produce help")
    ("input,i", po::value< std::string>(), "native input file")
    ("output,o", po::value<std::string>(), "output file")
    ("bench 0", "turn off")
    ("bench 1", "turn on")
    ("input-curve", po::value< std::string>(), "including curves")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  //Check whether there exists input or output file
  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 1;
  }
  std::string filename;
   // std::ifstram fin;
  if (vm.count("input")) {
    filename = vm["input"].as<std::string>();
    std::ifstream fin(filename.c_str());
    if (fin.is_open()) {
      fin >> arr;
    } else {
      std::cout << "Invalid input file" << std::endl;
      fin.close();
      return EXIT_SUCCESS;
    }
    fin.close();
  } else if (vm.count("input-curve")) {
    filename = vm["input-curve"].as<std::string>();
    std::ifstream fin(filename.c_str());
    //Check whether the inpur file is valid
    if (fin.is_open()) {
      while (!fin.eof()) {
        Segment s;
        Number_type a, b, c, d;
        // Read the coordinates
        fin >> a >> b >> c >> d;
        //Create the segment
        s = Segment(Point_2(a, b), Point_2(c, d));
        // add this curve into the Arrangement
        insert(arr, s);
      }
    } else {
      std::cout << "Invalid input file" << std::endl;
      fin.close();
      return EXIT_SUCCESS;
    }
    fin.close();
  } else {
    filename = "arr.dat";
    //Read the input file
    std::ifstream fin(filename.c_str());
    //Check whether the inpur file is valid
    if (fin.is_open()) {
      fin >> arr;
    } else {
      std::cout << "Invalid input file" << std::endl;
      fin.close();
      return EXIT_SUCCESS;
    }
    fin.close();
  }
  
  //Call the function that transfer the Arrangement to Linear_cell_complex
  if (vm.count("bench 1")) {
    boost::timer t;
    dart = CGAL::arr2lcc<Arrangement, LCC>(arr, lcc);
    std::cout<<"time is "<<t.elapsed()<<std::endl;
  } else {
    dart = CGAL::arr2lcc<Arrangement, LCC>(arr, lcc);
  }
  //Get the number of vertices and edges
  int num_ver = arr.number_of_vertices();
  int num_edge = arr.number_of_edges();
    
    /*
    std::cout<<"LCC characteristics: " << std::endl;
    lcc.display_characteristics(std::cout) <<std::endl;
    std::cout<<"valid=" << lcc.is_valid() << std::endl;
     */
     
  //Write file
  std::ofstream myfile;
  //Check whether the output file name is typed
  if (vm.count("output")) {
    myfile.open((vm["output"].as<std::string>()).c_str());
  }
  else myfile.open("lcc.dat");
  //Write the number of vertices and edges
  myfile << num_ver << "  " << num_edge << std::endl;
  std::vector<Point> vert;
  //Write the vertices
  for (LCC::Vertex_attribute_range::iterator it =
       lcc.vertex_attributes().begin(), itend=lcc.vertex_attributes().end();
       it!=itend; ++it) {
    Point temp = it->point();
    vert.push_back(temp);
    myfile << it->point();
    if (vert.size() != num_ver) myfile << "  ";
    else myfile << std::endl;
  }
           
  //Write the indexes of vertices according to the associate edges
  int count = 0;
  for (LCC::One_dart_per_cell_range<1>::iterator
       it = lcc.one_dart_per_cell<1>().begin(),
       itend = lcc.one_dart_per_cell<1>().end();
         it != itend; ++it) {
    ++count;
    Point p1 = LCC::point(it);
    Point p2 = LCC::point(it->other_extremity());
    int first = 0, second = 0;
    bool b1 = false, b2 = false;
    //Find the indexes
    for (int i = 0; i<vert.size(); ++i) {
      if (vert[i] == p1) {
        first = i;
        b1 = true;
      } else if (vert[i] == p2) {
        second = i;
        b2 = true;
      }
      if (b1 && b2) break;
    }
    //Print the indexes in the output file
    myfile << first <<' ' << second;
    if (count != num_edge) myfile << "  ";
  }
  myfile.close();
  return 0;
}
