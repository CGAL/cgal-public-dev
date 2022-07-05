// ----------------------------------------------------------------------------
// Includes
// ----------------------------------------------------------------------------

#include <iostream>
#include <stdlib.h>

//CGAL
#include <CGAL/Random.h>

//Mesh
#include <CGAL/Surface_mesh.h>

//file includes
#include "isr_test_types.h"
#include "isr_test_util_bbox.h"
#include "isr_test_util_reconstruction.h"
#include "isr_test_io_utils.h"

//boost
#include <boost/foreach.hpp>

// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------

int main(int argc, char **argv) //arguments : 1.input file name, 2.output .xyz file name, 3.lvl
{
  if (argc != 4)
    return EXIT_FAILURE;

  //reads input .xyz file
  std::string input_filename(argv[1]);
  PwnList input_pwnl;
  if(!get_point_set_from_file(input_filename, input_pwnl)) {
    std::cerr << "Unable to read file" << std::endl;
    return EXIT_FAILURE;
  }

  //initialization 
  PwnList modified_pwnl;
  size_t lvl = atoi(argv[3]);
  double phip = (lvl/10.0) * std::acos(0);
  std::cout << phip << std::endl;
  double thetap = (lvl/10.0) * std::acos(0);

  //alters intput pwnlist
  BOOST_FOREACH(Point_with_normal pwn, input_pwnl) {
    Point p = pwn.first;
    Vector n = pwn.second;

    //get phi and theta 
    double phi = std::acos(n[1]);
    double theta = std::acos(n[0] / (std::sin(phi)));

    //get new L,h,l
    double phirdm = phip * CGAL::get_default_random().uniform_real<double>(-1.0,1.0);
    double thetardm = thetap * CGAL::get_default_random().uniform_real<double>(-1.0,1.0);

    double new_L = std::sin(phi+phirdm) * std::cos(theta+thetardm);
    double new_h = std::cos(phi+phirdm);
    double new_l = std::sin(phi+phirdm) * std::sin(theta+thetardm);
    Vector new_n(new_L, new_h, new_l);

    modified_pwnl.push_back(std::make_pair(p, new_n));
  }

  //stores altered pwnlist into .xyz file
  std::string output_filename(argv[2]);  
  std::ofstream output_xyz_file(output_filename);
  BOOST_FOREACH(Point_with_normal pwn , modified_pwnl) {
    Point p = pwn.first;
    Vector n = pwn.second;
    output_xyz_file << p << " " << n << "\n";
  }

  return EXIT_SUCCESS;
}