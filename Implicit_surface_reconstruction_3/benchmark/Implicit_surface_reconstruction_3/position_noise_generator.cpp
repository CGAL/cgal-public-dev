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

  //alters intput pwnlist
  PwnList modified_pwnl;
  size_t lvl = atoi(argv[3]);
  double delta = (util_bb_diag(input_pwnl) * lvl / 1000);
  BOOST_FOREACH(Point_with_normal pwn, input_pwnl) {
    Point p = pwn.first;
    Vector n = pwn.second;

    double new_x = p.x() + delta * CGAL::get_default_random().uniform_real<double>(-1.0,1.0);
    double new_y = p.y() + delta * CGAL::get_default_random().uniform_real<double>(-1.0,1.0);
    double new_z = p.z() + delta * CGAL::get_default_random().uniform_real<double>(-1.0,1.0);
    Point new_p(new_x, new_y, new_z);

    modified_pwnl.push_back(std::make_pair(new_p, n));
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