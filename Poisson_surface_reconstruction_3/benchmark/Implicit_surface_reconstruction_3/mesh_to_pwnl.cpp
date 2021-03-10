// ----------------------------------------------------------------------------
// Includes
// ----------------------------------------------------------------------------

#include <iostream>
#include <stdlib.h>

//Mesh
#include <CGAL/Surface_mesh.h>

//file includes
#include "isr_test_types.h"
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

int main(int argc, char **argv) //arguments : 1. input file name, 2. output file name
{
  if (argc != 3)
    return EXIT_FAILURE;

  std::string input_filename(argv[1]);
  if (is_mesh_file(input_filename)) {
    //convert mesh into pwnl
    PwnList converted_pwnl;
    if(!get_point_set_from_file(input_filename, converted_pwnl)) {
      std::cerr << "Unable to read file" << std::endl;
      return EXIT_FAILURE;
    }

    //stores converted pwnl into .xyz file
    std::string output_filename(argv[2]);  
    std::ofstream output_xyz_file(output_filename);
    BOOST_FOREACH(Point_with_normal pwn , converted_pwnl) {
      Point p = pwn.first;
      Vector n = pwn.second;
      output_xyz_file << p << " " << n << "\n";    
    }
  }
  return EXIT_SUCCESS;
}