// ----------------------------------------------------------------------------
// Includes
// ----------------------------------------------------------------------------

#include <iostream>

//Mesh
#include <CGAL/Surface_mesh.h>
#include <CGAL/property_map.h>

//file includes
#include "include/isr_test_util_reconstruction.h"
#include "include/isr_test_types.h"
#include "include/isr_test_util_topo.h"
#include "include/isr_test_io_utils.h"

//boost
#include "boost/filesystem.hpp"
#include <boost/foreach.hpp>
#include <boost/property_map/property_map.hpp>

// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------

class TestHoleFilling
{
  public :

  TestHoleFilling() {} ;

  bool run(const Param &parameter, PwnList &input_pwn)
  {
    Mesh reconstructed_mesh;

    if (!surface_mesh_reconstruction(parameter,
                  input_pwn, reconstructed_mesh)) {
      std::cerr << "Error : Reconstruction failed" << std::endl;
      return false;
    }

    if (  (nb_boundaries(reconstructed_mesh) != 0) || //check #boundaries == 0
          (compute_genus(reconstructed_mesh) != 0) || //check genus == 0
          (nb_cc(reconstructed_mesh) != 1) )          //check #cc == 1
      return false;

    return true;
  }

};


int main()
{
  TestHoleFilling test_hole_filling;
  int accumulated_fatal_err = EXIT_SUCCESS ;
  std::cout << "|-------------------------------------------------------------------------|" << std::endl;
  std::cout << "|                            TEST : HOLE FILLING                          |" << std::endl;
  std::cout << "|-------------------------------------------------------------------------|" << std::endl << std::endl;

  boost::filesystem::path targetDir("./data/holes");
  boost::filesystem::recursive_directory_iterator iter(targetDir), eod;

  BOOST_FOREACH(boost::filesystem::path const& i, std::make_pair(iter, eod)) {

    //READS INPUT FILE
      PwnList pwnl;
      if(!get_point_set_from_file(i.string(), pwnl)) {
        std::cout << "Unable to read file" << std::endl;
        std::cout << "Test skipped for this file" << std::endl << std::endl;
        continue;
      }

      //TESTS
    if (is_regular_file(i)) {
      std::cout << "=============== Filename : " << i.string() << " ===============" << std::endl << std::endl;
      if (!test_all_param(test_hole_filling, pwnl)) 
        accumulated_fatal_err = EXIT_FAILURE;
      std::cout << "=========================================================================" << std::endl << std::endl;
    }      
  }

  return accumulated_fatal_err;
}