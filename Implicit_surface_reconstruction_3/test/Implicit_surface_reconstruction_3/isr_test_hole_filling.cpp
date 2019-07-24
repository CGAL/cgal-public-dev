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
#include "include/isr_test_util_nb_boundaries.h"

//boost
#include "boost/filesystem.hpp"
#include <boost/foreach.hpp>
#include <boost/property_map/property_map.hpp>

// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------

bool test_check_no_boundary(const std::string &input_file, const Param &parameter) 
{
  Mesh reconstructed_mesh;
  PwnList input_pwn;

  if (!mesh_reconstruction(input_file, parameter,
                input_pwn, reconstructed_mesh)) {
    std::cerr << "Error : Reconstruction failed" << std::endl;
    return false;
  }

  return (nb_boundaries(reconstructed_mesh) == 0);
}

bool test_check_no_boundary_all_params(const std::string &input_file)
{
  bool success = true;
  bool curr_par_success;
  Parameters plist;
  for (std::list<Param>::const_iterator param = plist.begin() ; param != plist.end() ; param++) {
    curr_par_success = true;
    std::cout << "///////////" << " " << *param << " "<< "///////////" << std::endl;
    if (!test_check_no_boundary(input_file, *param)) {
      success = false ;
      curr_par_success = false;
    }
    std::cout << "/////////////////////////// " << (curr_par_success ? "PASSED" : "FAILED") << " ///////////////////////////" << std::endl;
    std::cout << std::endl;
  }
  return (success);
}

int main()
{
  bool found_fail = false;
  std::cout << "|-------------------------------------------------------------------------|" << std::endl;
  std::cout << "|                  TEST : CHECK IF GENUS IS PRESERVED                     |" << std::endl;
  std::cout << "|-------------------------------------------------------------------------|" << std::endl << std::endl;

  boost::filesystem::path targetDir("./data/holes");
  boost::filesystem::recursive_directory_iterator iter(targetDir), eod;

  BOOST_FOREACH(boost::filesystem::path const& i, std::make_pair(iter, eod)) {
    if (is_regular_file(i)) {
      std::cout << "=============== Filename : " << i.string() << " ===============" << std::endl << std::endl;
      if (!test_check_no_boundary_all_params(i.string())) 
        found_fail = true;
      std::cout << "=========================================================================" << std::endl << std::endl;
    }      
  }

  int accumulated_fatal_err = found_fail ? EXIT_FAILURE : EXIT_SUCCESS ;
  return (accumulated_fatal_err);
}