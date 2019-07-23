// ----------------------------------------------------------------------------
// Includes
// ----------------------------------------------------------------------------
#include <iostream>
#include <filesystem>

#include "boost/filesystem.hpp"
#include <boost/foreach.hpp>

//files includes
#include "include/isr_test_util_bbox.h"
#include "include/isr_test_util_reconstruction.h"

//Mesh
#include <CGAL/Surface_mesh.h>
#include <CGAL/property_map.h>
#include <boost/property_map/property_map.hpp>

#include "include/isr_test_types.h"
#include "include/isr_test_util_nb_boundaries.h"


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
  Parameters plist;
  for (std::list<Param>::const_iterator param = plist.begin() ; param != plist.end() ; param++) {
    std::cout << *param << std::endl;
    if (!test_check_no_boundary(input_file, *param))
      success = false ;
    std::cout << (success ? "Passed" : "Failed") << std::endl ;
    std::cout << std::endl;
  }
  return (success);
}

int main()
{
  bool found_fail = false;
  std::cerr << "Test : Check if genus is preserved" << std::endl << std::endl;

  boost::filesystem::path targetDir("./data/holes");
  boost::filesystem::recursive_directory_iterator iter(targetDir), eod;

  BOOST_FOREACH(boost::filesystem::path const& i, std::make_pair(iter, eod)) {
    if (is_regular_file(i)) {
      std::cout << "Filename : " << i.string() << std::endl;
      if (!test_check_no_boundary_all_params(i.string())) 
        found_fail = true;
      std::cout << std::endl << std::endl;\
    }      
  }

  int accumulated_fatal_err = found_fail ? EXIT_FAILURE : EXIT_SUCCESS ;
  return (accumulated_fatal_err);
}