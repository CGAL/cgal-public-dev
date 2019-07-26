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
#include "include/isr_test_util_file_reading.h"

//boost
#include "boost/filesystem.hpp"
#include <boost/foreach.hpp>
#include <boost/property_map/property_map.hpp>


namespace PMP = CGAL::Polygon_mesh_processing;


// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------


long int test_check_genus(const std::string &input_file, const Param &parameter) 
{
  Mesh reconstructed_mesh;
  PwnList input_pwn;

  if (!mesh_reconstruction(input_file, parameter,
                input_pwn, reconstructed_mesh)) {
    std::cerr << "Error : Reconstruction failed" << std::endl;
    return (-1);
  }

  return (compute_genus(reconstructed_mesh));
}

bool test_check_genus_all_params(const std::string &input_file, const long int &in_gen)
{
  bool success = true;
  bool curr_par_success;
  Parameters plist;
  for (std::list<Param>::const_iterator param = plist.begin() ; param != plist.end() ; param++) {
    curr_par_success = true;
    std::cout << "///////////" << " " << *param << " "<< "///////////" << std::endl;
    long int out_gen = test_check_genus(input_file, *param);
    std::cout << "-> in_gen = " << in_gen << std::endl;
    if(out_gen >= 0)
      std::cout << "-> out_gen = " << out_gen << std::endl;
    else
      std::cout << "Unable to compute output genus because reconstruction failed" << std::endl;
    if (out_gen != in_gen) {
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
  int accumulated_fatal_err = EXIT_SUCCESS ;
  std::cerr << "|-------------------------------------------------------------------------|" << std::endl;
  std::cerr << "|                  TEST : CHECK IF GENUS IS PRESERVED                     |" << std::endl;
  std::cerr << "|-------------------------------------------------------------------------|" << std::endl << std::endl;

  boost::filesystem::path targetDir("./data/genus");
  boost::filesystem::recursive_directory_iterator iter(targetDir), eod;

  BOOST_FOREACH(boost::filesystem::path const& i, std::make_pair(iter, eod)) {
    if (is_regular_file(i)) {

      std::cout << "=============== Filename : " << i.string() << " ===============" << std::endl << std::endl;

      if(is_mesh(i.string())) //compute genus
      {
        Mesh input_m;
        if(!read_input_mesh_file(i.string(), input_m))
          return accumulated_fatal_err;
        size_t in_gen = compute_genus(input_m);
        if (!test_check_genus_all_params(i.string(), in_gen))
          accumulated_fatal_err = EXIT_FAILURE;  
      }
      else //get genus in file name if possible
      {
        std::string delimiter = "genus_" ;
        std::string str_gen = i.string();
        size_t pos = str_gen.find(delimiter);
        
        if (pos == std::string::npos) {
          std::cerr << "Impossible to find input genus" << std::endl;
          std::cerr << "Test skipped for this file" << std::endl;
          continue;
        }
        else {
          str_gen.erase(0, pos + delimiter.length());
          str_gen = str_gen.substr(0,str_gen.find("."));
          size_t in_gen = std::stoi(str_gen);
          if (!test_check_genus_all_params(i.string(), in_gen))
            accumulated_fatal_err = EXIT_FAILURE;
        }
      }
      
      std::cout << "=========================================================================" << std::endl << std::endl;

    }
  }

  return (accumulated_fatal_err);
}