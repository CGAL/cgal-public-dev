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


namespace PMP = CGAL::Polygon_mesh_processing;


// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------

class TestCheckGenus
{
  public :

  TestCheckGenus(size_t in_genus) : _in_genus(in_genus) {} ;

  bool run(const Param &parameter, PwnList &input_pwn)
  {
    Mesh reconstructed_mesh;

    if (!surface_mesh_reconstruction(parameter,
                  input_pwn, reconstructed_mesh)) {
      std::cerr << "Error : Reconstruction failed" << std::endl;
      return false;
    }

    size_t out_genus = compute_genus(reconstructed_mesh);
    std::cout << "-> in_genus = " << _in_genus << std::endl;
    std::cout << "-> out_genus = " << out_genus << std::endl;

    return (out_genus == _in_genus);
  }

  private :

  size_t _in_genus;

};

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

      //READS INPUT FILE
      PwnList pwnl;
      Mesh input_mesh;
      if(!get_point_set_from_file(i.string(), pwnl)) {
        std::cout << "Unable to read file" << std::endl;
        std::cout << "Test skipped for this file" << std::endl << std::endl;
        continue;
      }

      //TESTS
      if(is_mesh_file(i.string())) //compute genus
      {
        Mesh input_m;
        if(!read_mesh_file(i.string(), input_m))
          return accumulated_fatal_err;
        size_t in_gen = compute_genus(input_m);
        TestCheckGenus test_check_genus(in_gen);
        if (!test_all_param(test_check_genus, pwnl))
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
          TestCheckGenus test_check_genus(in_gen);
          if (!test_all_param(test_check_genus, pwnl))
            accumulated_fatal_err = EXIT_FAILURE;
        }
      }
      
      std::cout << "=========================================================================" << std::endl << std::endl;

    }
  }

  return (accumulated_fatal_err);
}