// ----------------------------------------------------------------------------
// Includes
// ----------------------------------------------------------------------------

#include <iostream>

//boost
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/property_map/property_map.hpp>

//Mesh
#include <CGAL/Surface_mesh.h>
#include <CGAL/property_map.h>

//file includes
#include "include/isr_test_util_reconstruction.h"
#include "include/isr_test_types.h"
#include "include/isr_test_io_utils.h"
#include "include/isr_test_util_topo.h"

// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------

class TestCheckCC
{
  public:
  
  TestCheckCC(size_t in_cc) : _in_cc(in_cc) {} ;

  bool run(const Param &parameter, PwnList &input_pwn)
  {
    Mesh reconstructed_mesh;

    if (!surface_mesh_reconstruction(parameter,
                input_pwn, reconstructed_mesh)) {
    std::cerr << "Error : Reconstruction failed" << std::endl;
    return false;
    }

    size_t out_cc = nb_cc(reconstructed_mesh);
    std::cout << "-> in_cc = " << _in_cc << std::endl;
    std::cout << "-> out_cc = " << out_cc << std::endl;

    return (out_cc == _in_cc);
  }

  private :

  size_t _in_cc;
  
};

int main()
{
  int accumulated_fatal_err = EXIT_SUCCESS ;
  std::cout << "|-------------------------------------------------------------------------|" << std::endl;
  std::cout << "|         TEST : CHECK IF NB OF CONNECTED COMPONENTS IS PRESERVED         |" << std::endl;
  std::cout << "|-------------------------------------------------------------------------|" << std::endl << std::endl;

  boost::filesystem::path targetDir("./data/con_comp");
  boost::filesystem::recursive_directory_iterator iter(targetDir), eod;

  BOOST_FOREACH(boost::filesystem::path const& i, std::make_pair(iter, eod)) {
    if (is_regular_file(i)) {
      std::cout << "=============== Filename : " << i.string() << " ===============" << std::endl << std::endl;

      //READS INPUT FILE
      PwnList pwnl;

      if(!get_point_set_from_file(i.string(), pwnl)) {
        std::cout << "Unable to read file" << std::endl;
        std::cout << "Test skipped for this file" << std::endl << std::endl;
        continue;
      }

      //TESTS
      if(is_mesh_file(i.string())) //compute #cc
      {
        Mesh input_m;
        if(!read_mesh_file(i.string(), input_m))
          return accumulated_fatal_err;
        size_t in_cc = nb_cc(input_m);
        TestCheckCC test_check_cc(in_cc);
        if (!test_all_param(test_check_cc, pwnl, true))
          accumulated_fatal_err = EXIT_FAILURE;
      }
      else //get #cc in file name if possible
      {
        std::string delimiter = "cc_" ;
        std::string str_cc = i.string();
        size_t pos = str_cc.find(delimiter);
        
        if (pos == std::string::npos) {
          std::cerr << "Impossible to find input nb of connected components" << std::endl;
          std::cerr << "Test skipped for this file" << std::endl << std::endl;
          continue;
        }
        else {
          str_cc.erase(0, pos + delimiter.length());
          str_cc = str_cc.substr(0,str_cc.find("."));
          size_t in_cc = std::stoi(str_cc);
          TestCheckCC test_check_cc(in_cc);
          if (!test_all_param(test_check_cc, pwnl, true))
            accumulated_fatal_err = EXIT_FAILURE;
        }
      }

      std::cout << "=========================================================================" << std::endl << std::endl;
    }      
  }

  return (accumulated_fatal_err);
}
