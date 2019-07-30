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
#include "include/isr_test_util_file_reading.h"
#include "include/isr_test_util_topo.h"

// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------


long int test_check_cc(const std::string &input_file, const Param &parameter, const int &i) 
{
  Mesh reconstructed_mesh;
  PwnList input_pwn;

  if (!mesh_reconstruction(input_file, parameter,
                input_pwn, reconstructed_mesh)) {
    std::cerr << "Error : Reconstruction failed" << std::endl;
    return (-1);
  }
      /*trucs de tests a enlever apres*/
      /*std::cout << "-> nb of input points : " << input_pwn.size() << std::endl;
      std::cout << "-> nb of mesh vertices : " << reconstructed_mesh.number_of_vertices() << std::endl;
      // saves reconstructed surface mesh
      std::string curr_outfile("test"+ std::to_string(i) +".off");
      std::ofstream out(curr_outfile);
      out << reconstructed_mesh;*/

  return (nb_cc(reconstructed_mesh));
}

bool test_check_cc_all_params(const std::string &input_file, const long int &in_cc)
{
  bool success = true;
  bool curr_par_success; // needs init curr_param_success
  Parameters plist(true);

  int j = 0; // to remove
  for (std::list<Param>::const_iterator param = plist.begin() ; param != plist.end() ; param++) { // for (const TestParameter &param : param_list)
    ++j;
    curr_par_success = true;
    std::cout << "///////////" << " " << *param << " " << "///////////" << std::endl;
    long int out_cc = test_check_cc(input_file, *param, j); // test_check_cc return true of false si le test passe, if false return false
    std::cout << "-> in_cc = " << in_cc << std::endl; // tout ça devrait être dans test_check_cc
    if(out_cc >= 0)
      std::cout << "-> out_cc = " << out_cc << std::endl;
    else
      std::cout << "Unable to compute output nb of connected components because reconstruction failed" << std::endl;
    if (out_cc != in_cc) {
      success = false ;
      curr_par_success = false;
    }

    std::cout << "/////////////////////////// " << (curr_par_success ? "PASSED" : "FAILED") << " ///////////////////////////" << std::endl;
    std::cout << std::endl;
  }
  return (success); // return true
}

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

      if(is_mesh(i.string())) //compute #cc
      {
        Mesh input_m;
        if(!read_input_mesh_file(i.string(), input_m))
          return accumulated_fatal_err;
        size_t in_cc = nb_cc(input_m);
        if (!test_check_cc_all_params(i.string(), in_cc))
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
          if (!test_check_cc_all_params(i.string(), in_cc))
            accumulated_fatal_err = EXIT_FAILURE;
        }
      }

      std::cout << "=========================================================================" << std::endl << std::endl;
    }      
  }

  return (accumulated_fatal_err);
}
