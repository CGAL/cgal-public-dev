// ----------------------------------------------------------------------------
// Includes
// ----------------------------------------------------------------------------

#include <iostream>
#include <stdlib.h>

//Mesh
#include <CGAL/Surface_mesh.h>
#include <CGAL/property_map.h>

//file includes
#include "isr_test_types.h"
#include "isr_test_util_reconstruction.h"
#include "isr_benchmark_dist_utils.h"

//boost
#include <boost/foreach.hpp>

//memory
#include <sys/resource.h>

//Time
#include <CGAL/Timer.h>

// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------

int main(int argc, char **argv) //arguments : 1.input xyz file name wo ext, 2.output mesh file 3.parameter
{
  //std::cerr displays are for debugging
  int who = RUSAGE_SELF;
  struct rusage usage;
  size_t curr_param = atoi(argv[3]);
  Param p0 = {false,  true, true, false, true, false}; //DR/P
  Param p1 = {false,  true, true, false, false, true};
  size_t p_len = 2;
  Param param_array[p_len] = {p0 , p1} ;

  if (argc != 4)
    return EXIT_FAILURE;

  //reads input .xyz file
  std::string input_xyz_filename(argv[1]);
  std::cerr << "INPUT XYZ FILE " << input_xyz_filename << std::endl;
  std::cerr << "INPUT PARAMETER " << curr_param << std::endl;
  PwnList input_pwnl;
  if(!get_point_set_from_file(input_xyz_filename, input_pwnl)) {
    std::cerr << "Unable to read file" << std::endl;
    return EXIT_FAILURE;
  }

  //starts timer
  CGAL::Timer task_timer; 
  task_timer.start();

  //reconstruts mesh
  Mesh reconstructed_mesh;
  if (!surface_mesh_reconstruction(param_array[curr_param],
              input_pwnl, reconstructed_mesh)) {
    std::cerr << "Error : Reconstruction failed" << std::endl;
    return EXIT_FAILURE;
  }

  //prints time + memory peak
  double time_taken = task_timer.time();
  task_timer.stop();
  int ret = getrusage(who, &usage);
  long mem_peak = usage.ru_maxrss;
  std::cout << "DAT"  << "_nb_points_" << input_pwnl.size()
                      << "_time_"  << time_taken
                      << "_mem_peak_" << mem_peak << "\n";

  //computes distances
  Measure_type_list mt_list;

  BOOST_FOREACH(Measure_type* mt, mt_list.list) {
      std::cout << "DAT"  << "_measure_type_" << mt->get_name_in_file()
                          << "_value_" << mt->run(reconstructed_mesh, input_pwnl) << "\n";
  }

  //stores mesh into .off file
  std::string out_off_file(argv[2]);
  std::cerr << "OUTPUT OFF FILE NAME " << out_off_file << std::endl;
  std::ofstream out_mesh(out_off_file);
  out_mesh << reconstructed_mesh;
  out_mesh.close();

  return EXIT_SUCCESS;
}