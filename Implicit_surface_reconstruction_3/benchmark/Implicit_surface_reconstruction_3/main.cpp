// ----------------------------------------------------------------------------
// Includes
// ----------------------------------------------------------------------------

#include <iostream>
#include <stdlib.h>

//CGAL
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

//Mesh
#include <CGAL/Surface_mesh.h>
#include <CGAL/property_map.h>

//file includes
#include "isr_test_util_reconstruction.h"
#include "isr_test_types.h"
#include "isr_benchmark_artefact_utils.h"
#include "isr_benchmark_dist_utils.h"

//boost
#include <boost/foreach.hpp>

//memory
#include <sys/resource.h>

// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------

int main(int argc, char **argv)
{
  int who = RUSAGE_SELF;
  struct rusage usage;
  int ret;

  if (argc != 3)
    return EXIT_FAILURE;

  //reads input file
  std::string curr_filename(argv[1]);
  PwnList pwnl;
  if(!get_point_set_from_file(curr_filename, pwnl)) {
    std::cout << "Unable to read file" << std::endl;
    std::cout << "Test skipped for this file" << std::endl << std::endl;
    return EXIT_FAILURE;
  }

  //parameters
  size_t curr_param = atoi(argv[2]);
  Param p0 = {false,  true, true, false, true, false}; //DR/P
  Param p1 = {false,  true, true, false, false, true};
  size_t p_len = 2;
  Param param_array[p_len] = {p0 , p1} ;

  //artefacts
  size_t a_len = 1;
  std::string art_array[a_len] = {"pos_noise"};

  //distances
  MeanDistPTM mean_dist_ptm;
  MeanDistPTM* ptr = &mean_dist_ptm; 
  size_t m_len = 1;
  Measure_type* meas_array[m_len] = {ptr};

  //starts timer
  CGAL::Timer task_timer; 
  task_timer.start();

  //reconstruts mesh
  Mesh reconstructed_mesh;
  if (!surface_mesh_reconstruction(param_array[curr_param],
              pwnl, reconstructed_mesh)) {
    std::cerr << "Error : Reconstruction failed" << std::endl;
    return EXIT_FAILURE;
  }

  //stores time
  double time_taken = task_timer.time();
  task_timer.stop();

  //stores memory peak
  ret = getrusage(who, &usage);
  long mem_peak = usage.ru_maxrss;

  size_t nb_pts = pwnl.size();
  std::cout << "DAT"  << "_file_" << "./dat_files/" << curr_param << "_mem_peak.dat" 
                      << "_xy_values_" << nb_pts << "\t" << mem_peak << "\n";
  std::cout << "DAT"  << "_file_" << "./dat_files/" << curr_param << "_time.dat" 
                      << "_xy_values_" << nb_pts << "\t" << time_taken << "\n";

  //for each artefact
  for (size_t art_idx = 0 ; art_idx < a_len ; ++art_idx) {
    //stores lvl 0
    std::cout << "DAT"  << "_file_" << "./dat_files/" << curr_param << "_" << art_array[art_idx] << "_mean_dist_ptm.dat" 
                        << "_xy_values_" << "0" << "\t" << mean_dist_ptm.run(reconstructed_mesh, pwnl) << "\n";
  
    //for each lvl
    for (size_t lvl = 1 ; lvl <=5 ; lvl++) {
      PwnList new_pwnl;
      position_noise_generator(pwnl, lvl, new_pwnl);

      Mesh new_mesh;
      if (!surface_mesh_reconstruction(param_array[curr_param],
            new_pwnl, new_mesh)) {
      std::cerr << "Error : Reconstruction failed" << std::endl;
      return EXIT_FAILURE;
      }

/*      MeanDistPTM mean_dist_ptm;
      std::cout << "DAT"  << "_file_" << "./dat_files/" << curr_param << "_" << art_array[0] << "_mean_dist_ptm.dat" 
                          << "_xy_values_" << lvl << "\t" << mean_dist_ptm.run(new_mesh, new_pwnl) << "\n";*/

      //for each measure type
      for (size_t mt_idx = 0 ; mt_idx < m_len ; ++mt_idx) {
      std::cout << "DAT"  << "_file_" << "./dat_files/" << curr_param << "_" << art_array[art_idx] << "_mean_dist_ptm.dat" 
                          << "_xy_values_" << lvl << "\t" << (meas_array[mt_idx])->run(new_mesh, new_pwnl) << "\n";
      }
    }
  }

  return EXIT_SUCCESS;
}