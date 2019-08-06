// ----------------------------------------------------------------------------
// Includes
// ----------------------------------------------------------------------------

#include <iostream>
#include <stdlib.h>

//CGAL
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>


#include <iostream>

//Mesh
#include <CGAL/Surface_mesh.h>
#include <CGAL/property_map.h>

//file includes
#include "isr_test_util_reconstruction.h"
#include "isr_test_types.h"

//boost
#include "boost/filesystem.hpp"
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

  //stores nb of points in memory dat file
  std::string mem_chart_file_name("./dat_files/mem_chart.dat");
  std::ofstream mem_chart(mem_chart_file_name, std::ios::app);
  int param = atoi(argv[2]);
  if (param == 0)
    mem_chart << pwnl.size() << "\t";
  mem_chart.close();

  //stores nb of points in time dat file
  std::string time_chart_file_name("./dat_files/time_chart.dat");
  std::ofstream time_chart(time_chart_file_name, std::ios::app);
  if (param == 0)
    time_chart << pwnl.size() << "\t";
  time_chart.close();

  //parameters
  Param p0 = {false,  true, true, false, true, false}; //DR/P
  Param p1 = {false,  true, true, false, false, true};
  Param param_array[2] = {p0 , p1} ;

  //starts timer
  CGAL::Timer task_timer; 
  task_timer.start();

  //reconstruts mesh
  Mesh reconstructed_mesh;
  if (!surface_mesh_reconstruction(param_array[param],
              pwnl, reconstructed_mesh)) {
    std::cerr << "Error : Reconstruction failed" << std::endl;
    return EXIT_FAILURE;
  }

  //stores time
  time_chart.open(time_chart_file_name, std::ios::app);
  time_chart << task_timer.time() << "\t" ;
  task_timer.stop();
  time_chart.close();

  //stores memory peak
  ret = getrusage(who, &usage); //this function is not compatible with Windows
  mem_chart.open(mem_chart_file_name, std::ios::app);
  mem_chart << usage.ru_maxrss << "\t" ;
  mem_chart.close();

  return EXIT_SUCCESS;
}