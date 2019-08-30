// ----------------------------------------------------------------------------
// Includes
// ----------------------------------------------------------------------------

#include <iostream>
#include <stdlib.h>

//CGAL
#include <CGAL/Random.h>
#include <CGAL/Origin.h>

//Mesh
#include <CGAL/Surface_mesh.h>

//file includes
#include "isr_test_types.h"
#include "isr_test_util_bbox.h"
#include "isr_test_util_reconstruction.h"
#include "isr_test_io_utils.h"

//boost
#include <boost/foreach.hpp>

// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------

typename Kernel::Point_3 origin = CGAL::ORIGIN;
// ----------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------

int main(int argc, char **argv) //arguments : 1.input file name, 2.output .xyz file name, 3.lvl
{
  if (argc != 4)
    return EXIT_FAILURE;

  //reads input file
  std::string input_filename(argv[1]);
  PwnList input_pwnl;
  if(!get_point_set_from_file(input_filename, input_pwnl)) {
    std::cerr << "Unable to read file" << std::endl;
    return EXIT_FAILURE;
  }

  //index points
  std::unordered_map<int, std::pair<Point_with_normal,bool>> pt_map; //bool = true if point is getting modified
  size_t i = 0;
  BOOST_FOREACH(Point_with_normal pwn, input_pwnl) {
    pt_map[i] = std::make_pair(pwn,false);
    i++; 
  }

  //alters intput pwnlist
  size_t lvl = atoi(argv[3]);
  size_t nb_clusters = lvl/2+1;
  const size_t b = 0;
  const size_t nb_pts = input_pwnl.size();

  size_t rdm_point_index = 0;
  //for each cluster, move the point with delta but by the point's normal
  while (nb_clusters) {
    rdm_point_index = CGAL::get_default_random().uniform_int(b,nb_pts-1);
    std::pair<Point_with_normal,bool> curr_ptb = pt_map[rdm_point_index];
    if (!curr_ptb.second) {
      Point_with_normal pwn = curr_ptb.first;
      Point p = pwn.first; 
      Vector n = pwn.second;
      double dist = CGAL::sqrt(CGAL::squared_distance(p, origin));
      double curr_delta = (lvl+4) * dist / 40.0 * CGAL::get_default_random().uniform_real<double>(0.75,1.0);
      p = p + n * curr_delta;
      pt_map[rdm_point_index] = std::make_pair(std::make_pair(p,n),true);
      nb_clusters--;
    }
  }

  //stores altered pwnlist into .xyz file
  std::string output_filename(argv[2]);  
  std::ofstream output_xyz_file(output_filename);
  for (i = 0 ; i  < nb_pts ; ++i) {
    std::pair<Point_with_normal,bool> ptb = pt_map[i];
    Point_with_normal pwn = ptb.first;
    Point p = pwn.first;
    Vector n = pwn.second;
    output_xyz_file << p << " " << n << "\n";
  }

  return EXIT_SUCCESS;
}