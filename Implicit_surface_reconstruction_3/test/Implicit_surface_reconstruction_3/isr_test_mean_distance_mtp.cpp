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
#include "include/isr_test_util_bbox.h"

//boost
#include "boost/filesystem.hpp"
#include <boost/foreach.hpp>
#include <boost/property_map/property_map.hpp>

//PMP
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>

//dD Tree
#include <CGAL/point_generators_2.h>


// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------

//dD_tree
typedef CGAL::Search_traits_3<Kernel> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree dD_Tree;

int threshold_mult = 10;


// ----------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------

bool test_mean_dist_mtp(const std::string &input_file, const Param &parameter) 
{
  Mesh reconstructed_mesh;
  PwnList input_pwn;
  if (!mesh_reconstruction(input_file, parameter,
                input_pwn, reconstructed_mesh)) {
    std::cerr << "Error : Reconstruction failed" << std::endl;
    return false;
  }
  double bbdiag = util_bb_diag(input_pwn);

  //sampling mesh
  std::list<Point> sample_points;
  CGAL::Polygon_mesh_processing::sample_triangle_mesh(reconstructed_mesh,
                                                      std::back_inserter(sample_points),
                                                      4000);

  //putting input points into dD_Tree
  typedef typename PwnList::value_type PwnList_t;
  boost::function<Point(const PwnList_t&)> pwn_it_to_point_it = boost::bind(&PwnList_t::first, _1);
  dD_Tree tree(boost::make_transform_iterator(input_pwn.begin(), pwn_it_to_point_it), 
                boost::make_transform_iterator(input_pwn.end(), pwn_it_to_point_it));

  //computation
  FT sum = 0;
  for (std::list<Point>::iterator it = sample_points.begin(); it != sample_points.end(); ++it) {
    Neighbor_search search(tree, *it, 1);
    sum += CGAL::sqrt(search.begin()->second);
  }

  double mean_dist = sum/(sample_points.size());
  std::cerr << "-> mean_dist_mtp = " << mean_dist << std::endl;
  return(mean_dist * threshold_mult < bbdiag);
}

bool test_mean_dist_mtp_all_params(const std::string &input_file)
{
  bool success = true;
  bool curr_par_success;
  Parameters plist;
  for (std::list<Param>::const_iterator param = plist.begin() ; param != plist.end() ; param++) {
    curr_par_success = true;
    std::cout << "///////////" << " " << *param << " "<< "///////////" << std::endl;
    if (!test_mean_dist_mtp(input_file, *param)) {
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
  std::cerr << "|            TEST : MEAN DISTANCE FROM MESH TO INPUT POINTS               |" << std::endl;
  std::cerr << "|-------------------------------------------------------------------------|" << std::endl << std::endl;

  boost::filesystem::path targetDir("./data/regular_data");
  boost::filesystem::recursive_directory_iterator iter(targetDir), eod;

  BOOST_FOREACH(boost::filesystem::path const& i, std::make_pair(iter, eod)) {
    if (is_regular_file(i)) {
      std::cout << "=============== Filename : " << i.string() << " ===============" << std::endl << std::endl; 
      if (!test_mean_dist_mtp_all_params(i.string())) 
        accumulated_fatal_err = EXIT_FAILURE;  
      std::cout << "=========================================================================" << std::endl << std::endl;
    }
  }

  return (accumulated_fatal_err);
}