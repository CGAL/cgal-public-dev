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
#include "include/isr_test_io_utils.h"

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

int threshold_ratio = 10;


// ----------------------------------------------------------------------------

class TestMeanDistMTP
{
  public:
  
  TestMeanDistMTP() {} ;

  bool run(const Param &parameter, PwnList &input_pwn)
  {
    Mesh reconstructed_mesh;

    if (!surface_mesh_reconstruction(parameter,
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
    return(mean_dist < bbdiag / threshold_ratio);
  }
};

int main()
{
  TestMeanDistMTP test_mean_dist_mtp;
  int accumulated_fatal_err = EXIT_SUCCESS ;
  std::cerr << "|-------------------------------------------------------------------------|" << std::endl;
  std::cerr << "|            TEST : MEAN DISTANCE FROM MESH TO INPUT POINTS               |" << std::endl;
  std::cerr << "|-------------------------------------------------------------------------|" << std::endl << std::endl;

  boost::filesystem::path targetDir("./data/regular_data");
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
      if (!test_all_param(test_mean_dist_mtp, pwnl)) 
        accumulated_fatal_err = EXIT_FAILURE; 

      std::cout << "=========================================================================" << std::endl << std::endl;
    }
  }

  return (accumulated_fatal_err);
}