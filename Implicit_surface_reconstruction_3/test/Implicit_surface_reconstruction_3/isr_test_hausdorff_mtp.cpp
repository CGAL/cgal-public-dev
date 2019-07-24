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

//Hausdorff & PMP
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>


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

bool test_hausdorff_mtp(const std::string &input_file, const Param &parameter) 
{
  Mesh reconstructed_mesh;
  PwnList input_pwn;

  if (!mesh_reconstruction(input_file, parameter,
                input_pwn, reconstructed_mesh)) {
    std::cerr << "Error : Reconstruction failed" << std::endl;
    return false;
  }

  double bbdiag = util_bb_diag(input_pwn);

	typedef typename PwnList::value_type PwnList_t;
  boost::function<Point(const PwnList_t&)> pwn_it_to_point_it = boost::bind(&PwnList_t::first, _1);
  double max_dist = PMP::approximate_max_distance_to_point_set(
                                            reconstructed_mesh,
                                            CGAL::make_range( boost::make_transform_iterator(input_pwn.cbegin(), pwn_it_to_point_it),
                                                               boost::make_transform_iterator(input_pwn.cend(), pwn_it_to_point_it)) ,
                                            4000 );

  std::cout << "-> haudorff_mtp = " << max_dist << std::endl;
  return ( max_dist * threshold_mult < bbdiag );
}

bool test_hausdorff_mtp_all_params(const std::string &input_file)
{
	bool success = true;
  bool curr_par_success;
	Parameters plist;
	for (std::list<Param>::const_iterator param = plist.begin() ; param != plist.end() ; param++) {
    curr_par_success = true;
    std::cout << "///////////" << " " << *param << " "<< "///////////" << std::endl;
		if (!test_hausdorff_mtp(input_file, *param)) {
			success = false ;
      curr_par_success = false;
    }
    std::cout << "/////////////////////////// " << (curr_par_success ? "PASSED" : "FAILED") << " ///////////////////////////" << std::endl;
		std::cout << std::endl;
	}
	return (success);
}

int	main()
{
	bool found_fail = false;
  std::cout << "|-------------------------------------------------------------------------|" << std::endl;
  std::cout << "|          TEST : HAUSDORFF DISTANCE FROM MESH TO INPUT POINTS            |" << std::endl;
  std::cout << "|-------------------------------------------------------------------------|" << std::endl << std::endl;

	boost::filesystem::path targetDir("./data/regular_data");
	boost::filesystem::recursive_directory_iterator iter(targetDir), eod;

	BOOST_FOREACH(boost::filesystem::path const& i, std::make_pair(iter, eod)) {
    if (is_regular_file(i)) {
    	std::cout << "=============== Filename : " << i.string() << " ===============" << std::endl << std::endl;
			if (!test_hausdorff_mtp_all_params(i.string())) 
				found_fail = true;
    	std::cout << "=========================================================================" << std::endl << std::endl;
  	}
	}

  int accumulated_fatal_err = found_fail ? EXIT_FAILURE : EXIT_SUCCESS ;
  return (accumulated_fatal_err);
}