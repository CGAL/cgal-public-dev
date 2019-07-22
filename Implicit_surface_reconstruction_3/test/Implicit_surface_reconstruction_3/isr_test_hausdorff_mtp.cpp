// ----------------------------------------------------------------------------
// Includes
// ----------------------------------------------------------------------------
#include <iostream>
#include <filesystem>

#include "boost/filesystem.hpp"
#include <boost/foreach.hpp>

//files includes
#include "include/isr_test_util_bbox.h"
#include "include/isr_test_util_reconstruction.h"

//Mesh
#include <CGAL/Surface_mesh.h>

//Hausdorff & PMP
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>

// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------
#include "include/isr_test_types.h"

//dD_tree
typedef CGAL::Search_traits_3<Kernel> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree dD_Tree;

int threshold = 10; /*changer le nom*/

// ----------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------

bool test_hausdorff_mtp_param(std::string input_file, const Param &parameter) 
{
  Mesh reconstructed_mesh;
  PwnList input_pwn;

  if (!reconstruction_param(reconstructed_mesh, input_pwn,
                parameter, input_file)) {
    std::cerr << "reconstruction failed" << std::endl;
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

  std::cerr << "haudorff_mtp = " << max_dist << std::endl;
  return ( max_dist * threshold < bbdiag );
}

bool test_param(std::string input_file)
{
	bool success = true;
	Parameters plist;
	for (std::list<Param>::const_iterator param = plist.begin() ; param != plist.end() ; param++) {
		std::cout << *param << std::endl;
		if (!test_hausdorff_mtp_param(input_file, *param))
			success = false ;
		std::cout << (success ? "Passed" : "Failed") << std::endl ;
		std::cout << std::endl;
	}
	return (success);
}

int	main()
{
	bool found_fail = false;
	std::cerr << "Test : Hausdorff distance from mesh to input points" << std::endl << std::endl;

	boost::filesystem::path targetDir("./data/");
	boost::filesystem::recursive_directory_iterator iter(targetDir), eod;

	BOOST_FOREACH(boost::filesystem::path const& i, std::make_pair(iter, eod)) {
    if (is_regular_file(i) && ((i.string()).find("big_data") == std::string::npos)) {
    	std::cout << "Filename : " << i.string() << std::endl;
			if (!test_param(i.string())) 
				found_fail = true;
    	std::cout << std::endl << std::endl;
  	}
	}

  int accumulated_fatal_err = found_fail ? EXIT_FAILURE : EXIT_SUCCESS ;
  return (accumulated_fatal_err);
}