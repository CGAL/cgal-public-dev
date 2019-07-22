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

//AABB_tree
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------
#include "include/isr_test_types.h"

//AABB Tree
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;

int threshold = 10; /*changer le nom*/

// ----------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------

bool test_mean_dist_ptm_param(std::string input_file, const Param &parameter)
{
	Mesh reconstructed_mesh;
	PwnList input_pwn;
	if (!reconstruction_param(reconstructed_mesh, input_pwn,
								parameter, input_file)) {
		std::cerr << "reconstruction failed" << std::endl;
		return false;
	}
	double bbdiag = util_bb_diag(input_pwn);

	//charging faces into AABB Tree
	Tree tree(faces(reconstructed_mesh).first, faces(reconstructed_mesh).second, reconstructed_mesh);

  //computation
  tree.accelerate_distance_queries();
  FT sum;
  FT sqd_dist;

  for (PwnList::const_iterator it = input_pwn.begin(); it != input_pwn.end(); ++it) {
    const Point& current_pt = it->first;
    sqd_dist = tree.squared_distance(current_pt);
    sum += CGAL::sqrt(sqd_dist);
  }
  FT mean_dist = sum / (input_pwn.size());

  std::cerr << "mean_dist_ptm = " << mean_dist << std::endl;
  return( mean_dist * threshold < bbdiag);
}

bool test_param(std::string input_file)
{
	bool success = true;
	Parameters plist;
	for (std::list<Param>::const_iterator param = plist.begin() ; param != plist.end() ; param++) {
		std::cout << *param << std::endl;
		if (!test_mean_dist_ptm_param(input_file, *param))
			success = false;
		std::cout << (success ? "Passed" : "Failed") << std::endl ;
		std::cout << std::endl;
	}
	return (success);
}

int	main()
{
	bool found_fail = false;
	std::cerr << "Test : Mean distance from input points to mesh" << std::endl << std::endl;

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

	std::cout << std::endl;

  int accumulated_fatal_err = found_fail ? EXIT_FAILURE : EXIT_SUCCESS ;
  return (accumulated_fatal_err);
}