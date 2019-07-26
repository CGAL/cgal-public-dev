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
#include "isr_test_util_file_reading.h"

//boost
#include "boost/filesystem.hpp"
#include <boost/foreach.hpp>
#include <boost/property_map/property_map.hpp>


// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------

//Mesh index
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef Mesh::Halfedge_index halfedge_descriptor;
typedef Mesh::Vertex_index Vertex_index;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

typedef CGAL::First_of_pair_property_map<Point_with_normal> Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;

//octree
typedef typename CGAL::OCTREE::Octree<Kernel, PwnList, Point_map, Normal_map> Octree;
typedef typename CGAL::OCTREE::Octree_node<Kernel, PwnList> Node;

// ----------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------


bool test_octree(const std::string &input_file, const std::pair<size_t, size_t> octree_params)
{
	bool success = true;
	PwnList pwnl;
	Mesh m;

	//FILE READING
	//input reading
	if (!read_file(input_file, pwnl, m)) {
		success = false;
		return success;
	}

	//requirements
	std::size_t nb_points = pwnl.size();

	if (nb_points == 0)
	{
		std::cerr << "Error: empty point set" << std::endl;
		success = false;
		return (success);
	}

	bool points_have_normals = (pwnl.begin()->second != CGAL::NULL_VECTOR);
	if ( ! points_have_normals )
	{
		std::cerr << "Input point set not supported: this reconstruction method requires unoriented normals" << std::endl;
		// this is not a bug => do not set success
	}


	//OCTREE CREATION
	//octree init
	PwnList octree_pwn;
	std::vector<Point> octree_steiner;
	Point_map pm;
	Normal_map nm;

	std::cout << "init octree...\n";
	Octree octree(pwnl, pm, nm);

	//octree refine
	size_t max_depth = octree_params.first;
	size_t max_pts_cell = octree_params.second;

	std::cout << "refine octree...\n";
	octree.refine(max_depth, max_pts_cell);

	//octree grading
	std::cout << "2:1 octree grading...\n";
	octree.grade();

	std::cout << "generate octree new points with normal...\n";
	octree.generate_points(std::back_inserter(octree_pwn), std::back_inserter(octree_steiner)); 


	//TESTS
	std::queue<Node *> leaf_nodes;
	Node *root = octree.root();
	octree.fill_leaf_queue(root, leaf_nodes);
	bool max_depth_exceeded = false;

	size_t nb_leaves = leaf_nodes.size();

	while (!leaf_nodes.empty())
	{
		Node *curr_node = leaf_nodes.front();
		if ((*curr_node).depth() > max_depth && !max_depth_exceeded) {
			success = false;
			max_depth_exceeded = true;
			std::cerr << "Error : Maximum depth exceeded" << std::endl;
		}

		if((*curr_node).num_points() > max_pts_cell && !max_depth_exceeded) {
			success = false;
			std::cerr << "Error : Maximum nb of points per cell exceeded" << std::endl;
			break;
		}
		leaf_nodes.pop();
	}

	if (!octree.debug_grading()) {
		success = false;
		std::cerr << "Error : Failed in grading octree" << std::endl;
	}

	if (octree.num_corner() + nb_leaves != octree_pwn.size() + octree_steiner.size()){
		success = false;
		std::cerr << "Error : Failed in generating the right amount of points";
	}

	std::cout << " -> " << (success ? "PASSED" : "FAILED") << std::endl;
	return (success);
}


int main()
{
	int accumulated_fatal_err = EXIT_SUCCESS ;
	std::cerr << "|-------------------------------------------------------------------------|" << std::endl;
	std::cerr << "|                              TEST : OCTREE                              |" << std::endl;
	std::cerr << "|-------------------------------------------------------------------------|" << std::endl << std::endl;

	boost::filesystem::path targetDir("./data/regular_data");
	boost::filesystem::recursive_directory_iterator iter(targetDir), eod;

	BOOST_FOREACH(boost::filesystem::path const& i, std::make_pair(iter, eod)) {
		if (is_regular_file(i)) {
			std::cout << "=============== Filename : " << i.string() << " ===============" << std::endl << std::endl;

			std::cout << " OCTREE PARAMS : max_depth = 9, max_nb_pts = 1" << std::endl;    
			if (!test_octree(i.string(), std::make_pair(9,1))) 
				accumulated_fatal_err = EXIT_FAILURE;
			std::cout << std::endl;

			std::cout << " OCTREE PARAMS : max_depth = 4, max_nb_pts = 1" << std::endl;    
			if (!test_octree(i.string(), std::make_pair(4,1))) 
				accumulated_fatal_err = EXIT_FAILURE;
			std::cout << std::endl;

			std::cout << "=========================================================================" << std::endl << std::endl;
		}
	}

	return (accumulated_fatal_err);
}