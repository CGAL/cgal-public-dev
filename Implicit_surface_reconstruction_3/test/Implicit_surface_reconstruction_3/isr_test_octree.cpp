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
#include "isr_test_io_utils.h"

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


bool test_octree(const std::pair<size_t, size_t> octree_params, PwnList &pwnl)
{
  bool success = true;
  
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

  size_t nb_leaves = leaf_nodes.size();

  while (!leaf_nodes.empty())
  {
    Node *curr_node = leaf_nodes.front();

    if (curr_node->depth() > max_depth) {
      success = false;
      std::cerr << "Error : Maximum depth exceeded" << std::endl;
      break;
    }

    if(curr_node->depth() < max_depth && curr_node->num_points() > max_pts_cell) {
      success = false;
      std::cerr << "Error : Maximum nb of points per cell exceeded" << std::endl;
      break;
    }
    
    leaf_nodes.pop();
  }

  if (!octree.debug_grading()) { // me: remove debug function in octree include, put here
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

      //READS INPUT FILE
      PwnList pwnl;
      if (!get_point_set_from_file(i.string(), pwnl)) {
        std::cout << "Unable to read file" << std::endl;
        std::cout << "Test skipped for this file" << std::endl << std::endl;
        continue;
      }

      //TESTS
      std::vector<size_t> max_nb_pts_vect = {1,2,3,5,8,10,15};
      for (size_t max_depth = 0 ; max_depth <= 11 ; ++max_depth) {
        for (size_t max_nb_pts : max_nb_pts_vect) {
          std::cout << " OCTREE PARAMS : max_depth = " << max_depth
                    << ", max_nb_pts = " << max_nb_pts << std::endl;
          if (!test_octree(std::make_pair(max_depth,max_nb_pts), pwnl)) 
            accumulated_fatal_err = EXIT_FAILURE;
          std::cout << std::endl;
        }
      }

      std::cout << "=========================================================================" << std::endl << std::endl;
    }
  }

  return (accumulated_fatal_err);
}
