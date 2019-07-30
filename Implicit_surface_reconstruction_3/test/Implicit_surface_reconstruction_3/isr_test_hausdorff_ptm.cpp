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

//AABB_tree
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>


// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------

//AABB Tree
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;

int threshold_mult = 10;


// ----------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------

bool test_hausdorff_ptm(const std::string &input_file, const Param &parameter) 
{
  Mesh reconstructed_mesh;
  PwnList input_pwn;
  if (!mesh_reconstruction(input_file, parameter,
                input_pwn, reconstructed_mesh)) {
  std::cerr << "Error : Reconstruction failed" << std::endl;
    return false;
  }
  double bbdiag = util_bb_diag(input_pwn);

  //charging faces into AABB Tree
  Tree tree(faces(reconstructed_mesh).first, faces(reconstructed_mesh).second, reconstructed_mesh);

  //computation
  tree.accelerate_distance_queries();
  FT sum;
  FT sqd_dist;
  FT max_sqd_dist = tree.squared_distance(input_pwn.begin()->first);

  for (PwnList::const_iterator it = input_pwn.begin(); it != input_pwn.end(); ++it) {
  const Point& current_pt = it->first;
  sqd_dist = tree.squared_distance(current_pt);
  max_sqd_dist = (sqd_dist > max_sqd_dist) ? sqd_dist : max_sqd_dist; // if(sqd_dist > max_sqd_dist) max_sqd_dist = sqd_dist
  }
  FT max_dist = CGAL::sqrt(max_sqd_dist);
  std::cout << "-> haudorff_ptm = " << max_dist << std::endl;
  return( max_dist * threshold_mult < bbdiag );
}

bool test_hausdorff_ptm_all_params(const std::string &input_file)
{
  bool success = true;
  bool curr_par_success;
  Parameters plist;
  for (std::list<Param>::const_iterator param = plist.begin() ; param != plist.end() ; param++) {
    curr_par_success = true;
    std::cout << "///////////" << " " << *param << " "<< "///////////" << std::endl;
    if (!test_hausdorff_ptm(input_file, *param)) {
      success = false ;
      curr_par_success = false;
    }
    std::cout << "/////////////////////////// " << (curr_par_success ? "PASSED" : "FAILED") << " ///////////////////////////" << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
  }
  return (success);
}

int main()
{
  int accumulated_fatal_err = EXIT_SUCCESS ;
  std::cout << "|-------------------------------------------------------------------------|" << std::endl;
  std::cout << "|          TEST : HAUSDORFF DISTANCE FROM INPUT POINTS TO MESH            |" << std::endl;
  std::cout << "|-------------------------------------------------------------------------|" << std::endl << std::endl;

  boost::filesystem::path targetDir("./data/regular_data");
  boost::filesystem::recursive_directory_iterator iter(targetDir), eod;

  BOOST_FOREACH(boost::filesystem::path const& i, std::make_pair(iter, eod)) {
    if (is_regular_file(i)) {
      std::cout << "=============== Filename : " << i.string() << " ===============" << std::endl << std::endl;
      if (!test_hausdorff_ptm_all_params(i.string())) 
        accumulated_fatal_err = EXIT_FAILURE;
      std::cout << "=========================================================================" << std::endl << std::endl;
    }
  }

  return (accumulated_fatal_err);
}
