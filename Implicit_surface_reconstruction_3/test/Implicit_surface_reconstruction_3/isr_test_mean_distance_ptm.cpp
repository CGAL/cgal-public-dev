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

int threshold_ratio = 10; 


// ----------------------------------------------------------------------------

class TestMeanDistPTM
{
  public :

  TestMeanDistPTM() {} ;

  bool run(const Param &parameter, PwnList &input_pwn)
  {
    Mesh reconstructed_mesh;

    if (!surface_mesh_reconstruction(parameter,
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

    for (PwnList::const_iterator it = input_pwn.begin(); it != input_pwn.end(); ++it) {
      const Point& current_pt = it->first;
      sqd_dist = tree.squared_distance(current_pt);
      sum += CGAL::sqrt(sqd_dist);
    }
    FT mean_dist = sum / (input_pwn.size());

    std::cout << "-> mean_dist_ptm = " << mean_dist << std::endl;
    return( mean_dist < bbdiag / threshold_ratio);    
  }
};

int main()
{
  TestMeanDistPTM test_mean_dist_ptm;
  int accumulated_fatal_err = EXIT_SUCCESS ;
  std::cerr << "|-------------------------------------------------------------------------|" << std::endl;
  std::cerr << "|            TEST : MEAN DISTANCE FROM INPUT POINTS TO MESH               |" << std::endl;
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
      if (!test_all_param(test_mean_dist_ptm, pwnl)) 
        accumulated_fatal_err = EXIT_FAILURE; 

      std::cout << "=========================================================================" << std::endl << std::endl;
    }
  }

  return (accumulated_fatal_err);
}