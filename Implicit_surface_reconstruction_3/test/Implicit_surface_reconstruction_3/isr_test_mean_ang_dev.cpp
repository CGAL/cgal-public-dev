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
#include "include/isr_test_normal_utils.h"

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

//dD_tree surface_mesh
typedef boost::graph_traits<Mesh>::vertex_descriptor dDT_Point;
typedef boost::graph_traits<Mesh>::vertices_size_type size_type;
typedef boost::property_map<Mesh,CGAL::vertex_point_t>::type Vertex_point_pmap;
typedef CGAL::Search_traits_adapter<dDT_Point,Vertex_point_pmap,TreeTraits> SurfaceMeshTreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<SurfaceMeshTreeTraits> SurfaceMeshNeighbor_search;
typedef SurfaceMeshNeighbor_search::Tree SurfaceMeshdD_Tree;
typedef SurfaceMeshdD_Tree::Splitter Splitter;
typedef SurfaceMeshNeighbor_search::Distance Distance;

//Mesh index
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef Mesh::Halfedge_index halfedge_descriptor;
typedef Mesh::Vertex_index Vertex_index;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

double threshold = 2.6; 

// ----------------------------------------------------------------------------

class TestMeanAngDev
{
  public :

  TestMeanAngDev() {} ;

  bool run(const Param &parameter, PwnList &input_pwn)
  {
    Mesh reconstructed_mesh;

    if (!surface_mesh_reconstruction(parameter,
                input_pwn, reconstructed_mesh)) {
    std::cerr << "Error : Reconstruction failed" << std::endl;
    return false;
    }
    double bbdiag = util_bb_diag(input_pwn); //a enlever

    double sum = 0.0;

    //compute each output mesh vertex's normal
    Mesh::Property_map<vertex_descriptor, Vector> vnormals_pm = 
      reconstructed_mesh.add_property_map<vertex_descriptor, Vector>("v:normals", CGAL::NULL_VECTOR).first;
    compute_area_weighted_vertex_normals(reconstructed_mesh, vnormals_pm);

    //putting mesh vertices into dD Tree
    Vertex_point_pmap vppmap = get(CGAL::vertex_point,reconstructed_mesh);
    SurfaceMeshdD_Tree tree(
              vertices(reconstructed_mesh).begin(),
              vertices(reconstructed_mesh).end(),
              Splitter(),
              SurfaceMeshTreeTraits(vppmap)
    );
    Distance tr_dist(vppmap);

    //for each input point, compute closest vertex of the mesh
    for(PwnList::const_iterator it = input_pwn.begin(); it != input_pwn.end(); ++it) {
      SurfaceMeshNeighbor_search search(tree, it->first, 1, 0, true, tr_dist);
      Vertex_index nearest_v = search.begin()->first;
      //compute deviation between input point normal and computed mesh vertex normal
      //add it up to the sum
      Vector in_normal = it->second;
      Vector out_normal = vnormals_pm[nearest_v];
      sum += std::acos( (CGAL::scalar_product(in_normal , out_normal)) / 
                          (CGAL::sqrt(in_normal.squared_length() * out_normal.squared_length())) );
    }

    std::cout << "-> ang_dev = " << sum/(input_pwn.size()) << std::endl;
    return( sum/(input_pwn.size()) < threshold ); /*changer ca, aucun sens math*/
    }
};

int main()
{
  TestMeanAngDev test_mean_ang_dev;
  int accumulated_fatal_err = EXIT_SUCCESS ;
  std::cerr << "|-------------------------------------------------------------------------|" << std::endl;
  std::cerr << "|      TEST : MEAN ANGLE DEVIATION BETWEEN INPUT AND OUTPUT NORMALS       |" << std::endl;
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
      if (!test_all_param(test_mean_ang_dev, pwnl)) 
        accumulated_fatal_err = EXIT_FAILURE;
      std::cout << "=========================================================================" << std::endl << std::endl;
    }
  }

  return (accumulated_fatal_err);
}
