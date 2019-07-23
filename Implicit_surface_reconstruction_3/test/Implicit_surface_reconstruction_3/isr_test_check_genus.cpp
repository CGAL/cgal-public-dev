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
#include <CGAL/property_map.h>
#include <boost/property_map/property_map.hpp>

namespace PMP = CGAL::Polygon_mesh_processing;

//topo
typedef boost::graph_traits<Mesh>::face_descriptor          face_descriptor;
typedef boost::graph_traits<Mesh>::faces_size_type          faces_size_type;
typedef Mesh::Property_map<face_descriptor, faces_size_type> FCCmap;

typedef Mesh::Halfedge_index halfedge_descriptor;

typedef boost::graph_traits<Mesh>::vertex_descriptor          vertex_descriptor;
typedef boost::graph_traits<Mesh>::vertices_size_type          vertex_size_type;

// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------
#include "include/isr_test_types.h"
#include "include/isr_test_util_nb_boundaries.h"
#include "include/isr_test_util_process_mesh_files.h"

//dD_tree
typedef CGAL::Search_traits_3<Kernel> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree dD_Tree;

int threshold = 10; /*changer le nom*/

// ----------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------

size_t compute_genus(Mesh &mesh) 
{
  size_t nb_vertices = mesh.number_of_vertices();
  size_t nb_edges = mesh.number_of_edges();
  size_t nb_faces = mesh.number_of_faces();
  FCCmap fccmap = mesh.add_property_map<face_descriptor, faces_size_type>("f:CC").first;
  faces_size_type nb_con_comp = PMP::connected_components(mesh,fccmap);
  size_t nb_bound = nb_boundaries(mesh);
  size_t out_gen = (nb_edges - nb_faces - nb_bound - nb_vertices + 2*nb_con_comp) / 2; //euler poincare
  return ( out_gen );
}

size_t test_check_genus_param(std::string input_file, const Param &parameter) 
{
  Mesh reconstructed_mesh;
  PwnList input_pwn;

  if (!reconstruction_param(reconstructed_mesh, input_pwn,
                parameter, input_file)) {
    std::cerr << "reconstruction failed" << std::endl;
    return (-1);
  }

  return (compute_genus(reconstructed_mesh));
}

bool test_param(std::string input_file, const size_t &in_gen)
{
  bool success = true;
  Parameters plist;
  for (std::list<Param>::const_iterator param = plist.begin() ; param != plist.end() ; param++) {
    std::cout << *param << std::endl;
    std::cout << "in_gen = "<< in_gen << std::endl;
    if (test_check_genus_param(input_file, *param) != in_gen)
      success = false ;
    std::cout << (success ? "Passed" : "Failed") << std::endl ;
    std::cout << std::endl;
  }
  return (success);
}

int main()
{
  bool found_fail = false;
  std::cerr << "Test : Check if genus is preserved" << std::endl << std::endl;

  boost::filesystem::path targetDir("./data/genus");
  boost::filesystem::recursive_directory_iterator iter(targetDir), eod;

  BOOST_FOREACH(boost::filesystem::path const& i, std::make_pair(iter, eod)) {
    if (is_regular_file(i) && ((i.string()).find("big_data") == std::string::npos)) {

      std::cout << "Filename : " << i.string() << std::endl;

      if(is_mesh(i.string())) //compute genus
      {
        Mesh input_m;
        if(!read_input_mesh_file(i.string(), input_m))
          return EXIT_FAILURE; /*a modifier*/
        size_t in_gen = compute_genus(input_m);
        if (!test_param(i.string(), in_gen))
            found_fail = true;
      }
      else //get genus in file name if possible
      {
        std::string delimiter = "genus_" ;
        std::string str_gen = i.string();
        size_t pos = str_gen.find(delimiter);
        
        if (pos == std::string::npos) {
          std::cerr << "Impossible to find input genus" << std::endl;
          std::cerr << "Test skipped for this file" << std::endl;
          continue;
        }
        else {
          str_gen.erase(0, pos + delimiter.length());
          str_gen = str_gen.substr(0,str_gen.find("."));
          size_t in_gen = std::stoi(str_gen);
          if (!test_param(i.string(), in_gen))
            found_fail = true;
        }
      }
      std::cout << std::endl << std::endl;
    }
  }

  int accumulated_fatal_err = found_fail ? EXIT_FAILURE : EXIT_SUCCESS ;
  return (accumulated_fatal_err);
}