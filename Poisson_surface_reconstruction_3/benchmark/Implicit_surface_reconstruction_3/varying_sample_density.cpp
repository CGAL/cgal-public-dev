// ----------------------------------------------------------------------------
// Includes
// ----------------------------------------------------------------------------

#include <iostream>
#include <stdlib.h>
#include <array> 

//Mesh
#include <CGAL/Surface_mesh.h>
#include <CGAL/property_map.h>


//file includes
#include "isr_test_types.h"
#include "isr_test_util_reconstruction.h"
#include "isr_test_io_utils.h"
#include "isr_benchmark_dist_utils.h"
#include "isr_test_normal_utils.h"

//boost
#include <boost/foreach.hpp>

#include <CGAL/subdivision_method_3.h>
#include <CGAL/Polygon_mesh_processing/distance.h>

// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------

//Mesh index
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
using namespace CGAL;
namespace params = CGAL::parameters;

// ----------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------

int main(int argc, char* argv[]) //arguments : 1.input mesh file name, 2.output .xyz file name, 3.lvl
{
  if (argc != 4)
    return EXIT_FAILURE;

  //read input mesh
  Mesh new_mesh;
  std::string in_file(argv[1]);
  std::ifstream in(in_file);
  if(in.fail()) {
    std::cerr << "Could not open input file." << in_file << std::endl;
    return EXIT_FAILURE;
  }
  in >> new_mesh;
  in.close();

  //subdivision
  size_t lvl = atoi(argv[3]);
  size_t nb_vert = new_mesh.vertices().size();
  Subdivision_method_3::Loop_subdivision(new_mesh, params::number_of_iterations(lvl/2+1));
  if(new_mesh.is_empty() || !new_mesh.is_valid())
  {
    std::cerr << "Error: mesh is not valid." << std::endl;
    return EXIT_FAILURE;
  }

  //index vertices
  std::unordered_map<int, vertex_descriptor> vmap;
  int i = 0;
  BOOST_FOREACH(vertex_descriptor v, new_mesh.vertices()) {
    vmap[i] = v;
    i++; 
  }

  //remove some vertices randomly (at the end we want the number of points in input multiplied by the lvl)
  Mesh::Property_map<vertex_descriptor, bool> vrm_pm = 
    new_mesh.add_property_map<vertex_descriptor, bool>("v:removed", false).first;
  size_t new_nb_vert = new_mesh.vertices().size();
  size_t curr_v_index = 0;
  size_t b = 0;
  do {
    curr_v_index = CGAL::get_default_random().uniform_int(b,new_nb_vert-1);
    if (!vrm_pm[vmap[curr_v_index]]) {
      vrm_pm[vmap[curr_v_index]] = true;
      i--;
    }
  } while (i != (nb_vert * std::pow(2,lvl)));

  //store pwn list into .xyz out file
  std::string out_file(argv[2]);
  std::ofstream out_xyz_file(out_file);
  Mesh::Property_map<vertex_descriptor, Vector> vnormals_pm = 
    new_mesh.add_property_map<vertex_descriptor, Vector>("v:normals", CGAL::NULL_VECTOR).first;
  compute_area_weighted_vertex_normals(new_mesh, vnormals_pm); 
  BOOST_FOREACH(vertex_descriptor v, new_mesh.vertices()) {
    if (!vrm_pm[v]) {
      const Point& p = new_mesh.point(v);
      Vector n = vnormals_pm[v];
      out_xyz_file << p << " " << n << "\n";
    }
  }
  out_xyz_file.close();

  return EXIT_SUCCESS;
}