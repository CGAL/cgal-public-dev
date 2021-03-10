// ----------------------------------------------------------------------------
// Includes
// ----------------------------------------------------------------------------

#include <iostream>
#include <stdlib.h>
#include <array> 
#include <queue> 

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
typedef Mesh::Halfedge_index halfedge_descriptor;
typedef Mesh::Vertex_index Vertex_index;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

using namespace CGAL;
namespace params = CGAL::parameters;

// ----------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------


int main(int argc, char* argv[]) //arguments : 1.input mesh file name, 2.output .xyz file name, 3.lvl
{
  if (argc != 4)
    return EXIT_FAILURE;

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
  Subdivision_method_3::Loop_subdivision(new_mesh, params::number_of_iterations(2));
  if(new_mesh.is_empty() || !new_mesh.is_valid()) // first is_empty() as it is faster, can use our is_valid()
  {
    std::cerr << "Error: mesh is not valid." << std::endl;
    return EXIT_FAILURE;
  }

  //store pwn list into .xyz out file
  std::string out_file(argv[2]);
  std::ofstream out_xyz_file(out_file);
  Mesh::Property_map<vertex_descriptor, Vector> vnormals_pm = 
    new_mesh.add_property_map<vertex_descriptor, Vector>("v:normals", CGAL::NULL_VECTOR).first;
  compute_area_weighted_vertex_normals(new_mesh, vnormals_pm);  
  size_t new_nb_faces = new_mesh.faces().size();

  //index faces
  std::unordered_map<int, face_descriptor> fmap;
  int j = 0;
  BOOST_FOREACH(face_descriptor f, new_mesh.faces()) {
    fmap[j] = f;
    j++; 
  }

  //generate clusters
  size_t b = 0;
  size_t nb_clusters = 2;
  std::list<vertex_descriptor> cluster_vertices_list;

  std::unordered_map<size_t, vertex_descriptor> cluster_vertices;

  for (size_t i = 0; i < nb_clusters ; ++i) {
    size_t face_rdm_cluster = CGAL::get_default_random().uniform_int(b,new_nb_faces-1);
    face_descriptor curr_face = fmap[face_rdm_cluster];
    CGAL::Vertex_around_face_iterator<Mesh> vbegin, vend;
    size_t j = 0;
    for(boost::tie(vbegin, vend) = vertices_around_face(new_mesh.halfedge(curr_face), new_mesh);
        vbegin != vend; 
        ++vbegin) {
      cluster_vertices[(3*i+j)] = *vbegin;
      ++j;
    }
  }

  //compute distance from cluster for each point (~BFS)
  Mesh::Property_map<vertex_descriptor, size_t> vvs_pm = 
    new_mesh.add_property_map<vertex_descriptor, size_t>("v:visited", 0).first; //0 : not visited, not in queue; 1 : not visited, in queue; 2 : visited
  Mesh::Property_map<vertex_descriptor, size_t> vdist_pm = 
    new_mesh.add_property_map<vertex_descriptor, size_t>("v:distance",  std::numeric_limits<size_t>::max()).first;
  std::queue<vertex_descriptor> v_to_visit;

    //initialize queue
  for (size_t i = 0; i < nb_clusters*3 ; ++i) {
    vertex_descriptor v = cluster_vertices[i];
    v_to_visit.push(v);
    vdist_pm[v] = 0;
    vvs_pm[v] = 1;
  }

    //search
  size_t dmax = 0;
  while (!v_to_visit.empty()) {
    vertex_descriptor v = v_to_visit.front();
    vvs_pm[v] = 2;
    halfedge_descriptor h = new_mesh.halfedge(v);
    vertex_descriptor vbegin = new_mesh.source(h);
    vertex_descriptor curr_v = vbegin;
    do {
      if(vvs_pm[curr_v] == 0) {
        v_to_visit.push(curr_v);
        vvs_pm[curr_v] = 1;
      }
      if(vdist_pm[curr_v] > vdist_pm[v] + 1)
        vdist_pm[curr_v] = vdist_pm[v] + 1;
      if(vdist_pm[curr_v] > dmax)
          dmax = vdist_pm[curr_v];
      h = new_mesh.opposite(new_mesh.next(h));
      curr_v = new_mesh.source(h);
    } while (curr_v != vbegin);
    std::cout << vdist_pm[v] << std::endl;
    v_to_visit.pop();
  }

  //remove points
  Mesh::Property_map<vertex_descriptor, bool> vrm_pm = 
    new_mesh.add_property_map<vertex_descriptor, bool>("v:removed", false).first;

  BOOST_FOREACH(vertex_descriptor v, new_mesh.vertices()) { 
    //compute probability to renive vertex : p(x) = ax^c with c depending on the lvl, a = 1/(dlim^c), and dlim = dmax+1 
    //the goal is to have a fixed point for dlim and a higher growth as lvl increases
    double dlim = dmax + 1.0;
    double c = lvl/10.0 + 0.9 ;
    double a = 1.0/(std::pow(dlim,c));
    double p = a*std::pow(vdist_pm[v],c);
    if (CGAL::get_default_random().uniform_01<double>() < p)
      vrm_pm[v] = true;
  }

  //store result in xyz file
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