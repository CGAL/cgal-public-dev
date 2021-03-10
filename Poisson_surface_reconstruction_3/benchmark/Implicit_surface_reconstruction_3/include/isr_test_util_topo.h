// rename isr_test_topology_utils.h
#ifndef ISR_TEST_UTIL_TOPO_H
#define ISR_TEST_UTIL_TOPO_H

//includes
#include "isr_test_types.h"
#include <CGAL/Polygon_mesh_processing/connected_components.h>

namespace PMP = CGAL::Polygon_mesh_processing;

//types
typedef Mesh::Halfedge_index halfedge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor          face_descriptor;
typedef boost::graph_traits<Mesh>::faces_size_type          faces_size_type;
typedef Mesh::Property_map<face_descriptor, faces_size_type> FCCmap;

size_t nb_boundaries(const Mesh &mesh)
{
    size_t nb = 0;
    std::set<halfedge_descriptor> he_set;

    BOOST_FOREACH(halfedge_descriptor hd, mesh.halfedges())
    {
      if(mesh.is_border(hd) && (he_set.find(hd) == he_set.end()))
      {
        nb++;
        halfedge_descriptor curr = hd;
        do
        {
          curr = mesh.next(curr);
          he_set.insert(curr);
        }
        while(curr != hd);
      }
    }
    return nb;
}

size_t nb_cc(Mesh &mesh) 
{
  FCCmap fccmap = mesh.add_property_map<face_descriptor, faces_size_type>("f:CC").first;
  faces_size_type nb_con_comp = PMP::connected_components(mesh,fccmap);
  return ( nb_con_comp );
}

size_t compute_genus(Mesh &mesh) 
{
  size_t nb_vertices = mesh.number_of_vertices();
  size_t nb_edges = mesh.number_of_edges();
  size_t nb_faces = mesh.number_of_faces();
  faces_size_type nb_con_comp = nb_cc(mesh);
  size_t nb_bound = nb_boundaries(mesh);
  size_t genus = (nb_edges - nb_faces - nb_bound - nb_vertices + 2*nb_con_comp) / 2; //euler poincare : V - E + F - B = 2 (C - G)
  return ( genus );
}

#endif //ISR_TEST_UTIL_TOPO_H
