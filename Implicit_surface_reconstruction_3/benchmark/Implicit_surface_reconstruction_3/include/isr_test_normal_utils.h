#ifndef ISR_TEST_NORMAL_UTILS_H
#define ISR_TEST_NORMAL_UTILS_H

// ----------------------------------------------------------------------------
// Includes
// ----------------------------------------------------------------------------

#include <iostream>

//Mesh
#include <CGAL/Surface_mesh.h>
#include <CGAL/property_map.h>

//file includes
#include "isr_test_types.h"
#include "isr_test_util_bbox.h"

//boost
#include <boost/foreach.hpp>
#include <boost/property_map/property_map.hpp>

//PMP
#include <CGAL/Polygon_mesh_processing/distance.h>

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


// ----------------------------------------------------------------------------

void compute_area_weighted_vertex_normals(Mesh &mesh, Mesh::Property_map<vertex_descriptor, Vector> &vnormals_pm)
{
  Mesh::Property_map<face_descriptor, FT> farea_pm = mesh.add_property_map<face_descriptor, FT>("f:area", 0.0).first;
  Mesh::Property_map<face_descriptor, Vector> fnormals_pm = mesh.add_property_map<face_descriptor, Vector>("f:normals", CGAL::NULL_VECTOR).first;

  //computing each face's normal
  CGAL::Polygon_mesh_processing::compute_face_normals(mesh,
        fnormals_pm,
        CGAL::Polygon_mesh_processing::parameters::vertex_point_map(mesh.points()).
        geom_traits(Kernel()));

  //computing each face's area
  BOOST_FOREACH(face_descriptor fd, mesh.faces()) {
    std::vector<Point> face_points; //face_points
    CGAL::Vertex_around_face_iterator<Mesh> vbegin, vend;
    for(boost::tie(vbegin, vend) = vertices_around_face(mesh.halfedge(fd), mesh);
      vbegin != vend; 
      ++vbegin) {
        face_points.push_back(mesh.point(*vbegin));
    }
    const Triangle t(face_points[0], face_points[1], face_points[2]);
    farea_pm[fd] = CGAL::sqrt(t.squared_area());
  }

  //computing every vertex's normal
  BOOST_FOREACH(vertex_descriptor vd, mesh.vertices()) {
    Vector n = CGAL::NULL_VECTOR;
    halfedge_descriptor hbegin = mesh.halfedge(vd);
    halfedge_descriptor curr_he = hbegin;
    face_descriptor curr_face;

    do {
      if (!(mesh.is_border(curr_he))) {
        curr_face = mesh.face(curr_he);
        n += farea_pm[curr_face] * fnormals_pm[curr_face];
      }
      curr_he = mesh.opposite(mesh.next(curr_he));
    }
    while (curr_he != hbegin);

    n = n/(CGAL::sqrt(n.squared_length()));
    vnormals_pm[vd] = n;
  }
}

#endif //ISR_TEST_NORMAL_UTILS_H