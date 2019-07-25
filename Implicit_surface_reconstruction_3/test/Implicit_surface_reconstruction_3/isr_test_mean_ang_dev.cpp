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

int threshold_mult = 10; 


// ----------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------

bool test_mean_ang_dev(const std::string &input_file, const Param &parameter)
{
  Mesh reconstructed_mesh;
  PwnList input_pwn;
  if (!mesh_reconstruction(input_file, parameter,
                input_pwn, reconstructed_mesh)) {
    std::cerr << "Error : Reconstruction failed" << std::endl;
    return false;
  }
  double bbdiag = util_bb_diag(input_pwn);

  double sum = 0.0;
  Mesh::Property_map<face_descriptor, FT> farea_pm = reconstructed_mesh.add_property_map<face_descriptor, FT>("f:area", 0.0).first;
  Mesh::Property_map<face_descriptor, Vector> fnormals_pm = reconstructed_mesh.add_property_map<face_descriptor, Vector>("f:normals", CGAL::NULL_VECTOR).first;
  Mesh::Property_map<vertex_descriptor, Vector> vnormals_pm = reconstructed_mesh.add_property_map<vertex_descriptor, Vector>("v:normals", CGAL::NULL_VECTOR).first;
  
  //computing each face's normal
  CGAL::Polygon_mesh_processing::compute_face_normals(reconstructed_mesh,
        fnormals_pm,
        CGAL::Polygon_mesh_processing::parameters::vertex_point_map(reconstructed_mesh.points()).
        geom_traits(Kernel()));

  //computing each face's area
  BOOST_FOREACH(face_descriptor fd, reconstructed_mesh.faces()) {
    std::vector<Point> fvertices;
    CGAL::Vertex_around_face_iterator<Mesh> vbegin, vend;
    for(boost::tie(vbegin, vend) = vertices_around_face(reconstructed_mesh.halfedge(fd), reconstructed_mesh);
      vbegin != vend; 
      ++vbegin){
        fvertices.insert(fvertices.end(), reconstructed_mesh.point(*vbegin));
    }
    const Triangle t(fvertices[0], fvertices[1], fvertices[2]);
    farea_pm[fd] = CGAL::sqrt(t.squared_area());
  }

  //computing every vertex's normal
  BOOST_FOREACH(vertex_descriptor vd, reconstructed_mesh.vertices()) {
    Vector n = CGAL::NULL_VECTOR;
    halfedge_descriptor hbegin = reconstructed_mesh.halfedge(vd);
    halfedge_descriptor curr_he = hbegin;
    face_descriptor curr_face;

    do
    {
      curr_face = reconstructed_mesh.face(curr_he);
      n += farea_pm[curr_face] * fnormals_pm[curr_face];
      curr_he = reconstructed_mesh.next_around_target(curr_he);
    }
    while (curr_he != hbegin);

    n = n/(CGAL::sqrt(n.squared_length()));
    vnormals_pm[vd] = n;
  }

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
  return( sum/(input_pwn.size()) * threshold_mult < bbdiag ); /*changer ca, aucun sens math*/
}

bool test_mean_ang_dev_all_params(const std::string &input_file)
{
  bool success = true;
  bool curr_par_success;
  Parameters plist;
  for (std::list<Param>::const_iterator param = plist.begin() ; param != plist.end() ; param++) {
    curr_par_success = true;
    std::cout << "///////////" << " " << *param << " "<< "///////////" << std::endl;
    if (!test_mean_ang_dev(input_file, *param)) {
      success = false ;
      curr_par_success = false;
    }
    std::cout << "/////////////////////////// " << (curr_par_success ? "PASSED" : "FAILED") << " ///////////////////////////" << std::endl;
    std::cout << std::endl;
  }
  return (success);
}

int main()
{
  int accumulated_fatal_err = EXIT_SUCCESS ;
  std::cerr << "|-------------------------------------------------------------------------|" << std::endl;
  std::cerr << "|      TEST : MEAN ANGLE DEVIATION BETWEEN INPUT AND OUTPUT NORMALS       |" << std::endl;
  std::cerr << "|-------------------------------------------------------------------------|" << std::endl << std::endl;

  boost::filesystem::path targetDir("./data/regular_data");
  boost::filesystem::recursive_directory_iterator iter(targetDir), eod;

  BOOST_FOREACH(boost::filesystem::path const& i, std::make_pair(iter, eod)) {
    if (is_regular_file(i)) {
      std::cout << "=============== Filename : " << i.string() << " ===============" << std::endl << std::endl;     
      if (!test_mean_ang_dev_all_params(i.string())) 
        accumulated_fatal_err = EXIT_FAILURE;
      std::cout << "=========================================================================" << std::endl << std::endl;
    }
  }

  return (accumulated_fatal_err);
}