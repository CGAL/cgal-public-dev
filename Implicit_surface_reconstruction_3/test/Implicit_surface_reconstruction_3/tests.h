#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Timer.h>
#include <CGAL/trace.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Implicit_reconstruction_function.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Polyhedron_VRML_1_ostream.h>

#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include "boost/program_options.hpp"

#include <CGAL/disable_warnings.h>

//include tests
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <list>
#include <cstdlib>
#include <fstream>
#include <math.h>

//Mesh
#include <CGAL/Surface_mesh.h>

//AABB_tree
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

//Hausdorff & PMP
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>

//dD Tree
#include <CGAL/point_generators_2.h>

//topological features
#include <boost/function_output_iterator.hpp>
#include <boost/property_map/property_map.hpp>

#include <CGAL/bounding_box.h>

namespace PMP = CGAL::Polygon_mesh_processing;



// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------

// kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

// Simple geometric types
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Surface_mesh<Point>                           Mesh;
typedef Mesh::Vertex_index                           Vertex_index;
typedef std::pair<Point, Vector> Point_with_normal;
typedef Kernel::Sphere_3 Sphere;
typedef Kernel::Triangle_3 Triangle;
typedef std::list<Point_with_normal> PointList;
typedef CGAL::First_of_pair_property_map<Point_with_normal> Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;

// Spectral implicit function
typedef CGAL::Implicit_reconstruction_function<Kernel, PointList, Normal_map> Implicit_reconstruction_function;

//typedef CGAL::Implicit_reconstruction_function<Kernel, PointList, CGAL::Identity_property_map<Point_with_normal> > Implicit_reconstruction_function;

// Surface mesher
typedef CGAL::Surface_mesh_default_triangulation_3 STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
typedef CGAL::Implicit_surface_3<Kernel, Implicit_reconstruction_function> Surface_3;


//includes rajoutes par Roxane

typedef Kernel::Segment_3 Segment;
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;

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

//topo
typedef boost::graph_traits<Mesh>::face_descriptor          face_descriptor;
typedef boost::graph_traits<Mesh>::faces_size_type          faces_size_type;
typedef Mesh::Property_map<face_descriptor, faces_size_type> FCCmap;

typedef Mesh::Halfedge_index halfedge_descriptor;

typedef boost::graph_traits<Mesh>::vertex_descriptor          vertex_descriptor;
typedef boost::graph_traits<Mesh>::vertices_size_type          vertex_size_type;

size_t nb_boundaries(const Mesh &m)
{
    size_t nb = 0;
    std::set<halfedge_descriptor> he_set;

    BOOST_FOREACH(halfedge_descriptor hd, m.halfedges())
    {
      if(m.is_border(hd) && (he_set.find(hd) == he_set.end()))
      {
        nb++;
        halfedge_descriptor curr = hd;
        do
        {
          curr = m.next(curr);
          he_set.insert(curr);
        }
        while(curr != hd);
      }
    }
    return nb;
}

std::pair<FT,FT> mean_dist_pt_mesh(const Mesh &m, const PointList &points) 
{
  //charging faces into AABB Tree
	Tree tree(faces(m).first, faces(m).second, m);

  //computation
  tree.accelerate_distance_queries();
  FT sum;
  FT max_sqd_dist = tree.squared_distance(points.begin()->first);
  FT sqd_dist;
  for (PointList::const_iterator it = points.begin(); it != points.end(); ++it) {
    const Point& current_pt = it->first;
    sqd_dist = tree.squared_distance(current_pt);
    sum += CGAL::sqrt(sqd_dist);
    max_sqd_dist = (sqd_dist > max_sqd_dist) ? sqd_dist : max_sqd_dist;
  }
  return( std::make_pair( (sum/(points.size())) , CGAL::sqrt(max_sqd_dist) ) );
}

FT mean_dist_mesh_pt(const Mesh &m, PointList &points) /*idem pour const*/
{
  //sampling mesh
  std::list<Point> sample_points;
  CGAL::Polygon_mesh_processing::sample_triangle_mesh(m,
                                                      std::back_inserter(sample_points),
                                                      4000);

  //putting input points into dD_Tree
  typedef typename PointList::value_type PointList_t;
  boost::function<Point(PointList_t&)> pwn_it_to_point_it = boost::bind(&PointList_t::first, _1);
  dD_Tree tree(boost::make_transform_iterator(points.begin(), pwn_it_to_point_it), 
                boost::make_transform_iterator(points.end(), pwn_it_to_point_it));

  //computation
  FT sum = 0;
  for (std::list<Point>::iterator it = sample_points.begin(); it != sample_points.end(); ++it) {
    Neighbor_search search(tree, *it, 1);
    sum += CGAL::sqrt(search.begin()->second);
  }
  return(sum/(sample_points.size()));
}

FT hausdorff_dist_mesh_pt(const Mesh &m, PointList &points) /*idem pour const*/
{
  typedef typename PointList::value_type PointList_t;
  boost::function<Point(PointList_t&)> pwn_it_to_point_it = boost::bind(&PointList_t::first, _1);
  double max_dist = PMP::approximate_max_distance_to_point_set(
                                            m,
                                            CGAL::make_range( boost::make_transform_iterator(points.begin(), pwn_it_to_point_it),
                                                               boost::make_transform_iterator(points.end(), pwn_it_to_point_it)) ,
                                            4000 );
  return (max_dist);
}

std::tuple<size_t,size_t,size_t,faces_size_type,size_t,size_t> topo(Mesh &m)
{
	size_t nb_vertices = m.number_of_vertices();
	size_t nb_edges = m.number_of_edges();
  size_t nb_faces = m.number_of_faces();
  FCCmap fccmap = m.add_property_map<face_descriptor, faces_size_type>("f:CC").first;
  faces_size_type nb_con_comp = PMP::connected_components(m,fccmap);
  size_t nb_bound = nb_boundaries(m);
  size_t genus = (nb_edges - nb_faces - nb_bound - nb_vertices + 2*nb_con_comp) / 2; //euler poincare

	return(std::make_tuple(nb_vertices , nb_edges, nb_faces, nb_con_comp, nb_bound, genus));
}

double mean_angle_dev( Mesh &m, const PointList &points) /*squared area -> sqrt*/
{
  double sum = 0.0;
  Mesh::Property_map<face_descriptor, FT> farea_pm = m.add_property_map<face_descriptor, FT>("f:area", 0.0).first;
  Mesh::Property_map<face_descriptor, Vector> fnormals_pm = m.add_property_map<face_descriptor, Vector>("f:normals", CGAL::NULL_VECTOR).first;
  Mesh::Property_map<vertex_descriptor, Vector> vnormals_pm = m.add_property_map<vertex_descriptor, Vector>("v:normals", CGAL::NULL_VECTOR).first;
  
  //computing each face's normal
  CGAL::Polygon_mesh_processing::compute_face_normals(m,
        fnormals_pm,
        CGAL::Polygon_mesh_processing::parameters::vertex_point_map(m.points()).
        geom_traits(Kernel()));

  //computing each face's area
  BOOST_FOREACH(face_descriptor fd, m.faces()) {
    std::vector<Point> fvertices;
    CGAL::Vertex_around_face_iterator<Mesh> vbegin, vend;
    for(boost::tie(vbegin, vend) = vertices_around_face(m.halfedge(fd), m);
      vbegin != vend; 
      ++vbegin){
        fvertices.insert(fvertices.end(), m.point(*vbegin));
    }
    const Triangle t(fvertices[0], fvertices[1], fvertices[2]);
    farea_pm[fd] = CGAL::sqrt(t.squared_area());
  }

  //computing every vertex's normal
  BOOST_FOREACH(vertex_descriptor vd, m.vertices()) {
    Vector n = CGAL::NULL_VECTOR;
    halfedge_descriptor hbegin = m.halfedge(vd);
    halfedge_descriptor curr_he = hbegin;
    face_descriptor curr_face;

    do
    {
      curr_face = m.face(curr_he); /*peut etre verifier si la face n'est pas nulle (bord)?*/
      n += farea_pm[curr_face] * fnormals_pm[curr_face];
      curr_he = m.next_around_target(curr_he);
    }
    while (curr_he != hbegin);

    n = n/(CGAL::sqrt(n.squared_length()));
    vnormals_pm[vd] = n;
  }

  //putting mesh vertices into dD Tree
  Vertex_point_pmap vppmap = get(CGAL::vertex_point,m);
  SurfaceMeshdD_Tree tree(
            vertices(m).begin(),
            vertices(m).end(),
            Splitter(),
            SurfaceMeshTreeTraits(vppmap)
  );
  Distance tr_dist(vppmap);

  //for each input point, compute closest vertex of the mesh
  for(PointList::const_iterator it = points.begin(); it != points.end(); ++it) {
    SurfaceMeshNeighbor_search search(tree, it->first, 1, 0, true, tr_dist);
    Vertex_index nearest_v = search.begin()->first;
    //compute deviation between input point normal and computed mesh vertex normal
    //add it up to the sum
    Vector in_normal = it->second;
    Vector out_normal = vnormals_pm[nearest_v];
    sum += std::acos( (CGAL::scalar_product(in_normal , out_normal)) / 
                        (CGAL::sqrt(in_normal.squared_length() * out_normal.squared_length())) );
  }

  return(sum/(points.size()));
}

double util_bb_diag(PointList pwnl)
{
  typedef typename PointList::value_type PointList_t;
  boost::function<Point(PointList_t&)> pwn_it_to_point_it = boost::bind(&PointList_t::first, _1);
  Kernel::Iso_cuboid_3 c3 = CGAL::bounding_box(boost::make_transform_iterator(pwnl.begin(), pwn_it_to_point_it), 
                                               boost::make_transform_iterator(pwnl.end(), pwn_it_to_point_it));
  double d_squared = (c3[7][0] - c3[0][0]) * (c3[7][0] - c3[0][0])
                    +(c3[7][1] - c3[0][1]) * (c3[7][1] - c3[0][1])
                    +(c3[7][2] - c3[0][2]) * (c3[7][2] - c3[0][2]) ;

  return(CGAL::sqrt(d_squared));
}

void run_tests(std::string file_input, PointList input_points) 
{
  //charging file into mesh
	Mesh reconstructed_mesh;
  std::ifstream in(file_input);
  in >> reconstructed_mesh;

  if(!in || !reconstructed_mesh.is_valid() || reconstructed_mesh.is_empty())
  {
    std::cerr << "Error: cannot read the ouput file containing the reconstructed mesh. " << file_input << std::endl;
    return;
  }

  //geometry
  std::pair<FT,FT> dist_pair = mean_dist_pt_mesh(reconstructed_mesh, input_points);
  FT l2_ptm = dist_pair.first;
  FT l2_mtp = mean_dist_mesh_pt(reconstructed_mesh, input_points);
  FT hausdorff_ptm = dist_pair.second;
  FT hausdorff_mtp = hausdorff_dist_mesh_pt(reconstructed_mesh, input_points);
  double mad = mean_angle_dev(reconstructed_mesh, input_points);

  //topology
  std::tuple<size_t,size_t,size_t,faces_size_type,size_t,size_t> topo_ft = topo(reconstructed_mesh);
  size_t v =  std::get<0>(topo_ft);
  size_t e =  std::get<1>(topo_ft);
  size_t f =  std::get<2>(topo_ft);
  size_t cc =  std::get<3>(topo_ft);
  size_t b =  std::get<4>(topo_ft);
  size_t g =  std::get<5>(topo_ft);

  double bb_diag = util_bb_diag(input_points);
  //display
  std::cerr << std::endl << "--------  TESTS  --------" << std::endl;

  std::cerr << std::endl << "1. Geometry" << std::endl;
  std::cerr << "  1.1. Mean distance" << std::endl;
  std::cerr << "    points -> mesh : d_ptm = " << l2_ptm << std::endl;
  std::cerr << "    mesh -> points : d_mtp = " << l2_mtp << std::endl;
  std::cerr << "  1.2. Hausdorff distance" << std::endl;
  std::cerr << "    points -> mesh : h_ptm = " << hausdorff_ptm << std::endl;
  std::cerr << "    mesh -> points : h_mtp = " << hausdorff_mtp << std::endl;
  std::cerr << "  1.3. Mean angle deviation between normals" << std::endl;
  std::cerr << "    theta = " << mad << std::endl;

  std::cerr << std::endl << "2. Topology" << std::endl;
/*  std::cerr << "    v = " << v << std::endl;
  std::cerr << "    e = " << e << std::endl;
  std::cerr << "    f = " << f << std::endl;*/
  std::cerr << "    nb of connected components = " << cc << std::endl;
  std::cerr << "    nb of nb_boundaries = " << b << std::endl;
  std::cerr << "    genus = " << g << std::endl;
  std::cerr << std::endl << "-------------------------" << std::endl;

  std::cerr << std::endl ;
}