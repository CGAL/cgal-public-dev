#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

// Graphs
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Regular_triangulation_2.h>
#include <CGAL/Regular_triangulation_2.h>

#include <CGAL/Polygon_mesh_processing/locate.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Dimension.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/property_map.h>
#include <CGAL/Random.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel                  EPICK;
typedef CGAL::Exact_predicates_exact_constructions_kernel                    EPECK;

namespace PMP = CGAL::Polygon_mesh_processing;

template<typename G>
void test_snappers(const G& g)
{
  typedef typename PMP::internal::Locate_types<G>::FT                        FT;

  typename PMP::internal::Locate_types<G>::Barycentric_coordinates coords = CGAL::make_array(FT(1e-11), FT(0.9999999999999999), FT(1e-10));
  typename PMP::internal::Locate_types<G>::Face_location loc = std::make_pair(*(faces(g).first), coords);

  // ---------------------------------------------------------------------------
  PMP::internal::snap_coordinates_to_border<G>(coords, 1e-10);
  PMP::internal::snap_coordinates_to_border<G>(coords);
  assert(coords[0] == 0. && coords[1] == 1. && coords[2] == 0.);

  // ---------------------------------------------------------------------------
  PMP::internal::snap_location_to_border<G>(loc, 1e-10);
  PMP::internal::snap_location_to_border<G>(loc);
  assert(PMP::is_on_face_border(loc, g));
}

template<typename G>
void test_constructions(const G& g, CGAL::Random& rnd)
{
  typedef typename boost::graph_traits<G>::vertex_descriptor                 vertex_descriptor;
  typedef typename boost::graph_traits<G>::halfedge_descriptor               halfedge_descriptor;
  typedef typename boost::graph_traits<G>::face_descriptor                   face_descriptor;
  typedef typename PMP::internal::Locate_types<G>::descriptor_variant        descriptor_variant;

  typedef typename boost::property_map_value<G, CGAL::vertex_point_t>::type  Point;
  typedef typename CGAL::Kernel_traits<Point>::type                          Kernel;
  typedef typename Kernel::FT                                                FT;

  using boost::choose_param;
  using boost::get_param;

  typedef typename boost::property_map<G, CGAL::vertex_point_t>::const_type  VPM;

  vertex_descriptor v = *(vertices(g).first);
  halfedge_descriptor h = halfedge(v, g);
  face_descriptor f = face(h, g);

  Point p;
  typename PMP::internal::Locate_types<G>::Face_location loc;

  VPM vpm = CGAL::get_const_property_map(boost::vertex_point, g);

  // ---------------------------------------------------------------------------
  // @todo
  PMP::barycentric_coordinates(p, p, p, p, Kernel());
  PMP::barycentric_coordinates(p, p, p, p);

  // ---------------------------------------------------------------------------
  loc = PMP::random_location_on_mesh(g, rnd);
  assert(loc.first != boost::graph_traits<G>::null_face());
  assert(loc.second[0] >= 0.0 && loc.second[0] <= 1.0 &&
         loc.second[1] >= 0.0 && loc.second[1] <= 1.0 &&
         loc.second[2] >= 0.0 && loc.second[2] <= 1.0);

  loc = PMP::random_location_on_face(f, g, rnd);
  assert(loc.first == f);
  assert(loc.second[0] >= 0.0 && loc.second[0] <= 1.0 &&
         loc.second[1] >= 0.0 && loc.second[1] <= 1.0 &&
         loc.second[2] >= 0.0 && loc.second[2] <= 1.0);

  loc = PMP::random_location_on_halfedge(h, g, rnd);
  assert(loc.first == face(h, g));
  assert(loc.second[0] >= 0.0 && loc.second[0] <= 1.0 &&
         loc.second[1] >= 0.0 && loc.second[1] <= 1.0 &&
         loc.second[2] >= 0.0 && loc.second[2] <= 1.0);
  int h_id = PMP::halfedge_index_in_face(h, g);
  assert(loc.second[(h_id+2)%3] == 0.0);

  // ---------------------------------------------------------------------------
  loc = std::make_pair(f, CGAL::make_array(FT(0.3), FT(0.4), FT(0.3)));
  descriptor_variant dv = PMP::get_descriptor_from_location(loc, g);
  assert(bool(face_descriptor* fd = boost::get<face_descriptor*>(&dv)));

//  loc = std::make_pair(f, CGAL::make_array(FT(0.5), FT(0.5), FT(0.0)));
//  dv = PMP::get_descriptor_from_location(loc, g);
//  assert(bool(boost::get<halfedge_descriptor*>(&dv)));

//  loc = std::make_pair(f, CGAL::make_array(FT(1.0), FT(0.0), FT(0.0)));
//  dv = PMP::get_descriptor_from_location(loc, g);
//  assert(bool(boost::get<vertex_descriptor*>(&dv)));
  // ---------------------------------------------------------------------------

  Point q = PMP::location_to_point(loc, g, CGAL::parameters::all_default());
  q = PMP::location_to_point(loc, g);
  assert(q == get(vpm, source(halfedge(f, g), g)));
}

template<typename G>
void test_helpers(const G& g)
{
  typedef typename boost::graph_traits<G>::vertex_descriptor                 vertex_descriptor;
  typedef typename boost::graph_traits<G>::halfedge_descriptor               halfedge_descriptor;
  typedef typename boost::graph_traits<G>::face_descriptor                   face_descriptor;

  typename PMP::internal::Locate_types<G>::Face_location loc;
  vertex_descriptor v;
  halfedge_descriptor h;
  face_descriptor f;
  std::set<face_descriptor> s;

  PMP::internal::incident_faces(loc, g, std::inserter(s, s.begin()));

  PMP::internal::common_halfedge(f, f, g);

  PMP::vertex_index_in_face(v, f, g);
  PMP::halfedge_index_in_face(h, g);
}

template<typename G>
void test_predicates(const G& g)
{
  typedef typename boost::graph_traits<G>::vertex_descriptor                 vertex_descriptor;
  typedef typename boost::graph_traits<G>::halfedge_descriptor               halfedge_descriptor;

  typename PMP::internal::Locate_types<G>::Face_location loc;
  vertex_descriptor v;
  halfedge_descriptor h;

  PMP::is_on_vertex(loc, v, g);
  PMP::is_on_halfedge(loc, h, g);
  PMP::is_in_face(loc, g);
  PMP::is_on_face_border(loc, g);
  PMP::is_on_mesh_border(loc, g);
}

template<typename G>
void test_locate_in_face(const G& g)
{
  typedef typename boost::property_map_value<G, CGAL::vertex_point_t>::type  Point;

  typedef typename boost::graph_traits<G>::vertex_descriptor                 vertex_descriptor;
  typedef typename boost::graph_traits<G>::halfedge_descriptor               halfedge_descriptor;
  typedef typename boost::graph_traits<G>::face_descriptor                   face_descriptor;

  typename PMP::internal::Locate_types<G>::Face_location loc;
  typename PMP::internal::Locate_types<G>::FT a = 0.1;
  Point p;
  vertex_descriptor v;
  halfedge_descriptor h;
  face_descriptor f;

  PMP::locate_in_face(v, g);
  PMP::locate_in_face(v, f, g);
  PMP::locate_in_face(h, a, g);
  PMP::locate_in_face(p, f, g, CGAL::parameters::all_default());
  PMP::locate_in_face(p, f, g);

  PMP::locate_in_adjacent_face(loc, f, g);

  PMP::locate_in_common_face(p, loc, loc, g);
  PMP::locate_in_common_face(p, loc, loc, g, 1e-10);
  PMP::locate_in_common_face(loc, loc, g);
}

template<typename G>
void test_locate_with_AABB_tree(const G& g)
{
  typedef typename boost::property_map_value<G, CGAL::vertex_point_t>::type  Point;
  typedef typename PMP::internal::Ray_type_selector<Point>::type             Ray;

  typedef typename boost::property_map<G, CGAL::vertex_point_t>::const_type  VertexPointMap;
  typedef PMP::internal::Point_to_Point_3_VPM<G, VertexPointMap>             VPM;

  typedef typename CGAL::Kernel_traits<Point>::type                          Kernel;
  typedef typename Kernel::Ray_3                                             Ray_3;

  // ---------------------------------------------------------------------------
  typedef CGAL::AABB_face_graph_triangle_primitive<G, VertexPointMap>        AABB_face_graph_primitive;
  typedef CGAL::AABB_traits<typename PMP::internal::Locate_types<G>::Kernel,
                            AABB_face_graph_primitive>                       AABB_face_graph_traits;

  CGAL::AABB_tree<AABB_face_graph_traits> tree_a;
  typename AABB_face_graph_traits::Point_3 p3_a;
  VertexPointMap vpm_a;
  // ---------------------------------------------------------------------------
  typedef CGAL::AABB_face_graph_triangle_primitive<G, VPM>                   AABB_face_graph_primitive_with_VPM;
  typedef CGAL::AABB_traits<typename PMP::internal::Locate_types<G>::Kernel,
                            AABB_face_graph_primitive_with_VPM>              AABB_face_graph_traits_with_VPM;

  CGAL::AABB_tree<AABB_face_graph_traits_with_VPM> tree_b;
  typename AABB_face_graph_traits_with_VPM::Point_3 p3_b;
  VPM vpm_b;
  // ---------------------------------------------------------------------------

  Point p;
  Ray r;
  Ray_3 r3;

  PMP::build_AABB_tree(g, tree_a);
  PMP::build_AABB_tree(g, tree_b, CGAL::parameters::vertex_point_map(vpm_b));

  PMP::locate_with_AABB_tree(p3_a, tree_a, g);
  PMP::locate_with_AABB_tree(p3_a, tree_a, g, CGAL::parameters::vertex_point_map(vpm_a));
  PMP::locate_with_AABB_tree(p3_b, tree_b, g, CGAL::parameters::vertex_point_map(vpm_b));

  PMP::locate(p, g);
  PMP::locate(p, g, CGAL::parameters::vertex_point_map(vpm_b));

  // ---------------------------------------------------------------------------

  PMP::locate_with_AABB_tree(r3, tree_a, g);
  PMP::locate_with_AABB_tree(r3, tree_b, g, CGAL::parameters::vertex_point_map(vpm_b));

  PMP::locate(r, g);
  PMP::locate(r, g, CGAL::parameters::vertex_point_map(vpm_b));
}

template<typename G>
void test_locate(const G & g, CGAL::Random& rnd)
{
  assert(num_vertices(g) != 0 && num_faces(g) != 0);

  test_snappers(g);
  test_constructions(g, rnd);
  test_helpers(g);
  test_predicates(g);
  test_locate_in_face(g);
  test_locate_with_AABB_tree(g);
}

template<typename K>
void test_polyhedron(const char* fname, CGAL::Random& rnd)
{
  typedef CGAL::Polyhedron_3<K>                               Polyhedron;

  std::cout << "Testing Polyhedron_3 " << fname << "..." << std::flush;
  std::ifstream input(fname);
  Polyhedron poly;
  if (!input || !(input >> poly)){
    std::cerr << "Error: can not read file.";
    return;
  }

  test_locate(poly, rnd);
}

template<typename K>
void test_surface_mesh(const char* fname, CGAL::Random& rnd)
{
  typedef typename K::Point_3                                 Point;
  typedef CGAL::Surface_mesh<Point>                           Mesh;

  std::cout << "Testing Surface_mesh " << fname << "..." << std::flush;
  std::ifstream input(fname);
  Mesh m;
  if (!input || !(input >> m)){
    std::cerr << "Error: can not read file.";
    return;
  }
  
  test_locate(m, rnd);
}

template<typename K>
void test_2D_mesh(const char* fname, CGAL::Random& rnd)
{
  typedef CGAL::Regular_triangulation_2<K>                    RT;
  RT tr;

  std::cout << "Testing regular_triangulation " << fname << "..." << std::flush;
  std::ifstream input(fname);
  CGAL::read_off(input, tr);

  test_locate(tr, rnd);

  // some additionnal tests of locate(), comparing it with tr.locate();
}

int main()
{
  CGAL::Random rnd(CGAL::get_default_random());
  std::cout << "seed: " << rnd.get_seed() << std::endl;

  test_surface_mesh<EPECK>("data/mech-holes-shark.off", rnd);
  test_surface_mesh<EPICK>("data/mech-holes-shark.off", rnd);

  test_polyhedron<EPECK>("data/mech-holes-shark.off", rnd);
  test_polyhedron<EPICK>("data/mech-holes-shark.off", rnd);

//  test_2D_mesh<EPECK>("data/two_tris_collinear.off", rnd);
//  test_2D_mesh<EPICK>("data/two_tris_collinear.off", rnd);

  return 0;
}

