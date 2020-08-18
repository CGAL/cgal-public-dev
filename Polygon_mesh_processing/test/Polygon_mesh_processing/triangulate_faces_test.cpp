#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <boost/foreach.hpp>
#include <CGAL/boost/graph/Dual.h>
#include <boost/graph/filtered_graph.hpp>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/centroid.h>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic;
typedef CGAL::Exact_predicates_exact_constructions_kernel Epec;

template <typename K>
bool
test_triangulate_faces()
{
  typedef typename K::Point_3                    Point;
  typedef CGAL::Surface_mesh<Point>          Surface_mesh;

  Surface_mesh mesh;
  std::ifstream input("data/cube_quad.off");

  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << "Not a valid off file." << std::endl;
    return false;
  }
  
  bool success =
  CGAL::Polygon_mesh_processing::triangulate_faces(mesh);
  assert(CGAL::is_triangle_mesh(mesh));

  return success;
}

template <typename K>
bool
test_triangulate_face_range()
{
  typedef typename K::Point_3                    Point;
  typedef CGAL::Surface_mesh<Point>          Surface_mesh;

  Surface_mesh mesh;
  std::ifstream input("data/cube_quad.off");

  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << "Not a valid off file." << std::endl;
    return false;
  }

  bool success =
    CGAL::Polygon_mesh_processing::triangulate_faces(faces(mesh), mesh);
  assert(CGAL::is_triangle_mesh(mesh));

  return success;
}

template <typename K>
bool
test_triangulate_face()
{
  typedef typename K::Point_3                    Point;
  typedef CGAL::Surface_mesh<Point>          Surface_mesh;

  Surface_mesh mesh;
  std::ifstream input("data/cube_quad.off");

  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << "Not a valid off file." << std::endl;
    return false;
  }

  unsigned int nb = 0;
  BOOST_FOREACH(typename boost::graph_traits<Surface_mesh>::face_descriptor fit, faces(mesh))
  {
    if (nb > 4)
      break;
    else if (next(next(halfedge(fit, mesh), mesh), mesh)
             !=   prev(halfedge(fit, mesh), mesh))
    {
      if(CGAL::Polygon_mesh_processing::triangulate_face(fit, mesh))
        ++nb;
      else assert(false);
    }
  }
  return true;
}

template <typename K>
bool
test_triangulate_triangle_face()
{
  typedef typename K::Point_3                    Point;
  typedef CGAL::Surface_mesh<Point>          Surface_mesh;

  Surface_mesh mesh;
  std::ifstream input("data/tetra3.off");

  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << "Not a valid off file." << std::endl;
    return false;
  }

  BOOST_FOREACH(typename boost::graph_traits<Surface_mesh>::face_descriptor fit, faces(mesh))
  {
    if(!CGAL::Polygon_mesh_processing::triangulate_face(fit, mesh))
      assert(false);
  }
  return true;
}

// todo: add this in Dual.h
template <class SurfaceMesh, class Point, class Primal_map>
struct Dual_vpm
{
  typedef typename boost::graph_traits<SurfaceMesh>::face_descriptor key_type;
  typedef Point reference;
  typedef Point value_type;
  typedef boost::readable_property_map_tag category;

  typedef typename boost::graph_traits<SurfaceMesh>::vertex_descriptor vertex_descriptor;

  Dual_vpm(const SurfaceMesh& primal,
           const Primal_map& primal_map)
    : primal_(primal)
    , primal_map_(primal_map) {}

  const Primal_map& primal_map() const
  {
    return primal_map_;
  }

  const SurfaceMesh& primal() const
  {
    return primal_;
  }

  friend
  Point get(Dual_vpm& map, key_type& f)
  {
    std::vector<Point> face_points;

    BOOST_FOREACH(vertex_descriptor v,
                  CGAL::vertices_around_face(halfedge(f, map.primal()), map.primal()))
    {
      face_points.push_back( get(map.primal_map(), v) );
    }

    // temp extra copy
    Point centroid = CGAL::centroid(face_points.begin(), face_points.end(),
                                    CGAL::Dimension_tag<0>());

    return centroid;
  }

  const SurfaceMesh& primal_;
  Primal_map primal_map_;
};

template <typename K>
bool
test_dual_with_various_faces()
{
  typedef typename K::Point_3                    Point;
  typedef CGAL::Surface_mesh<Point>          Surface_mesh;

  Surface_mesh mesh;
  std::ifstream input("data/elephant.off");

  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << "Not a valid off file." << std::endl;
    return false;
  }

  typedef typename boost::property_map<Surface_mesh, boost::vertex_point_t>::type Pmap;
  Pmap vpmap = get_property_map(boost::vertex_point, mesh);

  CGAL::Dual<Surface_mesh> dual(mesh);
  // copy dual to a sm
  Surface_mesh sm_dual;
  CGAL::copy_face_graph(dual, sm_dual,
                        CGAL::Emptyset_iterator(),
                        CGAL::Emptyset_iterator(),
                        CGAL::Emptyset_iterator(),
                        Dual_vpm<Surface_mesh, Point, Pmap>(mesh, vpmap));

  BOOST_FOREACH(typename boost::graph_traits<Surface_mesh>::face_descriptor fit, faces(sm_dual))
  {
    if(!CGAL::Polygon_mesh_processing::triangulate_face(fit, sm_dual))
      assert(false);
  }
  return true;
}





int main()
{
  assert(test_triangulate_faces<Epic>());
  assert(test_triangulate_face_range<Epic>());
  assert(test_triangulate_face<Epic>());
  assert(test_triangulate_triangle_face<Epic>());
  assert(test_dual_with_various_faces<Epic>());

  assert(test_triangulate_faces<Epec>());
  assert(test_triangulate_face_range<Epec>());
  assert(test_triangulate_face<Epec>());
  assert(test_triangulate_triangle_face<Epec>());
  assert(test_dual_with_various_faces<Epec>());

  return EXIT_SUCCESS;
}
