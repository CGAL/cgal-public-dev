
#include <fstream>
#include <iostream>

#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/Real_timer.h>
#include <CGAL/point_generators_3.h>

#include <boost/core/demangle.hpp>
#include <random>
#include <CGAL/Surface_mesh/Surface_mesh.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Timer.h>

static std::size_t C = 25;
static std::size_t T = 100000;


template<typename Surface_mesh>
CGAL::Bbox_3 bbox(typename boost::graph_traits<Surface_mesh>::face_descriptor fd,
                  const Surface_mesh &p) {
  typename boost::graph_traits<Surface_mesh>::halfedge_descriptor hd = halfedge(fd, p);
  CGAL::Bbox_3 res = p.point(source(hd, p)).bbox();
  res += p.point(target(hd, p)).bbox();
  res += p.point(target(next(hd, p), p)).bbox();
  return res;
}

template<typename K>
std::vector<CGAL::Ray_3<K>> generate_queries(std::size_t n) {
  typedef CGAL::Point_3<K> Point_3;
  typedef CGAL::Ray_3<K> Ray_3;

  // Generate some points
  CGAL::Random r(23);
  CGAL::Random_points_in_cube_3<Point_3, CGAL::Creator_uniform_3<typename K::FT, Point_3> > g(2.0, r);
  std::vector<Point_3> points;
  points.reserve(n * 2);
  std::copy_n(g, n * 2, std::back_inserter(points));

  // Combine those points into Rays
  std::vector<Ray_3> queries;
  queries.reserve(n);
  for (std::size_t i = 0, j = points.size() - 1; i < j; ++i, --j)
    queries.push_back(Ray_3(points[i], points[j]));

  return queries;
}

template<typename Traits, typename Surface_mesh, typename Queries>
double benchmark_traversal(Surface_mesh data, Queries queries) {
  typedef CGAL::AABB_tree<Traits> Tree;

//  typedef typename boost::graph_traits<Surface_mesh>::face_descriptor face_descriptor;
//  typedef typename Surface_mesh::template Property_map<face_descriptor, CGAL::Bbox_3> Bbox_pmap;
//
//  Bbox_pmap bbox_pmap = data.template add_property_map<face_descriptor, CGAL::Bbox_3>("f:bbox", CGAL::Bbox_3()).first;
//  for (face_descriptor fd : faces(data)) {
//    put(bbox_pmap, fd, bbox(fd, data));
//  }
//
//  Traits traits(bbox_pmap);
//  Tree tree(traits);
//  tree.insert(faces(data).first, faces(data).second, data);
//  tree.build();
//  data.remove_property_map(bbox_pmap);

  Tree tree(faces(data).first, faces(data).second, data);
  tree.build();

  typedef typename Tree::AABB_traits::template Intersection_and_primitive_id<decltype(queries[0])>::Type Result_type;
  std::vector<Result_type> v;
  v.reserve(queries.size());

  CGAL::Real_timer t;
  t.start();
  for (const auto &query : queries)
    tree.template all_intersections(query, std::back_inserter(v));
  t.stop();

  return t.time() / queries.size();
}

template<typename Traits, typename Surface_mesh>
double benchmark_construction(Surface_mesh data) {
  typedef CGAL::AABB_tree<Traits> Tree;

  typedef typename boost::graph_traits<Surface_mesh>::face_descriptor face_descriptor;
  typedef typename Surface_mesh::template Property_map<face_descriptor, CGAL::Bbox_3> Bbox_pmap;

  CGAL::Real_timer t;
  t.start();
  for (int i = 0; i < C; ++i) {

//    Bbox_pmap bbox_pmap = data.template add_property_map<face_descriptor, CGAL::Bbox_3>("f:bbox", CGAL::Bbox_3()).first;
//    for (face_descriptor fd : faces(data)) {
//      put(bbox_pmap, fd, bbox(fd, data));
//    }
//
//    Traits traits(bbox_pmap);
//    Tree tree(traits);
//    tree.insert(faces(data).first, faces(data).second, data);
//    tree.build();
//    data.remove_property_map(bbox_pmap);

    Tree tree(faces(data).first, faces(data).second, data);
    tree.build();

  }
  t.stop();

  return t.time() / C;
}

template<typename K>
void benchmark(std::string input_path) {

  typedef typename K::FT FT;
  typedef typename K::Point_3 Point_3;
  typedef CGAL::Bbox_3 Bbox_3;
  typedef CGAL::Surface_mesh<Point_3> Surface_mesh;
  typedef CGAL::Polyhedron_3<K> Polyhedron_3;
  typedef CGAL::Timer Timer;

  typedef CGAL::AABB_face_graph_triangle_primitive<Surface_mesh> Primitive;

  typedef typename boost::graph_traits<Surface_mesh>::face_descriptor face_descriptor;
  typedef typename Surface_mesh::template Property_map<face_descriptor, Bbox_3> Bbox_pmap;
  typedef CGAL::AABB_face_graph_triangle_primitive<Surface_mesh> Primitive;
  typedef CGAL::AABB_traits<K, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> Tree;

  std::ifstream in(input_path);
  Surface_mesh data;
  in >> data;

  auto queries = generate_queries<K>(T);

  std::cout << "{| class=\"wikitable\"\n";
  std::cout << "! Performance !! {version}"
            << "\n|-" << std::endl;
  std::cout << "| Construction Time || " << std::flush
            << benchmark_construction<Traits>(data) << " s "
            << "\n|-" << std::endl;
  std::cout << "| Traversal Time || " << std::flush
            << benchmark_traversal<Traits>(data, queries) << " s "
            << "\n|-" << std::endl;
  std::cout << "|}\n" << std::endl;
}

int main(int argc, char **argv) {

  // Determine our data source, with a default if no path is provided
  std::string input_path = argc > 1 ? argv[1] : "data/handle.off";

//  benchmark<CGAL::Simple_cartesian<float>>(input_path);
//  benchmark<CGAL::Simple_cartesian<double>>(input_path);
  benchmark<CGAL::Exact_predicates_inexact_constructions_kernel>(input_path);

}
