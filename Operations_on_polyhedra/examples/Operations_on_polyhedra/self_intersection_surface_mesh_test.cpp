#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/Self_intersection_polyhedron_3.h>
#include <CGAL/Timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;


int main(int, char** argv) {
  std::ifstream input(argv[1]);
  Mesh m;

  if ( !input || !(input >> m) ){
    std::cerr << "Error: can not read file.";
    return 1;
  }
  
  CGAL::Timer timer;
  timer.start();

  std::vector<std::pair<face_descriptor, face_descriptor> > intersected_tris;
  bool intersecting_1 = CGAL::self_intersect<K>(m, back_inserter(intersected_tris)).first;
  assert(intersecting_1 == !intersected_tris.empty());

  std::cerr << "Self-intersection test took " << timer.time() << " sec." << std::endl;
  std::cerr << intersected_tris.size() << " pair of triangles are intersecting." << std::endl;

  timer.reset();
  bool intersecting_2 = CGAL::self_intersect<K>(m);
  assert(intersecting_1 == intersecting_2);

  std::cerr << "Is self-intersection test took " << timer.time() << " sec." << std::endl;
  std::cerr << (intersecting_2 ? "There is a self-intersection." : "There is no self-intersection.") << std::endl;

  return 0;
}
