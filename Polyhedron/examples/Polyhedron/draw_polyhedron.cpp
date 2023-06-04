#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/draw_polyhedron.h>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel2;
typedef CGAL::Polyhedron_3<Kernel2>                       Polyhedron2;

int main(int argc, char* argv[])
{
  Polyhedron2 P;
  std::ifstream in1((argc>1)?argv[1]:CGAL::data_file_path("meshes/cross_quad.off"));
  in1 >> P;
  CGAL::draw(P);

  return EXIT_SUCCESS;
}
