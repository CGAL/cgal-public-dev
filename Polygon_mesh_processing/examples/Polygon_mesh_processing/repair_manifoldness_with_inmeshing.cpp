#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/repair_manifoldness.h>

typedef CGAL::Simple_cartesian<double>                                        Kernel;
typedef Kernel::Point_3                                                       Point_3;
typedef CGAL::Surface_mesh<Point_3>                                           Surface_mesh;


int main(int argc, char** argv)
{
  std::string mesh_file = "/home/felix/cgal-public-dev/Polygon_mesh_processing/test/Polygon_mesh_processing/data_repair/nm_vertices_simple.off";
  Surface_mesh sm;
  std::ifstream is(mesh_file);
  CGAL::read_off(is, sm);

  CGAL::Polygon_mesh_processing::treat_non_manifold_vertices(sm);

}