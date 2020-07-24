#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#define CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE
#include <CGAL/Polygon_mesh_processing/repair_manifoldness.h>

typedef CGAL::Simple_cartesian<double>                                        Kernel;
typedef Kernel::Point_3                                                       Point_3;
typedef CGAL::Surface_mesh<Point_3>                                           Surface_mesh;


int main(int argc, char** argv)
{
  std::string mesh_file = "/home/felix/Bureau/Geo_Facto/PSR/tests-repair/jeux-de-test/blobby/blobby_non_manifold.off";
  Surface_mesh sm;
  std::ifstream is(mesh_file);
  CGAL::read_off(is, sm);

  CGAL::Polygon_mesh_processing::treat_non_manifold_vertices(sm);

}