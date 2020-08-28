#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/repair_self_intersections.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>

typedef CGAL::Simple_cartesian<double>                                        Kernel;
typedef Kernel::Point_3                                                       Point_3;
typedef CGAL::Surface_mesh<Point_3>                                           Surface_mesh;


int main(int argc, char** argv)
{
  std::string mesh_file = "/home/felix/Bureau/Geo_Facto/PSR/tests-repair/jeux-de-test/self_intersections/F1_Splitter_2.off";
  Surface_mesh sm;
  std::ifstream is(mesh_file);
  CGAL::read_off(is, sm);

  CGAL::Polygon_mesh_processing::experimental::remove_self_intersections(sm, CGAL::parameters::apply_per_connected_component(true)
  .number_of_iterations(5));

  std::ofstream os ("/home/felix/Bureau/Geo_Facto/PSR/tests-repair/dumps/self_intersections/F1_Splitter_2_repaired.off");
  CGAL::write_off(os, sm);

}