#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
//  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/lion-head.off");
//  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/hand.off");
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/nefertiti.off");
//  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/cube.off");
  std::ifstream input(filename);

  Mesh mesh;
  if (!input || !(input >> mesh) || !CGAL::is_triangle_mesh(mesh)) {
    std::cerr << "Not a valid input file." << std::endl;
    return 1;
  }

  std::cout << "Start remeshing of " << filename
    << " (" << num_faces(mesh) << " faces)..." << std::endl;

  const double tol = 0.001;
  const std::pair<double, double> edge_min_max{0.001, 0.5};
  PMP::Adaptive_sizing_field<Mesh> sizing_field(tol, edge_min_max, mesh);
  unsigned int nb_iter = 5;

  PMP::isotropic_remeshing(
      faces(mesh),
      sizing_field,
      mesh,
      PMP::parameters::number_of_iterations(nb_iter)
      );

  CGAL::IO::write_polygon_mesh("out.off", mesh, CGAL::parameters::stream_precision(17));

  std::cout << "Remeshing done." << std::endl;

  return 0;
}
