#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/splat_surface_reconstruction.h>

#include <vector>
#include <fstream>
#include <string>
#include <filesystem>

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef std::pair<Point, Vector> Pwn;
typedef CGAL::Surface_mesh<Point> Polyhedron;

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("points_3/kitten.xyz");
  std::vector<Pwn> points;

  if(!CGAL::IO::read_points(filename, std::back_inserter(points),
                            CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>())
                                             .normal_map(CGAL::Second_of_pair_property_map<Pwn>())))
  {
    std::cerr << "Error: cannot read input file!" << std::endl;
    return EXIT_FAILURE;
  }

  Polyhedron output_mesh;

  double average_spacing = CGAL::compute_average_spacing<CGAL::Parallel_if_available_tag>
    (points, 6, CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>()));

  if (CGAL::splat_surface_reconstruction
      (points.begin(), points.end(),
       CGAL::First_of_pair_property_map<Pwn>(),
       CGAL::Second_of_pair_property_map<Pwn>(),
       output_mesh, average_spacing))
    {
        std::string fname =  std::filesystem::path(filename).stem().string() + ".off";
        CGAL::IO::write_polygon_mesh(fname, output_mesh);
    }
  else
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
