#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
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
  std::vector<Point> points;
  std::vector<Vector> normals;


  struct Appender {
    std::vector<Point>& points;
    std::vector<Vector>& normals;

    Appender(std::vector<Point>& points, std::vector<Vector>& normals)
      : points(points), normals(normals)
    {}

    void operator () (const Pwn& pwn)
    {
      points.push_back(pwn.first);
      normals.push_back(pwn.second);
    }

  };

  if(!CGAL::IO::read_points<Pwn>(filename, boost::make_function_output_iterator(Appender(points, normals)),
                                 CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>())
                                                  .normal_map(CGAL::Second_of_pair_property_map<Pwn>())))
  {
    std::cerr << "Error: cannot read input file!" << std::endl;
    return EXIT_FAILURE;
  }

  double average_spacing = CGAL::compute_average_spacing<CGAL::Parallel_if_available_tag>(points, 6);


  if(normals[0] == CGAL::NULL_VECTOR)
  {
    std::cerr << "Estimating normals..." << std::endl;
    std::vector<std::size_t> indices(points.size());
    std::iota(indices.begin(), indices.end(), 0);
    CGAL::jet_estimate_normals<CGAL::Sequential_tag>
      (indices, 6, CGAL::parameters::point_map(CGAL::make_random_access_property_map(points))
                                    .normal_map(CGAL::make_random_access_property_map(normals)));
    CGAL::mst_orient_normals(indices, 6,
                             CGAL::parameters::point_map(CGAL::make_random_access_property_map(points))
                                              .normal_map(CGAL::make_random_access_property_map(normals)));
  }

  std::cout << "First 10 points with normals:" << std::endl;
   for (std::size_t i = 0; i < 10; ++ i)
     std::cout << points[i] << " with normal " << normals[i] << std::endl;
  for (std::size_t i = 0; i < 10; ++ i)
    std::cerr << points[i] << " with normal " << normals[i] << std::endl;


  Polyhedron output_mesh;

  if (CGAL::splat_surface_reconstruction(points, normals,
                                         output_mesh, average_spacing))
    {
        std::string fname =  std::filesystem::path(filename).stem().string() + ".off";
        CGAL::IO::write_polygon_mesh(fname, output_mesh);
    }
  else
    return EXIT_FAILURE;


  return EXIT_SUCCESS;
}
