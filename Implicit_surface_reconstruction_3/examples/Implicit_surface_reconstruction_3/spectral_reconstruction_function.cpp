#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/spectral_surface_reconstruction.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/property_map.h>

#include <vector>
#include <fstream>

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef std::pair<Point, Vector> Pwn;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

int main(void)
{
  std::vector<Pwn> points;
  if (!CGAL::IO::read_points(
           "../data/kitten.xyz",
           std::back_inserter(points),
           CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>()).
           normal_map(CGAL::Second_of_pair_property_map<Pwn>())))
    {
      std::cerr << "Error: cannot read file data/kitten.xyz" << std::endl;
      return EXIT_FAILURE;
    }

  Polyhedron output_mesh;
  double data_fitting = 100, laplacian=1, bilaplacian=1;
  bool use_octree = true, use_marching_tets = true;

  if (CGAL::spectral_surface_reconstruction_delaunay_new(
      points,
      CGAL::First_of_pair_property_map<Pwn>(),
      CGAL::Second_of_pair_property_map<Pwn>(),
      output_mesh, 
      data_fitting,
      laplacian,
      bilaplacian,
      !use_octree, 
      !use_marching_tets
      ))
    {
        std::ofstream out("kitten_spectral-Delaunay-C2T3-20-100-0.025.off");
        out << output_mesh;
    }
  else
    return EXIT_FAILURE;

  if (CGAL::spectral_surface_reconstruction_delaunay_new(
      points,
      CGAL::First_of_pair_property_map<Pwn>(),
      CGAL::Second_of_pair_property_map<Pwn>(),
      output_mesh, 
      data_fitting,
      laplacian,
      bilaplacian,
      use_octree, 
      use_marching_tets
  ))
  {
      std::ofstream out("kitten_spectral-octree-mt.off");
      out << output_mesh;
  }
  else
      return EXIT_FAILURE;


  return EXIT_SUCCESS;
}
