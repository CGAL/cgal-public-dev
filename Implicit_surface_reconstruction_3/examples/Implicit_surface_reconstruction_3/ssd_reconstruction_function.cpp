#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/ssd_surface_reconstruction.h>
#include <CGAL/IO/read_points.h>

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

  if(!CGAL::IO::read_points("../data/kitten.xyz", std::back_inserter(points),
                        CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>())
                                         .normal_map(CGAL::Second_of_pair_property_map<Pwn>())))
  {
    std::cerr << "Error: cannot read input file!" << std::endl;
    return EXIT_FAILURE;
  }

  Polyhedron output_mesh_1, output_mesh_2;

  double data_fitting = 10, laplacian=1, hessian=1e-3;
  bool use_octree = true, use_marching_tets = true;

  // Delaunay refinement with C2t3
  if (CGAL::ssd_surface_reconstruction_delaunay(
      points,
      CGAL::First_of_pair_property_map<Pwn>(),
      CGAL::Second_of_pair_property_map<Pwn>(),
      output_mesh_1, 
      data_fitting, 
      laplacian, 
      hessian, 
      !use_octree,
      !use_marching_tets))
    {
        std::ofstream out("kitten_ssd-Delaunay-C2T3-20-100-0.025.off");
        out << output_mesh_1;
    }
  else
    return EXIT_FAILURE;
  
  // octree with marching tets
  if (CGAL::ssd_surface_reconstruction_delaunay(
      points,
      CGAL::First_of_pair_property_map<Pwn>(),
      CGAL::Second_of_pair_property_map<Pwn>(),
      output_mesh_2, 
      data_fitting, 
      laplacian, 
      hessian, 
      use_octree, 
      use_marching_tets))
  {
      std::ofstream out("kitten_ssd-octree-mt.off");
      out << output_mesh_2;
  }
  else
      return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
