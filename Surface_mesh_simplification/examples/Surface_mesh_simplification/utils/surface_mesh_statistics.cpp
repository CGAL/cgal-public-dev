#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <chrono>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;

int main(int argc, char** argv)
{
  Surface_mesh surface_mesh;

  std::ifstream is(argv[1]);
  is >> surface_mesh;
  if(!CGAL::is_triangle_mesh(surface_mesh)){
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << argv[1] << " "
            << surface_mesh.number_of_edges() << " edges.\n";


  return EXIT_SUCCESS;
}
