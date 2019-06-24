#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <iostream>
#include <fstream>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <vector>

#if defined(CGAL_LINKED_WITH_TBB)
#define TAG CGAL::Parallel_tag
#else
#define TAG CGAL::Sequential_tag
#endif


typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

double gh_distance(Surface_mesh& tm1, Surface_mesh& tm2) {
  double error = 0;
  long long int count = 0;

  std::vector<Point_3> samples;
  PMP::sample_triangle_mesh(
    tm1,
    std::back_inserter(samples),
    PMP::parameters::number_of_points_per_area_unit(1)
  );
  for(const Point_3& pt: samples) {
    const std::vector<Point_3> vec = { pt };
    double dist = PMP::max_distance_to_triangle_mesh<TAG>(vec, tm2);
    dist = dist * dist;
    error = ((double)count / (count+1)) * error + dist / (count+1);
    count++;
  }

  samples.clear();
  PMP::sample_triangle_mesh(
    tm2,
    std::back_inserter(samples),
    PMP::parameters::number_of_points_per_area_unit(1)
  );
  for(const Point_3& pt: samples) {
    const std::vector<Point_3> vec = { pt };
    double dist = PMP::max_distance_to_triangle_mesh<TAG>(vec, tm1);
    dist = dist * dist;
    error = ((double)count / (count+1)) * error + dist / (count+1);
    count++;
  }

  return error;
}

int main(int argc, char *argv[]) {
  if(argc < 3) {
    std::cerr << "Usage: " << argv[0] << " file1.off file2.off" << std::endl;
  }

  std::ifstream file1(argv[1]);
  Surface_mesh tm1;
  file1 >> tm1;
  if(!CGAL::is_triangle_mesh(tm1)) {
    std::cerr << argv[1] << ": input geometry is not triangulated." << std::endl;
    return 0;
  }

  std::ifstream file2(argv[2]);
  Surface_mesh tm2;
  file2 >> tm2;
  if(!CGAL::is_triangle_mesh(tm2)) {
    std::cerr << argv[2] << ": input geometry is not triangulated." << std::endl;
    return 0;
  }

  file1.close();
  file2.close();

  std::cout << "Symmetric Hausdorff distance:"
            << PMP::approximate_symmetric_Hausdorff_distance<TAG>(
                tm1,
                tm2,
                PMP::parameters::number_of_points_per_area_unit(10),
                PMP::parameters::number_of_points_per_area_unit(10)
              ) << std::endl;

  std::cout << "GH Distance:" << gh_distance(tm1, tm2) << std::endl;
  return 0;
}
