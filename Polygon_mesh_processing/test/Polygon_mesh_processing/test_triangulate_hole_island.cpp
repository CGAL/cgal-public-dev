#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole_island.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel  Epic;
typedef Epic::Point_3 Point_3;

template<typename PointRange>
using Domain = CGAL::internal::Domain<PointRange>;


void read_polyline_boundary_and_holes(const char* file_name,
                                      std::vector<Point_3>& points_b,
                                      std::vector<Point_3>& points_h)
{
  std::ifstream stream(file_name);
  if(!stream) {assert(false);}

  for(int i =0; i < 2; ++i) {
    int count;
    if(!(stream >> count)) { assert(false); }
    while(count-- > 0) {
      Point_3 p;
      if(!(stream >> p)) { assert(false); }
      i == 0 ? points_b.push_back(p) : points_h.push_back(p);
    }
  }
}

void read_polyline_one_line(const char* file_name, std::vector<Point_3>& points) {
  std::ifstream stream(file_name);
  if(!stream) { assert(false); }

  int count;
  if(!(stream >> count)) { assert(false); }
  while(count-- > 0) {
    Point_3 p;
    if(!(stream >> p)) { assert(false); }
    points.push_back(p);
  }
}

template <typename PointRange>
void test_split_domain(PointRange& boundary)
{
  std::cout << "test_split_domain" << std::endl;
  // e_D (i, k)
  const int i = 1;
  const int k = 2;
  // trird vertex - on the boundary
  const int v = 4;

  PointRange boundary1;
  Domain<PointRange> D(boundary);
  Domain<PointRange> D1(boundary1);
  Domain<PointRange> D2(boundary1);

  CGAL::internal::split_domain(D, D1, D2, i, v, k);

  std::cout << "left  : \n";
  CGAL::internal::print(D1.boundary);
  std::cout << "right: \n";
  CGAL::internal::print(D2.boundary);

}

template <typename PointRange>
void test_join_domains(PointRange& boundary, PointRange& hole)
{
  std::cout << "--test_join_domains--" << std::endl;
  // e_D (i, k)
  const int i = 1;
  const int k = 2;
  // trird vertex - index of hole vertices
  const int v = 1;

  Domain<PointRange> domain(boundary);
  domain.add_hole(hole);
  PointRange b_vertices;
  Domain<PointRange> new_domain(b_vertices);

  CGAL::internal::join_domain(domain, new_domain, i, v, k);

}

template <typename PointRange>
void test_permutations(PointRange& boundary, PointRange& hole)
{
  std::cout << "--test_permutations--" << std::endl;

  Domain<PointRange> domain(boundary);
  // add 3 holes
  domain.add_hole(hole);
  domain.add_hole(hole);
  domain.add_hole(hole);

  using Phi = CGAL::internal::Phi;
  Phi partitions;
  auto holes = domain.holes;

  CGAL::internal::do_permutations(holes, partitions);

  // all possible combinations: divide the number of holes to 2 sets.
  assert(partitions.size() == pow(2, domain.holes.size()));


}

template <typename PointRange>
void triangulate_hole_island(PointRange& boundary, PointRange& hole)
{
  std::cout << "--test_triangulate_hole_island--" << std::endl;

  Domain<PointRange> domain(boundary);

  domain.add_hole(hole);

  // access edge (1, 2)
  const int i = 1;
  const int k = 2;
  std::size_t count;

  // weight calculator
  /*
  typedef CGAL::internal::Weight_min_max_dihedral_and_area      Weight;
  typedef CGAL::internal::Weight_calculator<Weight,
                CGAL::internal::Is_not_degenerate_triangle>  WC;

  CGAL::internal::processDomain(domain, i, k, count, WC());
  */

  CGAL::internal::processDomain(domain, i, k, count);

  std::cout << "Possible triangles tested: " << count << std::endl;
}




int main()
{

  std::vector<std::string> input_file = {"data/bighole.polylines.txt"};

  const char* file_name = input_file[0].c_str();
  std::vector<Point_3> points_b; // this will contain n and +1 repeated point
  std::vector<Point_3> points_h; // points on the holes

  //read_polyline_boundary_and_holes(file_name, points_b, points_h);

  read_polyline_one_line(file_name, points_b);

  // low level functions
  // test_split_domain(points_b);
  // test_join_domains(points_b, points_h);
  // test_permutations(points_b, points_h);

  triangulate_hole_island(points_b, points_h);










  return 0;
}
