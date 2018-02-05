#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole_island.h>
#include <CGAL/Polyhedron_3.h>


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
void test_split_domain(PointRange boundary)
{
  std::cout << "test_split_domain" << std::endl;
  // e_D (i, k)
  const int i = 1;
  const int k = 2;
  // trird vertex - on the boundary
  const int v = 4;

  boundary.pop_back(); // take out last(=first)

  std::vector<int> g_indices(boundary.size());
  for(std::size_t i = 0; i < boundary.size(); ++i)
    g_indices[i] = i;

  Domain<PointRange> D(g_indices);
  Domain<PointRange> D1;
  Domain<PointRange> D2;

  CGAL::internal::split_domain(D, D1, D2, i, v, k);

  assert(D1.b_ids.size() == 4);
  assert(D2.b_ids.size() == 3);

}

template <typename PointRange>
void test_join_domains(PointRange boundary, PointRange hole)
{
  std::cout << "test_join_domains" << std::endl;
  // e_D (i, k)
  const int i = 1;
  const int k = 2;

  boundary.pop_back();
  hole.pop_back();


  std::vector<int> b_indices;
  for(std::size_t i = 0; i < boundary.size(); ++i)
    b_indices.push_back(i);

  std::vector<int> h_ids;
  std::size_t n_b =  b_indices.size();
  for(std::size_t i = b_indices.size(); i < n_b + hole.size(); ++i)
    h_ids.push_back(i);

  // trird vertex - index of hole vertices
  const int v = h_ids.front();


  Domain<PointRange> domain(b_indices); //todo: overload the constructor
  domain.add_hole(h_ids);
  Domain<PointRange> new_domain;

  CGAL::internal::join_domain(domain, new_domain, i, v, k);

  // todo: hardcode points to assert
}


template <typename PointRange>
void test_permutations(PointRange boundary, PointRange hole)
{
  std::cout << "test_permutations" << std::endl;

  boundary.pop_back();
  hole.pop_back();

  Domain<PointRange> domain(boundary); // boundary points not really needed here, just construct the object for now
  // add 3 holes

  std::vector<int> b_indices;
  for(std::size_t i = 0; i < boundary.size(); ++i)
    b_indices.push_back(i);

  std::vector<int> h_ids;
  std::size_t n_b =  b_indices.size();
  for(std::size_t i = n_b; i < n_b + hole.size(); ++i)
    h_ids.push_back(i);

  domain.add_hole(h_ids);
  domain.add_hole(h_ids);
  domain.add_hole(h_ids);

  using Phi = CGAL::internal::Phi;
  Phi partitions;
  auto holes = domain.holes_list;

  CGAL::internal::do_permutations(holes, partitions);

  // all possible combinations: divide the number of holes to 2 sets.
  assert(partitions.size() == pow(2, domain.holes_list.size()));

}



void test_single_triangle(const char* file_name)
{
  std::cout << std::endl << "--- test_single_triangle ---" << std::endl;
  std::vector<Point_3> points_b;
  std::vector<Point_3> points_h;
  read_polyline_one_line(file_name, points_b);

  CGAL::Polyhedron_3<Epic> mesh;

  std::size_t count =
  CGAL::Polygon_mesh_processing::triangulate_hole_islands(points_b, points_h, mesh);

  std::cout << "Possible triangles tested: " << count << std::endl;

  std::ofstream out("data/single_triangle.off");
  out << mesh;
  out.close();

  assert(count == 1);

}

void test_hexagon(const char* file_name)
{
  std::cout << std::endl << "--- test_hole_without_island ---" << std::endl;
  std::vector<Point_3> points_b;
  std::vector<Point_3> points_h;
  read_polyline_one_line(file_name, points_b);

  test_split_domain(points_b);

  CGAL::Polyhedron_3<Epic> mesh;

  std::size_t count =
  CGAL::Polygon_mesh_processing::triangulate_hole_islands(points_b, points_h, mesh);

  std::cout << "Possible triangles tested: " << count << std::endl;

  std::ofstream out("data/hexagon.off");
  out << mesh;
  out.close();

  assert(count == 40); // 40 without memoization, 20 with.
}

void test_quad(const char* file_name)
{
  std::cout << std::endl << "--- test_single_triangle ---" << std::endl;
  std::vector<Point_3> points_b;
  std::vector<Point_3> points_h;
  read_polyline_one_line(file_name, points_b);

  CGAL::Polyhedron_3<Epic> mesh;

  std::size_t count =
  CGAL::Polygon_mesh_processing::triangulate_hole_islands(points_b, points_h, mesh);

  std::cout << "Possible triangles tested: " << count << std::endl;

  std::ofstream out("data/quad.off");
  out << mesh;
  out.close();

  assert(count == 4);
}

void test_triangle_with_triangle_island(const char* file_name)
{
  std::cout << std::endl << "--- test_triangle_with_triangle_island ---" << std::endl;
  std::vector<Point_3> points_b;
  std::vector<Point_3> points_h;
  read_polyline_boundary_and_holes(file_name, points_b, points_h);

  test_permutations(points_b, points_h);
  test_join_domains(points_b, points_h);

  CGAL::Polyhedron_3<Epic> mesh;

  std::size_t count =
  CGAL::Polygon_mesh_processing::triangulate_hole_islands(points_b, points_h, mesh);

  std::cout << "Possible triangles tested: " << count << std::endl;

  std::ofstream out("data/triangle_island.off");
  out << mesh;
  out.close();

  assert(count == 292);
}

int main()
{

  std::vector<std::string> input_file = {"data/triangle.polylines.txt",
                                         "data/bighole.polylines.txt",
                                         "data/triangle-island2.polylines.txt",
                                         "data/quad.polylines.txt"};

  const char* file_name0 = input_file[0].c_str();
  const char* file_name1 = input_file[1].c_str();
  const char* file_name2 = input_file[2].c_str();
  const char* file_name3 = input_file[3].c_str();

  test_single_triangle(file_name0);
  test_quad(file_name3);
  test_hexagon(file_name1);
  test_triangle_with_triangle_island(file_name2);






  return 0;
}
