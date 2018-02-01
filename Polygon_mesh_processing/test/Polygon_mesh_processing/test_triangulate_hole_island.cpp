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
void test_split_domain(PointRange boundary)
{
  std::cout << "--test_split_domain--" << std::endl;
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
  std::cout << "--test_join_domains--" << std::endl;
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
  std::cout << "--test_permutations--" << std::endl;

  boundary.pop_back();
  hole.pop_back();

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
void triangulate_hole_island_i(PointRange boundary, PointRange hole)
{
  std::cout << "--triangulate_hole_island(boundary, hole)--" << std::endl;

  // remove the last(=first) stupid point. It's only a source of evil.
  boundary.pop_back();
  hole.pop_back();


  std::vector<int> b_indices;
  for(std::size_t i = 0; i < boundary.size(); ++i)
    b_indices.push_back(i);

  std::vector<int> h_ids;
  std::size_t n_b =  b_indices.size();
  for(std::size_t i = n_b; i < n_b + hole.size(); ++i)
    h_ids.push_back(i);

  Domain<PointRange> domain(b_indices);
  domain.add_hole(h_ids); // to assert if hole is empty - to add an empty hole maybe

  // access edge (1, 2)
  const int i = 1;
  const int k = 2;
  std::size_t count = 0;

  // weight calculator - temp
  typedef CGAL::internal::Weight_min_max_dihedral_and_area      Weight;
  typedef CGAL::internal::Weight_calculator<Weight,
                CGAL::internal::Is_not_degenerate_triangle>  WC;

  // look up maps - to make them template
  typedef std::pair<int, int> Edge;
  typedef WC::Weight Weight;
  std::map<Edge, Weight> w_map;
  std::map<Edge, int> lambda_map;


  CGAL::internal::processDomain(domain, i, k, count, w_map, lambda_map);

  std::cout << "Possible triangles tested: " << count << std::endl;

  assert(count == 105);
}

template <typename PointRange>
void triangulate_hole_island_i(PointRange boundary)
{
  std::cout << "--triangulate_hole_island(boundary)--" << std::endl;

  // remove the last(=first) stupid point. It's only a source of evil.
  boundary.pop_back();

  std::vector<int> b_indices(boundary.size());
  for(std::size_t i = 0; i < boundary.size(); ++i)
    b_indices[i] = i;

  Domain<PointRange> domain(b_indices);

  // access edge (1, 2)
  const int i = 1;
  const int k = 2;
  std::size_t count = 0;

  // weight calculator
  typedef CGAL::internal::Weight_min_max_dihedral_and_area      Weight;
  typedef CGAL::internal::Weight_calculator<Weight,
                CGAL::internal::Is_not_degenerate_triangle>  WC;

  // look up maps - to make them template
  typedef std::pair<int, int> Edge;
  typedef WC::Weight Weight;
  std::map<Edge, Weight> w_map;
  std::map<Edge, int> lambda_map;


  CGAL::internal::processDomain(domain, i, k, count, w_map, lambda_map);

  std::cout << "Possible triangles tested: " << count << std::endl;

  assert(count == 18);
}


void test_hole_without_island(const char* file_name)
{
  std::cout << std::endl << "--- test_hole_without_island ---" << std::endl;
  std::vector<Point_3> points_b; // this will contain n and +1 repeated point
  std::vector<Point_3> points_h; // points on the holes (contains +1)

  read_polyline_one_line(file_name, points_b);

  test_split_domain(points_b);
  triangulate_hole_island_i(points_b);
}


void test_triangle_with_triangle_island(const char* file_name)
{
  std::cout << std::endl << "--- test_triangle_with_triangle_island ---" << std::endl;
  std::vector<Point_3> points_b; // this will contain n and +1 repeated point
  std::vector<Point_3> points_h; // points on the holes (contains +1)

  read_polyline_boundary_and_holes(file_name, points_b, points_h);

  test_permutations(points_b, points_h);
  test_join_domains(points_b, points_h);
  triangulate_hole_island_i(points_b, points_h);
}

int main()
{

  std::vector<std::string> input_file = {"data/bighole.polylines.txt",
                                         "data/triangle-island2.polylines.txt"};

  const char* file_name1 = input_file[0].c_str();
  const char* file_name2 = input_file[1].c_str();

  test_hole_without_island(file_name1);
  test_triangle_with_triangle_island(file_name2);






  return 0;
}
