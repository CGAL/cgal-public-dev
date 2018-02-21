#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole_island.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

/// extra code necessary for cgal's hole filling

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef typename K::Point_3 Point_3;
typedef CGAL::Polyhedron_3<K> Polyhedron;

template<class HDS, class K>
class Polyhedron_builder : public CGAL::Modifier_base<HDS> {
  typedef typename K::Point_3 Point_3;
public:
  Polyhedron_builder(std::vector<boost::tuple<int, int, int> >* triangles,
    std::vector<Point_3>* polyline)
    : triangles(triangles), polyline(polyline)
  { }

  void operator()(HDS& hds) {
    CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
    B.begin_surface(polyline->size() -1, triangles->size());

    for(typename std::vector<Point_3>::iterator it = polyline->begin();
      it != --polyline->end(); ++it) {
        B.add_vertex(*it);
    }

    for(typename std::vector<boost::tuple<int, int, int> >::iterator it = triangles->begin();
      it != triangles->end(); ++it) {
        B.begin_facet();
        B.add_vertex_to_facet(it->get<0>());
        B.add_vertex_to_facet(it->get<1>());
        B.add_vertex_to_facet(it->get<2>());
        B.end_facet();
    }

    B.end_surface();
  }

private:
  std::vector<boost::tuple<int, int, int> >* triangles;
  std::vector<Point_3>* polyline;
};

void check_triangles(std::vector<Point_3>& points, std::vector<boost::tuple<int, int, int> >& tris) {
  if(points.size() - 3 != tris.size()) {
    std::cerr << "  Error: there should be n-2 triangles in generated patch." << std::endl;
    assert(false);
  }

  const int max_index = static_cast<int>(points.size())-1;
  for(std::vector<boost::tuple<int, int, int> >::iterator it = tris.begin(); it != tris.end(); ++it) {
    if(it->get<0>() == it->get<1>() ||
      it->get<0>() == it->get<2>() ||
      it->get<1>() == it->get<2>() )
    {
      std::cerr << "Error: indices of triangles should be all different." << std::endl;
      assert(false);
    }

    if(it->get<0>() >= max_index ||
      it->get<1>() >= max_index ||
      it->get<2>() >= max_index )
    {
      std::cerr << "  Error: max possible index check failed." << std::endl;
      assert(false);
    }
  }
}

void check_constructed_polyhedron(const char* file_name,
  std::vector<boost::tuple<int, int, int> >* triangles,
  std::vector<Point_3>* polyline,
  const bool save_poly)
{
  Polyhedron poly;
  Polyhedron_builder<typename Polyhedron::HalfedgeDS,K> patch_builder(triangles, polyline);
  poly.delegate(patch_builder);

  if(!poly.is_valid()) {
    std::cerr << "  Error: constructed patch does not constitute a valid polyhedron." << std::endl;
    assert(false);
  }

  if (!save_poly)
    return;

  std::string out_file_name;
  out_file_name.append(file_name).append(".off");
  std::ofstream out(out_file_name.c_str());
  out << poly; out.close();
}

////////////////

typedef CGAL::Exact_predicates_inexact_constructions_kernel  Epic;
//typedef Epic::Point_3 Point_3;

template<typename PointRange>
using Domain = CGAL::internal::Domain<PointRange>;

void read_boundary_and_holes(const std::string& file_name,
                             std::vector<Point_3>& points_b,
                             std::vector<std::vector<Point_3>>& points_h,
                             const std::size_t& number_of_islands)
{
  std::ifstream stream(file_name);
  if(!stream) {assert(false);}

  // import boundary
  int count;
  if(!(stream >> count)) { assert(false); }
  while(count-- > 0) {
    Point_3 p;
    if(!(stream >> p)) { assert(false); }
    points_b.push_back(p);
  }

  // import islands
  points_h.resize(number_of_islands);
  for(std::size_t i=0; i < number_of_islands; ++i)
  {
    int count;
    if(!(stream >> count)) { assert(false); }
    while(count-- > 0) {
      Point_3 p;
      if(!(stream >> p)) { assert(false); }
      points_h[i].push_back(p); // fix this
    }

  }

}



void read_polyline_boundary_and_holes(const std::string& file_name,
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

void read_polyline_one_line(const std::string& file_name, std::vector<Point_3>& points) {
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

  boundary.pop_back(); // take out last(=first)

  // e_D (i, k)
  const int i = 0;
  const int k = boundary.size() - 1;
  // trird vertex - on the boundary
  //const int v = 4;
  std::vector<int> inds = {4};
  typename std::vector<int>::iterator v_it = inds.begin();


  std::vector<int> g_indices(boundary.size());
  for(std::size_t i = 0; i < boundary.size(); ++i)
    g_indices[i] = i;

  Domain<PointRange> D(g_indices);
  Domain<PointRange> D1;
  Domain<PointRange> D2;

  CGAL::internal::split_domain_case_2(D, D1, D2, i, v_it, k);

  assert(D1.b_ids.size() == 5);
  assert(D2.b_ids.size() == 2);

}

template <typename PointRange>
void test_join_domains(PointRange boundary, PointRange hole)
{
  std::cout << "test_join_domains" << std::endl;

  boundary.pop_back();
  hole.pop_back();

  // e_D (i, k)
  const int i = 0;
  const int k = boundary.size() - 1;


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
  Domain<PointRange> new_domain1;
  Domain<PointRange> new_domain2;

  CGAL::internal::split_domain_case_1(domain, new_domain1, new_domain2, i, v, k);

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
  auto islands = domain.islands_list;

  CGAL::internal::do_permutations(islands, partitions);

  // all possible combinations: divide the number of holes to 2 sets.
  assert(partitions.size() == pow(2, domain.islands_list.size()));
}

void test_single_triangle(const std::string& file_name)
{
  std::cout << std::endl << "--- test_single_triangle ---" << std::endl;
  std::vector<Point_3> points_b;
  std::vector<Point_3> points_h;
  read_polyline_one_line(file_name, points_b);

  CGAL::Polyhedron_3<Epic> mesh;

  std::size_t count =
  CGAL::Polygon_mesh_processing::triangulate_hole_islands(points_b, points_h, mesh);

  std::cout << "Possible triangles tested: " << count << std::endl;

  std::ofstream out(file_name + ".off");
  out << mesh;
  out.close();

  assert(count == 1);
}

void test_hexagon(const std::string& file_name)
{
  std::cout << std::endl << "--- test_hexagon ---" << std::endl;
  std::vector<Point_3> points_b;
  std::vector<Point_3> points_h;
  read_polyline_one_line(file_name, points_b);

  //test_split_domain(points_b); // todo: test this or get rid of.

  CGAL::Polyhedron_3<Epic> mesh;

  std::size_t count =
  CGAL::Polygon_mesh_processing::triangulate_hole_islands(points_b, points_h, mesh);

  std::cout << "Possible triangles tested: " << count << std::endl;
  assert(count == 40); // 40 without memoization, 20 with.

  std::ofstream out(file_name + ".off");
  out << mesh;
  out.close();
}

void test_quad(const std::string& file_name)
{
  std::cout << std::endl << "--- test_quad ---" << std::endl;
  std::vector<Point_3> points_b;
  std::vector<Point_3> points_h;
  read_polyline_one_line(file_name, points_b);

  CGAL::Polyhedron_3<Epic> mesh;

  std::size_t count =
  CGAL::Polygon_mesh_processing::triangulate_hole_islands(points_b, points_h, mesh);

  std::cout << "Possible triangles tested: " << count << std::endl;
  assert(count == 4);

  std::ofstream out(file_name + ".off");
  out << mesh;
  out.close();
}

void test_triangle_with_triangle_island(const std::string& file_name)
{
  std::cout << std::endl << "--- test_triangle_with_triangle_island ---" << std::endl;
  std::vector<Point_3> points_b;
  std::vector<Point_3> points_h;
  read_polyline_boundary_and_holes(file_name, points_b, points_h);

  test_permutations(points_b, points_h);
  //test_join_domains(points_b, points_h);

  CGAL::Polyhedron_3<Epic> mesh;

  std::size_t count =
  CGAL::Polygon_mesh_processing::triangulate_hole_islands(points_b, points_h, mesh);

  std::cout << "Possible triangles tested: " << count << std::endl;

  std::ofstream out(file_name + ".off");
  out << mesh;
  out.close();
}


void test_square_triangle(const std::string& file_name)
{
  std::cout << std::endl << "--- test_square_triangle ---" << std::endl;
  std::vector<Point_3> points_b;
  std::vector<Point_3> points_h;
  read_polyline_boundary_and_holes(file_name, points_b, points_h);

  CGAL::Polyhedron_3<Epic> mesh;

  std::size_t count =
  CGAL::Polygon_mesh_processing::triangulate_hole_islands(points_b, points_h, mesh);

  std::cout << "Possible triangles tested: " << count << std::endl;

  std::ofstream out(file_name + ".off");
  out << mesh;
  out.close();
}


void test_non_convex(const std::string& file_name)
{
  std::cout << std::endl << "--- test_non_convex ---" << std::endl;
  std::vector<Point_3> points_b;
  std::vector<Point_3> points_h;
  read_polyline_one_line(file_name, points_b);

  CGAL::Polyhedron_3<Epic> mesh;

  std::size_t count =
  CGAL::Polygon_mesh_processing::triangulate_hole_islands(points_b, points_h, mesh);

  std::cout << "Possible triangles tested: " << count << std::endl;

  std::ofstream out(file_name + ".off");
  out << mesh;
  out.close();
}

void test_triangle_quad(const std::string& file_name)
{
  std::cout << std::endl << "--- test_triangle_quad ---" << std::endl;
  std::vector<Point_3> points_b;
  std::vector<Point_3> points_h;
  read_polyline_boundary_and_holes(file_name, points_b, points_h);

  CGAL::Polyhedron_3<Epic> mesh;

  std::size_t count =
  CGAL::Polygon_mesh_processing::triangulate_hole_islands(points_b, points_h, mesh);

  std::cout << "Possible triangles tested: " << count << std::endl;

  std::ofstream out(file_name + ".off");
  out << mesh;
  out.close();
}

void test_quad_in_quad(const std::string& file_name)
{
  std::cout << std::endl << "--- test_quad_in_quad ---" << std::endl;
  std::vector<Point_3> points_b;
  std::vector<Point_3> points_h;
  read_polyline_boundary_and_holes(file_name, points_b, points_h);

  CGAL::Polyhedron_3<Epic> mesh;

  std::size_t count =
  CGAL::Polygon_mesh_processing::triangulate_hole_islands(points_b, points_h, mesh);

  std::cout << "Possible triangles tested: " << count << std::endl;
  // assert

  std::ofstream out(file_name + ".off");
  out << mesh;
  out.close();
}

void test_quad_quad_non_convex(const std::string& file_name)
{
  std::cout << std::endl << "--- test_quad_quad_non_convex ---" << std::endl;
  std::vector<Point_3> points_b;
  std::vector<Point_3> points_h;
  read_polyline_boundary_and_holes(file_name, points_b, points_h);

  CGAL::Polyhedron_3<Epic> mesh;

  std::size_t count =
  CGAL::Polygon_mesh_processing::triangulate_hole_islands(points_b, points_h, mesh);

  std::cout << "Possible triangles tested: " << count << std::endl;

  std::ofstream out(file_name + ".off");
  out << mesh;
  out.close();
}


void test_non_convex_non_convex(const std::string& file_name)
{
  std::cout << std::endl << "--- test_non_convex_non_convex ---" << std::endl;
  std::vector<Point_3> points_b;
  std::vector<Point_3> points_h;
  read_polyline_boundary_and_holes(file_name, points_b, points_h);

  CGAL::Polyhedron_3<Epic> mesh;

  std::size_t count =
  CGAL::Polygon_mesh_processing::triangulate_hole_islands(points_b, points_h, mesh);

  std::cout << "Possible triangles tested: " << count << std::endl;
  // assert count

  std::ofstream out(file_name + ".off");
  out << mesh;
  out.close();
}

void test_triangles_zaxis(const std::string& file_name)
{
  std::cout << std::endl << "--- test_triangles_zaxis ---" << std::endl;
  std::vector<Point_3> points_b;
  std::vector<Point_3> points_h;
  read_polyline_boundary_and_holes(file_name, points_b, points_h);

  CGAL::Polyhedron_3<Epic> mesh;

  std::size_t count =
  CGAL::Polygon_mesh_processing::triangulate_hole_islands(points_b, points_h, mesh);

  std::cout << "Possible triangles tested: " << count << std::endl;
  assert(count == 630);

  std::ofstream out(file_name + ".off");
  out << mesh;
  out.close();
}

void test_triangles_planes_cross(const std::string& file_name)
{
  std::cout << std::endl << "--- test_triangles_zaxis ---" << std::endl;
  std::vector<Point_3> points_b;
  std::vector<Point_3> points_h;
  read_polyline_boundary_and_holes(file_name, points_b, points_h);

  CGAL::Polyhedron_3<Epic> mesh;

  std::size_t count =
  CGAL::Polygon_mesh_processing::triangulate_hole_islands(points_b, points_h, mesh);

  std::cout << "Possible triangles tested: " << count << std::endl;
  //assert(count == 630);

  std::ofstream out(file_name + ".off");
  out << mesh;
  out.close();
}

void test_triangles_planes_cross_opposite(const std::string& file_name)
{
  std::cout << std::endl << "--- test_triangles_zaxis_opposite ---" << std::endl;
  std::vector<Point_3> points_b;
  std::vector<Point_3> points_h;
  read_polyline_boundary_and_holes(file_name, points_b, points_h);

  CGAL::Polyhedron_3<Epic> mesh;

  std::size_t count =
  CGAL::Polygon_mesh_processing::triangulate_hole_islands(points_b, points_h, mesh);

  std::cout << "Possible triangles tested: " << count << std::endl;
  //assert(count == 630);

  //std::ofstream out("data/triangles_cross_opposite.off");
  std::ofstream out(file_name + ".off");
  out << mesh;
  out.close();
}

void test_both_algorithms(const std::string& file_name)
{
  std::cout << std::endl << "--- test_both_algorithms ---" << std::endl;
  // cgal's hole filling

  std::vector<Point_3> points; // this will contain n and +1 repeated point
  read_polyline_one_line(file_name, points);

  std::vector<boost::tuple<int, int, int> > tris;
  CGAL::Polygon_mesh_processing::triangulate_hole_polyline(
    points, std::back_inserter(tris),
    CGAL::Polygon_mesh_processing::parameters::use_delaunay_triangulation(false));

  check_triangles(points, tris);
  const char* char_ptr = file_name.c_str();
  check_constructed_polyhedron(char_ptr, &tris, &points, true);
  std::cerr << "  Done!" << std::endl;


  // recursive algorithm
  std::vector<Point_3> points_b;
  std::vector<Point_3> points_h;
  read_polyline_one_line(file_name, points_b);

  CGAL::Polyhedron_3<Epic> mesh;

  std::size_t count =
  CGAL::Polygon_mesh_processing::triangulate_hole_islands(points_b, points_h, mesh);

  std::cout << "Possible triangles tested: " << count << std::endl;

  // to optimize the output file name
  std::ofstream out(file_name + "-recursive" + ".off");
  out << mesh;
  out.close();
}

void test_two_triangle_islands(const std::string& file_name)
{
  std::cout << std::endl << "--- test_two_islands_triangle ---" << std::endl;
  std::vector<Point_3> points_b;
  std::vector<std::vector<Point_3>> points_h;
  const std::size_t number_of_islands = 2;
  read_boundary_and_holes(file_name, points_b, points_h, number_of_islands);

  CGAL::Polyhedron_3<Epic> mesh;

  std::size_t count =
  CGAL::Polygon_mesh_processing::triangulate_hole_islands(points_b, points_h, mesh);

  std::cout << "Possible triangles tested: " << count << std::endl;
  //assert(count == 630);

  //std::ofstream out("data/triangles_cross_opposite.off");
  std::ofstream out(file_name + ".off");
  out << mesh;
  out.close();
}




template<class Point>
bool load_polylines(std::ifstream& input,
                    std::vector<std::vector<Point>>& points)
{
  int counter = 0;
  std::size_t n;
  while(input >> n) {
    ++counter;
    std::vector<Point> new_polyline;
    points.push_back(new_polyline);
    std::vector<Point>&polyline = points.back();
    polyline.reserve(n);
    while(n--){
      Point p;
      input >> p;
      polyline.push_back(p);
      if(!input.good()) return 0;
    }
    std::string line_remainder;
    std::getline(input, line_remainder);

    if(input.bad() || input.fail()) return 0;
    }

  return 1;
}

void test_hole_filling_with_islands(const std::string& filename)
{
  std::cout << "\n--- testing " + filename + " ---\n";

  std::ifstream input(filename);
  std::vector<std::vector<Point_3>> points;

  if(!input || !load_polylines(input, points))
  {
    std::cerr << "Error loading file.\n";
    return;
  }

  std::vector<Point_3> b_points = points[0];
  std::vector<std::vector<Point_3>> islands(points.begin() + 1, points.end());
  std::cout << "Number of islands: " << islands.size() << std::endl;

  CGAL::Polyhedron_3<Epic> mesh;
  std::size_t count =
  CGAL::Polygon_mesh_processing::triangulate_hole_islands(b_points, islands, mesh);
  std::cout << "Possible triangles tested: " << count << std::endl;

  std::ofstream out(filename + ".off");
  out << mesh;
  out.close();
}


void run_unit_tests()
{
  std::vector<std::string> tests =
  {
    // 2D holes
    /*"data/triangle.polylines.txt",
    "data/quad.polylines.txt",
    "data/hexagon.polylines.txt",
    "data/non-convex.polylines.txt",
    "data/hexagon.polylines.txt",
    // 2D holes with island
    "data/triangle-island.polylines.txt",
    "data/square_triangle.polylines.txt",
    "data/triangle_quad.polylines.txt",
    "data/quad_in_quad.polylines.txt",
    "data/quad_quad_non_convex.polylines.txt",*/
    "data/triangles_cross.polylines.txt",
    // 3D tests
    /*"data/triangles-zaxis.polylines.txt",
    "data/triangles_cross.polylines.txt",
    "data/triangles_cross_opposite.polylines.txt"
    // 2 islands
    "data/two_islands_triangles.polylines.txt"*/
  };

  for(std::string& filename : tests)
  {
    test_hole_filling_with_islands(filename);
  }

}


int main(int argc, char* argv[])
{
  /*

  const char* filename =
      (argc > 1) ? argv[1] : "data/two_islands_triangles.polylines.txt";

  std::ifstream input(filename);
  std::vector<std::vector<Point_3>> points;

  if(!input || !load_polylines(input, points))
  {
    std::cerr << "Error loading file.\n";
    return 1;
  }

  std::vector<Point_3> b_points = points[0];
  std::vector<std::vector<Point_3>> islands(points.begin() + 1, points.end());
  std::cout << "Number of islands: " << islands.size() << std::endl;

  CGAL::Polyhedron_3<Epic> mesh;
  std::size_t count =
  CGAL::Polygon_mesh_processing::triangulate_hole_islands(b_points, islands, mesh);
  std::cout << "Possible triangles tested: " << count << std::endl;

  std::ofstream out(std::string(filename) + ".off");
  out << mesh;
  out.close();

  */


  run_unit_tests();


  return 0;
}



