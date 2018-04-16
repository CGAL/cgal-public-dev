#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole_island.h>
#include <CGAL/Polyhedron_3.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;

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

void test_hole_filling_with_islands(const std::string& filename, const bool& use_DT, const bool& correct_orientation)
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

  CGAL::Polyhedron_3<K> mesh;
  CGAL::Polygon_mesh_processing::triangulate_hole_with_islands(b_points, islands, mesh, use_DT, correct_orientation);

  std::ofstream out(filename + ".off");
  out << mesh;
  out.close();
}

void run_unit_tests()
{
  std::vector<std::string> tests =
  {
    // 2D holes
    "data/triangle.polylines.txt",
    "data/quad.polylines.txt",
    "data/hexagon.polylines.txt",
    "data/non_convex.polylines.txt",
    // 2D holes with island
    "data/triangle_island.polylines.txt",
    "data/square_triangle.polylines.txt",
    "data/triangle_quad.polylines.txt",
    "data/quad_in_quad.polylines.txt",
    "data/quad_quad_non_convex.polylines.txt",
    // 3D tests
    "data/triangles_cross.polylines.txt",
    "data/triangles_zaxis.polylines.txt",
    "data/triangles_cross_opposite.polylines.txt",
    // 2 islands
    "data/two_islands_triangles.polylines.txt",
    "data/two_crossing_islands.polylines.txt",
    "data/two_in_quad.polylines.txt",
    "data/pentagon_two_islands.polylines.txt",
    // 3 islands
    "data/three_islands_incorrectly_oriented.polylines.txt",
  };

  std::vector<std::string> tests_correct_orientation =
  {
    // 3 islands
    "data/three_islands_correct_orientation.polylines.txt",
    "data/three_in_hexagon.polylines.txt",
    "data/three_non_convex_heptagon.polylines.txt",
    "data/three_various_islands.polylines.txt",
    // 4 islands
    "data/four_islands.polylines.txt",
    // elephant
    "data/elephant_one_island.cgal",
    "data/elephant_two_islands.cgal",
    "data/elephant_three_islands.cgal"
  };

  bool use_DT = true;
  bool correct_orientation = false;
  for(std::string& filename : tests)
  {
    test_hole_filling_with_islands(filename, use_DT, correct_orientation);
  }

  correct_orientation = true;
  for(std::string& filename : tests_correct_orientation)
  {
    test_hole_filling_with_islands(filename, use_DT, correct_orientation);
  }
}

// parser
char* getCmdOption(char** begin, char** end, const std::string& option)
{
  char** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end)
  {
      return *itr;
  }
  return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
  return std::find(begin, end, option) != end;
}

// example usage:                                          (defaults: with DT and correctly oriented)
// ./test_triangulate_hole_island -all -f filename
// ./test_triangulate_hole_island -all -both -f filename   (all search space and both orientations)
// ./test_triangulate_hole_island                          (runs all unit tests)
// ./test_triangulate_hole_island -f filename              (runs with DT, correct island orientation assumed)

int main(int argc, char* argv[])
{
  bool use_DT = true;
  bool correct_orientation = true;

  if(cmdOptionExists(argv, argv+argc, "-all"))
  {
    use_DT = false;
  }

  if(cmdOptionExists(argv, argv+argc, "-both"))
  {
    correct_orientation = false;
  }

  if(argc > 1)
  {
    //const char* filename = argv[1];
    const char * filename = getCmdOption(argv, argv + argc, "-f");

    std::ifstream input(filename);
    std::vector<std::vector<Point_3>> points;

    if(!input || !load_polylines(input, points))
    {
      std::cerr << "Error loading file.\n";
      return 1;
    }

    std::vector<Point_3> b_points = points[0];
    std::vector<std::vector<Point_3>> islands(points.begin() + 1, points.end());
    std::cout << "Number of islands: " << islands.size() << "\n";
    std::cout << "Triangulating...\n";

    CGAL::Polyhedron_3<K> mesh;
    CGAL::Polygon_mesh_processing::triangulate_hole_with_islands(b_points, islands, mesh, use_DT, correct_orientation);

    std::ofstream out(std::string(filename) + ".off");
    out << mesh;
    out.close();
  }
  else
  {
    run_unit_tests();
  }

  return 0;
}



