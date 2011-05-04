// PolygonApproximation
//
// the program receives 3 arguments:
// [1] the name of file with input polygon points
// [2] the radius of the second polygon (1 by default)
// [3] the number of edges in the second polygon (6 by default)
//
// It computes MinkowskySum of two known polygons,
// the input polygon P, and the n-gon created according
// to the input parameters Br.
//
// Then it "reverse-engineers" the best approximation
// (for some criteria of the "best") of the input polygon P',
// and shows the result in QT application.


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <boost/format.hpp>

#include <CGAL/config.h>

// number types
#include <CGAL/Gmpq.h>

// kernels
#include <CGAL/Cartesian.h>

// polygon
#include <CGAL/Polygon_2.h>

#include <CGAL/minkowski_sum_2.h>

#include <CGAL/Circle_approximation_2.h>

typedef CGAL::Gmpq              NT;
typedef CGAL::Cartesian<NT>          Kernel;
typedef CGAL::Polygon_2<Kernel>        Polygon_2;
typedef Circle_approximation_2<Polygon_2>   Circle_app_2;

typedef Circle_app_2::Polygon_with_holes_2  Polygon_with_holes_2;
typedef Circle_app_2::Point_2        Point_2;
typedef Circle_app_2::Segment_2        Segment_2;


void read_input(const std::string& i_file_name, Polygon_2& o_polygon);
void create_kgon(const NT& i_radius, const int& i_edge_num, Polygon_2& o_polygon);

int main(int argc, char** argv)
{
  // check arguments
  if(argc != 2 && argc != 3 && argc != 4)
  {
    std::cerr
      << "Usage: PolygonApproximation InputPolygonFile "
      << "[Radius [NumberOfEdges]]" << std::endl;
    return 1;
  };

  // the program receives up to 3 arguments:
  // [1] the name of file with input polygon points
  // [2] the radius of the second polygon (0.1 by default)
  // [3] the number of edges in the second polygon (4 by default)
  std::string input_file_name(argv[1]);
  NT radius(1, 1);
  int edge_num = 4;

  if (argc > 2)
  {
    std::istringstream str_radius(argv[2]);
    str_radius >> radius;
    std::clog << "radius: " << radius << std::endl;
  }

  if (argc > 3)
  {
    edge_num = atoi(argv[3]);
    std::clog << "edge_num: " << edge_num << std::endl;
  }

  // Read the input polygon
  Polygon_2 polygon_P;
  read_input(input_file_name, polygon_P);

  // Create the "disk" poligon,
  // i.e. n-gon of given radius
  // with unique slopes
  Circle_app_2 circle_app(polygon_P);

  Polygon_2 polygon_Br;
  //  create_kgon(radius, edge_num, polygon_Br);
  //  circle_app.random_kgon(polygon_Br, edge_num, radius);
  circle_app.regular_kgon(polygon_Br, edge_num, radius);

  // Compute the Minkowski sum.
  Polygon_with_holes_2 polygon_PBr = minkowski_sum_2 (polygon_P, polygon_Br);

  CGAL::set_pretty_mode(std::cout);
  std::cout << "created the polygon P+Br:" << std::endl;
  std::cout << polygon_PBr << std::endl;
  std::cout << std::endl;

  // Reverse-engineer the original polygon
  // find all offset edges and move them back.

  // maintain a set of all slopes from Br polygon
  // to do so.
  circle_app.log_circle_edges(polygon_PBr);

  return 0;
}


void read_input(const std::string& i_file_name, Polygon_2& o_polygon)
{
  o_polygon.clear();

  std::ifstream input_file;
  input_file.open(i_file_name.c_str());

  input_file >> o_polygon;

  input_file.close();

  CGAL::set_pretty_mode(std::cout);
  std::cout << "read the polygon P:" << std::endl;
  std::cout << o_polygon << std::endl;
  std::cout << std::endl;
}

void create_kgon(const NT& i_radius, const int& i_edge_num, Polygon_2& o_polygon)
{
  o_polygon.clear();

  // creating by default square polygon
  Point_2 points[] = {
    Point_2(0.0, i_radius),
    Point_2(-i_radius, 0.0),
    Point_2(0.0, -i_radius),
    Point_2(i_radius, 0.0)
  };

//  o_polygon = Polygon_2(points, points+4);

//  o_polygon.push_back(Point_2(0.0, i_radius));
//  o_polygon.push_back(Point_2(i_radius, 0.0));
//  o_polygon.push_back(Point_2(0.0, -i_radius));
//  o_polygon.push_back(Point_2(-i_radius, 0.0));

  o_polygon.push_back(Point_2(0.0, i_radius));
  o_polygon.push_back(Point_2(-i_radius, 0.0));
  o_polygon.push_back(Point_2(0.0, -i_radius));
  o_polygon.push_back(Point_2(i_radius, 0.0));

  CGAL::set_pretty_mode(std::cout);
  std::cout << "created the polygon Br:" << std::endl;
  std::cout << o_polygon << std::endl;
  std::cout << std::endl;

/*  // check if the polygon is simple.
  std::cout << "The polygon is " <<
    (o_polygon.is_simple() ? "" : "not ") << "simple." << std::endl;

  // check if the polygon is convex
  std::cout << "The polygon is " <<
    (o_polygon.is_convex() ? "" : "not ") << "convex." << std::endl;

  // check if the polygon is clockwise
  std::cout << "The polygon is " <<
    ((o_polygon.orientation() == CGAL::CLOCKWISE) ? "" : "not ") <<
    "clockwise." << std::endl;
*/
}
