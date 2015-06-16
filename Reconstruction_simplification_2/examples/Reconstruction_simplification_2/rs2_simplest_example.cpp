// Example for Reconstruction_simplification_2, with no mass
// attributes for the input points

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Reconstruction_simplification_2.h>

#include <fstream>
#include <iostream>
#include <string>
#include <iterator>
#include <list>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                                          Point;
typedef K::Segment_2                                        Segment;

typedef K::FT                                               FT;

typedef CGAL::Reconstruction_simplification_2<K>            Rs_2;


void load_xy_file(const std::string& fileName, std::list<Point>& points)
{
  std::ifstream ifs(fileName.c_str());

  Point point;
  while (ifs >> point)
    points.push_back(point);
  
  ifs.close();
}

int main ()
{
  std::list<Point> points;
  load_xy_file("data/stair.xy", points);
  
  Rs_2 rs2(points);

  rs2.run(100); // 100 steps
  
  return 0;
}
