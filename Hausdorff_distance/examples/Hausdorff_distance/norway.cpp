
#include <CGAL/Simple_cartesian.h>
#include <CGAL/bounding_box.h>

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <set>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point_2;
typedef K::Vector_2 Vector_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;

int main(int argc, char* argv[])
{
  std::set<Point_2> points;

  int n;
  std::ifstream in("norway.edg");
  in >> n;
  Point_2 p;
  while(in >> p){
    points.insert(p);
  }
  
  Iso_rectangle_2 ir = CGAL::bounding_box(points.begin(), points.end());

  double w = ir.xmax() - ir.xmin();
  double h = ir.ymax() - ir.ymin();
  Vector_2 v(w/1000.0, h/1000.0);
  Vector_2 v2(w,h);

  std::ofstream N1("N1.pts");
  for(int i=0; i<10; i++){
    BOOST_FOREACH(Point_2 p, points){
      N1 << p + i*v2 << std::endl;
    }
  }  
  std::ofstream N2("N2.pts");
  for(int i=0; i<10; i++){
    BOOST_FOREACH(Point_2 p, points){
      N2 << p + i*v2 + v << std::endl;
    }
  }
  return 0;
}
