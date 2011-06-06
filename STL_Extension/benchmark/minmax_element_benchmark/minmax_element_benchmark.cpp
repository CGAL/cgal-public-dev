#include <cstdio>
#include <utility>
#include <iterator>
#include <algorithm>
#include <vector>
#include <list>
#include <iostream>
#include <cassert>

#include <boost/timer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/minmax_element.hpp>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/algorithm.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_selection.h>

int format(const char* lib, const char* container, int elements, double time) {
  return std::printf(" %s || %s || %d || %.4f M items/sec \n", lib, container, elements, time);
}

template<typename Iterator>
void test(Iterator begin, Iterator end, const char* container, int repeats) {
  int n = std::distance(begin, end);

  boost::timer t;
  std::pair<Iterator,Iterator> res1;
  std::pair<Iterator,Iterator> res2;
  std::pair<Iterator,Iterator> res3;
  double time;

  //CGAL::min_max_element
  t.restart();
  for (int i = 0; i < repeats; ++i) {
    res1 = CGAL::min_max_element(begin, end);
  }
  time = (double)n*repeats/t.elapsed()/1.0E6;
  format("CGAL", container, n, time);
  
  //boost::minmax_element
  t.restart();
  for (int i = 0; i < repeats; ++i) {
    res2 = boost::minmax_element(begin, end);
  }
  time = (double)n*repeats/t.elapsed()/1.0E6;
  format("boost", container, n, time);

  //std::minmax_element
#ifndef CGAL_CFG_NO_CPP0X_MINMAX_ELEMENT
  t.restart();
  for (int i = 0; i < repeats; ++i) {
    res3 = std::minmax_element(begin, end);
  }
  time = (double)n*repeats/t.elapsed()/1.0E6;
  format("stdlib", container, n, time);
#endif

  if(!(*(res1.first) == *(res2.first) && *(res1.second) == *(res2.second)
                                      && 
       *(res1.first) == *(res3.first) && *(res1.second) == *(res3.second)
                                      &&
       *(res2.first) == *(res3.first) && *(res2.second) == *(res3.second)))
    std::cerr << "Different results";
}

typedef CGAL::Simple_cartesian<double>         R;
typedef R::Point_2                             Point;
typedef CGAL::Creator_uniform_2<double,Point>  Creator;

int main(int argc, char** argv)
{
  int n = 20000000;
  int repeats = 10;

  if(argc > 1)
    n = boost::lexical_cast<int>(argv[1]);

  if(argc > 2)
    repeats = boost::lexical_cast<int>(argv[2]);

  //generate n random points
  std::vector<Point> points;
  points.reserve(n);

  CGAL::Random_points_in_disc_2<Point,Creator> g( 1000.0);
  CGAL::copy_n( g, n, std::back_inserter(points));

  std::vector<Point> vector(points.begin(), points.end());
  test(vector.begin(), vector.end(), "vector", repeats);
  vector.clear();

  std::list<Point> list(points.begin(), points.end());
  test(list.begin(), list.end(), "list", repeats);

  return 0;
}
