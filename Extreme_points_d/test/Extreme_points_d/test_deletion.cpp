#include <getopt.h>
#include <vector>
#include <fstream>
#include <string>
#include <iterator>
#include <iostream>
#include <algorithm>
#include <cassert>

#include <CGAL/config.h>
#include <CGAL/basic.h>
#include <CGAL/Extreme_points_d.h>
#include <CGAL/Extreme_points_options_d.h>
#include <CGAL/Extreme_points_traits_d.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian_d.h>

#include <CGAL/Convex_hull_d.h>


#include <CGAL/point_generators_d.h>

#include <time.h>
#include <algorithm>
#include <math.h>

typedef CGAL::Cartesian_d<double>               Kernel;
typedef Kernel::Point_d                         Point_d;
typedef CGAL::Extreme_points_traits_d<Point_d>  EP_Traits_d;
typedef Kernel::Less_lexicographically_d        Less_lexicographically_d;

void test(int d, int n) {
  srand (time(NULL));
  CGAL::Extreme_points_options_d op;
  op.set_algorithm(CGAL::EP_CHOOSE_APPROPRIATE);
  op.set_deletion(true);

  std::set<Point_d, Less_lexicographically_d> points;
  std::vector<Point_d> extreme_points;
  double coord;
  for (int i=0; i<n; i++) {
    std::vector<double> np;
    coord = rand();
    for (int j=0; j<d; j++) {
      np.push_back(pow(coord, j+1));
      //std::cout << np[j] << " ";
    }
    //std::cout << std::endl;
    points.insert(Point_d(d,np.begin(),np.end()));
  }

  CGAL::Extreme_points_d<EP_Traits_d> ep(d,op);
  ep.insert(points.begin(),points.end());
  ep.extreme_points(std::back_inserter(extreme_points));

  assert(extreme_points.size() == points.size());
  for (int i=0; i<extreme_points.size(); i++) {
    assert(points.find(extreme_points[i]) != points.end());
  }

  std::set<Point_d>::iterator it;
  while (points.size() != 10) {
    ep.remove(*points.begin());
    points.erase(points.begin());
    extreme_points.clear();
    ep.extreme_points(std::back_inserter(extreme_points));

    assert(extreme_points.size() == points.size());
    for (int i=0; i<extreme_points.size(); i++) {
      assert(points.find(extreme_points[i]) != points.end());
    }
  }
  while (points.size() != 0) {
    ep.remove(*points.begin());
    points.erase(points.begin());
    extreme_points.clear();
    ep.extreme_points(std::back_inserter(extreme_points));

    assert(extreme_points.size() == points.size());
    for (int i=0; i<extreme_points.size(); i++) {
      assert(points.find(extreme_points[i]) != points.end());
    }
  }

 /* while (points.size() < n) {
    points.push_back(*rnd2);
    rnd2++;
  }

  //Extreme_points_d
  t1 = clock();
  CGAL::Extreme_points_d<EP_Traits_d> ep(d,op);
  ep.insert(points.begin(),points.end());
  ep.extreme_points(std::back_inserter(extreme_points));
  t2 = clock();*/
}

int main(int argc, char **argv) {
  test(5, 100);
  return 0;
}
