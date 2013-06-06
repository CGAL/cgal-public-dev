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
#include <CGAL/Random.h>
#include <CGAL/point_generators_d.h>

#include <time.h>
#include <algorithm>

typedef CGAL::Cartesian_d<double>               Kernel;
typedef Kernel::Point_d                         Point_d;
typedef CGAL::Extreme_points_traits_d<Point_d>  EP_Traits_d;
typedef Kernel::Less_lexicographically_d        Less_lexicographically;

void test(int d, int n) {
  clock_t t1,t2,t3,t4;
  srand (time(NULL));
  CGAL::Extreme_points_options_d op;
  op.set_algorithm(CGAL::EP_CHOOSE_APPROPRIATE);
  op.set_deletion(false);
  std::vector<Point_d> points;
  std::vector<Point_d> extreme_points,extreme_points_ch;
  std::vector<double> point1;
  double dim1,dimm1;
  for (int i=0; i<n; i++) {
    dim1=rand()*1000;
    dimm1=dim1;
    for (int j=0; j<d; j++) {
      point1.push_back(dim1);
      dim1*=dimm1;
    }
    points.push_back(Point_d(d,point1.begin(),point1.end()));
    point1.clear();
  }

  //Extreme_points_d
  t1 = clock();
  CGAL::Extreme_points_d<EP_Traits_d> ep(d,op);
  ep.insert(points.begin(),points.end());
  ep.extreme_points(std::back_inserter(extreme_points));
  t2 = clock();
  t3 = clock();
  CGAL::Convex_hull_d<Kernel> chull(d);
  for (int i=0; i<n; ++i)
    chull.insert(points[i]);
    for (CGAL::Convex_hull_d<Kernel>::Vertex_iterator it=chull.vertices_begin(); it!=chull.vertices_end(); ++it)
      extreme_points_ch.push_back(it->point());
  t4 = clock();

  std::cout << d << " " << n << " " << (t2-t1) << " " << (t4-t3);
  if (t2-t1 > t4-t3) {
    std::cout << " CHull" << std::endl;
  } else {
    std::cout << " Extreme" << std::endl;
  }
}

int main(int argc, char **argv) {
  for (int d=2; d<8; d++) {
    for (int n=100; n<200; n+=100) {
      test(d,n);
    }
  }
    return 0;
}
