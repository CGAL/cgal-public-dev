#include <CGAL/Extreme_points_d.h>
#include <CGAL/Extreme_points_options_d.h>
#include <CGAL/Extreme_points_traits_d.h>
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

//Testing bulk deletion with moment curve points

void test(int d, int n) {
  srand (time(NULL));
  CGAL::Extreme_points_options_d op;
  op.set_algorithm(CGAL::EP_CHOOSE_APPROPRIATE);
  op.set_deletion(true);

  std::vector<Point_d> extreme_points;
  std::set<Point_d, Less_lexicographically_d> points;
  double coord;
  for (int i=0; i<n; i++) {
    std::vector<double> np;
    coord = rand();
    for (int j=0; j<d; j++) {
      np.push_back(pow(coord, j+1));
    }
    points.insert(Point_d(d,np.begin(),np.end()));
  }

  CGAL::Extreme_points_d<EP_Traits_d> ep(d,op);
  ep.insert(points.begin(),points.end());
  ep.extreme_points(std::back_inserter(extreme_points));

  //this should hold because points are on a moment curve
  assert(extreme_points.size() == points.size());
  for (int i=0; i<extreme_points.size(); i++) {
    assert(points.find(extreme_points[i]) != points.end());
  }

  std::set<Point_d, Less_lexicographically_d>::iterator it,itt;
  while (points.size() != 0) {
    it=points.begin();
    itt=points.begin();
    itt++; itt++;
    ep.remove(it, itt);
    points.erase(points.begin());
    points.erase(points.begin());
    extreme_points.clear();
    ep.extreme_points(std::back_inserter(extreme_points));

    assert(extreme_points.size() == points.size());
    for (int i=0; i<extreme_points.size(); i++) {
      assert(points.find(extreme_points[i]) != points.end());
    }
  }

  std::cout<<"test finished successfully!"<<std::endl;
}

int main(int argc, char **argv) {
  test(5, 100);
  return 0;
}
