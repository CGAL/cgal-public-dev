#include <CGAL/Extreme_points_d.h>
#include <CGAL/Extreme_points_options_d.h>
#include <CGAL/Extreme_points_traits_d.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Convex_hull_d.h>
#include <CGAL/point_generators_d.h>

typedef CGAL::Cartesian_d<double>               Kernel;
typedef Kernel::Point_d                         Point_d;
typedef CGAL::Extreme_points_traits_d<Point_d>  EP_Traits_d;
typedef Kernel::Less_lexicographically_d        Less_lexicographically;

clock_t t1,t2,t3;

void test(int d, int n, double ratio) {
  CGAL::Random_points_on_sphere_d<Point_d> rnd(d,100.);
  CGAL::Random_points_on_sphere_d<Point_d> rnd2(d,1000.);
  CGAL::Extreme_points_options_d op;
  op.set_algorithm(CGAL::EP_CHOOSE_APPROPRIATE);
  op.set_deletion(false);
  std::vector<Point_d> points;
  std::vector<Point_d> extreme_points,extreme_points_ch;
  for (int i=0; i<(int) ratio*n; i++) {
    points.push_back(*rnd);
    rnd++;
  }
  while (points.size() < n) {
    points.push_back(*rnd2);
    rnd2++;
  }

  //Extreme_points_d
  t1 = clock();
  CGAL::Extreme_points_d<EP_Traits_d> ep(d,op);
  ep.insert(points.begin(),points.end());
  ep.extreme_points(std::back_inserter(extreme_points));
  t2 = clock();

  //Convex_hull_d
  CGAL::Convex_hull_d<Kernel> chull(d);
  for (int i=0; i<n; ++i)
    chull.insert(points[i]);
  for (CGAL::Convex_hull_d<Kernel>::Vertex_iterator it=chull.vertices_begin(); it!=chull.vertices_end(); ++it)
    extreme_points_ch.push_back(it->point());
  t3 = clock();

  assert(extreme_points.size() == extreme_points_ch.size());
  //Output
  std::cout << std::fixed << std::setprecision(2) << ratio << "\t" << d << "\t" << n;
  if (t2-t1 > t3-t2) {
    std::cout << "\tConvex_Hull_d" << std::endl;
  } else {
    std::cout << "\tExtreme_points_d" << std::endl;
  }
}

int main(int argc, char **argv) {
  for (int d=2; d<8; d++) {
    for (int n=100; n<200; n+=100) {
      test(d,n,1.0);
      test(d,n,0.5);
      test(d,n,0.1);
    }
  }
  return 0;
}
