#include <vector>
#include <iterator>
#include <iostream>
#include <CGAL/Extreme_points_d.h>
#include <CGAL/Extreme_points_options_d.h>
#include <CGAL/Extreme_points_traits_d.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Convex_hull_d.h>
#include <CGAL/point_generators_d.h>
#include <time.h>

/* Compares Convex_Hull_d and Extreme_points_d for different dimensions and numbers of points.
   Points are random on a moment curve.
*/

typedef CGAL::Cartesian_d<double>               Kernel;
typedef Kernel::Point_d                         Point_d;
typedef CGAL::Extreme_points_traits_d<Point_d>  EP_Traits_d;
typedef Kernel::Less_lexicographically_d        Less_lexicographically;

int main(int argc, char **argv) {
  //std::cout.width(30);
  std::cout << "Dim   Points   Faster package" << std::endl;
  clock_t t1,t2,t3,t4;
  srand (time(NULL));
  std::vector<Point_d> points;
  std::vector<Point_d> extreme_points,extreme_points_ch;
  std::vector<double> point;
  /*Testing for dimensions 2 to 6*/
  for (int d=2; d<=6; d++) {
    /*Testing for a 100 points*/
    for (int n=100; n<200; n+=100) {
      CGAL::Extreme_points_options_d op;
      op.set_algorithm(CGAL::EP_CHOOSE_APPROPRIATE);
      op.set_deletion(false);

      //Creating the points on the moment curve
      double coord;
      for (int i=0; i<n; i++) {
        coord=rand()%100000 + 1;
        for (int j=0; j<d; j++) {
          point.push_back(coord);
          coord*=coord;
        }
        points.push_back(Point_d(d,point.begin(),point.end()));
        point.clear();
      }

      //Extreme_points_d
      t1 = clock();
      CGAL::Extreme_points_d<EP_Traits_d> ep(d,op);
      ep.insert(points.begin(),points.end());
      ep.extreme_points(std::back_inserter(extreme_points));
      t2 = clock();

      //Convex_Hull_d
      CGAL::Convex_hull_d<Kernel> chull(d);
      for (int i=0; i<n; ++i)
        chull.insert(points[i]);
      for (CGAL::Convex_hull_d<Kernel>::Vertex_iterator it=chull.vertices_begin(); 
           it!=chull.vertices_end(); ++it)
        extreme_points_ch.push_back(it->point());
      t3 = clock();

      //Output
      std::cout << "  " << d << "      " << n;
      if (t2-t1 > t3-t2) {
        std::cout << "   Convex_Hull_d" << std::endl;
      } else {
        std::cout << "   Extreme_points_d" << std::endl;
      }
      points.clear(); extreme_points.clear(); extreme_points_ch.clear();
    }
  }

  return 0;
}
