#include <CGAL/point_generators_d.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Extreme_points_d.h>

typedef CGAL::Cartesian_d<double>  Kernel_d;
typedef Kernel_d::Point_d      Point_d;

int main() {
  const int D = 10;
  const int N = 500;
  
  // D dimensional points with coordinates in the range [-10, 10]
  CGAL::Random_points_in_cube_d<Point_d> gen (D, 10.);
  
  // generate N points randomly in the D dimensional iso box
  // and copy them to a vector
  std::vector<Point_d> points;
  CGAL::copy_n( gen, N, std::back_inserter(points) );

  std::cout<<"Running the Dula-Helgason algorithm on a point set of "
       <<N<<" points randomly chosen out of a "<<D<<"-dimensional "
       <<"iso box"<<std::endl;
  
  // compute the extreme points
  std::vector<Point_d> extreme_points;
  CGAL::extreme_points_d_dula_helgason(points.begin(), points.end(),
        std::back_inserter(extreme_points));
  
  std::cout<<"Found "<<extreme_points.size()
           <<" extreme points:"<<std::endl;

  for (std::vector<Point_d>::iterator it=extreme_points.begin();
     it!=extreme_points.end();
     it++) {
    std::cout<<*it<<std::endl;
  }

  return 0;
}
