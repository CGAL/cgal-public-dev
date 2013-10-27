#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Extreme_points_d.h>
#include <CGAL/Extreme_points_traits_d.h>
#include <vector>
#include <iostream>
#include <cassert>

#include <CGAL/convex_hull_2.h>
#include <cmath>

typedef CGAL::Cartesian<double>                 Kernel_2;
typedef Kernel_2::Point_2                       Point_2;
typedef CGAL::Extreme_points_traits_d<Point_2>  EP_Traits_2;

bool compare_results(std::vector<Point_2> m, std::vector<Point_2> n) {
  std::cout << m.size() << "  " << n.size() << std::endl;
  bool equal;
  for (int i=0; i<m.size(); i++) {

    for (std::vector<Point_2>::iterator it=m.begin();
      it!=m.end();
      it++) {
      std::cout<<*it<<std::endl;
    }
    std::cout << std::endl;
    for (std::vector<Point_2>::iterator it=n.begin();
      it!=n.end();
      it++) {
      std::cout<<*it<<std::endl;
    }
    std::cout << std::endl;

    for (int j=0; j<n.size(); j++) {
      equal = true;
      for (int k=0; k<2; k++) {
        if (abs(m[i][k]-n[j][k]) > 1e-3){
          equal=false;
          std::cout << "FALSE " << m[i] << " " << n[j] << std::endl;
        }
      }
      if (equal) {
	std::cout << "TRUE " << m[i] << " " << n[j] << std::endl;
        m.erase(m.begin()+i, m.begin()+i+1);
        n.erase(n.begin()+j, n.begin()+j+1);
        break;
      }
    }
  }

  std::cout << m.size() << "  " << n.size() << std::endl;
  if ((m.size() == 0) && (n.size() == 0))
    return true;
  return false;
}

int main() {
    const int N = 100;     // number of points for every batch
    const int BATCHES = 1; // number of batches
    
    CGAL::Extreme_points_d<EP_Traits_2> ep(2);
    
    // Generator for D-dimensional points with coordinates
    // in the range [-10, 10]
    CGAL::Random_points_in_square_2<Point_2> gen (10.);
    
    for (int i=0;i<BATCHES;++i) {
//        std::cout<<"Inserting the points:"<<std::endl;
        
        // generate N points randomly from the D dimensional iso box
        std::vector<Point_2> points;
        for (int j=0;j<N;++j) {
            Point_2 p = *gen++;
        //    std::cout<<p<<std::endl;
            points.push_back(p);
        }
        
        // add these points to the current point set maintained by ep
        ep.insert(points.begin(), points.end());
        
        // compute the extreme points
        //std::cout<<"\nExtreme points of the current set: "<<std::endl;
        std::vector<Point_2> extreme_points;
        ep.extreme_points(std::back_inserter(extreme_points));
     
        /*for (std::vector<Point_2>::iterator it=extreme_points.begin();
             it!=extreme_points.end();
             it++) {
            std::cout<<*it<<std::endl;
        }*/
        std::vector<Point_2> extreme_points2;
        //std::ostream_iterator< Point_2 >  out( std::cout, "\n" );
        CGAL::convex_hull_2( points.begin(), points.end(), std::back_inserter(extreme_points2) );
        /*for (std::vector<Point_2>::iterator it=extreme_points2.begin();
             it!=extreme_points2.end();
             it++) {
            std::cout<<*it<<std::endl;
        }*/
	std::cout<<std::endl;
      std::cout << compare_results(extreme_points,extreme_points2) << std::endl;
    }
    
    return 0;
}
