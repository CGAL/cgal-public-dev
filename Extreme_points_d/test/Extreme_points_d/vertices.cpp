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
#include <CGAL/Point_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian.h>

#include <CGAL/convex_hull_2.h>
#include <CGAL/point_generators_2.h>

#include <time.h>

#include <CGAL/convex_hull_2.h>

typedef CGAL::Cartesian<double>                 Kernel;
typedef Kernel::Point_2                         Point_2;
typedef CGAL::Extreme_points_traits_d<Point_2>  EP_Traits_2;
typedef Kernel::Less_xy_2                       Less_lexicographically;

void test1() {
    std::vector<Point_2> points;
    std::vector<Point_2> extreme_points;
    clock_t start, end,t1,t2,t3,t4;
   
    CGAL::Random_points_on_circle_2<Point_2> rnd(100.);
    for (int i=0; i<100; ++i) {
      points.push_back(*rnd);
      rnd++;
    }

    CGAL::Extreme_points_options_d op;
    
    extreme_points.clear();
    op.set_algorithm(CGAL::EP_SIMPLE);
    start = clock();
    CGAL::Extreme_points_d<EP_Traits_2> ep3(2,op);
    ep3.insert(points.begin(),points.end());
    ep3.extreme_points(std::back_inserter(extreme_points));
    end = clock();
    std::cout << "Test2 - extreme << points - SIMPLE " << (end-start)/CLOCKS_PER_SEC << std::endl;
    
    int vertices=0;
    for (int i=0; i<points.size(); i++) {
      for (int j=i; j<points.size(); j++) {
        if (ep3.classify(Point_2((points[i][0] + points[j][0])/2,(points[i][1] + points[j][1])/2))==CGAL::ON_BOUNDARY)
          vertices++;
      }
    }

    extreme_points.clear();
    start = clock();
    CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(extreme_points));
    end = clock();
    std::cout << "Test1 - extreme ~ points - CONVEX HULL " << (end-start)/CLOCKS_PER_SEC << std::endl;
    

    std::cout << "Vertices: " << vertices << std::endl;

    std::cout << "number of extreme points: " << extreme_points.size() << std::endl;  
    std::cout << "number of points: " << points.size() << std::endl;  
}

int main(int argc, char **argv) {
    test1();
    return 0;
}
