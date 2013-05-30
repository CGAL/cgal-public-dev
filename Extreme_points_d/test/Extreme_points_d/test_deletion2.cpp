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
    for (int i=0; i<2000; ++i) {
      points.push_back(*rnd);
      rnd++;
    }

    CGAL::Extreme_points_options_d op;
    op.set_algorithm(CGAL::EP_SIMPLE);
    op.set_deletion(true);
    
    //2000 POINTS ON A CIRCLE'S CIRCUMFERENCE
    start = clock();
    CGAL::Extreme_points_d<EP_Traits_2> ep(2,op);
    ep.insert(points.begin(),points.end());
    ep.get_extreme_points(std::back_inserter(extreme_points));
    end = clock();
    std::cout << "Test1 - extreme ~ points - SIMPLE " << (end-start)/CLOCKS_PER_SEC << std::endl;
    std::cout << "CV size: " << extreme_points.size() << std::endl;
    //std::cout << ep.all_points.size() << std::endl;

    t1=end-start;

    points.push_back(Point_2(-1000.,-1000.));
    points.push_back(Point_2(1000.,1000.));
    points.push_back(Point_2(-1000.,1000.));
    points.push_back(Point_2(1000.,-1000.));
    points.push_back(Point_2(1.,1.));
    extreme_points.clear();

    //5 NEW POINTS, 4 OF THEM CONSIST THE NEW CV
    //2000 POINTS ON A CIRCLE'S CIRCUMFERENCE
    start = clock();
    ep.insert(points.begin(),points.end());
    ep.get_extreme_points(std::back_inserter(extreme_points));
    end = clock();
    std::cout << "Test2 - extreme << points - SIMPLE " << (end-start)/CLOCKS_PER_SEC << std::endl;
    std::cout << "CV size: " << extreme_points.size() << std::endl;
    assert(extreme_points.size() == 4);
    //for (int i=0; i<extreme_points.size(); i++)
    //  std::cout << extreme_points[i] << std::endl;
 
    //REMOVE A POINT FROM THE CV
    //5 NEW POINTS, 4 OF THEM CONSIST THE NEW CV
    //2000 POINTS ON A CIRCLE'S CIRCUMFERENCE
    ep.remove(Point_2(-1000.,-1000.));
    extreme_points.clear();
    ep.get_extreme_points(std::back_inserter(extreme_points));
    std::cout << "Test2 - extreme << points - SIMPLE " << (end-start)/CLOCKS_PER_SEC << std::endl;
    std::cout << "CV size: " << extreme_points.size() << std::endl;
    assert(extreme_points.size() > 4);
    //for (int i=0; i<extreme_points.size(); i++)
    //  std::cout << extreme_points[i] << std::endl;
    
    //REMOVE EXTERNAL POINT
    ep.remove(Point_2(-10000.,-10000.));
    extreme_points.clear();
    ep.get_extreme_points(std::back_inserter(extreme_points));
    std::cout << "Test3 - extreme << points - SIMPLE " << (end-start)/CLOCKS_PER_SEC << std::endl;
    std::cout << "CV size: " << extreme_points.size() << std::endl;
 
    //REMOVE INTERNAL
    ep.remove(Point_2(1.,1.));
    extreme_points.clear();
    ep.get_extreme_points(std::back_inserter(extreme_points));
    std::cout << "Test4 - extreme << points - SIMPLE " << (end-start)/CLOCKS_PER_SEC << std::endl;
    std::cout << "CV size: " << extreme_points.size() << std::endl;
 
    std::cout<<"test1 finished successfully!"<<std::endl;
}

int main(int argc, char **argv) {
    test1();
    return 0;
}
