#include <CGAL/Extreme_points_d.h>
#include <CGAL/Extreme_points_options_d.h>
#include <CGAL/Extreme_points_traits_d.h>
#include <CGAL/Point_2.h>
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
    start = clock();
    CGAL::Extreme_points_d<EP_Traits_2> ep(2,op);
    ep.insert(points.begin(),points.end());
    ep.extreme_points(std::back_inserter(extreme_points));
    end = clock();
    std::cout << "Test1 - extreme ~ points - SIMPLE " << (end-start)/CLOCKS_PER_SEC << std::endl;

    t1=end-start;

    extreme_points.clear();
    op.set_algorithm(CGAL::EP_DULA_HELGASON);
    start = clock();
    CGAL::Extreme_points_d<EP_Traits_2> ep2(2,op);
    ep2.insert(points.begin(),points.end());
    ep2.extreme_points(std::back_inserter(extreme_points));
    end = clock();
    std::cout << "Test1 - extreme ~ points - DULA HELGASON " << (end-start)/CLOCKS_PER_SEC << std::endl;
    
    t2=end-start;


    extreme_points.clear();
    start = clock();
    CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(extreme_points));
    end = clock();
    std::cout << "Test1 - extreme ~ points - CONVEX HULL " << (end-start)/CLOCKS_PER_SEC << std::endl;

    
    std::cout << "number of extreme points: " << extreme_points.size() << std::endl;  
    std::cout << "number of points: " << points.size() << std::endl;  

    assert(t1<t2);

    points.push_back(Point_2(-1000.,-1000.));
    points.push_back(Point_2(1000.,1000.));
    points.push_back(Point_2(-1000.,1000.));
    points.push_back(Point_2(1000.,-1000.));
    
    extreme_points.clear();
    op.set_algorithm(CGAL::EP_SIMPLE);
    start = clock();
    CGAL::Extreme_points_d<EP_Traits_2> ep3(2,op);
    ep3.insert(points.begin(),points.end());
    ep3.extreme_points(std::back_inserter(extreme_points));
    end = clock();
    std::cout << "Test2 - extreme << points - SIMPLE " << (end-start)/CLOCKS_PER_SEC << std::endl;
    
    t3=end-start;
    
    extreme_points.clear();
    op.set_algorithm(CGAL::EP_DULA_HELGASON);
    start = clock();
    CGAL::Extreme_points_d<EP_Traits_2> ep4(2,op);
    ep4.insert(points.begin(),points.end());
    ep4.extreme_points(std::back_inserter(extreme_points));
    end = clock();
    std::cout << "Test1 - extreme ~ points - DULA HELGASON " << (end-start)/CLOCKS_PER_SEC << std::endl;
    
    t4=end-start;

    extreme_points.clear();
    start = clock();
    CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(extreme_points));
    end = clock();
    std::cout << "Test1 - extreme ~ points - CONVEX HULL " << (end-start)/CLOCKS_PER_SEC << std::endl;
    
    std::cout << "number of extreme points: " << extreme_points.size() << std::endl;  
    std::cout << "number of points: " << points.size() << std::endl;  
    
    assert(t3>t4);
 
    std::cout<<"test1 finished successfully!"<<std::endl;
}

int main(int argc, char **argv) {
    test1();
    return 0;
}
