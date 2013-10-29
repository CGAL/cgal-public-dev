#include <CGAL/Extreme_points_d.h>
#include <CGAL/Extreme_points_options_d.h>
#include <CGAL/Extreme_points_traits_d.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Convex_hull_d.h>
#include <CGAL/Random.h>
#include <time.h>
#include <algorithm>

typedef CGAL::Cartesian_d<double>               Kernel;
typedef Kernel::Point_d                         Point_d;
typedef CGAL::Extreme_points_traits_d<Point_d>  EP_Traits_d;
typedef Kernel::Less_lexicographically_d        Less_lexicographically;

int d=6;

void test1() {
    std::vector<Point_d> points;
    std::vector<Point_d> extreme_points,extreme_points2;
    clock_t start,end,t1,t2,t3,t4;

    //random points inside [-100,100]^d
    CGAL::Random rnd(12345);
    double p[d];
    for (int i=0; i<500; ++i) {
      for (int j=0; j<d; ++j) {
        p[j] = rnd.get_double(-100.,100.);
      }
      points.push_back(Point_d(d,p,p+d));
    }

    //EP_SIMPLE
    CGAL::Extreme_points_options_d op;
    op.set_algorithm(CGAL::EP_SIMPLE);
    start = clock();
    CGAL::Extreme_points_d<EP_Traits_d> ep(d,op);
    ep.insert(points.begin(),points.end());
    ep.extreme_points(std::back_inserter(extreme_points));
    end = clock();
    std::cout << "Test1 - extreme ~ points - SIMPLE " << (end-start)/CLOCKS_PER_SEC << std::endl;

    //t1=end-start;

    extreme_points.clear();
    op.set_algorithm(CGAL::EP_DULA_HELGASON);
    start = clock();
    CGAL::Extreme_points_d<EP_Traits_d> ep2(d,op);
    ep2.insert(points.begin(),points.end());
    ep2.extreme_points(std::back_inserter(extreme_points));
    end = clock();
    std::cout << "Test1 - extreme ~ points - DULA HELGASON " << (end-start)/CLOCKS_PER_SEC << std::endl;
    
    //t2=end-start;

    start = clock();
    CGAL::Convex_hull_d<Kernel> chull(d);
    for (int i=0; i<points.size(); ++i)
      chull.insert(points[i]);
    for (CGAL::Convex_hull_d<Kernel>::Vertex_iterator it=chull.vertices_begin(); it!=chull.vertices_end(); ++it)
      extreme_points2.push_back(it->point());
    end = clock();
    std::cout << "Test1 - extreme ~ points - CONVEX HULL " << (end-start)/CLOCKS_PER_SEC << std::endl;

/*    assert(extreme_points.size()==extreme_points2.size());
    std::sort(extreme_points.begin(),extreme_points.end(), Less_lexicographically());
    std::sort(extreme_points2.begin(),extreme_points2.end(), Less_lexicographically());
    assert(std::equal(extreme_points.begin(),extreme_points.end(),extreme_points2.begin()));
*/
    std::cout << "number of extreme points: " << extreme_points.size() << std::endl;  
    std::cout << "number of points: " << points.size() << std::endl;  

    //Add more points, see what happens
/*    for (int i=0; i<50000; ++i) {
      //std::cout << *rnd << std::endl;
      points.push_back(*rnd);
      rnd++;
    }*/
    /*
    op.set_algorithm(CGAL::EP_SIMPLE);
    start = clock();
    ep.insert(points.begin(),points.end());
    ep.extreme_points(std::back_inserter(extreme_points));
    end = clock();
    std::cout << "Test1 - extreme ~ points - SIMPLE " << (end-start)/CLOCKS_PER_SEC << std::endl;

    t1=end-start;

    extreme_points.clear();
    op.set_algorithm(CGAL::EP_DULA_HELGASON);
    start = clock();
    ep2.insert(points.begin(),points.end());
    ep2.extreme_points(std::back_inserter(extreme_points));
    end = clock();
    std::cout << "Test1 - extreme ~ points - DULA HELGASON " << (end-start)/CLOCKS_PER_SEC << std::endl;
    
    t2=end-start;
    */
/*    extreme_points2.clear();
    start = clock();
    CGAL::convex_hull_3(points.begin(), points.end(), poly);
    for (CGAL::Polyhedron_3<Kernel>::Vertex_iterator it=poly.vertices_begin(); it!=poly.vertices_end(); ++it)
      extreme_points2.push_back(it->point());
    end = clock();
    std::cout << "Test1 - extreme ~ points - CONVEX HULL " << (end-start)/CLOCKS_PER_SEC << std::endl;

    //assert(extreme_points.size()==extreme_points2.size());
    //std::sort(extreme_points.begin(),extreme_points.end(), Less_lexicographically());
    //std::sort(extreme_points2.begin(),extreme_points2.end(), Less_lexicographically());
    //assert(std::equal(extreme_points.begin(),extreme_points.end(),extreme_points2.begin()));

    std::cout << "number of extreme points: " << extreme_points.size() << std::endl;  
    std::cout << "number of points: " << points.size() << std::endl;  
 

    points.push_back(Point_3(-1000.,-1000.,-1000.));
    points.push_back(Point_3(1000.,-1000.,-1000.));
    points.push_back(Point_3(-1000.,1000.,-1000.));
    points.push_back(Point_3(-1000.,-1000.,1000.));
    points.push_back(Point_3(1000.,1000.,-1000.));
    points.push_back(Point_3(-1000.,1000.,1000.));
    points.push_back(Point_3(-1000.,-1000.,1000.));
    points.push_back(Point_3(1000.,1000.,1000.));
    
    extreme_points.clear();
    op.set_algorithm(CGAL::EP_SIMPLE);
    start = clock();
    CGAL::Extreme_points_d<EP_Traits_d> ep3(3,op);
    ep3.insert(points.begin(),points.end());
    ep3.extreme_points(std::back_inserter(extreme_points));
    end = clock();
    std::cout << "Test2 - extreme << points - SIMPLE " << (end-start)/CLOCKS_PER_SEC << std::endl;
    
    t3=end-start;
    
    extreme_points.clear();
    op.set_algorithm(CGAL::EP_DULA_HELGASON);
    start = clock();
    CGAL::Extreme_points_d<EP_Traits_d> ep4(3,op);
    ep4.insert(points.begin(),points.end());
    ep4.extreme_points(std::back_inserter(extreme_points));
    end = clock();
    std::cout << "Test1 - extreme << points - DULA HELGASON " << (end-start)/CLOCKS_PER_SEC << std::endl;
    
    t4=end-start;

    extreme_points2.clear();
    start = clock();
    CGAL::convex_hull_3(points.begin(), points.end(), poly);
    for (CGAL::Polyhedron_3<Kernel>::Vertex_iterator it=poly.vertices_begin(); it!=poly.vertices_end(); ++it)
      extreme_points2.push_back(it->point());
    end = clock();
    std::cout << "Test1 - extreme << points - CONVEX HULL " << (end-start)/CLOCKS_PER_SEC << std::endl;

    assert(extreme_points.size()==extreme_points2.size());
    std::sort(extreme_points.begin(),extreme_points.end(), Less_lexicographically());
    std::sort(extreme_points2.begin(),extreme_points2.end(), Less_lexicographically());
    assert(std::equal(extreme_points.begin(),extreme_points.end(),extreme_points2.begin()));
    
    std::cout << "number of extreme points: " << extreme_points.size() << std::endl;  
    std::cout << "number of points: " << points.size() << std::endl;  
    
    assert(t3>t4);
 
//    assert(op.get_last_used_algorithm() == NULL);

    std::cout<<"test1 finished successfully!"<<std::endl;
*/
}

int main(int argc, char **argv) {
    test1();
    return 0;
}
