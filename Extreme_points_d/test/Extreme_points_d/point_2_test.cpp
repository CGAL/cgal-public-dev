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
#include <CGAL/Extreme_points_traits_d.h>
#include <CGAL/Point_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian.h>

typedef CGAL::Cartesian<double>                 Kernel;
typedef Kernel::Point_2                         Point_2;
// typedef CGAL::Point_2<Kernel>                   Point_2;
typedef CGAL::Extreme_points_traits_d<Point_2>  EP_Traits_2;
typedef Kernel::Less_xy_2                       Less_lexicographically;

void test1() {
	std::vector<Point_2> points;
	std::vector<Point_2> extreme_points;
    
    points.push_back(Point_2(0.,0.));
    points.push_back(Point_2(1.,0.));
    points.push_back(Point_2(0.,1.));
    points.push_back(Point_2(1.,1.));
    points.push_back(Point_2(0.,.5));
    points.push_back(Point_2(.5,0.));
    points.push_back(Point_2(.5,.5));
    points.push_back(Point_2(.5,1.1));
    points.push_back(Point_2(1.1,.1));
    
    
    const int number_of_extreme_points = 6;
    
    CGAL::extreme_points_d_dula_helgason(points.begin(), points.end(), std::back_inserter(extreme_points));
	
	std::cout<<"found "<<extreme_points.size()<<" extreme points"<<std::endl;
	
	std::vector<Point_2> extreme_points2;
    CGAL::extreme_points_d_simple(points.begin(), points.end(), std::back_inserter(extreme_points2));

    std::vector<Point_2> extreme_points3;
    CGAL::extreme_points_d(points.begin(), points.end(), std::back_inserter(extreme_points3));
    
	assert(extreme_points.size()==number_of_extreme_points);
    assert(extreme_points2.size()==number_of_extreme_points);
    assert(extreme_points3.size()==number_of_extreme_points);
    
    for (int i=0;i<extreme_points.size();++i)
        std::cout<<extreme_points[i]<<std::endl;
    
    sort(extreme_points.begin(),extreme_points.end(), Less_lexicographically());
    sort(extreme_points2.begin(),extreme_points2.end(), Less_lexicographically());
    sort(extreme_points3.begin(),extreme_points3.end(), Less_lexicographically());
    
	// check that the different implementations produce the same output
	assert(std::equal(extreme_points.begin(),extreme_points.end(),extreme_points2.begin()));
    assert(std::equal(extreme_points.begin(),extreme_points.end(),extreme_points3.begin()));
    
    // testing dynamic class
    CGAL::Extreme_points_d<EP_Traits_2> ep(2);
    ep.insert(points.begin(),points.end());
    std::vector<Point_2> extreme_points4;
    ep.get_extreme_points(std::back_inserter(extreme_points4));
    
    // check that the different implementations produce the same output
    assert(extreme_points4.size()==number_of_extreme_points);
    sort(extreme_points4.begin(),extreme_points4.end(), Less_lexicographically());
    assert(std::equal(extreme_points.begin(),extreme_points.end(),extreme_points4.begin()));
    
    // testing classification functions..
    assert(ep.classify(Point_2(0.,0.))==CGAL::EXTREME_POINT);
    assert(ep.classify(Point_2(1.,0.))==CGAL::EXTREME_POINT);
    assert(ep.classify(Point_2(0.,1.))==CGAL::EXTREME_POINT);
    assert(ep.classify(Point_2(1.,1.))==CGAL::EXTREME_POINT);
    assert(ep.classify(Point_2(0.,.5))==CGAL::INTERNAL_POINT);
    assert(ep.classify(Point_2(.5,0.))==CGAL::INTERNAL_POINT);
    assert(ep.classify(Point_2(.5,.5))==CGAL::INTERNAL_POINT);
    assert(ep.classify(Point_2(.5,1.1))==CGAL::EXTREME_POINT);
    assert(ep.classify(Point_2(1.1,.1))==CGAL::EXTREME_POINT);
    assert(ep.classify(Point_2(-1.,-1.))==CGAL::EXTERNAL_POINT);
    assert(ep.classify(Point_2(2.,.5))==CGAL::EXTERNAL_POINT);
    
    assert(ep.classify(Point_2(0.,0.),true)==CGAL::EXTREME_POINT);
    assert(ep.classify(Point_2(1.,0.),true)==CGAL::EXTREME_POINT);
    assert(ep.classify(Point_2(0.,1.),true)==CGAL::EXTREME_POINT);
    assert(ep.classify(Point_2(1.,1.),true)==CGAL::EXTREME_POINT);
    assert(ep.classify(Point_2(0.,.5),true)==CGAL::INTERNAL_POINT);
    assert(ep.classify(Point_2(.5,0.),true)==CGAL::INTERNAL_POINT);
    assert(ep.classify(Point_2(.5,.5),true)==CGAL::INTERNAL_POINT);
    assert(ep.classify(Point_2(.5,1.1),true)==CGAL::EXTREME_POINT);
    assert(ep.classify(Point_2(1.1,.1),true)==CGAL::EXTREME_POINT);
    
    
    
	std::cout<<"test1 finished successfully!"<<std::endl;
}

void test2() {
    std::vector<Point_2> points;
    std::vector<Point_2> extreme_points;
    
    points.push_back(Point_2(-0.786408,-0.262136));
    points.push_back(Point_2(0.699029,-0.791262));
    points.push_back(Point_2(0.631068,0.524272));
    points.push_back(Point_2(0.218447,-0.15534));
    
    
    CGAL::extreme_points_d_dula_helgason(points.begin(), points.end(), std::back_inserter(extreme_points));
    
    std::cout<<"found "<<extreme_points.size()<<" extreme points"<<std::endl;
    
    std::vector<Point_2> extreme_points2;
    CGAL::extreme_points_d_simple(points.begin(), points.end(), std::back_inserter(extreme_points2));
    
    std::vector<Point_2> extreme_points3;
    CGAL::extreme_points_d(points.begin(), points.end(), std::back_inserter(extreme_points3));

    assert(extreme_points.size()==3);
    assert(extreme_points2.size()==3);
    assert(extreme_points3.size()==3);
    
    for (int i=0;i<extreme_points.size();++i)
        std::cout<<extreme_points[i]<<std::endl;
    
    sort(extreme_points.begin(),extreme_points.end(), Less_lexicographically());
    sort(extreme_points2.begin(),extreme_points2.end(), Less_lexicographically());
    sort(extreme_points3.begin(),extreme_points3.end(), Less_lexicographically());
    
    // check that the different implementations produce the same output
    assert(std::equal(extreme_points.begin(),extreme_points.end(),extreme_points2.begin()));
    assert(std::equal(extreme_points.begin(),extreme_points.end(),extreme_points3.begin()));
    
    std::cout<<"test2 finished successfully!"<<std::endl;
}

int main(int argc, char **argv) {
    test1();
    test2();
    return 0;
}