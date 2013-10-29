#include <CGAL/Extreme_points_d.h>
#include <CGAL/Extreme_points_traits_d.h>
#include <CGAL/Point_3.h>
#include <CGAL/Cartesian.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polyhedron_3.h>

typedef CGAL::Cartesian<double>                 Kernel;
typedef Kernel::Point_3                         Point_3;
typedef CGAL::Extreme_points_traits_d<Point_3>  EP_Traits_3;
typedef Kernel::Less_xyz_3                      Less_lexicographically;

void test1() {
    std::vector<Point_3> points;
    std::vector<Point_3> extreme_points;
    
    Point_3 p1(0.,0.,0.); points.push_back(p1);
    Point_3 p2(0.,0.,1.); points.push_back(p2);
    Point_3 p3(0.,1.,0.); points.push_back(p3);
    Point_3 p4(0.,1.,1.); points.push_back(p4);
    Point_3 p5(1.,0.,0.); points.push_back(p5);
    Point_3 p6(1.,0.,1.); points.push_back(p6);
    Point_3 p7(1.,1.,0.); points.push_back(p7);
    Point_3 p8(1.,1.,1.); points.push_back(p8);
    Point_3 p9(.5,.5,2.); points.push_back(p9);
    Point_3 p10(.5,.5,1.5); points.push_back(p10);
    Point_3 p11(.5,.5,.5); points.push_back(p11);
    Point_3 p12(1.,.5,.5); points.push_back(p12);
    Point_3 p13(0.,0.,.5); points.push_back(p13);
    Point_3 p14(0.,0.,.5); points.push_back(p14); // duplicate, not extreme
    Point_3 p15(0.,0.,0.); points.push_back(p15); // duplicate, extreme

    const int number_of_extreme_points = 9;

    CGAL::extreme_points_d_dula_helgason(points.begin(), points.end(), std::back_inserter(extreme_points));

    std::cout<<"found "<<extreme_points.size()<<" extreme points"<<std::endl;

    std::vector<Point_3> extreme_points2;
    CGAL::extreme_points_d_simple(points.begin(), points.end(), std::back_inserter(extreme_points2));
	
    std::vector<Point_3> extreme_points3;
    CGAL::extreme_points_d(points.begin(), points.end(), std::back_inserter(extreme_points3));
 
    CGAL::Polyhedron_3<Kernel> polyh;
    CGAL::convex_hull_3(points.begin(), points.end(), polyh);
    std::vector<Point_3> extreme_points5;
    for (CGAL::Polyhedron_3<Kernel>::Vertex_iterator it=polyh.vertices_begin(); it!=polyh.vertices_end(); ++it)
      extreme_points5.push_back(it->point());

    assert(extreme_points.size()==number_of_extreme_points);
    assert(extreme_points2.size()==number_of_extreme_points);
    assert(extreme_points3.size()==number_of_extreme_points);
    assert(extreme_points5.size()==number_of_extreme_points);
    
    sort(extreme_points.begin(),extreme_points.end(), Less_lexicographically());
    sort(extreme_points2.begin(),extreme_points2.end(), Less_lexicographically());
    sort(extreme_points3.begin(),extreme_points3.end(), Less_lexicographically());
    sort(extreme_points5.begin(),extreme_points5.end(), Less_lexicographically());
    
    for (int i=0;i<extreme_points.size();++i)
        std::cout<<extreme_points[i]<<std::endl;
    
    // check that the different implementations produce the same output
    assert(std::equal(extreme_points.begin(),extreme_points.end(),extreme_points2.begin()));
    assert(std::equal(extreme_points.begin(),extreme_points.end(),extreme_points3.begin()));
    assert(std::equal(extreme_points.begin(),extreme_points.end(),extreme_points5.begin()));
    
    // testing dynamic class
    CGAL::Extreme_points_d<EP_Traits_3> ep(3);
    ep.insert(points.begin(),points.end());
    std::vector<Point_3> extreme_points4;
    ep.extreme_points(std::back_inserter(extreme_points4));
    
    // check that the different implementations produce the same output
    assert(extreme_points4.size()==number_of_extreme_points);
    sort(extreme_points4.begin(),extreme_points4.end(), Less_lexicographically());
    assert(std::equal(extreme_points.begin(),extreme_points.end(),extreme_points4.begin()));
    
    // testing classification functions..
    assert(ep.classify(p1)==CGAL::ON_BOUNDARY);
    assert(ep.classify(p2)==CGAL::ON_BOUNDARY);
    assert(ep.classify(p3)==CGAL::ON_BOUNDARY);
    assert(ep.classify(p4)==CGAL::ON_BOUNDARY);
    assert(ep.classify(p5)==CGAL::ON_BOUNDARY);
    assert(ep.classify(p6)==CGAL::ON_BOUNDARY);
    assert(ep.classify(p7)==CGAL::ON_BOUNDARY);
    assert(ep.classify(p8)==CGAL::ON_BOUNDARY);
    assert(ep.classify(p9)==CGAL::ON_BOUNDARY);
    assert(ep.classify(p15)==CGAL::ON_BOUNDARY);
    assert(ep.classify(p10)==CGAL::ON_UNBOUNDED_SIDE);
    assert(ep.classify(p11)==CGAL::ON_UNBOUNDED_SIDE);
    assert(ep.classify(p12)==CGAL::ON_UNBOUNDED_SIDE);
    assert(ep.classify(p13)==CGAL::ON_UNBOUNDED_SIDE);
    assert(ep.classify(p14)==CGAL::ON_UNBOUNDED_SIDE);

    assert(ep.classify(Point_3(.1,.1,.1))==CGAL::ON_UNBOUNDED_SIDE);
    assert(ep.classify(Point_3(1.,0.,.5))==CGAL::ON_UNBOUNDED_SIDE);
    assert(ep.classify(Point_3(1.,.3,0.))==CGAL::ON_UNBOUNDED_SIDE);
    assert(ep.classify(Point_3(.1,.1,.1))==CGAL::ON_UNBOUNDED_SIDE);
    
    assert(ep.classify(Point_3(.5,.5,2.1))==CGAL::ON_BOUNDED_SIDE);
    assert(ep.classify(Point_3(.5,.5,-.1))==CGAL::ON_BOUNDED_SIDE);
    assert(ep.classify(Point_3(-.5,0.,0.))==CGAL::ON_BOUNDED_SIDE);
    
    std::cout<<"test1 finished successfully!"<<std::endl;
}

int main(int argc, char **argv) {
    test1();
    return 0;
}
