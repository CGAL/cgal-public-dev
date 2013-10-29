#include <CGAL/Extreme_points_d.h>
#include <CGAL/Extreme_points_traits_d.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Convex_hull_d.h>

int d = 4;

typedef CGAL::Cartesian_d<double>                 Kernel;
typedef Kernel::Point_d                         Point_d;
// typedef CGAL::Point_2<Kernel>                   Point_2;
typedef CGAL::Extreme_points_traits_d<Point_d>  EP_Traits_d;
typedef Kernel::Less_lexicographically_d        Less_lexicographically;

void test1() {
	std::vector<Point_d> points;
	std::vector<Point_d> extreme_points;

	double coord1[4]={0.,0.,0.,0.};
	points.push_back(Point_d(d,coord1,coord1+d));
	double coord2[4]={1.,0.,0.,0.};
	points.push_back(Point_d(d,coord2,coord2+d));
	double coord3[4]={0.,1.,0.,0.};
	points.push_back(Point_d(d,coord3,coord3+d));
	double coord4[4]={0.,0.,1.,0.};
	points.push_back(Point_d(d,coord4,coord4+d));
	double coord5[4]={0.,0.,0.,1.};
	points.push_back(Point_d(d,coord5,coord5+d));
	double coord6[4]={1.,1.,0.,0.};
	points.push_back(Point_d(d,coord6,coord6+d));
	double coord7[4]={0.,1.,1.,0.};
	points.push_back(Point_d(d,coord7,coord7+d));
	double coord8[4]={0.,0.,1.,1.};
	points.push_back(Point_d(d,coord8,coord8+d));
	double coord9[4]={1.,0.,1.,0.};
	points.push_back(Point_d(d,coord9,coord9+d));
	double coord10[4]={0.,1.,0.,1.};
	points.push_back(Point_d(d,coord10,coord10+d));
	double coord11[4]={1.,0.,0.,1.};
	points.push_back(Point_d(d,coord11,coord11+d));
	double coord12[4]={1.,1.,1.,0.};
	points.push_back(Point_d(d,coord12,coord12+d));
	double coord13[4]={1.,0.,1.,1.};
	points.push_back(Point_d(d,coord13,coord13+d));
	double coord14[4]={1.,1.,0.,1.};
	points.push_back(Point_d(d,coord14,coord14+d));
	double coord15[4]={0.,1.,1.,1.};
	points.push_back(Point_d(d,coord15,coord15+d));
	double coord16[4]={1.,1.,1.,1.};
	points.push_back(Point_d(d,coord16,coord16+d));

    double coord[4] = {0.,0.,0.,0.};
    for (int i=0; i<d; ++i) {
      coord[i]=0.5;
      points.push_back(Point_d(d,coord,coord+d));
    }
    for (int i=0; i<d; ++i) coord[i]=0.0;
    for (int i=d-1; i>=0; --i) {
      coord[i]=0.5;
      points.push_back(Point_d(d,coord,coord+d));
    }

    
    const int number_of_extreme_points = 16;
    
    CGAL::extreme_points_d_dula_helgason(points.begin(), points.end(), std::back_inserter(extreme_points));
	
	std::cout<<"found "<<extreme_points.size()<<" extreme points"<<std::endl;
	
	std::vector<Point_d> extreme_points2;
    CGAL::extreme_points_d_simple(points.begin(), points.end(), std::back_inserter(extreme_points2));

    std::vector<Point_d> extreme_points3;
    CGAL::extreme_points_d(points.begin(), points.end(), std::back_inserter(extreme_points3));
    
    CGAL::Convex_hull_d<Kernel> chull(d);
    for (int i=0; i<points.size(); ++i)
      chull.insert(points[i]);
    std::vector<Point_d> extreme_points5;
    for (CGAL::Convex_hull_d<Kernel>::Vertex_iterator it=chull.vertices_begin(); it!=chull.vertices_end(); ++it)
      extreme_points5.push_back(it->point());

    assert(extreme_points.size()==number_of_extreme_points);
    assert(extreme_points2.size()==number_of_extreme_points);
    assert(extreme_points3.size()==number_of_extreme_points);
    assert(extreme_points5.size()==number_of_extreme_points);
    
    //for (int i=0;i<extreme_points.size();++i)
    //    std::cout<<extreme_points[i]<<std::endl;
    
    sort(extreme_points.begin(),extreme_points.end(), Less_lexicographically());
    sort(extreme_points2.begin(),extreme_points2.end(), Less_lexicographically());
    sort(extreme_points3.begin(),extreme_points3.end(), Less_lexicographically());
    sort(extreme_points5.begin(),extreme_points5.end(), Less_lexicographically());
    
	// check that the different implementations produce the same output
	assert(std::equal(extreme_points.begin(),extreme_points.end(),extreme_points2.begin()));
    assert(std::equal(extreme_points.begin(),extreme_points.end(),extreme_points3.begin()));
    assert(std::equal(extreme_points.begin(),extreme_points.end(),extreme_points5.begin()));
    
    // testing dynamic class
    CGAL::Extreme_points_d<EP_Traits_d> ep(d);
    ep.insert(points.begin(),points.end());
    std::vector<Point_d> extreme_points4;
    ep.extreme_points(std::back_inserter(extreme_points4));
    
    // check that the different implementations produce the same output
    assert(extreme_points4.size()==number_of_extreme_points);
    sort(extreme_points4.begin(),extreme_points4.end(), Less_lexicographically());
    assert(std::equal(extreme_points.begin(),extreme_points.end(),extreme_points4.begin()));
    
    // testing classification functions..
    assert(ep.classify(Point_d(d,coord1,coord1+d))==CGAL::ON_BOUNDARY);
    assert(ep.classify(Point_d(d,coord2,coord2+d))==CGAL::ON_BOUNDARY);
    assert(ep.classify(Point_d(d,coord3,coord3+d))==CGAL::ON_BOUNDARY);
    assert(ep.classify(Point_d(d,coord4,coord4+d))==CGAL::ON_BOUNDARY);
    assert(ep.classify(Point_d(d,coord5,coord5+d))==CGAL::ON_BOUNDARY);
    assert(ep.classify(Point_d(d,coord6,coord6+d))==CGAL::ON_BOUNDARY);
    assert(ep.classify(Point_d(d,coord7,coord7+d))==CGAL::ON_BOUNDARY);
    assert(ep.classify(Point_d(d,coord8,coord8+d))==CGAL::ON_BOUNDARY);
    assert(ep.classify(Point_d(d,coord9,coord9+d))==CGAL::ON_BOUNDARY);
    assert(ep.classify(Point_d(d,coord10,coord10+d))==CGAL::ON_BOUNDARY);
    assert(ep.classify(Point_d(d,coord11,coord11+d))==CGAL::ON_BOUNDARY);
    assert(ep.classify(Point_d(d,coord12,coord12+d))==CGAL::ON_BOUNDARY);
    assert(ep.classify(Point_d(d,coord13,coord13+d))==CGAL::ON_BOUNDARY);
    assert(ep.classify(Point_d(d,coord14,coord14+d))==CGAL::ON_BOUNDARY);
    assert(ep.classify(Point_d(d,coord15,coord15+d))==CGAL::ON_BOUNDARY);
    assert(ep.classify(Point_d(d,coord16,coord16+d))==CGAL::ON_BOUNDARY);


    assert(ep.classify(Point_d(d,coord,coord+d))==CGAL::ON_UNBOUNDED_SIDE);
    for (int i=0; i<d; ++i)
      coord[i]=10.;
    assert(ep.classify(Point_d(d,coord,coord+d))==CGAL::ON_BOUNDED_SIDE);

    std::cout<<"test1 finished successfully!"<<std::endl;
}

int main(int argc, char **argv) {
    test1();
    return 0;
}
