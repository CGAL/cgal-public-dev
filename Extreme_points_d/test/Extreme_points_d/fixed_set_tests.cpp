#include <vector>
#include <fstream>
#include <string>
#include <iterator>
#include <iostream>
#include <algorithm>
#include <cassert>

#include <CGAL/config.h>
#include <CGAL/Timer.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Extreme_points_d.h>

typedef CGAL::Cartesian_d<double> Kernel_d;
typedef Kernel_d::Point_d Point_d;

template <class InputStream, class OutputIterator>
void read_input(InputStream &in, OutputIterator out, int &n, int &d) {
    // read *.con file
    std::string s;
    do {
        getline(in,s);
    } while (s[0]=='%'); // skipping comments
    std::stringstream ss(s);
    ss>>d>>n;
    
    std::vector<Point_d> points(n);
    for (int i=0;i<n;++i) {
        std::vector<double> p(d);
        for (int j=0;j<d;++j)
            in>>p[j];
        *out++=Point_d(d,p.begin(),p.end());
    }
}

int main(int argc, char **argv) {
    std::ifstream fin("./data/05by02500at01.con");
    std::vector<Point_d> points;
    int n,d,p=1;
    
    read_input(fin,std::back_inserter(points),n,d);
    std::cout<<"input read"<<std::endl;
    
    std::vector<Point_d> extreme_points;
    extreme_points_d_dula_helgason(points.begin(), points.end(),std::back_inserter(extreme_points));
	
    std::cout<<"extreme_points_dula_helgason found "<<extreme_points.size()<<" extreme points"<<std::endl;
    assert(extreme_points.size()==n*p/100);
	
    std::vector<Point_d> extreme_points2;
    extreme_points_d_simple(points.begin(), points.end(),std::back_inserter(extreme_points2));
    
    std::cout<<"extreme_points_d_simple found "<<extreme_points2.size()<<" extreme points"<<std::endl;
    assert(extreme_points2.size()==n*p/100);
    
    std::vector<Point_d> extreme_points3;
    extreme_points_d(points.begin(), points.end(),std::back_inserter(extreme_points3));
    
    std::cout<<"extreme_points_d found "<<extreme_points3.size()<<" extreme points"<<std::endl;
    assert(extreme_points3.size()==n*p/100);
    
    sort(extreme_points.begin(),extreme_points.end(), Kernel_d::Less_lexicographically_d());
    sort(extreme_points2.begin(),extreme_points2.end(), Kernel_d::Less_lexicographically_d());
    sort(extreme_points3.begin(),extreme_points3.end(), Kernel_d::Less_lexicographically_d());
    
    // check that the different implementations produce the same output
    assert(std::equal(extreme_points.begin(),extreme_points.end(),extreme_points2.begin()));
    assert(std::equal(extreme_points.begin(),extreme_points.end(),extreme_points3.begin()));
    
    std::cout<<"Checking behaviour of duplicated points.. "<<std::flush;
    
    std::vector<Point_d> points2(points);
    points2.insert(points2.end(),points.begin(),points.end()); // every point twice
    
    std::vector<Point_d> extreme_points4;
    extreme_points_d_dula_helgason(points2.begin(), points2.end(),std::back_inserter(extreme_points4));
    assert(extreme_points4.size()==n*p/100);
    
    std::vector<Point_d> extreme_points5;
    extreme_points_d_simple(points2.begin(), points2.end(),std::back_inserter(extreme_points5));
    assert(extreme_points5.size()==n*p/100);

    std::vector<Point_d> extreme_points6;
    extreme_points_d(points2.begin(), points2.end(),std::back_inserter(extreme_points6));
    assert(extreme_points6.size()==n*p/100);
    
    sort(extreme_points4.begin(),extreme_points4.end(), Kernel_d::Less_lexicographically_d());
    sort(extreme_points5.begin(),extreme_points5.end(), Kernel_d::Less_lexicographically_d());
    sort(extreme_points6.begin(),extreme_points6.end(), Kernel_d::Less_lexicographically_d());
    
    assert(std::equal(extreme_points.begin(),extreme_points.end(),extreme_points4.begin()));
    assert(std::equal(extreme_points.begin(),extreme_points.end(),extreme_points5.begin()));
    assert(std::equal(extreme_points.begin(),extreme_points.end(),extreme_points6.begin()));
    
    std::cout<<"OK"<<std::endl;
    
    std::cout<<"Finished successfully!"<<std::endl;
    return 0;
}
