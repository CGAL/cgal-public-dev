// (c) 2007-2008 George Tzoumas <geortz@gmail.com>

#include <iostream>

#include <CGAL/basic.h>
#include <CGAL/Timer.h>
#include <CGAL/Gmpq.h>

#define VERBOSE 1
//#define VCAPXNOCACHE

#include <CGAL/Gmpfi.h>
//#include <CGAL/Boost_interval_Gmpfr.h>

#include <CGAL/Ellipse_traits.h>
#include <CGAL/Ellipse_2.h>


#include <CGAL/Voronoi_circle_apx.h>

typedef CGAL::Gmpq QQ;
typedef CGAL::Gmpfi INT;
typedef CGAL::Ellipse_traits<QQ, INT> ET;
typedef CGAL::Ellipse_2<ET> Ellipse_2;
typedef CGAL::Ellipse_triplet<ET> Ellipse_triplet;
typedef CGAL::VORELL::Voronoi_circle_apx<ET> Voronoi_circle_apx;

Ellipse_2 e1, e2, e3;

using namespace std;

void test_vorcircle()
{
    Ellipse_triplet triplet(e1,e2,e3);
    Voronoi_circle_apx vc(triplet);

    cout << "num vc = " << triplet.get_num_voronoi_circles() << endl;
    cout << vc << endl;
    ET::IT x, y, r;
    vc.get_coords(x, y, r);
    cout << x << y << r << endl;
    
    CGAL::Gmpfr::set_default_precision(2*CGAL::Gmpfr::get_default_precision());
//    CGAL::Gmpfi::set_default_precision(CGAL::Gmpfr::get_default_precision());

    // beware of caching! -- TODO
    Voronoi_circle_apx vc2(vc);
//    Voronoi_circle_apx vc2(triplet);
    cout << vc2 << endl; // this should not be refined, if caching enabled
}

int main(int argc, char **argv)
{
    std::cin >> e1 >> e2 >> e3;

    CGAL::set_pretty_mode (cout);
    CGAL::set_pretty_mode (cerr);
    test_vorcircle();

    return 0;
}

