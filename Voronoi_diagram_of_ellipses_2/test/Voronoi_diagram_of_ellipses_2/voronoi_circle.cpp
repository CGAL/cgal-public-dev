// (c) 2007-2009 George Tzoumas <geortz@gmail.com>

#define VERBOSE 2

#include <iostream>

#include <CGAL/basic.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpfi.h>

#include <CGAL/Ellipse_traits.h>
#include <CGAL/Ellipse_2.h>
#include <CGAL/Ellipse_triplet.h>
//#include <CGAL/Voronoi_circle.h>
#include <CGAL/Voronoi_circle_exact.h>



//typedef CGAL::mpfr_interval INT;
typedef CGAL::Gmpfi INT;
typedef CGAL::Gmpq QQ;
typedef CGAL::Ellipse_traits<QQ, INT> ET;
typedef CGAL::Ellipse_2<ET> Ellipse_2;

Ellipse_2 e1, e2, e3;

using namespace std;

// TODO:
//        check polz / polq (what to use)
//        could be faster with NTL?


void compute_vc()
{
    CGAL::Ellipse_triplet<ET> triplet(e1, e2, e3);
    CGAL::VORELL::Voronoi_circle_exact<ET> vc(triplet);
    CGAL::VORELL::Medial_axis_location<ET> maloc(e1, e2, triplet.get_bt12().internal_range());
//    cout << ET::to_interval(vc.circle_apx()()) << endl;
}

int main(int argc, char **argv)
{
    std::cin >> e1 >> e2 >> e3;

    CGAL::set_pretty_mode (cerr);
    compute_vc();

    return 0;
}

