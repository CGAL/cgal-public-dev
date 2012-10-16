// (c) 2007-2009 George Tzoumas <geortz@gmail.com>

#include <iostream>

#include <CGAL/basic.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpfi.h>

#include <CGAL/Ellipse_traits.h>
#include <CGAL/Ellipse_2.h>

#include <CGAL/Bitangent.h>
#include <CGAL/Visible_arc.h>
#include <CGAL/visibility_apx.h>

#include <CGAL/Timer.h>

typedef CGAL::Gmpq QQ;
typedef CGAL::Gmpfi INT;
typedef CGAL::Ellipse_traits<QQ, INT> ET;
typedef CGAL::Ellipse_2<ET> Ellipse_2;
typedef CGAL::VORELL::Range<ET::Root> Range;
typedef CGAL::VORELL::Range<ET::BT> RangeX;

Ellipse_2 e1, e2, e3;

using namespace std;

//void bench_visibility()
//{
//    CGAL::VORELL::Bitangent<ET> bt12(e1, e2);
//    CGAL::VORELL::Bitangent<ET> bt13(e1, e3);
//    CGAL::VORELL::Range<ET::Root> vis0 = CGAL::VORELL::Visible_arc<ET>(bt12, bt13);
//    cout << "visible arc of e1 from {e2,e3}:" << vis0 << endl;
//    CGAL::VORELL::Range<ET::Root> vis, vis2;
//    QQ xp = e3.get_x_coord(1);
//    QQ yp = e3.get_y_coord(1);
//    CGAL::Timer timer;
//    while (timer.time() < 3) {
//        timer.start();
//        vis = CGAL::VORELL::Visible_arc<ET>(e1, xp, yp);
//        timer.stop();
//    }
//    cout << vis << std::endl;
//    cout << timer.intervals()/timer.time() << " instantiations/sec" << endl;
//    timer.reset();
//    while (timer.time() < 3) {
//        timer.start();
//        vis2 = CGAL::VORELL::Visible_arc<ET>(e1, e3, 1);
//        timer.stop();
//    }
//    cout << vis2 << std::endl;
//    cout << timer.intervals()/timer.time() << " instantiations/sec" << endl;
//}

void test_visibility()
{
    CGAL::VORELL::Bitangent<ET> bt12(e1, e2);
    CGAL::VORELL::Bitangent<ET> bt13(e1, e3);
    Range vis0 = CGAL::VORELL::Visible_arc<ET>(bt12, bt13);
    cout << "visible arc of e1 from {e2,e3}: " << vis0 << endl;
    Range vis, vis2;
    QQ xp = e3.boundary_x(1);
    QQ yp = e3.boundary_y(1);
    vis = CGAL::VORELL::Visible_arc<ET>(e1, xp, yp);
    cout << "visible arc of e1 from e3(1)  : " << vis << endl;
    RangeX apoll = CGAL::VORELL::Apollonius_arc_apx<ET>(e1, e3, 1);
    cout << "apollonius arc of e1 from e3(1) = " << apoll << endl;
}

int main(int argc, char **argv)
{
    std::cin >> e1 >> e2 >> e3;

//    CGAL::set_pretty_mode (cerr);
    test_visibility();

    return 0;
}

