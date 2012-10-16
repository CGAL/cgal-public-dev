// (c) 2007-2009 George Tzoumas <geortz@gmail.com>

#include <iostream>

#include <CGAL/basic.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpfi.h>

#include <CGAL/Ellipse_traits.h>
#include <CGAL/Ellipse_2.h>

#include <CGAL/Distance_from_bitangent.h>

#include <CGAL/Timer.h>

typedef CGAL::Gmpq QQ;
typedef CGAL::Gmpfi INT;
typedef CGAL::Ellipse_traits<QQ, INT> ET;
typedef CGAL::Ellipse_2<ET> Ellipse_2;

Ellipse_2 e1, e2, e3;

using namespace std;

void test_k2()
{
    QQ x1 = -8, y1 = -8;
    QQ x2 = 6, y2 = 1;
    QQ xstep = (x2-x1) / 50;
    QQ ystep = (x2-x1) / 50;
    QQ x = x1;
    QQ y = y1;
    CGAL::Timer timer;
    CGAL::Comparison_result r;
    while (x < x2) {
        timer.start();
        e3.translate(x, y);
        r = CGAL::Distance_from_bitangent<ET>()(e1, e2, e3);
        timer.stop();
        cout << "Distance_from_bitangent(e1,e2,e3(" << x << '=' << 
                CGAL::to_double(x) << ',' << y << '=' << 
                CGAL::to_double(x) << ")) = " << r << endl;
        x += xstep;
        y += ystep;
    }
    cerr << "total time = " << timer.time() << endl;
    cerr << "total iters = " << timer.intervals() << endl;
    cerr << "avg time = " << timer.time()/timer.intervals() << endl;
//    for (int i = 0; i < 1; i++) {
//        timer.start();
//        CGAL::SideOfBisector<Ellipse_2>()(e1, e2, x, y);
//        timer.stop();
//    }
}

int main(int argc, char **argv)
{
    std::cin >> e1 >> e2 >> e3;

    test_k2();

    return 0;
}

