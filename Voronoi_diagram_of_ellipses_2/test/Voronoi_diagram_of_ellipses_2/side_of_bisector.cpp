// (c) 2007-2009 George Tzoumas <geortz@gmail.com>

#include <iostream>

#include <CGAL/basic.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpfi.h>

#include <CGAL/Ellipse_traits.h>
#include <CGAL/Ellipse_2.h>

#include <CGAL/Side_of_bisector.h>

#include <CGAL/Timer.h>

typedef CGAL::Gmpfi INT;
typedef CGAL::Gmpq QQ;
typedef CGAL::Ellipse_traits<QQ, INT> ET;
typedef CGAL::Ellipse_2<ET> Ellipse_2;

Ellipse_2 e1, e2;

using namespace std;

void test_k1()
{
    QQ x1 = -5, y1 = -5;
    QQ x2 = 5, y2 = 5;
    QQ xstep = (x2-x1) / 50;
    QQ ystep = (x2-x1) / 50;
    QQ x = x1;
    QQ y = y1;
    CGAL::Timer timer;
    CGAL::Comparison_result r;
    while (x < x2) {
        timer.start();
        r = CGAL::Side_of_bisector<ET>()(e1, e2, x, y);
        timer.stop();
        cout << "Side_of_bisector(e1,e2,(" << x << '=' << CGAL::to_double(x) << 
                ',' << y << '=' << CGAL::to_double(x) << ")) = " << r << endl;
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
    std::cin >> e1 >> e2;

    test_k1();

    return 0;
}

