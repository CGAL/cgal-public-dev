// (c) 2007-2008 George Tzoumas <geortz@gmail.com>

#include <iostream>

#include <CGAL/basic.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpfi.h>

#include <CGAL/Ellipse_traits.h>
#include <CGAL/Ellipse_2.h>

#include <CGAL/Timer.h>

//#define VERBOSE

#include <CGAL/In_circle.h>

#include <cstring>

typedef CGAL::Gmpq QQ;
typedef CGAL::Gmpfi INT;
typedef CGAL::Ellipse_traits<QQ, INT> ET;
typedef CGAL::Ellipse_2<ET> Ellipse_2;

Ellipse_2 e1, e2, e3, e4;

using namespace std;

void test_k3()
{
    QQ x1 = -12, y1 = 12;
    QQ x2 = 0, y2 = 0;
//    QQ x1 = QQ(-276,25), y1 = QQ(276,25);
//    QQ x2 = x1, y2 = y1;
    QQ xstep = (x2-x1) / 50;
    QQ ystep = (y2-y1) / 50;
    QQ x = x1;
    QQ y = y1;
    CGAL::Timer timer;
    CGAL::Comparison_result r;
    while (x <= x2) {
        timer.start();
        e4.translate(x, y);
        r = CGAL::In_circle<ET>()(e1, e2, e3, e4);
        timer.stop();
        cout << "InCircle(e1,e2,e3,e4(" << x << '=' 
             << CGAL::to_double(x) << ',' << y << '=' 
             << CGAL::to_double(y) << ")) = " << r << endl;
        x += xstep;
        y += ystep;
        if (x1 == x2 && y1 == y2) break;
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
    std::cin >> e1 >> e2 >> e3 >> e4;

    
//    CGAL::set_pretty_mode (cerr);
//    test_k3();
    if (argc > 1) {
        if (strcmp(argv[1], "slide") == 0) {
            test_k3();
            return 0;
        } else CGAL::Gmpfr::set_default_precision(atoi(argv[1]));
    }
    cout << CGAL::In_circle<ET>()(e1, e2, e3, e4) << endl;

    return 0;
}

