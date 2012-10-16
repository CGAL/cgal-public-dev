// (c) 2009 George Tzoumas <geortz@gmail.com>

#include <iostream>

#include <CGAL/basic.h>

#include <CGAL/Gmpq.h>
#include <CGAL/Gmpfi.h>

#include <CGAL/Ellipse_traits.h>
#include <CGAL/Ellipse_2.h>

#include <CGAL/Ellipse_bisector_2.h>

#include <CGAL/Timer.h>

typedef CGAL::Gmpfi INT;
typedef CGAL::Gmpq QQ;
typedef CGAL::Ellipse_traits<QQ, INT> ET;
typedef CGAL::Ellipse_2<ET> Ellipse_2;

Ellipse_2 e1, e2, e3, e4;

using namespace std;

void test_bisector(int res)
{
    
    CGAL::Timer tm;
    tm.start();
    CGAL::Ellipse_bisector_2<ET> bs12(e1, e2, e3, e4, res);
    tm.stop();
    cerr << "bisector points computed in " << tm.time() << " seconds" << endl;
    e1.draw(cout);
    e2.draw(cout);
    e3.draw(cout);
    e4.draw(cout);
    cout << bs12;
}

int main(int argc, char **argv)
{
    std::cin >> e1 >> e2 >> e3 >> e4;

//    CGAL::set_pretty_mode (cout);
    int res = 128;
    if (argc > 1) res = atoi(argv[1]);
    e1.set_resolution(res/2);
    e2.set_resolution(res/2);
    e3.set_resolution(res/2);
    e4.set_resolution(res/2);
    test_bisector(res);

    return 0;
}

