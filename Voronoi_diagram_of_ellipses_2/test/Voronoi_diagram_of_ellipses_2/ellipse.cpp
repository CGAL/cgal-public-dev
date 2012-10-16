#include <iostream>
#include <cstring>

#include <CGAL/basic.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpfi.h>

#include <CGAL/Ellipse_traits.h>
#include <CGAL/Ellipse_2.h>

typedef CGAL::Gmpfi INT;
typedef CGAL::Gmpz ZZ;
typedef CGAL::Gmpq QQ;
typedef CGAL::Ellipse_traits<QQ, INT> ET;
typedef CGAL::Ellipse_2<ET> Ellipse_2;

using namespace std;

template< class Stream >
void sample(Stream& w, Ellipse_2& e, int res = 16, int poly = false)
{
    QQ xc = e.x_center();
    QQ yc = e.y_center();

    QQ tstep = QQ(4,res);

    if (poly) w << res << std::endl;

    for (QQ t = -1; t < 1; t += tstep) {
        w << e.boundary_x(t) << ' ' << e.boundary_y(t);
        w << std::endl;
    }
    // symmetric part
    for (QQ t = -1; t < 1; t += tstep) {
        w << (2*e.x_center() - e.boundary_x(t)) << ' ' << (2*e.y_center() - e.boundary_y(t));
        w << std::endl;
    }
    if (poly) w << std::endl;
}


void test_ellipse(Ellipse_2& e)
{
    cout << "(a,b,w,xc,yc) = (" << e.parametric_coefficient(0) << ',' 
         << e.parametric_coefficient(1) << ',' 
         << e.parametric_coefficient(2) << ',' 
         << e.parametric_coefficient(3) << ',' 
         << e.parametric_coefficient(4) << ')' << endl;
    cout << "equation = ";
    e.print_equation(cout);
    cout << endl;
    cout << "is_circle = " << e.is_circle() << endl;
    cout << "'bottom' (x,y) = (" << e.boundary_x(-1) << ',' 
         << e.boundary_y(-1) << ')' << endl;
    cout << "'right' (x,y) = (" << e.boundary_x(0) << ',' 
         << e.boundary_y(0) << ')' << endl;
    cout << "'top' (x,y) = (" << e.boundary_x(1) << ',' 
         << e.boundary_y(1) << ')' << endl;
    cout << "'left' (x,y) = (" << e.boundary_x_inf() << ','
         << e.boundary_y_inf() << ')' << endl;
}

int bitsize(const QQ &a)
{
    ZZ num(0), denom(0);
    num = CGAL::abs(a.numerator());
    denom = CGAL::abs(a.denominator());
    ZZ max = CGAL::max(num, denom);
    return max.bit_size();
}

int bitsize(const Ellipse_2 &a)
{
    int bmax = 0;
    for (int i = 0; i < 5; i++) {
        bmax = CGAL::max(bmax, bitsize(a.parametric_coefficient(i)));
    }
    return bmax;
}

QQ perturb(const QQ &a, int p)
{
    ZZ ten(10);
    QQ delta = QQ(1)/QQ(CGAL::ipower(ten, p));
    if (lrand48()%2 == 1) delta = -delta;
    return a + delta;
}

Ellipse_2 perturb_ellipse(const Ellipse_2 &a, int p)
{
    QQ pc[5];
    for (int i = 0; i < 5; i++) {
        pc[i] = perturb(a.parametric_coefficient(i), p);
    }
    return Ellipse_2(pc[0], pc[1], pc[2], pc[3], pc[4]);
}

int main(int argc, char **argv)
{
    if (argc == 1) {
        cout << "usage:" << argv[0] << " info | plot [res] | [pt]bitsize | perturb [pow] | sample [res]" << endl;
        return 0;
    }

    if (strcmp(argv[1], "info") == 0) {
        Ellipse_2 e;
        std::cin >> e;
        test_ellipse(e);
        cout << "moving to origin..." << endl;
        e.translate(0,0);
        test_ellipse(e);
//        Ellipse_2 e
    } else if (strcmp(argv[1], "plot") == 0) {
        Ellipse_2 e;
        std::cin >> e;
        int res = 48;
        if (argc > 2) res = atoi(argv[2]);
        e.set_resolution(res);
        cout << e;
    } else if (strcmp(argv[1], "bitsize") == 0) {
        int count = 0;
        int bc = 0;
        int bmin = 1<<30;
        int bmax = 0;
        while (cin) {
            Ellipse_2 e;
            cin >> e;
            if (cin) {
                count++;
                int b = bitsize(e);
                bc += b;
                bmin = CGAL::min(b, bmin);
                bmax = CGAL::max(b, bmax);
            }
        }
        double bs = 0.0;
        if (count) bs = bc*1.0 / count;
        cout << "ellipses: " << count << " bitsize: " << bs
             << " min: " << bmin << " max: " << bmax << endl;
    } else if (strcmp(argv[1], "ptbitsize") == 0) {
        int count = 0;
        int bc = 0;
        int bmin = 1<<30;
        int bmax = 0;
        while (cin) {
            QQ x;
            cin >> x;
            if (cin) {
                count++;
                int b = bitsize(x);
                bc += b;
                bmin = CGAL::min(b, bmin);
                bmax = CGAL::max(b, bmax);
            }
        }
        double bs = 0.0;
        if (count) bs = bc*1.0 / count;
        cout << "pts: " << count << " bitsize: " << bs
             << " min: " << bmin << " max: " << bmax << endl;
    } else if (strcmp(argv[1], "sample") == 0) {
        int res = 16;
        if (argc > 2) res = atoi(argv[2]);
        while (cin) {
            Ellipse_2 e;
            cin >> e;
            if (cin) {
                sample(cout, e, res);
            }
        }
    } else if (strcmp(argv[1], "polysample") == 0) {
        int res = 16;
        if (argc > 2) res = atoi(argv[2]);
        while (cin) {
            Ellipse_2 e;
            cin >> e;
            if (cin) {
                sample(cout, e, res, true);
            }
        }
    } else if (strcmp(argv[1], "perturb") == 0) {
        int per = 2;
        if (argc > 2) per = atoi(argv[2]);
        while (cin) {
            Ellipse_2 e;
            cin >> e;
            if (cin) {
                cout << perturb_ellipse(e, per) << endl;
            }
        }
    }
    return 0;
}

