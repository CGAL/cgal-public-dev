// (c) 2008-2009 George Tzoumas <geortz@gmail.com>

#include <iostream>
#include <CGAL/basic.h>

#include <CGAL/Gmp_arithmetic_kernel.h>
#include <CGAL/Gmpfr.h>
//#include <CGAL/Boost_interval_Gmpfr.h>
#include <CGAL/Gmpfi.h>
#include <CGAL/Polynomial.h>

#include <CGAL/Interval_newton_solver.h>

#include <CGAL/Timer.h>

typedef CGAL::Gmpz ZZ;
//typedef CGAL::Boost_interval_Gmpfr INT;
typedef CGAL::Gmpfi INT;
typedef CGAL::Polynomial< INT > POLY;

// -1.309079643, -0.7847023334, 1.206303566

CGAL::Timer timer;

void demo(double a, double b) {
    std::vector<INT> p;
    p.reserve(5);
    p.push_back(INT(74));
    p.push_back(INT(6));
    p.push_back(INT(-92));
    p.push_back(INT(75));
    p.push_back(INT(23));
    p.push_back(INT(-50));
    
//    p.push_back(74);
//    p.push_back(6);
//    p.push_back(-92);
//    p.push_back(75);
//    p.push_back(23);
//    p.push_back(-50);
    
    POLY poly(p.begin(), p.end());

//    std::cerr << poly << std::endl;
//    std::cerr << poly.degree() << std::endl;
    INT range(a,b);
    CGAL::Interval_newton_solver<POLY, INT> 
            Solver(poly, CGAL::differentiate(poly), range);
    timer.start();
    INT sol = Solver();
    timer.stop();
    
    if (CGAL::Gmpfr::get_default_precision() < 256)
        std::cout << sol << std::endl;
    std::cout << CGAL::get_significant_bits(sol) << std::endl;
    std::cerr << "status = " << Solver.get_status() << ", iters = "
              << Solver.get_iters() << std::endl;
}

int main(int argc, char **argv)
{
    CGAL::set_pretty_mode (std::cerr);
    CGAL::set_pretty_mode (std::cout);

    if (argc > 1) CGAL::Gmpfr::set_default_precision(atoi(argv[1]));
    double a = 1.0, b = 2.0;
    if (argc > 2) a = atof(argv[2]);
    if (argc > 3) b = atof(argv[3]);
//    std::cerr << a << ' ' << b << std::endl;
//    std::cerr << "benchmarking newton for 2 sec ...";
    demo(a,b);
//    while (timer.time() < 2) demo(a, b);
    std::cerr << std::endl;
//    std::cerr << "timeprec = " << timer.precision() << std::endl;
//    std::cerr << "reps = " << timer.intervals() << std::endl;
//    std::cerr << "time = " << timer.time() << " time/rep = " << timer.time()/timer.intervals() << std::endl;
    std::cerr << CGAL::Gmpfr::get_default_precision() << '\t'
              << timer.time()/timer.intervals() << std::endl;
    
    return 0;
}

