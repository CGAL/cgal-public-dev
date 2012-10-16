// (c) 2008 George Tzoumas <geortz@gmail.com>

#include <iostream>

#include <CGAL/basic.h>

#include <CGAL/Gmp_arithmetic_kernel.h>
#include <CGAL/Gmpfi.h>
#include <CGAL/Polynomial.h>

#include <CGAL/Interval_trinomial_solver.h>

typedef CGAL::Gmpfi INT;
typedef CGAL::Polynomial< INT > POLY;

int main(int argc, char **argv)
{
    CGAL::Polynomial_traits_d<POLY>::Construct_polynomial construct;
    CGAL::set_pretty_mode (std::cerr);
    CGAL::set_pretty_mode (std::cout);

    double a = 1, b = 6, c = 6;
    
    if (argc > 1) a = atof(argv[1]);
    if (argc > 2) b = atof(argv[2]);
    if (argc > 3) c = atof(argv[3]);
    if (argc > 4) mpfr_set_default_prec(atoi(argv[4]));

//    std::cerr << CGAL::Gmpfr(0) << std::endl;
//    std::cerr << CGAL::Gmpfr(-0.000000125) << std::endl;
//    std::cerr << CGAL::Gmpfr(-0.00000125) << std::endl;
//    std::cerr << CGAL::Gmpfr(-0.0125) << std::endl;
//    std::cerr << CGAL::Gmpfr(-0.125) << std::endl;
//    std::cerr << CGAL::Gmpfr(-1.25) << std::endl;
//    std::cerr << CGAL::Gmpfr(12.5) << std::endl;
//    std::cerr << CGAL::Gmpfr(125) << std::endl;
//    std::cerr << CGAL::Gmpfr(1250) << std::endl;
//    std::cerr << CGAL::Gmpfr(12500000) << std::endl;
//    std::cerr << CGAL::Gmpfr(125000000) << std::endl;

//    INT range = CGAL::Interval_traits<INT>::Construct()(5.2,5.6);
//    std::cerr << range << std::endl;

    std::vector<INT> p;
    
    p.push_back(INT(c));
    p.push_back(INT(b));
    p.push_back(INT(a));

    POLY poly = construct(p.begin(), p.end());
    std::cerr << poly << std::endl;
    CGAL::Interval_trinomial_solver<POLY, INT> Solver(poly);

    std::pair<INT,INT> sol = Solver();
    
    std::cerr << "status = " << Solver.get_status() << std::endl;
    std::cout << sol.first << width(sol.first) 
              << " -> " << CGAL::singleton(sol.first) << std::endl;
    std::cout << sol.second << width(sol.second) 
              << " -> " << CGAL::singleton(sol.second) << std::endl;
    
    return 0;
}

