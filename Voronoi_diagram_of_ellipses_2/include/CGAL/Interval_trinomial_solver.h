//    (c) 2008-2009 National and Kapodistrian University of Athens
//    (c) 2009-2011 INRIA Nancy
//    (c) 2011-2012 National and Kapodistrian University of Athens
//
//    Author:
//        George M. Tzoumas <geortz@gmail.com>
//
//  THIS SOURCE CODE IS CONSIDERED EXPERIMENTAL AND PROVIDED WITHOUT ANY WARRANTY
//  YOU ARE FREE TO USE THIS CODE FOR NON-COMMERCIAL PURPOSES
//  BUT PROPER CREDITS MUST BE GIVEN
//  MORE DETAILS ON THE LICENCE WHEN OFFICIALLY RELEASED

#ifndef CGAL_INTERVAL_TRINOMIAL_SOLVER_H
#define CGAL_INTERVAL_TRINOMIAL_SOVLER_H

#include <CGAL/basic.h>
#include <CGAL/Polynomial.h>

namespace CGAL {

template<class IPOLY, class INT>
class Interval_trinomial_solver {
    IPOLY poly;
    int status;
    
public:
    Interval_trinomial_solver(const IPOLY& poly_): poly(poly_), status(-1) {  }
    std::pair<INT,INT> operator()();
    int get_status() const { return status; }
};

// status:
//     -3 = not a trinomial (lcoeff contains zero)
//     -2 = cannot separate roots, possibly double root
//     -1 = uncertain discriminant sign ==> roots may or may not exist
//      0 = no real roots
//      1 = 1 double root
//      2 = 2 real roots
//TODO: use coercion traits?
template<class IPOLY, class INT> 
std::pair<INT,INT> Interval_trinomial_solver<IPOLY,INT>::operator()()
{
    typedef typename CGAL::Interval_traits<INT>::Bound BT;
    typename CGAL::Polynomial_traits_d<IPOLY>::Get_coefficient get_coeff;
    
//    assert ( degree(p) == 2 );
    CGAL_precondition ( !zero_in(get_coeff(poly, 2)) );

    INT disc = get_coeff(poly, 1)*get_coeff(poly, 1) - 
               INT(4)*get_coeff(poly, 2)*get_coeff(poly, 0);
    if (CGAL::upper(disc) < BT(0)) {
        status = 0;
        return std::pair<INT,INT>();
    }
    if (CGAL::lower(disc) < BT(0)) {
        status = -1;
        disc = INT(BT(0), CGAL::upper(disc));
#if VERBOSE > 0
        std::cerr << "WARNING: TrinomialSolver: discriminant strictly contains zero\n";
#endif
    } else if (CGAL::lower(disc) > BT(0)) {
        status = 2;
    } else if (CGAL::upper(disc) == BT(0)) {
        CGAL_assertion (CGAL::singleton(disc));
        status = 1;
        INT r = (INT(-1)*get_coeff(poly, 1)) / (INT(2)*get_coeff(poly, 2));
        return std::pair<INT,INT>(r, r);
    } else status = -2;
    
    INT sqd = sqrt(disc);
#if VERBOSE > 1
    std::cerr << "d, sqd = " << disc << ' ' << sqd << std::endl;
#endif

    INT r1 = (INT(-1)*get_coeff(poly, 1) - sqd) / (INT(2)*get_coeff(poly, 2));
    INT r2 = (sqd - get_coeff(poly, 1)) / (INT(2)*get_coeff(poly, 2));

    if (status == -2 || CGAL::overlap(r1, r2)) {
#if VERBOSE > 0
        std::cerr << "WARNING: TrinomialSolver: cannot separate roots -- possibly double root\n";
#endif
        if (status == -2) 
            std::cerr << "WARNING: TrinomialSolver: intervals overlap and I knew it!\n";
        if (status == 2) status = -2;
        return std::pair<INT,INT>(r1,r2);
    }
    if (r1 < r2) return std::pair<INT,INT>(r1,r2);
    else return std::pair<INT,INT>(r2,r1);
}

} //namespace CGAL
#endif
