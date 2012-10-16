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

#ifndef CGAL_INTERVAL_NEWTON_SOLVER_H
#define CGAL_INTERVAL_NEWTON_SOLVER_H

#include <CGAL/basic.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Bigfloat_interval_traits.h>

namespace CGAL {

template<class IPOLY, class INT>
class Interval_newton_solver {
    typedef typename Interval_traits<INT>::Bound BT;
    typename Interval_traits<INT>::Construct construct;
    typename Polynomial_traits_d<IPOLY>::Evaluate evaluate;
    
    IPOLY poly, dpoly;
    INT approx;
    int iters, status;
    
    INT newton_iteration() {
        BT x = median(approx);
        iters++;
        return INT(x) - evaluate(poly, INT(x))/evaluate(dpoly, approx);
    } 
    
    INT bisect_iteration();
    
public:

    Interval_newton_solver(const IPOLY& poly_, const IPOLY& dpoly_, 
        const INT& init_approx): poly(poly_), dpoly(dpoly_),  approx(init_approx) {
        status = -1; iters = 0;
#if VERBOSE > 1
        std::cerr << "approx = " << approx << std::endl;
#endif
    }

    void restart(const INT& init_approx) { 
        approx = init_approx; status = -1; iters = 0;
    }
    INT operator()(BT delta = BT(0));
    
    int get_status() const { return status; }
    int get_iters() const { return iters; }
};


template<class IPOLY, class INT>
INT Interval_newton_solver<IPOLY, INT>::bisect_iteration() 
{
    typedef typename Interval_traits<INT>::Bound BT;
    typename Polynomial_traits_d<IPOLY>::Evaluate evaluate;

    BT x = CGAL::median(approx);
    INT lval, rval, mval;

    iters++;
    lval = evaluate(poly, INT(lower(approx)));
    rval = evaluate(poly, INT(upper(approx)));
    mval = evaluate(poly, INT(x));
#if VERBOSE > 1
        std::cerr << "approx = " << approx << std::endl;
        std::cerr << "lval = " << lval << " rval =  " << rval << " mval = " << mval << " x = " << x << std::endl;
#endif
    // BUG: check singleton here?
    if (CGAL::zero_in(lval) || CGAL::zero_in(rval)) {
        status = -3;
        return approx;
    }
    if (CGAL::sign(lval) == CGAL::sign(rval)) { // check for different sign at endpoints; TODO: just evaluate in approx?
        status = -1;
        return approx;
    }

    status = 0;
#if VERBOSE > 1
        std::cerr << "checking zeros" << std::endl;
#endif
    if (CGAL::singleton(mval) && CGAL::sign(mval) == ZERO)
        return INT(x);
    if (CGAL::singleton(lval) && CGAL::sign(lval) == ZERO)
        return CGAL::lower(approx);
    if (CGAL::singleton(rval) && CGAL::sign(rval) == ZERO)
        return CGAL::upper(approx);
#if VERBOSE > 1
        std::cerr << "done checking zeros" << std::endl;
#endif
    if (CGAL::zero_in(mval)) return approx;

//    if (CGAL::sign(lval) == CGAL::sign(mval)) return INT(x,CGAL::upper(approx));
//    else return INT(CGAL::lower(approx),x);
    if (CGAL::sign(lval) == CGAL::sign(mval)) return construct(x, CGAL::upper(approx));
    else return construct(CGAL::lower(approx), x);
}

// status:
//     -3 = precision not enough
//     -2 = no root exists (certainly)
//     -1 = no root found (cannot guarantee existence)
//      0 = some root found (cannot guarantee uniqueness)
//      1 = unique root found (certainly)
// TODO:: caution of coeffs including zero !
template<class IPOLY, class INT>
INT Interval_newton_solver<IPOLY, INT>::operator()(
        typename Interval_traits<INT>::Bound delta)
{
    typedef typename Interval_traits<INT>::Bound BT;
    typename Polynomial_traits_d<IPOLY>::Evaluate evaluate;

    INT slope, testapprox, newapprox;
    BT oldw = CGAL::width(approx);
    BT neww;
    bool z, init_z;
    INT lval, rval;

    lval = evaluate(poly, INT(CGAL::lower(approx)));
    rval = evaluate(poly, INT(CGAL::upper(approx)));
#if VERBOSE > 1
        std::cerr << "approx = " << approx << std::endl;
        std::cerr << "lval = " << lval << " rval =  " << rval << std::endl;
#endif
    if (CGAL::zero_in(lval) || CGAL::zero_in(rval)) {
        status = -3;
        return approx;
    }

    status = 0;
    iters = 0;
    slope = evaluate(dpoly, approx);
#if VERBOSE > 1
    std::cerr << "slope = " << slope << CGAL::width(slope) << std::endl;
#endif
    init_z = z = CGAL::zero_in(slope);
#if VERBOSE > 1
    std::cerr << "sign(lval) = " << CGAL::sign(lval) << " sign(rval) = " << CGAL::sign(rval) << std::endl;
#endif
    if (CGAL::sign(lval) == CGAL::sign(rval)) {
        if (z) status = -1; else status = -2;
        return approx;
    }

#if VERBOSE > 1
        int i = 1;
#endif

    int old_prec = -1;
    while (1) {
        if (z) {
            approx = bisect_iteration();
            if (status < 0) return approx;
            neww = CGAL::width(approx);
#if VERBOSE > 1
        std::cerr << '(' << i << ") Bisect: " << approx << ' ' << neww << std::endl;
#endif
            if (neww < delta) break;

            // check convergence
            CGAL_assertion (neww <= oldw);
            if (neww == oldw) break;

            oldw = neww;
        } else {
#if VERBOSE > 1
        std::cerr << "slope = " << slope << std::endl;
#endif
            testapprox = newton_iteration();
            newapprox = CGAL::intersection(testapprox, approx);

            // check proper convergence
//            CGAL_assertion ( ! empty(newapprox) );

#if VERBOSE > 1
        std::cerr <<"new =" << newapprox << " test = " << testapprox << " old = " << approx << std::endl;
#endif

            // new approx should contain the root
            CGAL_assertion ( CGAL::zero_in(evaluate(poly, newapprox)) );

            int new_prec = CGAL::get_significant_bits(newapprox);
            if (new_prec > 1) {
                if (new_prec > old_prec) old_prec = new_prec; else break;
            }
            neww = CGAL::width(newapprox);
#if VERBOSE > 1
        std::cerr << '(' << i << ") Newton: " << newapprox << ' ' << neww << std::endl;
#endif
            if (neww >= oldw) break;
//                if (BT(2)*neww >= oldw) break;
            approx = newapprox;
            if (neww < delta) break;
            oldw = neww;
        }
        slope = evaluate(dpoly, approx);
#if VERBOSE > 1
        std::cerr << "slope = " << slope << CGAL::width(slope) << std::endl;
#endif
        z = zero_in(slope);
#if VERBOSE > 1
        i++;
#endif
    }

    if (!init_z) status = 1;
    return approx;
}

} //namespace CGAL
#endif
