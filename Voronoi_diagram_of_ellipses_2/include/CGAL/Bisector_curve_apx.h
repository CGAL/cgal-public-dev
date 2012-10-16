//    (c) 2007-2009 National and Kapodistrian University of Athens
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


#ifndef CGAL_BISECTOR_CURVE_APX_H
#define CGAL_BISECTOR_CURVE_APX_H

#include<iostream>
#include<algorithm>

#include<CGAL/basic.h>
#include<CGAL/vorell.h>

#include<CGAL/Range.h>
#include<CGAL/Ellipse_2.h>
#include<CGAL/Ellipse_triplet.h>
#include<CGAL/visibility_apx.h>

#include <CGAL/Interval_newton_solver.h>

#include <CGAL/polynomial_utils.h>

namespace CGAL {

namespace VORELL {
    
template<class ET>
class Bisector_curve_apx {
    typedef typename ET::IT IT;
    typedef typename ET::BT BT;
    typedef typename ET::QQ QQ;
    typedef typename VORELL::Apollonius_arc_apx<ET> Apollonius_arc_apx;
    typedef typename VORELL::Bitangent<ET> Bitangent;
    typedef typename ET::ARX ARX;
    typedef typename ET::bpolz_t bpolz_t;
    typedef typename ET::bpolq_t bpolq_t;
    typedef typename ET::bpoli_t bpoli_t;
    typedef typename ET::upoli_t upoli_t;
    typedef typename ET::upolz_t upolz_t;
    typedef typename ET::Root Root;

    typename CGAL::Interval_traits<IT>::Construct make_interval;
    
    bpoli_t bx12;
    bpolq_t bs12;
    Bitangent bt12, bt21;
    int status;
    Ellipse_2<ET> e1, e2;
    ARX inrange1, inrange2;

    IT solve(IT t1, bool have_guess_interval, ARX guess_);
    IT solve_CC(IT t1);

    void init(const Ellipse_2<ET>& e1_, const Ellipse_2<ET>& e2_) {
        e1 = e1_; e2 = e2_;
        bpolz_t bs12z = CGAL::VORELL::Generated_code<ET>().bisector(
                       ELL_PARAM_COEFFS(e2), ELL_PARAM_COEFFS(e1)
                   );
        bx12 = typename 
               CGAL::Coercion_traits<bpoli_t, bpolz_t>::Cast()(bs12z);
        bs12 = typename 
               CGAL::Coercion_traits<bpolq_t, bpolz_t>::Cast()(bs12z);
        status = -100; // uninitialized
        if (bt12.relative_position() == VORELL::PSEUDO_CIRCLES) {
            inrange1 = ET::to_arx(bt12.internal_range(), false, 53); // TODO: adaptive prec.
            inrange2 = ET::to_arx(bt21.internal_range(), false, 53);
#if VERBOSE > 2
            std::cerr << "inrange1 = " << inrange1 << std::endl;
            std::cerr << "inrange2 = " << inrange2 << std::endl;
#endif
        }
    }

    // less_than
    static bool compare_circle_radii(const IT& xt, const IT& yt,
                                     const IT& xc1, const IT& yc1,
                                     const IT& xc2, const IT& yc2) {
        if (!CGAL::zero_in(xc1-xt)) {
            if (CGAL::zero_in(xc2-xc1)) return false;
            if (xc1 > xt) return xc1 < xc2; else return xc1 > xc2;
        } else if (!CGAL::zero_in(yc1-yt)) {
            if (CGAL::zero_in(yc2-yc1)) return false;
            if (yc1 > yt) return yc1 < yc2; else return yc1 > yc2;
//        } else if (!CGAL::zero_in(xc2-xt)) {
//            if (CGAL::zero_in(xc2-xc1)) return false;
//            if (xc2 > xt) return xc1 < xc2; else return xc1 > xc2;
//        } else if (!CGAL::zero_in(yc2-yt)) {
//            if (CGAL::zero_in(yc2-yc1)) return false;
//            if (yc2 > yt) return yc1 < yc2; else return yc1 > yc2;
        } else {
            CGAL_assertion(false); return true;
//            return false;
        }
    }

    struct compare_bitangent_circle_radii {
        IT t;
        QQ tq;
        IT xt, yt;
        const Bisector_curve_apx *bp;
        compare_bitangent_circle_radii(const IT& t_, const QQ& tq_, const Ellipse_2<ET>& e1,
                           const Bisector_curve_apx<ET>* bp_): t(t_), tq(tq_), bp(bp_) {
          // TODO: user lower()/upper()
            xt = ET::to_interval(e1.boundary_x(tq));
            yt = ET::to_interval(e1.boundary_y(tq));
        }

        bool operator()(const IT& r1, const IT& r2) const {
            IT xc1, yc1, xc2, yc2;
            int s = bp->get_coords(t, r1, xc1, yc1);
            CGAL_assertion( s == 0 );
            s = bp->get_coords(t, r2, xc2, yc2);
            CGAL_assertion( s == 0 );
            return compare_circle_radii(xt, yt, xc1, yc1, xc2, yc2);
        }
    };

    void process_bitangent_circles(const IT& t1, const ARX& apoll, IT& r2i);
    void process_internal_bitangent_circles(
                                   const IT& t1, const ARX& apoll, IT& r2i);
public:
    Bisector_curve_apx() { }

    Bisector_curve_apx(const Ellipse_2<ET>& e1_, const Ellipse_2<ET>& e2_) {
        bt12 = Bitangent(e1_, e2_);
        bt21 = Bitangent(e2_, e1_);
        init(e1_, e2_);
    }

    Bisector_curve_apx(const Ellipse_2<ET>& e1_, const Ellipse_2<ET>& e2_,
                       const Bitangent& bt12_, const Bitangent& bt21_) {
        bt12 = bt12_;
        bt21 = bt21_;
        init(e1_, e2_);
    }
    
    int get_status() const { return status; }
    
    inline IT operator()(IT t1, bool have_guess_interval = false, 
                         ARX guess_ = ARX()) {
        if (e1.is_circle() && e2.is_circle()) return solve_CC(t1);
        return solve(t1, have_guess_interval, guess_);
    }
    
    // ret = 0:  OK
    //      -2:  denom. is zero 
    int get_coords(const IT t, const IT r, IT& x, IT& y) const;
    
    const Ellipse_2<ET>& get_e1() const { return e1; }
    const Ellipse_2<ET>& get_e2() const { return e2; }
};

template<class ET>
typename ET::IT Bisector_curve_apx<ET>::solve_CC(IT t1) 
{
    IT dy = ET::to_interval(e2.y_center() - e1.y_center());
    IT dx = ET::to_interval(e2.x_center() - e1.x_center());
    IT dr = ET::to_interval(e2.major_axis() - e1.major_axis());

    IT num = (-dy)*t1 - (dx + dr);
    IT den = (dr - dx)*t1 + dy;
    if (CGAL::zero_in(den)) {
        if (CGAL::zero_in(num)) {
            status = -1;
#if VERBOSE > 0 
            std::cerr << "bad denom\n";
            std::cerr << "t = " << t1 << " dy = " << dy << 
                                         " dx = " << dx << " dr = " << 
                   dr << " num = " << num << " den = " << den << std::endl;
            std::cerr << "dr - dx = " << (dr-dx) << std::endl;
#endif
            return IT(10000);
        } else {
#if VERBOSE > 0 
            std::cerr << "i-point CONFLICT!\n";
            std::cerr << "t = " << t1 << " dy = " << dy 
                      << " dx = " << dx << " dr = " << dr 
                      << " num = " << num << " den = " << den << std::endl;
            std::cerr << "dr - dx = " << (dr-dx) << std::endl;
#endif
            status = 0;
            return num/upper(den);
        }
    } else {
        status = 0;
        return num/den;
    }
}

//#ifdef VERBOSE
//#define OLDVERB VERBOSE
//#undef VERBOSE
//#define VERBOSE 3
//#endif
// fallback, when Newton fails
// compute ALL bitangent circles using the algebraic method
// TODO: solve in rational interval?
template<class ET>
void Bisector_curve_apx<ET>::process_bitangent_circles(
                     const IT& t1, const ARX& apoll, IT& r2i) 
{ 
    QQ tq = ET::to_rational(CGAL::median(t1));
#if VERBOSE > 2
    std::cerr << "PROCESS BITANG. CIRCLE t = " << CGAL::median(t1) << std::endl;
#endif
    upolz_t polyq = VORELL::primpart(
                      typename ET::BPTQ::Evaluate()(bs12, tq));
    // TODO: list could be used instead of vector
    std::vector<Root> sol;
    std::vector<IT> candidates;
    sol = ET::Real_roots(polyq);
    typename std::vector<typename ET::Root>::iterator it;

#if VERBOSE > 2
    std::cerr << "r-sols: ";
#endif
    for (it = sol.begin(); it != sol.end(); it++) {
#if VERBOSE > 2
        std::cerr << CGAL::to_double(*it);
#endif
        IT tmpr = ET::to_interval(*it, 53); // TODO: adaptive prec.
        bool qualifies = apoll.contains(lower(tmpr)) &&  // TODO: relaxed (||) inclusion?
                         apoll.contains(upper(tmpr));
        if (qualifies) {
            candidates.push_back(tmpr);
        }
#if VERBOSE > 2
        if (qualifies) std::cerr << "* ";  else std::cerr << ' ';
#endif
    }
#if VERBOSE > 2
    std::cerr << std::endl;
#endif
    if (candidates.size() == 0) {
        status = -105;
        return;
    }
    if (candidates.size() == 1) {
        status = 1;
//        std::cerr << "front = " << candidates.front() << std::endl;
        r2i = candidates.front();
    } else {
        // pick smallest bitangent circle
        r2i = *std::min_element(candidates.begin(), 
                 candidates.end(), 
                 compare_bitangent_circle_radii(t1, tq, e1, this));
#if VERBOSE > 2
        std::cerr << "picked" << r2i << std::endl;
#endif
        status = 1;
    }
}

//#ifdef OLDVERB
//#undef VERBOSE
//#define VERBOSE OLDVERB
//#undef OLDVERB
//#endif


template<class ET>
void Bisector_curve_apx<ET>::process_internal_bitangent_circles(
                                     const IT& t1, const ARX& apoll, IT& r2i) 
{ 
    QQ tq = ET::to_rational(CGAL::median(t1));
    // compare against self-bitangent circle (medial axis)
    BT t1s;
    // find symmetric point (wrt medial axis)
    if (e1.major_axis() >= e1.minor_axis()) {
        t1s = -CGAL::median(t1);
    } else {
        if (zero_in(t1)) t1s = CGAL::median(t1);
        t1s = BT(1)/CGAL::median(t1);
    }
    if (inrange1.contains(t1s)) {
        IT xt = ET::to_interval(e1.boundary_x(tq));  
        IT yt = ET::to_interval(e1.boundary_y(tq));
        IT xcm = ET::to_interval(e1.medial_x(tq));  
        IT ycm = ET::to_interval(e1.medial_y(tq));
        IT xcr, ycr;
        int s = get_coords(t1, r2i, xcr, ycr);
        //CGAL_assertion( s == 0 );
        if (s != 0) {
#if VERBOSE > 2
            std::cerr << "failed at point t1 = " << t1
                        << " r2i = " << r2i << std::endl;
            std::cerr << "xt = " << xt << " yt = " << yt << std::endl;
            std::cerr << "xcm = " << xcm << " ycm = " << ycm << std::endl;
#endif
            status = -105;
        } else {
            if (compare_circle_radii(xt, yt, xcm, ycm, xcr, ycr)) {
#if VERBOSE > 2
                std::cerr << "rejected " << t1 << "," << r2i << std::endl;
#endif
                status = -104;
            }
        }
    }
}

// TODO: exploit internal bitangent range
//       improve Visible_arc/Apollonius_arc using all BT arcs
template<class ET>
typename ET::IT Bisector_curve_apx<ET>::solve(IT t1, bool have_guess_interval, 
                                     ARX guess_)
{
#if VERBOSE > 2
    std::cerr << "t1 = " << t1 << std::endl;
#endif
    upoli_t poly = typename ET::BPTI::Evaluate()(bx12, t1);
#if VERBOSE > 2
    std::cerr << "poly = " << poly << std::endl;
    std::cerr << "cauchy bound = " << 
            CGAL::VORELL::cauchy_bound<upoli_t, IT>(poly) << std::endl;
#endif
    IT r2i, guess;
    status = -100; // uninitialized
    if (have_guess_interval) {
        guess = make_interval(guess_.left(), guess_.right());
#if VERBOSE > 2
        std::cerr << "guess = " << guess_ << std::endl;
#endif
        CGAL::Interval_newton_solver<upoli_t, IT> 
                Solver(poly, ::CGAL::differentiate(poly), guess);
        r2i = Solver();
        status = Solver.get_status();
        // TODO: fallback, but t1 should not be interval (=wide)
        // ==> valid for 1st bounce of the subdivision
        // if (status != 1) process_bitangent_circles(t1, guess_, r2i);
    } else {
        ARX apoll;
        if (!inrange1.is_empty()) {
            if (CGAL::in(inrange1.left(),t1)) return inrange1.left();
            if (CGAL::in(inrange1.right(),t1)) return inrange1.right();
        }
        bool inside = inrange1.contains(CGAL::median(t1));
        if (inside) apoll = inrange2;
        else apoll = Apollonius_arc_apx (e2, e1, t1);
//        if (guess_.is_empty()) apoll = Apollonius_arc_apx (e2, e1, t1);
//        else apoll = guess_;
        if (apoll.is_empty()) {
            status = -102;
            return r2i;
        }
        //CGAL_assertion( !apoll.is_empty() );
#if VERBOSE > 2
        std::cerr << "apoll = " << apoll << std::endl;
#endif

        BT cb = CGAL::VORELL::cauchy_bound<upoli_t, IT>(poly);
        if (apoll.is_finite()) {
            ARX cauchy(-cb, cb);
            ARX apoll2 = apoll.intersection(cauchy);
            guess = make_interval(apoll2.left(), apoll2.right());
#if VERBOSE > 2
            std::cerr << "apoll (cauchy) = " << guess << std::endl;
#endif
            CGAL::Interval_newton_solver<upoli_t, IT> 
                    Solver(poly, differentiate(poly), guess);
            r2i = Solver();
            status = Solver.get_status();
            // TODO: fallback, but t1 should not be interval (=wide)
            // ==> valid for 1st bounce of the subdivision
            // if (status != 1) process_bitangent_circles(t1, apoll, r2i);
            if (inside) {
                if (status != 1) process_bitangent_circles(t1, apoll, r2i);
                process_internal_bitangent_circles(t1, apoll, r2i);
            }
        } else {
            if (cb < apoll.left() && (-cb) > apoll.right() ) {
                status = -103;
                return r2i;
            }
            //CGAL_assertion (cb > apoll.left() || (-cb) < apoll.right() );
            if (cb > apoll.left() ) {
                guess = make_interval(apoll.left(), cb);
#if VERBOSE > 2
                std::cerr << "apoll (cauchy) = " << guess << std::endl;
#endif
                CGAL::Interval_newton_solver<upoli_t, IT> 
                        Solver(poly, differentiate(poly), guess);
                r2i = Solver();
                status = Solver.get_status();
            }
            // TODO: status < 1?
            if (status < 0 && (-cb) < apoll.right()) {
                guess = make_interval(-cb, apoll.right());
#if VERBOSE > 2
                std::cerr << "apoll (cauchy) = " << guess << std::endl;
#endif
                CGAL::Interval_newton_solver<upoli_t, IT> 
                        Solver(poly, differentiate(poly), guess);
                r2i = Solver();
                status = Solver.get_status();
            }
            // TODO: fallback, but t1 should not be interval (=wide)
            // ==> valid for 1st bounce of the subdivision
            //if (status != 1) process_bitangent_circles(t1, 
                        //apoll, r2i);
            if (inside) {
                if (status != 1) process_bitangent_circles(t1, apoll, r2i);
                process_internal_bitangent_circles(t1, apoll, r2i);
            }
        }
    }
#if VERBOSE > 1
    std::cerr << "r2i = " << r2i << CGAL::width(r2i) << " ("
              << status << ')' << std::endl;
#endif
    return r2i;
}


// cstatus = 0:  OK
//          -2:  denom. is zero 
//          -3:  wrong part of bisector computed
template<class ET>
int Bisector_curve_apx<ET>::get_coords(const IT t, 
        const IT r, IT& x, IT& y) const
{
    IT t15, t11, t94, t118, t10, t122, t12, t8, t97, t101, t66, t9, t77, 
       t110, t128, t90, t6, t72, t47, t123, t85, t23, t126, t89, t124,
       t96, t109, t75, t105, t88, t104, t7, t92, t51, t86, t84, t79,
       t71, t69, t68, t67, t65, t64, t59, t58, t57, t42, t20, t1;

    IT a1, b1, w1, xc1, yc1;
    IT a2, b2, w2, xc2, yc2;

    a1  = ET::to_interval(e1.major_axis()); 
    b1  = ET::to_interval(e1.minor_axis()); 
    w1  = ET::to_interval(e1.rotation()); 
    xc1 = ET::to_interval(e1.x_center()); 
    yc1 = ET::to_interval(e1.y_center());
    a2  = ET::to_interval(e2.major_axis()); 
    b2  = ET::to_interval(e2.minor_axis()); 
    w2  = ET::to_interval(e2.rotation()); 
    xc2 = ET::to_interval(e2.x_center()); 
    yc2 = ET::to_interval(e2.y_center());

    t15 = t*t;
    t11 = t15*t15;
    t94 = -t11+IT(1);
    t118 = b1*t94;
    t10 = t15*t;
    t122 = (b1*b1-a1*a1)*(-t10+t);
    t12 = r*r;
    t8 = t12*t12;
    t97 = t8-IT(1);
    t101 = b2*t97;
    t66 = -yc2+yc1;
    t9 = r*t12;
    t77 = t9+r;
    t110 = t77*w2;
    t128 = t66*t110;
    t90 = xc2-xc1;
    t6 = w2*w2;
    t72 = -t6+IT(1);
    t47 = (t8-t8*t6-t72)*yc2;
    t123 = (-t9+r)*(-b2*b2+a2*a2);
    t85 = IT(-2)*t6;
    t23 = (IT(-2)+t85)*t123;
    t126 = -t23+(t47-t90*w2*(IT(2)*t8-IT(2)))*b2;
    t89 = t77*a2;
    t124 = w2*t101;
    t96 = -t10-t;
    t109 = t96*a1;
    t75 = IT(2)+t85;
    t105 = t90*t75;
    t88 = xc2*t77;
    t104 = t75*t88;
    t7 = w1*w1;
    t92 = t75*t7;
    t51 = t92+t75;
    t86 = t96*t77;
    t84 = IT(2)*xc1;
    t79 = IT(-4)+IT(4)*t6;
    t71 = a2*w2;
    t69 = w2*xc2;
    t68 = yc2*w2;
    t67 = w1*yc1;
    t65 = IT(-4)*t68;
    t64 = w1*t84;
    t59 = IT(-2)*xc2+IT(4)*t67+t84;
    t58 = t64-t66;
    t57 = IT(-2)*t69-t66;
    t42 = IT(-4)*t67+IT(2)*t90;
    t20 = (IT(4)*t68+t105)*t89+t126;
    t1 = ((w2-t7*w2+t6*w1-w1)*t101*t118+(-t89*t118+t101*t109)*
          (IT(4)*w2*w1+t7*t6-t7+t72)+
          (w1*t79-w2*(IT(-4)+IT(4)*t7))*a2*a1*t86);

    if (CGAL::zero_in(t1)) {
        // std::cerr << "t1 = " << t1 << std::endl;
        return -2;
    }

    IT tx1 = (IT(-1)*(IT(8)*t7-IT(8))*t71*t86-t96*(t92-t75)*t101)*xc1-
            t96*(((IT(8)-IT(8)*t6)*t88-IT(16)*t128)*a2+
            ((IT(8)-IT(8)*t8)*t69+t97*t66*t79)*b2+
            (IT(8)*t6+IT(8))*t123)*w1;

    IT tx2 = -(-t23+(t104-IT(4)*t128)*a2-(-t66*t6-t57)*t101)*t7-
              t23+(-t58*t6+t64+t57)*t101+(t104+IT(4)*t58*t110)*a2;

    x = (tx1*a1-tx2*t118+(-t77*t71*(-IT(8)-IT(8)*t7)+t51*t101)*t122)/
        (IT(2)*t1);

    y = ((((-t7+t7*t11+t94)*t124+(t7-IT(1))*t72*t94*t89)*yc1+
            ((t65-t105)*t89+t20*t11-t126)*w1)*b1-
            ((-t47+(t8*t42+t59)*w2)*b2+
            (t65+t6*t42+t59)*t89+t20*t7+t23)*t109-
            ((-IT(2)-IT(2)*t7)*t124+t89*t51)*t122)/t1;

//    if (CGAL::sign(e1.evaluate_equation(ET::to_rational(CGAL::median(x)),
//                                        ET::to_rational(CGAL::median(y)))) !=
//        CGAL::sign(e2.evaluate_equation(ET::to_rational(CGAL::median(x)),
//                                        ET::to_rational(CGAL::median(y))))) return -3;
    return 0;
}

} // VORELL

} //namespace CGAL
#endif
