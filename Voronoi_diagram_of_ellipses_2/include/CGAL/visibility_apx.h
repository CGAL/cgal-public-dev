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

#ifndef CGAL_VISIBILITY_APX_H
#define CGAL_VISIBILITY_APX_H

#include<iostream>
#include<CGAL/basic.h>
#include<CGAL/vorell.h>

#include<CGAL/Range.h>
#include<CGAL/Ellipse_2.h>

#include<CGAL/Interval_trinomial_solver.h>

#include<CGAL/Voronoi_diagram_of_ellipses_2/Polynomial_factory.h>

namespace CGAL {
namespace VORELL {

template<class ET>
class Visible_arc_apx: public Range<typename ET::BT> {
    typedef typename ET::QQ QQ;
    typedef typename ET::BT BT;
    typedef typename ET::IT IT;
    typedef typename ET::upoli_t upoli_t;
    typedef typename ET::upolz_t upolz_t;
    typedef typename ET::bpoli_t bpoli_t;
    typedef typename ET::bpolz_t bpolz_t;

    typedef typename VORELL::Generated_code<ET> gencode;
    typedef typename VORELL::Ellipse_univariate_polynomial_factory<ET> polfact1;
    typedef typename VORELL::Ellipse_bivariate_polynomial_factory<ET> polfact2;
    typedef typename Coercion_traits<bpoli_t, bpolz_t>::Cast BPI;
    typedef typename Coercion_traits<upoli_t, upolz_t>::Cast UPI;

public:

    Visible_arc_apx(const Ellipse_2<ET>& e, const Ellipse_2<ET>& e0, const IT& t) {
        std::pair<IT,IT> rsol;
        typename ET::PTI::Evaluate evaluate1;
        typename ET::BPTI::Evaluate evaluate2;
        typename ET::BPTI::Swap swap;

        bpoli_t polar_tr = BPI()(gencode().polar(ELL_PARAM_COEFFS(e0),
                                                    ELL_PARAM_COEFFS(e)));
//        bpoli_t polar_tr = BPI()(polfact2(e0,e).polar());
        upoli_t tanp = evaluate2(swap(polar_tr, 0, 1), t);

        typename CGAL::Interval_trinomial_solver<upoli_t, IT> solver(tanp);

        rsol = solver();

        if (solver.get_status() != 2) return;

        IT midr = -tanp[1] / (IT(2)*tanp[2]);

        // ***** sign of tangent @ midr ****
        upoli_t tan_center = UPI()(gencode().tan_poly_xy(ELL_PARAM_COEFFS(e),
                                                 e.x_center(), e.y_center()));
//        upoli_t tan_center = UPI()(polfact1(e).tan_poly_xy());
        Sign sm = sign(evaluate1(tanp, midr)) *
                  sign(evaluate1(tan_center, midr));

        // outwards rounding
        if (sm == NEGATIVE) {
            *static_cast<Range<BT> *>(this) =
                    Range<BT>(lower(rsol.first),
                                upper(rsol.second));
        } else {
        //  from x,y visible arc of e contains i-point
            *static_cast<Range<BT> *>(this) =
                    Range<BT>(lower(rsol.second),
                                upper(rsol.first));
        }
    }

};

template<class ET>
class Apollonius_arc_apx: public Range<typename ET::BT> {
    typedef typename ET::QQ QQ;
    typedef typename ET::BT BT;
    typedef typename ET::IT IT;
    typedef typename ET::upoli_t upoli_t;
    typedef typename ET::upolz_t upolz_t;
    typedef typename ET::bpoli_t bpoli_t;
    typedef typename ET::bpolz_t bpolz_t;
    typedef typename VORELL::Ellipse_univariate_polynomial_factory<ET> polfact1;
    typedef typename VORELL::Ellipse_bivariate_polynomial_factory<ET> polfact2;

public:
    
    Apollonius_arc_apx(const Ellipse_2<ET>& e, const Ellipse_2<ET>& e0, 
                       const IT& t) {
        std::pair<IT,IT> rsol;
        typename ET::PTI::Evaluate evaluate1;
        typename ET::BPTI::Evaluate evaluate2;
        typename ET::BPTI::Swap swap;
        
        bpoli_t cut_tr =
                typename Coercion_traits<bpoli_t, bpolz_t>::Cast()(
                    Generated_code<ET>().tan_poly_cut(
                        ELL_PARAM_COEFFS(e0), ELL_PARAM_COEFFS(e)));
//        bpoli_t cut_tr =
//                typename Coercion_traits<bpoli_t, bpolz_t>::Cast()(
//                    polfact2(e0,e).tan_poly_cut());

        upoli_t tcut = evaluate2(swap(cut_tr, 0, 1), t);

        typename CGAL::Interval_trinomial_solver<upoli_t, IT> solver(tcut);
        
        rsol = solver();
        
        Visible_arc_apx<ET> vis = Visible_arc_apx<ET>(e, e0, t);
        
        if (solver.get_status() == 0) {
            *(static_cast<Range<BT> *>(this)) = vis;
            return;
        }
        if (solver.get_status() != 2) return;

        IT Ept, Qpt;
        // TODO: Qpt is an interval, check degenerate cases

        // TODO: rsol is an interval... is it contained in whole?
        CGAL_assertion( vis.contains(lower(rsol.first)) && 
                        vis.contains(upper(rsol.first)) ||
                        vis.contains(lower(rsol.second)) && 
                        vis.contains(upper(rsol.second)) );
        if (vis.contains(lower(rsol.first))) Qpt = rsol.first; 
        else Qpt = rsol.second;
    //    std::cerr << "Qpt = " << Qpt << std::endl;

        // s1 is always positive (sign of center(e1) at tangent e1(t))
        Sign s2 = sign(evaluate1(tcut, IT(vis.left())));
    //    std::cerr << "s1 s2 = " << s1 << ' ' << s2 << std::endl;

        if (s2 == POSITIVE) Ept = vis.right(); else Ept = vis.left();
        // endpoint chosen

        // outwards rounding
        // TODO: can be optimized ?
        if (vis.is_finite()) {          // no ipt
            if (Qpt > Ept) *static_cast<Range<BT> *>(this) = 
                Range<BT>(lower(Ept), upper(Qpt));
            else *static_cast<Range<BT> *>(this) = 
                Range<BT>(lower(Qpt), upper(Ept));
        } else {
            if (Qpt < vis.right()) {
                if (Ept == vis.left()) *static_cast<Range<BT> *>(this) = 
                    Range<BT>(lower(Ept), upper(Qpt));
                else *static_cast<Range<BT> *>(this) = 
                    Range<BT>(lower(Qpt), upper(Ept));
            } else {
                if (Ept == vis.right()) *static_cast<Range<BT> *>(this) = 
                    Range<BT>(lower(Qpt), upper(Ept));
                else *static_cast<Range<BT> *>(this) = 
                    Range<BT>(lower(Ept), upper(Qpt));
            }
        }
    }
};


} // namespace
} //namespace CGAL

#endif
