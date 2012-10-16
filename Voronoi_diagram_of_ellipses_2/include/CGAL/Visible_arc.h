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

#ifndef CGAL_VISIBLE_ARC_H
#define CGAL_VISIBLE_ARC_H

#include<iostream>
#include<CGAL/basic.h>
#include<CGAL/vorell.h>

#include<CGAL/Range.h>
#include<CGAL/Ellipse_2.h>
#include<CGAL/Bitangent.h>

namespace CGAL {
namespace VORELL {

template<class ET>
class Visible_arc: public Range<typename ET::AK::Algebraic_real_1> {
    typedef typename ET::QQ QQ;
    typedef typename ET::AK AK;
    typedef typename ET::upolz_t upolz_t;
    typedef typename ET::upolq_t upolq_t;
    typedef typename ET::Root Root;
    typedef typename CGAL::VORELL::Generated_code<ET> gencode;

    // ccw with a point P(x,y), 
    // and the endpoints of the polar arc from P to the ellipse
    QQ ccw_point_polar_arc(const QQ& x, const QQ& y, const QQ& a, const QQ& b, 
                            const QQ& w, const QQ& xc, const QQ& yc, 
                            const upolz_t &tanpoly) {
	typename ET::PTZ::Get_coefficient coeff;
        QQ t1,t2,t3,t4,t5,t6,t7, res;
        t6 = -1*yc+y;
        t2 = 2*w;
        t7 = t6*t2;
        t5 = x-xc;
        t4 = -1*a-t5;
        t3 = -1*a+t5;
        t1 = w*w;
        res = (t5*t2+t6*t1-t6)*QQ(coeff(tanpoly,1))*a+
                ((t7+t4*t1+t3)*QQ(coeff(tanpoly,2))+
                (-1*t7+t3*t1+t4)*QQ(coeff(tanpoly,0)))*b;
        return res;
    }

public:
    Visible_arc(Bitangent<ET> &bt12, Bitangent<ET> &bt13) {
        Range<Root> t12 = bt12.CH_range();
        //std::cerr << "range12 = " << t12 << std::endl;

        Range<Root> t13 = bt13.CH_range();
        //std::cerr << "range13 = " << t13 << std::endl;

        *static_cast<Range<Root> *>(this) = t12.intersection(t13);
        //std::cerr << "range = " << *static_cast<Range<Root> *>(this) << std::endl;
    }
    
    // this is slower, TODO: remove?
    Visible_arc(const Ellipse_2<ET>& e, const Ellipse_2<ET>& e0, const QQ& t) {
        std::vector<Root> rsol;
        rsol.reserve(4);
        typename ET::PTQ::Evaluate evaluate1;
        typename ET::BPTQ::Evaluate evaluate2;
        typename ET::BPTQ::Swap swap;
        
//        upolz_t tanp = gencode().tan_poly_xy(ELL_PARAM_COEFFS(e), x, y);
//        typename ET::bpolq_t polar_tr = typename CGAL::Coercion_traits<typename ET::bpolq_t, typename ET::bpolz_t>::Cast()
//                                        (gencode().polar(ELL_PARAM_COEFFS(e0), ELL_PARAM_COEFFS(e)));
        typename ET::bpolq_t polar_tr =
                gencode().polar_q(ELL_PARAM_COEFFS(e0), ELL_PARAM_COEFFS(e));
//        std::cerr << polar_tr << std::endl;
        upolq_t tanp = evaluate2(swap(polar_tr, 0, 1), t);
//        std::cerr << tanp << std::endl;
        
        typename AK::Solve_1()(CGAL::VORELL::primpart(tanp), 
                                std::back_inserter(rsol));
        if (rsol.size() != 2) return;

        QQ midr = -tanp[1] / (2*tanp[2]);
        
        // ***** sign of tangent @ midr ****
        upolq_t tan_center = gencode().tan_poly_xy_q(
                                    ELL_PARAM_COEFFS(e), e.xc(), e.yc());
        Sign sm = CGAL::sign(evaluate1(tanp, midr)) *
                  CGAL::sign(evaluate1(tan_center, midr));

        if (sm == CGAL::NEGATIVE) {
            *static_cast<Range<Root> *>(this) = Range<Root>(rsol[0], rsol[1]);
        } else {
    //        from x,y visible arc of e contains i-point
            *static_cast<Range<Root> *>(this) = 
                    Range<Root>(rsol[1], rsol[0]);
        }
    }

    Visible_arc(const Ellipse_2<ET>& e, const QQ& x, const QQ& y) {
        std::vector<Root> rsol;
        rsol.reserve(4);

        upolz_t tanp = gencode().tan_poly_xy(ELL_PARAM_COEFFS(e), x, y);
        
        typename AK::Solve_1()(tanp, true, std::back_inserter(rsol));
        if (rsol.size() != 2) return;

        QQ midr = -tanp[1] / (2*tanp[2]);

        // ***** line r1 r2 (Determinant sign) (CCWs) ******
        Sign Lrr = CGAL::sign(
                    ccw_point_polar_arc(x, y, ELL_PARAM_COEFFS(e), tanp)) *
                   CGAL::sign(
                    ccw_point_polar_arc(e.boundary_x(midr), e.boundary_y(midr),
                                            ELL_PARAM_COEFFS(e), tanp));

        if (Lrr == CGAL::POSITIVE) {
            *static_cast<Range<Root> *>(this) = Range<Root>(rsol[0], rsol[1]);
        } else {
    //        from x,y visible arc of e contains i-point
            *static_cast<Range<Root> *>(this) = 
                    Range<Root>(rsol[1], rsol[0]);
        }
    }

};


} // namespace
} //namespace CGAL

#endif
