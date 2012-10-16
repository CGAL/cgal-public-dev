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

#ifndef CGAL_MEDIAL_AXIS_LOCATION_H
#define CGAL_MEDIAL_AXIS_LOCATION_H

#include<iostream>
#include<CGAL/basic.h>
#include<CGAL/vorell.h>

#include<CGAL/Range.h>
#include<CGAL/Ellipse_traits.h>

// #include<CGAL/Timer.h>

namespace CGAL {

namespace VORELL {

template<class ET>
class Medial_axis_location {
    typedef typename ET::Root Root;
    typedef typename VORELL::Range<Root> Range;
    typedef typename ET::IT IT;

    Ellipse_2<ET> e1, e2;
    Range trange;
    std::vector<typename ET::Root> sol;
    typename ET::upolz_t res; // resultant deg.16
    Range ma_range;

public:
    Medial_axis_location() { }
    Medial_axis_location(const Ellipse_2<ET> &e1_, const Ellipse_2<ET> &e2_, const Range &trange_):
                         e1(e1_), e2(e2_), trange(trange_) {

        // TODO: optimize by taking into account symmetry in trange;
        typename ET::QPTQ::Move move;
        typename ET::BPTZ::Move move2;
        typename ET::QPTQ::Shift shift;
        typename ET::TPTQ::Shift shift3;

        typename ET::upolq_t dt1;
        typename ET::qpolq_t xt, yt, xr, yr, x, y, dt, dr;
        typename ET::qpolz_t Mtr, Nr;

        typename CGAL::Coercion_traits<typename ET::qpolq_t,
                                        typename ET::upolq_t>::Cast QP;
        typename CGAL::Coercion_traits<typename ET::tpolq_t,
                                        typename ET::upolq_t>::Cast TP;
        typename CGAL::Coercion_traits<typename ET::bpolz_t,
                                        typename ET::upolz_t>::Cast BPZ;
        typename CGAL::VORELL::Generated_code<ET> GC;

        if (trange.is_finite() &&
            sign(trange.left()) == sign(trange.right())) return;

        xt = QP( GC.numerxt(ELL_PARAM_COEFFS(e1)) );
        yt = QP( GC.numeryt(ELL_PARAM_COEFFS(e1)) );
        dt1 = GC.denomt(ELL_PARAM_COEFFS(e1));
        dt = QP( dt1 );

        xr = move(QP( GC.numerxt(ELL_PARAM_COEFFS(e2)) ), 0, 3);
        yr = move(QP( GC.numeryt(ELL_PARAM_COEFFS(e2)) ), 0, 3);
        dr = move(QP( GC.denomt(ELL_PARAM_COEFFS(e2)) ), 0, 3);

        x = shift(typename ET::qpolq_t(1),1,1);
        y = shift(typename ET::qpolq_t(1),1,2);

        typename ET::qpolq_t t3 = -2*dt;
        Mtr = CGAL::VORELL::primpart( (-xr*xr-yr*yr+2*(x*xr+yr*y)*dr)*dt*dt+
                                      ((y*t3+yt)*yt+(x*t3+xt)*xt)*dr*dr );

        typename ET::upolq_t nxr, nyr, ncr;
        GC.ell_normal(ELL_PARAM_COEFFS(e2), nxr, nyr, ncr);
        Nr = CGAL::VORELL::primpart ( move(QP( nxr ), 0, 3)*x +
                                      move(QP( nyr ), 0, 3)*y +
                                      move(QP( ncr ), 0, 3) );

        typename ET::tpolz_t res1 = typename ET::QPTZ::Resultant()(Mtr, Nr);

        GC.ell_normal(ELL_PARAM_COEFFS(e1), nxr, nyr, ncr);
        typename ET::upolz_t nyt = CGAL::VORELL::primpart(nyr);
        typename ET::tpolz_t Nt =
                CGAL::VORELL::primpart (
                    TP( nxr ) * shift3(typename ET::tpolq_t(1),1,1) +
                    TP( nyr ) * shift3(typename ET::tpolq_t(1),1,2) +
                    TP( ncr ) );
        typename ET::bpolz_t res1y = typename ET::TPTZ::Resultant()(res1, Nt);

        typename ET::upolz_t dt2, dt4, dt6;
        dt2 = CGAL::VORELL::primpart(dt1);
        dt2 = dt2*dt2; dt4 = dt2 * dt2; dt6 = dt4 * dt2;
        res1y = typename ET::BPTZ::Swap()(res1y, 0, 1);
        typename ET::bpolz_t res1yq =
                typename ET::BPTZ::Pseudo_division_quotient()(
                    res1y, move2(BPZ(dt6),0,1));

        res1yq = typename ET::BPTZ::Swap()(res1yq, 0, 1);

        //std::cerr << res1yq << std::endl;

        typename ET::bpolz_t med_x = GC.medial_x(ELL_PARAM_COEFFS(e1));
        //std::cerr << med_x << std::endl;
        res = typename ET::BPTZ::Resultant()(res1yq, med_x);

        typename ET::upolz_t nyt2, nyt6;
        nyt2 = nyt*nyt; // 2
        nyt6 = nyt2*nyt2; // 4
        nyt6 = nyt6*nyt2; // 6

        res = typename ET::PTZ::Pseudo_division_quotient()(res, nyt6);
        res = typename ET::PTZ::Canonicalize()(res);
        //std::cerr << res << std::endl;

        //std::cerr << typename ET::PTZ::Degree()(res) << std::endl;

        std::vector<Root> sol2 = ET::Real_roots(res);
        typename std::vector<Root>::iterator it;

        for (it = sol2.begin(); it != sol2.end(); it++) {
            if (trange.contains(*it)) sol.push_back(*it);
        }

        if (sol.size() < 2) return;

        int l = 0;
        int r = sol.size() - 1;

        // TODO: fix this with exact
        if (trange.is_finite()) {
            while (l < r) {
                IT lint = - ET::to_interval(sol[l]);
                IT rint = ET::to_interval(sol[r]);
                if (CGAL::overlap(lint, rint)) {
                    ma_range = Range(sol[l], sol[r]);
                    break;
                } else if (lint > rint) l++;
                else r--;
            }
        } else {
            IT lint = - ET::to_interval(sol[l]);
            IT rint = ET::to_interval(sol[r]);
            if (CGAL::overlap(lint, rint)) {
                ma_range = Range(sol[r], sol[l]);
             } else l = r; // sentinel
        }
        CGAL_assertion( l < r );

#if VERBOSE > 1
        std::cerr << "MA = [" << to_double(ma_range.left()) << ',' << to_double(ma_range.right()) << ']' << std::endl;

        std::cerr << "t-sols: ";
        for (it = sol.begin(); it != sol.end(); it++) {
            std::cerr << CGAL::to_double(*it);
            std::cerr << ' ';
        }
        std::cerr << std::endl;
#endif

    }

    Range operator()() const { return ma_range; }

    const typename ET::upolz_t& resultant() const {
        return res;
    }

};


} // VORELL

} //namespace CGAL
#endif
