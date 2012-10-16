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

#ifndef CGAL_VORONOI_CIRCLE_EXACT_H
#define CGAL_VORONOI_CIRCLE_EXACT_H

#include<iostream>
#include<CGAL/basic.h>
#include<CGAL/vorell.h>

#include<CGAL/Range.h>
#include<CGAL/Ellipse_traits.h>
#include<CGAL/Ellipse_triplet.h>
#include<CGAL/Voronoi_circle_apx.h>

// #include<CGAL/Timer.h>

namespace CGAL {

namespace VORELL {

template<class ET>
class Voronoi_circle_exact {
    
    Ellipse_triplet<ET> triplet;
    std::vector<typename ET::Root> sol;
    typename ET::upolz_t res; // resultant deg.184
    
    // bad idea with iterator: somehow got invalidated!
    //typename std::vector<typename ET::Root>::iterator t_sol;
    
    typename ET::Root t_sol;
    
public:    
    Voronoi_circle_exact() { }
    Voronoi_circle_exact(const Ellipse_triplet<ET>& triplet_): 
                         triplet(triplet_) { 
//        CGAL::Timer timer;
//        timer.start();

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

        xt = QP( GC.numerxt(ELL_PARAM_COEFFS(triplet.get_e1())) );
        yt = QP( GC.numeryt(ELL_PARAM_COEFFS(triplet.get_e1())) );
        dt1 = GC.denomt(ELL_PARAM_COEFFS(triplet.get_e1()));
        dt = QP( dt1 );

        xr = move(QP( GC.numerxt(ELL_PARAM_COEFFS(triplet.get_e2())) ), 0, 3);
        yr = move(QP( GC.numeryt(ELL_PARAM_COEFFS(triplet.get_e2())) ), 0, 3);
        dr = move(QP( GC.denomt(ELL_PARAM_COEFFS(triplet.get_e2())) ), 0, 3);

        x = shift(typename ET::qpolq_t(1),1,1);
        y = shift(typename ET::qpolq_t(1),1,2);

        typename ET::qpolq_t t3 = -2*dt;
        Mtr = CGAL::VORELL::primpart( (-xr*xr-yr*yr+2*(x*xr+yr*y)*dr)*dt*dt+
                                      ((y*t3+yt)*yt+(x*t3+xt)*xt)*dr*dr );

        typename ET::upolq_t nxr, nyr, ncr;
        GC.ell_normal(ELL_PARAM_COEFFS(triplet.get_e2()), nxr, nyr, ncr);
        Nr = CGAL::VORELL::primpart ( move(QP( nxr ), 0, 3)*x + 
                                      move(QP( nyr ), 0, 3)*y + 
                                      move(QP( ncr ), 0, 3) );

    //    cout << Mtr << endl;
    //    cout << Nr << endl;
        typename ET::tpolz_t res1 = typename ET::QPTZ::Resultant()(Mtr, Nr);

        GC.ell_normal(ELL_PARAM_COEFFS(triplet.get_e1()), nxr, nyr, ncr);
        typename ET::upolz_t nyt = CGAL::VORELL::primpart(nyr);
        typename ET::tpolz_t Nt = 
                CGAL::VORELL::primpart ( 
                    TP( nxr ) * shift3(typename ET::tpolq_t(1),1,1) +
                    TP( nyr ) * shift3(typename ET::tpolq_t(1),1,2) +
                    TP( ncr ) );
        typename ET::bpolz_t res1y = typename ET::TPTZ::Resultant()(res1, Nt);
    //    cout << res1 << endl;
    //    cout << res1y << endl;
    //    cout << ET::TPTZ::Degree()(res1y, 0) << endl;
    //    cout << ET::TPTZ::Degree()(res1y, 1) << endl;

        typename ET::upolz_t dt2, dt4, dt6;
        dt2 = CGAL::VORELL::primpart(dt1);
        dt2 = dt2*dt2; dt4 = dt2 * dt2; dt6 = dt4 * dt2;
    //    cout << CGAL::VORELL::primpart(dt6) << endl;
        res1y = typename ET::BPTZ::Swap()(res1y, 0, 1);
        typename ET::bpolz_t res1yq = 
                typename ET::BPTZ::Pseudo_division_quotient()(
                    res1y, move2(BPZ(dt6),0,1));
    //    cout << res1yq << endl;
    //    cout << ET::TPTZ::Degree()(res1yq, 0) << endl;
    //    cout << ET::TPTZ::Degree()(res1yq, 1) << endl;

        // Mts
        xr = move(QP( GC.numerxt(ELL_PARAM_COEFFS(triplet.get_e3())) ), 0, 3);
        yr = move(QP( GC.numeryt(ELL_PARAM_COEFFS(triplet.get_e3())) ), 0, 3);
        dr = move(QP( GC.denomt(ELL_PARAM_COEFFS(triplet.get_e3())) ), 0, 3);
        Mtr = CGAL::VORELL::primpart( (-xr*xr-yr*yr+2*(x*xr+yr*y)*dr)*dt*dt+
                                      ((y*t3+yt)*yt+(x*t3+xt)*xt)*dr*dr );
        GC.ell_normal(ELL_PARAM_COEFFS(triplet.get_e3()), nxr, nyr, ncr);
        Nr = CGAL::VORELL::primpart ( move(QP( nxr ), 0, 3)*x + 
                                      move(QP( nyr ), 0, 3)*y + 
                                      move(QP( ncr ), 0, 3) );
        res1 = typename ET::QPTZ::Resultant()(Mtr, Nr);
        res1y = typename ET::TPTZ::Resultant()(res1, Nt);
        res1y = typename ET::BPTZ::Swap()(res1y, 0, 1);
        typename ET::bpolz_t res2yq = 
                typename ET::BPTZ::Pseudo_division_quotient()(
                    res1y, move2(BPZ(dt6),0,1));
    //    cout << ET::TPTZ::Degree()(res2yq, 0) << endl;
    //    cout << ET::TPTZ::Degree()(res2yq, 1) << endl;

        res1yq = typename ET::BPTZ::Swap()(res1yq, 0, 1);
        res2yq = typename ET::BPTZ::Swap()(res2yq, 0, 1);
        
        res = typename ET::BPTZ::Resultant()(res1yq, res2yq);
//        timer.stop();
        // std::cerr << "resultant computation @ " << timer.time() << std::endl;
//        timer.start();

        typename ET::upolz_t nyt4, nyt32;
        nyt4 = typename ET::PTZ::Pseudo_division_quotient()(
                nyt, CGAL::VORELL::primpart(dt1));
        nyt4 = nyt4*nyt4; // ^2
        nyt4 = nyt4*nyt4; nyt32 = nyt4*nyt4;
        nyt32 = nyt32*nyt32; nyt32 = nyt32*nyt32;
        nyt32 = nyt32*nyt4; // ^36

        typename ET::upolz_t dt40 = dt6*dt6; // ^12
        dt40 = dt40*dt40*dt40*dt4;
        res = typename ET::PTZ::Pseudo_division_quotient()(res, dt40);
        res = typename ET::PTZ::Pseudo_division_quotient()(res, nyt32);
        res = typename ET::PTZ::Canonicalize()(res);
//        timer.stop();
        // std::cerr << "factors & primpart @ " << timer.time() << std::endl;
        // std::cout << res << std::endl;

        //std::cerr << ET::PTZ::Degree()(res) << std::endl;

//        timer.start();
        sol = ET::Real_roots(res);
//        timer.stop();
        // std::cerr << "solve @ " << timer.time() << std::endl;
//        std::cerr << "t-range = " << triplet.t_range() << std::endl;


#if VERBOSE > 1
        typename std::vector<typename ET::Root>::iterator it;
        
        std::cerr << "t-sols: ";
        for (it = sol.begin(); it != sol.end(); it++) {
                std::cerr << CGAL::to_double(*it);
            if (triplet.t_range().contains(*it)) std::cerr << "* ";
            else std::cerr << ' ';
        }
        std::cerr << std::endl;
#endif
        
//        timer.start();
        Voronoi_circle_apx<ET> vc(triplet);
        t_sol = vc.isolate_external_circle(sol.begin(), sol.end());
//        timer.stop();

#if VERBOSE > 1
        std::cerr << "external = " << CGAL::to_double(t_sol) << std::endl;
        //std::cerr << "isolate @ " << timer.time() << std::endl;
#endif
    }
    
    const typename ET::upolz_t& resultant() const {
        return res;
    }
    
    const typename ET::Root& operator()() const {
        return t_sol;
    }
    
//    friend std::ostream& operator<<(std::ostream& o, const VoronoiCircleApx<ET>& v) {
//        return (o << "VC(" << v.get_t() << ',' << v.get_r() << ',' << v.get_s() << ')');
//    }

};

} // VORELL

} //namespace CGAL
#endif
