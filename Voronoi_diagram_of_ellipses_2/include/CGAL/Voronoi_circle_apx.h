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

#ifndef CGAL_VORONOI_CIRCLE_APX_H
#define CGAL_VORONOI_CIRCLE_APX_H

#include<iostream>
#include<algorithm>

#include<CGAL/basic.h>
#include<CGAL/vorell.h>

#include<CGAL/Range.h>
#include<CGAL/Ellipse_2.h>
#include<CGAL/Ellipse_triplet.h>
#include<CGAL/visibility_apx.h>

#include <CGAL/Interval_newton_solver.h>
#include <CGAL/Bisector_curve_apx.h>

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#ifdef VORELL_VCAPX_PROFILE
#include <CGAL/Timer.h>
#endif

#include <CGAL/Object_cache.h>

namespace CGAL {

namespace VORELL {
    
template<class ET>
class Voronoi_circle_apx {
    typedef VORELL::Range<typename ET::Root> Range;
    typedef typename ET::ARX ARX;
    typedef typename ET::IT IT;
    typedef typename ET::BT BT;
    typedef typename VORELL::Bisector_curve_apx<ET> Bisector_curve_apx;
    typedef typename ET::bpolz_t bpolz_t;
    typedef typename ET::bpoli_t bpoli_t;
    typedef typename ET::upoli_t upoli_t;
    typedef typename std::vector<typename ET::Root>::iterator rit;
    typename CGAL::Interval_traits<IT>::Construct make_interval;

    Ellipse_triplet<ET> triplet;
    
    int iters;
    Bisector_curve_apx bx12, bx23, bx31;
    ARX tt, rr, ss;
    ARX ch12, ch21, ch32, ch13;
    bool normalized;

#ifndef VORELL_VCAPX_NOCACHE
    typedef typename Ellipse_triplet<ET>::Key_type CKey;
    typedef ::boost::tuple<ARX, ARX, ARX> CValue;
    typedef Value_cache<CKey, CValue> Cache;
    Cache cache;
    bool cached;
#endif

    struct Ordering_interval_root {
        bool operator()(const typename ET::Root& a, const IT &x) {
            return CGAL::upper(ET::to_interval(a)) < CGAL::lower(x);
        }
        bool operator()(const IT& x, const typename ET::Root& a) {
            return CGAL::upper(x) < CGAL::lower(ET::to_interval(a));
        }
    };
    
#if VERBOSE > 0
    void info_all() {
        std::cerr << "RANGES:\n";
        std::cerr << "trange = [" << tt.left() << ',' << tt.right() << ']';
        if (tt.is_finite()) std::cerr << tt.width();
        std::cerr << std::endl;
        std::cerr << "rrange = [" << rr.left() << ',' << rr.right() << ']';
        if (rr.is_finite()) std::cerr << rr.width();
        std::cerr << std::endl;
        std::cerr << "srange = [" << ss.left() << ',' << ss.right() << ']';
        if (ss.is_finite()) std::cerr << ss.width();
        std::cerr << std::endl;
        std::cerr << "iters = " << iters << std::endl;
//        if (trange.is_finite() && rrange.is_finite())
//            std::cerr << "btr = " << bx12.polyx(trange.interval())(rrange.interval()) << std::endl;
//        if (rrange.is_finite() && srange.is_finite())
//            std::cerr << "brs = " << bx23.polyx(rrange.interval())(srange.interval()) << std::endl;
//        if (srange.is_finite() && trange.is_finite())
//            std::cerr << "bst = " << bx31.polyx(srange.interval())(trange.interval()) << std::endl;
    }
#endif

    inline int handle_pair_opposite(const IT &x, IT &y, Bisector_curve_apx &bx,
            const ARX &ch, const ARX &yy, const BT t1) {
#if VERBOSE > 1
        std::cerr << "hpo: x = " << x << std::endl;
#endif
        y = bx(x, normalized, yy);
        if (bx.get_status() < 0) {
#if VERBOSE > 1
            std::cerr << "hpo: BAD STATUS = " << bx.get_status() << std::endl;
#endif
            return -1;
        }
#if VERBOSE > 1
        std::cerr << "hpo: y = " << y << std::endl;
#endif
        if (ch.contains(yy.left()) && 
                ch.order(upper(y), yy.left()) < 0) { // y too left
            tt = ARX(tt.left(), t1); 
            // move t1 to the left (moves opposite)
            return 0;
        } else if (ch.contains(yy.right()) && 
                ch.order(lower(y), yy.right()) > 0) { // y too right
            tt = ARX(t1, tt.right()); // move t1 to the right
            return 0;
        }
        return 1;
    }

    inline int handle_pair_same(const IT &x, IT &y, Bisector_curve_apx &bx, 
            const ARX &ch, const ARX &yy, const BT t1) {
#if VERBOSE > 1
        std::cerr << "hpo: x = " << x << std::endl;
#endif
        if (normalized) y = bx(x, yy.is_finite(), yy);
        else y = bx(x);
        if (bx.get_status() < 0) {
#if VERBOSE > 1
            std::cerr << "hpo: BAD STATUS = " << bx.get_status() << std::endl;
#endif
            return -1;
        }
#if VERBOSE > 1
        std::cerr << "hpo: y = " << y << std::endl;
#endif
        if (ch.contains(yy.left()) && 
                ch.order(upper(y), yy.left()) < 0) { // y too left
            tt = ARX(t1, tt.right()); 
            // move t1 to the right (moves the same)
            return 0;
        } else if (ch.contains(yy.right()) && 
            ch.order(lower(y), yy.right()) > 0) { // y too right
            tt = ARX(tt.left(), t1); // move t1 to the left
            return 0;
        }
        return 1;
    }
    
    inline int update_interval(BT x1l, BT x1u, BT x2l, BT x2u, ARX& xx) {
        if (xx.is_finite()) {
            BT w = xx.width();
            if (x2l > x1u) xx = ARX(x1l, x2u);
            else if (x2u < x1l) xx = ARX(x2l, x1u);
            else xx = ARX(x2l, x2u).hull(ARX(x1l, x1u));
            // TODO: consider intersection with old?
            if (xx.width() >= w) return 0;
        } else {
            if (x2l > x1u || x2u < x1l) {
                if (xx.order(x1u, x2l) == -1) xx = ARX(x1l, x2u);
                else xx = ARX(x2l, x1u);
            } else xx = ARX(x2l, x2u).hull(ARX(x1l, x2u));
            // TODO: consider intersection with old?
        }
        return 1;
    }
    
    // returns 1 on sucessful refinement, 0 if precision limit reached
    // TODO: have to make it faster by not computing apoll.arc every time
    // we can use refined r2r1  instead!!

//    enum Refinement_state { STARTING = 0, REFINING = 1, REFINED = 2 };
//    Refinement_state refinement_state;
    int prec_achieved;

    int refine_t() {
        BT t1;
        IT r2, s1, t2;
        IT r1, s2;
        int v;

#ifndef VORELL_VCAPX_NOCACHE
        if (cached) {
#if VERBOSE > 0
            std::cerr << "vc cache hit!" << std::endl;
#endif
            return 0;
        }
#endif
        bool wanna_stop = false;
        while (true) {
            // TODO: fix this binary search (with optimal normalization (spm))

            t1 = ET::midpoint(tt);
#if VERBOSE > 1
            std::cerr << "midpoint t1 = " << t1 << std::endl;
#endif
            wanna_stop = t1 == tt.left();
            if (wanna_stop) return 0;
//            if (wanna_stop && refinement_state == STARTING) return 0;

            iters++;

            v = handle_pair_opposite(t1, r2, bx12, ch21, rr, t1);
            if (v == 0) continue; else if (v < 0) return 0;

            v = handle_pair_same    (r2, s1, bx23, ch32, ss, t1);
            if (v == 0) continue; else if (v < 0) return 0;

            v = handle_pair_opposite(s1, t2, bx31, ch13, tt, t1);
            if (v == 0) continue; else if (v < 0) return 0;

//            if (!wanna_stop && refinement_state == REFINING) break;
            
            v = handle_pair_same    (t2, r1, bx12, ch21, rr, t1);
            if (v == 0) continue; else if (v < 0) return 0;

            v = handle_pair_opposite(r1, s2, bx23, ch32, ss, t1);
            if (v == 0) continue; else if (v < 0) return 0;

            break;
        }
//        if (refinement_state == STARTING) refinement_state = REFINING;


        // t1 r2 s1 t2 computed successfully
        if (!update_interval(t1, t1, lower(t2), upper(t2), tt))
            return 0;

//        if (refinement_state != REFINED) {
//            if (!update_interval(t1, t1, lower(t2), upper(t2), tt)) {
//                wanna_stop = true;
//                if (refinement_state == REFINING) refinement_state = REFINED;
//            }
//            if (tt.is_finite()) {
//                int new_prec = CGAL::get_significant_bits(get_t());
//                if (new_prec > prec_achieved) prec_achieved = new_prec;
//                else refinement_state = REFINED;
//            }
//            if (!wanna_stop && refinement_state == REFINING || refinement_state == REFINED) return 1;
//        }

        // r1 s2 computed successfully
        if (!update_interval(lower(r1), upper(r1), 
                            lower(r2), upper(r2), rr)) return 0;
        if (!update_interval(lower(s1), upper(s1), 
                            lower(s2), upper(s2), ss)) return 0;


        if (!normalized) {
            normalized = tt.is_finite() && rr.is_finite() && ss.is_finite();
//            if (normalized) std::cerr << "normalized @ iters = " << iters << std::endl;
        }

        if (tt.is_finite()) {
            int new_prec = CGAL::get_significant_bits(get_t());
            if (new_prec > 1) {
                if (new_prec > prec_achieved) prec_achieved = new_prec;
                else return 0;
            }
        }
        return 1;
    }

    IT boundary_x(const IT& t) const
    {
        IT t1, t2, t3, t4;
        IT a, b, w, xc, yc;
        
        a = ET::to_interval(triplet.get_e1().major_axis()); 
        b = ET::to_interval(triplet.get_e1().minor_axis()); 
        w = ET::to_interval(triplet.get_e1().rotation()); 
        xc = ET::to_interval(triplet.get_e1().x_center()); 
        yc = ET::to_interval(triplet.get_e1().y_center());
        // maple gen code
        t4 = -a+xc;
        t3 = a+xc;
        t2 = t*t;
        t1 = w*w;
        return (-IT(4)*b*w*t+t4*t2+(t3*t2+t4)*t1+t3)/((IT(1)+t1)*(IT(1)+t2));
    }
    
    IT boundary_y(const IT& t) const
    {
        IT t1, t2, t3, t4, t5;
        IT a, b, w, xc, yc;
        
        a = ET::to_interval(triplet.get_e1().major_axis()); 
        b = ET::to_interval(triplet.get_e1().minor_axis()); 
        w = ET::to_interval(triplet.get_e1().rotation()); 
        xc = ET::to_interval(triplet.get_e1().x_center()); 
        yc = ET::to_interval(triplet.get_e1().y_center());
        // maple gen code
        t2 = t*t;
        t5 = IT(1)+t2;
        t4 = b*t;
        t3 = t5*yc;
        t1 = w*w;
        return -(-IT(2)*t4+(IT(2)*t2-IT(2))*w*a+(IT(2)*t4-t3)*t1-t3)/((IT(1)+t1)*t5);
    }
    
    void init() {
#ifdef VORELL_VCAPX_PROFILE
        tm_init.start();
#endif
        iters = 0;
        prec_achieved = -1;
//        refinement_state = STARTING;
        bx12 = Bisector_curve_apx(triplet.get_e1(), triplet.get_e2());
        bx23 = Bisector_curve_apx(triplet.get_e2(), triplet.get_e3());
        bx31 = Bisector_curve_apx(triplet.get_e3(), triplet.get_e1());

        // TODO: adaptive precision; also when used in int. bitang. bx12
//        tt = ET::to_arx(triplet.t_range());
//        rr = ET::to_arx(triplet.r_range());
//        ss = ET::to_arx(triplet.s_range());
//        ch12 = ET::to_arx(triplet.get_bt12().CH_range());
//        ch21 = ET::to_arx(triplet.get_bt21().CH_range());
//        ch32 = ET::to_arx(triplet.get_bt32().CH_range());
//        ch13 = ET::to_arx(triplet.get_bt13().CH_range());
        tt = ET::to_arx(triplet.t_range(), false, 53);
        rr = ET::to_arx(triplet.r_range(), false, 53);
        ss = ET::to_arx(triplet.s_range(), false, 53);
        ch12 = ET::to_arx(triplet.get_bt12().CH_range(), false, 53);
        ch21 = ET::to_arx(triplet.get_bt21().CH_range(), false, 53);
        ch32 = ET::to_arx(triplet.get_bt32().CH_range(), false, 53);
        ch13 = ET::to_arx(triplet.get_bt13().CH_range(), false, 53);

        normalized = false;
        
#if VERBOSE > 2
        std::cerr <<"e1 = " << triplet.get_e1() << std::endl;
        std::cerr <<"e2 = " << triplet.get_e2() << std::endl;
        std::cerr <<"e3 = " << triplet.get_e3() << std::endl;
        std::cerr <<"tt = " << tt << std::endl;
        std::cerr <<"rr = " << rr << std::endl;
        std::cerr <<"ss = " << ss << std::endl;
        std::cerr <<"ch12 = " << ch12 << std::endl;
        std::cerr <<"ch21 = " << ch21 << std::endl;
        std::cerr <<"ch32 = " << ch32 << std::endl;
        std::cerr <<"ch13 = " << ch13 << std::endl;
#endif

        CGAL_assertion( !tt.is_empty() );
        CGAL_assertion( !ss.is_empty() );
        CGAL_assertion( !ss.is_empty() );
#ifdef VORELL_VCAPX_PROFILE
        tm_init.stop();
#endif

#ifndef VORELL_VCAPX_NOCACHE
        read_cache();

        if (!cached) {
#endif
            (*this)(); // operator() on self
            CGAL_assertion( normalized );


#ifndef VORELL_VCAPX_NOCACHE
            write_cache();
        }
#endif
#ifdef VORELL_VCAPX_PROFILE
#endif
    }

#ifndef VORELL_VCAPX_NOCACHE
    void read_cache() {
        CKey k = triplet.key();
        cached = false;
        typename Cache::Found_value fv;
        fv = cache.read(k);
        if (fv.first) {
            tt = fv.second.get<0>();
            rr = fv.second.get<1>();
            ss = fv.second.get<2>();
            normalized = true;
//                normalized = it->second.get<3>();
            cached = true;
        } else {
        // check for CCW symmetries
            k = ::boost::make_tuple(triplet.get_e2().get_id(), triplet.get_e3().get_id(), triplet.get_e1().get_id());
            fv = cache.read(k);
            if (fv.first) {
                tt = fv.second.get<2>();
                rr = fv.second.get<0>();
                ss = fv.second.get<1>();
                normalized = true;
//                    normalized = it->second.get<3>();
                cached = true;
            } else {
                k = ::boost::make_tuple(triplet.get_e3().get_id(), triplet.get_e1().get_id(), triplet.get_e2().get_id());
                fv = cache.read(k);
                if (fv.first) {
                    tt = fv.second.get<1>();
                    rr = fv.second.get<2>();
                    ss = fv.second.get<0>();
                    normalized = true;
//                        normalized = it->second.get<3>();
                    cached = true;
                }
            }
        }
    }

    void write_cache() {
        CKey k = triplet.key();
        CValue v = ::boost::make_tuple(tt, rr, ss);
        cache.write(k,v);
        cached = true;
    }
#endif

    void operator()() {
#ifdef VORELL_VCAPX_PROFILE
        tm_refine.start();
#endif
        while (refine_t() > 0) {
#if VERBOSE > 0
            std::cerr << "trange = [" << to_double(tt.left()) << ',' << to_double(tt.right()) << ']';
            if (tt.is_finite()) std::cerr << "@" << CGAL::get_significant_bits(get_t());
            std::cerr << std::endl;
#endif
            ;
        }
#ifdef VORELL_VCAPX_PROFILE
        tm_refine.stop();
#endif
#if VERBOSE > 0
        info_all();
#endif
    }

public:    

#ifdef VORELL_VCAPX_PROFILE
    CGAL::Timer tm_init, tm_refine;
#endif

    Voronoi_circle_apx() { }
    Voronoi_circle_apx(const Ellipse_triplet<ET>& triplet_): triplet(triplet_) {
        init();
    }
    
    Voronoi_circle_apx(const Ellipse_2<ET>& e1, const Ellipse_2<ET>& e2, 
                       const Ellipse_2<ET>& e3): triplet(e1,e2,e3) {
        init();
    }

#if 0
    // TODO: check correctness, check usability (do we need this?)
    Voronoi_circle_apx(const Voronoi_circle_apx<ET>& vc_): triplet(vc_.triplet) {
        // is this proper?
//        *this = Voronoi_circle_apx(vc_.triplet);
        
        iters = vc_.iters;
        prec_achieved = vc_.prec_achieved;
        bx12 = Bisector_curve_apx(triplet.get_e1(), triplet.get_e2());
        bx23 = Bisector_curve_apx(triplet.get_e2(), triplet.get_e3());
        bx31 = Bisector_curve_apx(triplet.get_e3(), triplet.get_e1());

        ch12 = ET::to_arx(triplet.get_bt12().CH_range());
        ch21 = ET::to_arx(triplet.get_bt21().CH_range());
        ch32 = ET::to_arx(triplet.get_bt32().CH_range());
        ch13 = ET::to_arx(triplet.get_bt13().CH_range());

        tt = ARX(BT(1.0)*vc_.t_range().left(), BT(1)*vc_.t_range().right());
        rr = ARX(BT(1.0)*vc_.r_range().left(), BT(1)*vc_.r_range().right());
        ss = ARX(BT(1.0)*vc_.s_range().left(), BT(1)*vc_.s_range().right());
        normalized = vc_.normalized;

         (*this)(); // operator() on self
         CGAL_assertion( normalized );
    }
#endif

    ARX t_range() const { return tt; }
    ARX r_range() const { return rr; }
    ARX s_range() const { return ss; }
    
    IT get_t() const { return make_interval(tt.left(), tt.right()); }
    IT get_r() const { return make_interval(rr.left(), rr.right()); }
    IT get_s() const { return make_interval(ss.left(), ss.right()); }
    
    bool is_normalized() const { return normalized; }
    int  iterations() const { return iters; }
    
    void get_coords(IT& x, IT& y, IT& r) const { 
        int status = bx12.get_coords(get_t(), get_r(), x, y);
        IT ex = boundary_x(get_t());
        IT ey = boundary_y(get_t());
        r = CGAL::sqrt((ex-x)*(ex-x)+(ey-y)*(ey-y));
    }

    friend std::ostream& operator<<(std::ostream& o, 
                                        const Voronoi_circle_apx<ET>& v) {
        return (o << "VC(" << v.get_t() << ',' << v.get_r() << ','
                << v.get_s() << ")@(" << CGAL::get_significant_bits(v.get_t()) << ")");
    }

    // a.k.a. compare_with
    Sign order_on_boundary(Voronoi_circle_apx<ET>& v) const {
        CGAL_assertion( triplet.get_e1() == v.get_triplet().get_e1() );
        ARX ch;
        if (triplet.get_e2() == v.get_triplet().get_e2() ||
            triplet.get_e2() == v.get_triplet().get_e3()) ch = ch12;
        else if (triplet.get_e3() == v.get_triplet().get_e2() ||
                 triplet.get_e3() == v.get_triplet().get_e3()) ch = ch13;
        else {
            CGAL_assertion (false);
            ch = ch12;
        }
#if 0
        bool ref1 = true, ref2 = true;
        while (ref1 || ref2) {
            if (ref1) ref1 = (refine_t() > 0);
            if (ref2) ref2 = (v.refine_t() > 0);
            if ((!is_normalized()) || (!v.is_normalized())) continue;
            // if (tt.is_infinite() && v.t_range().is_infinite()) continue; 
#if VERBOSE > 0
            std::cerr << "comparing " << tt << " and " << v.t_range() << " in " << ch << std::endl;
#endif
            if (!tt.subset_of(ch)) continue;
//            if (!v.t_range().subset_of(ch12)) continue;
            if (tt.intersection(v.t_range()).is_empty())
                return static_cast<Sign>(ch.order(tt.left(), v.t_range().left()));
        }
#endif
        if (tt.intersection(v.t_range()).is_empty()) {
            return static_cast<Sign>(ch.order(tt.left(), v.t_range().left()));
        } else {
            std::cerr << "possible INCIRCLE degeneracy" << std::endl;
            return CGAL::ZERO;
        }
    }
    
    const typename ET::Root& isolate_external_circle(rit rbegin, rit rend) const {
//        (*this)();
        // TODO: iterate with increased precision when needed
        std::pair<rit,rit> er = std::equal_range(rbegin, rend, 
                        make_interval(tt.left(), tt.right()), Ordering_interval_root());
        rbegin = er.first;
        rend = er.second;
        if (er.second == er.first + 1) return *er.first;
        if (er.second == er.first) {
            std::cerr << "could not select external circle" << std::endl;
            return *er.first;
        }
        CGAL_postcondition( false );
    }

    const Ellipse_triplet<ET>& get_triplet() const { return triplet; }
};

} // VORELL


} //namespace CGAL
#endif
