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

#ifndef CGAL_ELLIPSE_BISECTOR_2_H
#define CGAL_ELLIPSE_BISECTOR_2_H

#include <CGAL/enum.h>

#include <CGAL/Ellipse_2.h>
#include <CGAL/Bisector_curve_apx.h>
#include <CGAL/Voronoi_circle_apx.h>

#include <CGAL/Object_cache.h>

#define ELLIPSE_BISECTOR_DEFAULT_RES 32

namespace CGAL {

template < class ET >
class Ellipse_bisector_2 {
    typedef Ellipse_2<ET> Site_2;
    typedef typename ET::Kernel_scr::Point_2 Point_2;
    typedef typename ET::Kernel_scr::Segment_2 Segment_2;
    typedef typename ET::BT BT;
    typedef typename ET::IT IT;
    typedef typename ET::ARX ARX;
    

    typedef ::std::pair<int, int> CKey2;
    typedef ::boost::tuple<int, int, int> CKey3;
    typedef ::boost::tuple<int, int, int, int> CKey4;
    typedef Ellipse_bisector_2<ET> CValue;
    typedef Value_cache<CKey2, CValue> Cache2;
    typedef Value_cache<CKey3, CValue> Cache3;
    typedef Value_cache<CKey4, CValue> Cache4;
    Cache2 cache2;
    Cache3 cache3;
    Cache4 cache4;

    int res; // draw resolution

    VORELL::Bisector_curve_apx<ET> bx12;
    ARX t_arx;
    BT last_trace;
    
    bool add_point(const BT t) {
        IT x, y;
        IT ri = bx12(t);
#if VERBOSE > 1
            std::cerr << "t = " << t << " r = " << ri << std::endl;
            std::cerr << bx12.get_status() << std::endl;
#endif
//        CGAL_assertion( CGAL::get_significant_bits(ri) > 20 );
        if (bx12.get_status() < 0) {
#if VERBOSE > 2
            std::cerr << "bad arc point (" << t << ") " << bx12.get_status() << std::endl;
#endif
            return false;
        }
        if (bx12.get_coords(t, ri, x, y) == 0) {
            xx.push_back(CGAL::to_double(x));
            yy.push_back(CGAL::to_double(y));
            last_trace = t;
#if VERBOSE > 1
            std::cerr << "x = " << x << " y = " << y << std::endl;
#endif
            return true;
        } // else std::cerr << "invalid bisector point (" << t << ") " << cs << ")\n";
        return false;
    }

    void generate_points(const ARX t_arx, const int npt) {
#if VERBOSE > 2
        std::cerr << "generate_points(" << t_arx << ',' << npt << ")\n";
#endif
        if (npt > 2) {
            generate_points(ARX(t_arx.left(), ET::midpoint(t_arx)), npt/2);
            generate_points(ARX(ET::midpoint(t_arx), t_arx.right()), npt - npt/2);
        } else {
#if VERBOSE > 2
            std::cerr << "add " << to_double(t_arx.left()) << std::endl;
            std::cerr << "add " << to_double(ET::midpoint(t_arx)) << std::endl;
#endif
            add_point(t_arx.left());
            add_point(ET::midpoint(t_arx));
        }
        if (npt == res) {
            if (!add_point(t_arx.right())) {
                BT alt = ET::midpoint(ARX(last_trace,t_arx.right()));
#if VERBOSE > 2
                std::cerr << "alternative: " << to_double(alt) << std::endl;
#endif
                add_point(alt);
            }
#if VERBOSE > 0
            std::cerr << "Ellipse_bisector_2: generated " << xx.size()
                      << " points" << std::endl;
#endif
        }
    }

    template< class Stream> void point_separator(Stream &W) const { }
    template< class Stream> void object_separator(Stream &W) const { }
    
    void point_separator(std::ostream &W) const { W << std::endl; }
    void object_separator(std::ostream &W) const { point_separator(W); }

public:
    std::vector<double> xx, yy;

    Ellipse_bisector_2(): res(0) { }

    Ellipse_bisector_2(const Site_2 &p, const Site_2 &q, 
            const int res_ = ELLIPSE_BISECTOR_DEFAULT_RES) {
#if VERBOSE > 0
        std::cerr << "BisectorCurve FOR (" << p.get_id()
                  << "," << q.get_id() << ")\n";
#endif
        CKey2 k = ::std::make_pair(p.get_id(), q.get_id());
        bool hit = false;
        typename Cache2::Found_value fv = cache2.read(k);
        if (fv.first) {
            *this = fv.second;
            hit = true;
        }
        if (hit) return;
        res = res_;
        bx12 = VORELL::Bisector_curve_apx<ET>(p, q);
        VORELL::Bitangent<ET> bt12(p, q);
        t_arx = ET::to_arx(bt12.CH_range(), true);
        generate_points(t_arx, res);
        cache2.write(k, *this);

#if VERBOSE > 0
        std::cerr << "BisectorCurve(" << t_arx << ") FOR (" << p.get_id()
                  << "," << q.get_id() << ")\n";
#endif
    }

    // ray
    Ellipse_bisector_2(const Site_2 &p, const Site_2 &q, const Site_2 &r,
                        const int res_ = ELLIPSE_BISECTOR_DEFAULT_RES) {
#if VERBOSE > 0
        std::cerr << "BisectorCurve/RAY FOR (" << p.get_id()
                  << "," << q.get_id() << "," << r.get_id() << ")\n";
#endif
        CKey3 k = ::boost::make_tuple(p.get_id(), q.get_id(), r.get_id());
        bool hit = false;
        typename Cache3::Found_value fv = cache3.read(k);
        if (fv.first) {
            *this = fv.second;
            hit = true;
        }
        if (hit) return;
        res = res_;
        bx12 = VORELL::Bisector_curve_apx<ET>(p, q);
        VORELL::Bitangent<ET> bt12(p, q);
        t_arx = ET::to_arx(bt12.CH_range(), true);

        VORELL::Voronoi_circle_apx<ET> vpqr(p, q, r);
//        vpqr();
        t_arx = typename ET::ARX(t_arx.left(), vpqr.t_range().right());
        generate_points(t_arx, res);
        cache3.write(k, *this);

#if VERBOSE > 0
        std::cerr << "BisectorCurve(" << t_arx << ")/RAY FOR (" << p.get_id()
                  << "," << q.get_id() << "," << r.get_id() << ")\n";
#endif
    }

    // segment
    Ellipse_bisector_2(const Site_2 &p, const Site_2 &q, 
                          const Site_2 &r, const Site_2 &s, 
                            const int res_ = ELLIPSE_BISECTOR_DEFAULT_RES) {
#if VERBOSE > 0
        std::cerr << "BisectorCurve/SEGMENT FOR (" << p.get_id()
                  << "," << q.get_id() << "," << r.get_id() << ","
                  << s.get_id() << ")\n";
#endif
        CKey4 k = ::boost::make_tuple(p.get_id(), q.get_id(), r.get_id(), s.get_id());
        bool hit = false;
        typename Cache4::Found_value fv = cache4.read(k);
        if (fv.first) {
            *this = fv.second;
            hit = true;
        }
        if (hit) return;
        res = res_;
        bx12 = VORELL::Bisector_curve_apx<ET>(p, q);

        VORELL::Voronoi_circle_apx<ET> vpqr(p, q, r);
//        vpqr();
        VORELL::Voronoi_circle_apx<ET> vpsq(p, s, q);
//        vpsq();
        t_arx = typename ET::ARX(vpsq.t_range().left(), vpqr.t_range().right());
        generate_points(t_arx, res);
        cache4.write(k, *this);

#if VERBOSE > 0
        std::cerr << "BisectorCurve(" << t_arx << ")/SEGMENT FOR (" << p.get_id()
                  << "," << q.get_id() << "," << r.get_id() << ","
                  << s.get_id() << ")\n";
#endif
    }

    template< class Stream > void draw(Stream &w) const {
        for (unsigned int i = 0; i < xx.size(); i++) {
            Point_2 p(xx[i], yy[i]);
            w << p;
            point_separator(w);
        }
//        if (xx.size() > 1) for (unsigned int i = 0; i < xx.size()-1; i++) {
//            Point_2 p1, p2;
//            p1 = Point_2(xx[i], yy[i]);
//            p2 = Point_2(xx[i+1], yy[i+1]);
//            w << Segment_2(p1, p2) << std::endl;
//        }
        object_separator(w);
    }

    int get_resolution() const { return res; }
    
    const VORELL::Bisector_curve_apx<ET>& get_bisector_curve_apx() const { 
        return bx12;
    }
    // int set_resolution(int res_) { return res = res_; }
};

template< class Stream, class ET >
inline Stream& operator<<(Stream& s, const Ellipse_bisector_2<ET> &b) {
    b.draw(s);
    return s;
}

} //namespace CGAL

#undef ELLIPSE_BISECTOR_DEFAULT_RES

#endif
