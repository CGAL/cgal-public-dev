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

#ifndef CGAL_ELLIPSE_TRIPLET_H
#define CGAL_ELLIPSE_TRIPLET_H

#include<iostream>
#include<CGAL/basic.h>
#include<CGAL/vorell.h>

#include<CGAL/Range.h>
#include<CGAL/Ellipse_2.h>
#include<CGAL/Voronoi_diagram_of_ellipses_2/Generated_code.h>

#include<CGAL/Bitangent.h>
#include<CGAL/Visible_arc.h>
#include<CGAL/Distance_from_bitangent.h>
#include<CGAL/Medial_axis_location.h>

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include<CGAL/Object_cache.h>

namespace CGAL {

template<class ET>
class Ellipse_triplet {
    typedef typename ET::QQ QQ;
    typedef typename ET::AK AK;
    typedef typename ET::upolz_t upolz_t;
    typedef typename ET::upolq_t upolq_t;
    typedef typename AK::Algebraic_real_1 Root;
    typedef typename AK::Sign_at_1 SignAt;
    typedef typename VORELL::Bitangent<ET> Bitangent;
    typedef typename VORELL::Range<Root> Range;

public:
    typedef ::boost::tuple<int, int, int> Key_type;
    inline Key_type key() const {
        return ::boost::make_tuple(e1.get_id(), e2.get_id(), e3.get_id());
    }

private:

    typedef int CValue;
//    typedef ::std::map<Key_type, CValue> MA_Cache;
    typedef Value_cache<Key_type, CValue> MA_Cache;
    MA_Cache ma_cache;

    Ellipse_2<ET> e1, e2, e3;
    Bitangent bt12, bt21, bt13, bt31, bt23, bt32;
    Range trange, rrange, srange;
    int num_voronoi_circles;
//    bool third_intersects_CH_interior;
    bool shadow_region_connected;
    int cover_index; // 0 = no_cover
    int shadow_region_type;

    bool ccw_maybe_internal;
    
    bool computed_bitangents;
    bool computed_ranges;
    bool computed_nvc;

    typedef Object_cache<Ellipse_triplet> Cache;
    Cache cache;
    
    void compute_bitangents() {
        if (computed_bitangents) return;

        bt12 = Bitangent(e1, e2); 
        bt21 = Bitangent(e2, e1);
        
        bt13 = Bitangent(e1, e3); 
        bt31 = Bitangent(e3, e1);
        
        bt23 = Bitangent(e2, e3); 
        bt32 = Bitangent(e3, e2);

        computed_bitangents = true;
    }

    // TODO: experimental PC stuff
    // compute candidate internal CW circle
    bool compute_candidate_internal(Range &tt, Range &rr, Range &ss,
                                    Bitangent &bt12, Bitangent &bt21,
                                    Bitangent &bt23, Bitangent &bt32,
                                    Bitangent &bt31, Bitangent &bt13) {
        if (bt12.relative_position() == VORELL::PSEUDO_CIRCLES &&
            bt23.relative_position() == VORELL::PSEUDO_CIRCLES &&
            bt31.relative_position() == VORELL::PSEUDO_CIRCLES) {
            // if all pairs intersect, check for possible INTERNAL circle
            tt = bt12.internal_range().intersection(bt13.internal_range());
#if VERBOSE > 0
            std::cerr << "INT:trange = " << tt << std::endl;
#endif
            if (!tt.is_empty()) rr = bt23.internal_range().intersection(bt21.internal_range());
#if VERBOSE > 0
            std::cerr << "INT:rrange = " << rr << std::endl;
#endif
            if (!tt.is_empty() && !rr.is_empty())
                ss = bt31.internal_range().intersection(bt32.internal_range());
#if VERBOSE > 0
            std::cerr << "INT:srange = " << ss << std::endl;
#endif
            if ((!tt.is_empty()) && (!rr.is_empty()) && (!ss.is_empty())) {
                if (tt.left() == bt12.internal() ||
                    rr.left() == bt23.internal() ||
                    ss.left() == bt31.internal()) {
                    return true; // internal, unless covered
                }
            }
        }
#if VERBOSE > 0
        std::cerr << "NOT INTERNAL" << std::endl;
#endif
        return false;
    }

    void compute_ranges() {
        if (computed_ranges) return;
        bool maybe_internal_ccw = false;
        Range tt, rr, ss;

        compute_bitangents();
        int nvc = get_num_voronoi_circles();

        computed_ranges = true;
        ccw_maybe_internal = false;

        if (nvc == 1) {
            maybe_internal_ccw = compute_candidate_internal(tt, ss, rr, bt13, bt31, bt32, bt23, bt21, bt12);
            if (maybe_internal_ccw) {
                trange = tt; rrange = rr; srange = ss;
                ccw_maybe_internal = true;
                return;
            }
        } else if (nvc == 2) {
            // TODO: ATTENTION direction cw/ccw reverses when both circles internal ??? !!!
            maybe_internal_ccw = compute_candidate_internal(tt, ss, rr, bt13, bt31, bt32, bt23, bt21, bt12);
            if (maybe_internal_ccw) {
                Range tt2, rr2, ss2;
                bool maybe_internal_cw = compute_candidate_internal(tt2, rr2, ss2, bt12, bt21, bt23, bt32, bt31, bt13);
                if (maybe_internal_cw) {
                    trange = tt2; rrange = rr2; srange = ss2;
                    ccw_maybe_internal = true;
                    return;
                }
                trange = tt; rrange = rr; srange = ss;
                ccw_maybe_internal = true;
                return;
            }
        }
        // TODO: needed to check for nvc = 0? do we need the ranges then?

        trange = VORELL::Visible_arc<ET>(bt12, bt13);
        rrange = VORELL::Visible_arc<ET>(bt23, bt21);
        srange = VORELL::Visible_arc<ET>(bt31, bt32);

#if VERBOSE > 0
        std::cerr << "BEFORE PRUNING" << std::endl;
        std::cerr << "trange = " << trange << " [" << e1.get_id() << "]" << std::endl;
        std::cerr << "rrange = " << rrange << " [" << e2.get_id() << "]" << std::endl;
        std::cerr << "srange = " << srange << " [" << e3.get_id() << "]" << std::endl;
#endif

        if (bt12.relative_position() == VORELL::PSEUDO_CIRCLES) {
            trange = trange.prune(bt12.internal_range());
            rrange = rrange.prune(bt21.internal_range());
        }
        if (bt23.relative_position() == VORELL::PSEUDO_CIRCLES) {
            rrange = rrange.prune(bt23.internal_range());
            srange = srange.prune(bt32.internal_range());
        }
        if (bt31.relative_position() == VORELL::PSEUDO_CIRCLES) {
            srange = srange.prune(bt31.internal_range());
            trange = trange.prune(bt13.internal_range());
        }


        if (bt12.relative_position() == VORELL::PSEUDO_CIRCLES) {
            // TODO: experimental PC stuff
            int nvc = get_num_voronoi_circles();
            if (nvc == 1) {
                if (e3.boundary_relpos(e1, bt12.internal()) == ON_BOUNDED_SIDE) {
                    trange = trange.intersection(bt12.first_range());
                } else {
                    trange = trange.intersection(bt12.last_range());
                }
            } else if (nvc == 2) {
                if (shadow_region_connected) {
                    if (e3.boundary_relpos(e1, bt12.internal()) == ON_BOUNDED_SIDE) {
                        trange = trange.intersection(bt12.first_range());
                    } else {
                        // TODO: handle separately all touching cases, when 1-2=PC 2-3=PC, but 3-1 not etc...
                        if (bt13.relative_position() == VORELL::PSEUDO_CIRCLES) {
                            if (bt13.internal_range().subset_of(bt12.first_range())) {
                                trange = trange.intersection(bt12.first_range());
                            } else if (bt13.internal_range().subset_of(bt12.last_range())) {
                                trange = trange.intersection(bt12.last_range());
                            } else {
                                CGAL_assertion(false);
                            }
                        } else if (bt23.relative_position() == VORELL::PSEUDO_CIRCLES) {
                            if (bt23.internal_range().subset_of(bt21.last_range())) {
                                trange = trange.intersection(bt12.first_range());
                            } else if (bt23.internal_range().subset_of(bt21.first_range())) {
                                trange = trange.intersection(bt12.last_range());
                            } else {
                                CGAL_assertion(false);
                            }
                        } else { // e3 not intersecting e1,e2
                            if (bt13.first_range().strict_subset_of(bt12.first_range())) {
                                trange = trange.intersection(bt12.first_range());
                            } else if (bt13.last_range().strict_subset_of(bt12.last_range())) {
                                trange = trange.intersection(bt12.last_range());
                            } else {
                                CGAL_assertion(false);
                            }
                        }
                    }
                } else {
                    trange = trange.intersection(bt12.last_range());
                }
            }

        }

#if VERBOSE > 0
        std::cerr << "AFTER (POSSIBLE) PRUNING" << std::endl;
        std::cerr << "trange = " << trange << std::endl;
        std::cerr << "rrange = " << rrange << std::endl;
        std::cerr << "srange = " << srange << std::endl;
#endif

    }
    
    bool is_medial_axis_covered(Bitangent &bt12, Bitangent &bt13) const {
        int r1 = 0, r2 = 0;
        int pos = bt12.relative_position();
        if (pos == VORELL::HIDDEN) r1 = 1;
        else if (pos == VORELL::SEPARATED || pos == VORELL::HIDING) r1 = -1;
        if (r1 > 0) return true;
        pos = bt13.relative_position();
        if (pos == VORELL::HIDDEN) r2 = 1;
        else if (pos == VORELL::SEPARATED || pos == VORELL::HIDING) r2 = -1;
        if (r2 > 0) return true;
        if (r1 < 0 || r2 < 0) return false;
        // pseudo-circles
        VORELL::Medial_axis_location<ET> ma12(bt12.get_e1(), bt12.get_e2(), bt12.internal_range());
        if (ma12().is_empty()) return false;
        VORELL::Medial_axis_location<ET> ma13(bt13.get_e1(), bt13.get_e2(), bt13.internal_range());
        if (ma13().is_empty()) return false;

        return ma12().complement().intersection(ma13().complement()).is_empty();
    }

    int get_nvc() {
        int b1, b2, b3, k;

        VORELL::Relative_position p12, p13, p23;
        p12 = bt12.relative_position();
        p13 = bt13.relative_position();
        p23 = bt23.relative_position();
        cover_index = 0; // obviously false if following test holds :)
        if (p12 <= VORELL::HIDDEN || p13 <= VORELL::HIDDEN ||
            p23 <= VORELL::HIDDEN) return 0;
        // CCW
        // t,r,s
        b1 = (Distance_from_bitangent<ET>(bt21, bt23).evaluate() == NEGATIVE);
        // r,s,t
        b2 = (Distance_from_bitangent<ET>(bt32, bt31).evaluate() == NEGATIVE);
        // s,t,r
        b3 = (Distance_from_bitangent<ET>(bt13, bt12).evaluate() == NEGATIVE);
        k = b1 + b2 + b3;
        if (k == 3) return 1;
        if (k == 2) {
            if (!read_ma_cache()) {
                bool c1, c2, c3;
                c1 = is_medial_axis_covered(bt12, bt13);
                c2 = is_medial_axis_covered(bt21, bt23);
                c3 = is_medial_axis_covered(bt31, bt32);
                if (c1) cover_index = 1;
                else if (c2) cover_index = 2;
                else if (c3) cover_index = 3;
                write_ma_cache();
            } else {
#if VERBOSE > 1
                std::cerr << "MA-loc cache hit!" << std::endl;
#endif
            }
#if VERBOSE > 1
            std::cerr << "MA-loc for triplet (" << e1.get_id() << ',' << e2.get_id() << ',' << e3.get_id() << ')' << std::endl;
#endif
            if (cover_index > 0) {
#if VERBOSE > 1
                std::cerr << "MA-loc: covering detected between "
                          <<  e1.get_id() << ',' << e2.get_id() << ',' << e3.get_id() << " index = " << cover_index << std::endl;
#endif
                return 0;
            }
            return 2;
        }
        if (k == 1) return 0;
        return -1;
    }

//     return value:
//        0 = no cw/ccw circle exists
//        1 = only 1 (ccw) circle exists
//       -1 = only 1 (cw) circle exists
//        2 = both cw & ccw circles exist

    void compute_nvc() {
        if (computed_nvc) return;
        
        //bool l1, l2;

        compute_bitangents();
        
        num_voronoi_circles = get_nvc();

        if (num_voronoi_circles == 0) {
            if (cover_index == 1) {
                shadow_region_connected = true;
                shadow_region_type = ENTIRE_EDGE;
            } else if (cover_index > 1) {
                shadow_region_connected = false;
                shadow_region_type = EMPTY;
            } else {
                shadow_region_connected = (Distance_from_bitangent<ET>(bt21, bt23).evaluate() == NEGATIVE);
                if (shadow_region_connected) shadow_region_type = ENTIRE_EDGE;
                else shadow_region_type = EMPTY;
            }
        } else if (num_voronoi_circles == 2) {
            // TODO: verify != NEGATIVE or == POSITIVE
            shadow_region_connected = (Distance_from_bitangent<ET>(bt21, bt23).evaluate() != NEGATIVE);
            if (shadow_region_connected) shadow_region_type = INTERIOR;
            else shadow_region_type = TWO_VERTICES;
        } else {
            shadow_region_connected = true;
            if (num_voronoi_circles == 1) shadow_region_type = LEFT_VERTEX;
            else shadow_region_type = RIGHT_VERTEX;
        }
    }

    bool read_ma_cache() {
        Key_type k = key();
        typename MA_Cache::Found_value fv;
        fv = ma_cache.read(k);
        if (fv.first) {
            cover_index = fv.second;
            return true;
        }
        // check for CCW symmetries
        k = ::boost::make_tuple(e2.get_id(), e3.get_id(), e1.get_id());
        fv = ma_cache.read(k);
        if (fv.first) {
            cover_index = fv.second;
            if (cover_index > 0) cover_index = cover_index % 3 + 1;
            return true;
        }
        k = ::boost::make_tuple(e3.get_id(), e1.get_id(), e2.get_id());
        fv = ma_cache.read(k);
        if (fv.first) {
            cover_index = fv.second;
            if (cover_index > 0) {
                cover_index = cover_index % 3 + 1;
                cover_index = cover_index % 3 + 1;
            }
            return true;
        }
        return false;
    }

    void write_ma_cache() {
        Key_type k = key();
        CValue v = cover_index;
        ma_cache.write(k,v);
    }

public:
    // shadow region type
    enum { EMPTY, ENTIRE_EDGE, RIGHT_VERTEX, LEFT_VERTEX, INTERIOR, TWO_VERTICES };
    // (), (0,1), (0,a), (b,1), (a,b), (0,a) + (b,1)

    Ellipse_triplet() { }

    Ellipse_triplet(const Ellipse_2<ET>& e1_, const Ellipse_2<ET>& e2_, 
                    const Ellipse_2<ET>& e3_): e1(e1_), e2(e2_), e3(e3_)
    {
        const Ellipse_triplet *h = cache.read(*this);
        if (h) {
            *this = *h;
            return;
        }
        computed_bitangents = false;
        computed_nvc = false;
        computed_ranges = false;

        compute_bitangents();
        compute_nvc();
        compute_ranges();

        cache.write(*this);
    }
    
    int get_num_voronoi_circles() const {
        return num_voronoi_circles;
    }
    
    bool is_shadow_region_connected() const {
        return shadow_region_connected;
    }

    int get_shadow_region_type() const {
        return shadow_region_type;
    }

    int get_cover_index() const {
        return cover_index;
    }

    bool is_ccw_possibly_internal() const {
        return ccw_maybe_internal;
    }
    
    const Range& t_range() const {
        return trange;
    }
    
    const Range& r_range() const {
        return rrange;
    }
    
    const Range& s_range() const {
        return srange;
    }
    
    const VORELL::Bitangent<ET>& get_bt12() const {
        return bt12;
    }
    
    const VORELL::Bitangent<ET>& get_bt21() const {
        return bt21;
    }
    
    const VORELL::Bitangent<ET>& get_bt32() const {
        return bt32;
    }
    
    const VORELL::Bitangent<ET>& get_bt13() const {
        return bt13;
    }

    const Ellipse_2<ET>& get_e1() const { return e1; }
    const Ellipse_2<ET>& get_e2() const { return e2; }
    const Ellipse_2<ET>& get_e3() const { return e3; }

};

} //namespace CGAL
#endif
