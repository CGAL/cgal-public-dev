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

#ifndef CGAL_IN_CIRCLE_H
#define CGAL_IN_CIRCLE_H

#include<iostream>
#include<CGAL/basic.h>
#include<CGAL/vorell.h>

#include<CGAL/Ellipse_2.h>
#include<CGAL/Ellipse_triplet.h>
#include<CGAL/Range.h>
#include<CGAL/Bitangent.h>
#include<CGAL/Visible_arc.h>
#include<CGAL/Voronoi_circle.h>

#include <CGAL/vorell_profile.h>

namespace CGAL {


template<class ET>
class In_circle {
    typedef VORELL::Range<typename ET::Root> Range;
    typedef typename ET::IT IT;
    typedef typename ET::tpoli_t tpoli_t;
    typedef typename ET::tpolz_t tpolz_t;
    typedef typename ET::bpoli_t bpoli_t;
    typedef typename ET::bpolz_t bpolz_t;

public:
    
    // POSITIVE = no conflict
    Sign operator()(const Ellipse_2<ET>& e1, const Ellipse_2<ET>& e2, 
                        const Ellipse_2<ET>& e3, const Ellipse_2<ET>& e4) {
        PROF_DECL;
        PROF_START;
        Ellipse_triplet<ET> triplet(e1,e2,e3);

#if VERBOSE > 0
        std::cerr << "IN_CIRCLE(" << e1.get_id() << ',' << e2.get_id() << ',' << e3.get_id() << ',' << e4.get_id() << ')' << std::endl;
#endif

        // existence
#if VERBOSE > 0
        std::cerr << "NVC1 = " << triplet.get_num_voronoi_circles() << std::endl;
#endif

        if (triplet.get_num_voronoi_circles() < 1) return CGAL::NEGATIVE;

#if VERBOSE > 0
        std::cerr << "trange123 = " << triplet.t_range() << std::endl;
#endif
        
        PROF("InCircle: initial ranges etc.: ");
        Voronoi_circle<ET> vc(triplet);
        
        Ellipse_triplet<ET> triplet2(e1,e2,e4);

        int nvc2 = triplet2.get_num_voronoi_circles();
        // TODO: check again [checked 3 times]
        //bool inside2 = triplet2.is_third_intersecting_CH_interior();
        bool src2 = triplet2.is_shadow_region_connected();

#if VERBOSE > 0
        std::cerr << "SRC2(inside) = " << src2 << " NVC2 = " << nvc2 << std::endl;
#endif

        if (nvc2 == 0) {
            if (src2) return CGAL::NEGATIVE;
            else return CGAL::POSITIVE;
        } else if (nvc2 == 1) {
            Voronoi_circle<ET> vc2(triplet2);
            Sign r = vc2.order_on_boundary(vc);
            //std::cerr << "ORDER ON BOUNDARY = " << r << std::endl;
            return r;
        } else if (nvc2 == -1) {
            triplet2 = Ellipse_triplet<ET>(e1,e4,e2);
            Voronoi_circle<ET> vc2(triplet2);
            Sign r = vc2.order_on_boundary(vc);
            //std::cerr << "ORDER ON BOUNDARY = " << r << std::endl;
            return -r;
        } else { // nvc2 == 2
            Voronoi_circle<ET> vc2(triplet2);
            Sign r1 = vc2.order_on_boundary(vc);
            if (r1 == CGAL::ZERO) return r1;
            if (!src2 && r1 == CGAL::NEGATIVE) return r1;
            triplet2 = Ellipse_triplet<ET>(e1,e4,e2);
            vc2 = Voronoi_circle<ET>(triplet2);
            Sign r2 = vc2.order_on_boundary(vc);
            if (r2 == CGAL::ZERO) return r2;
            if (src2) return r1 * r2;
            else return -r2;
        }
    }
};

} //namespace CGAL
#endif
