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

#ifndef CGAL_VORONOI_CIRCLE_H
#define CGAL_VORONOI_CIRCLE_H

#include<iostream>
#include<CGAL/basic.h>
#include<CGAL/vorell.h>

#include<CGAL/Range.h>
#include<CGAL/Ellipse_traits.h>
#include<CGAL/Ellipse_triplet.h>
#include<CGAL/Voronoi_circle_apx.h>
#include<CGAL/Voronoi_circle_exact.h>

// #include<CGAL/Timer.h>

namespace CGAL {

template<class ET>
class Voronoi_circle {
    
    Ellipse_triplet<ET> triplet;
    
    VORELL::Voronoi_circle_apx<ET> c_apx;
    VORELL::Voronoi_circle_exact<ET> c_exact;
    
    bool apx_init;
    bool exact_init;
    
    void init_apx() {
        if (!apx_init) {
            c_apx = VORELL::Voronoi_circle_apx<ET>(triplet);
            apx_init = true;
        }
    }
     
    void init_exact() {
        if (!exact_init) {
            // TODO: pass already computed circle_apx
            c_exact = VORELL::Voronoi_circle_exact<ET>(triplet);
            exact_init = true;
        }
    }

    VORELL::Voronoi_circle_apx<ET>& circle_apx() { 
        init_apx(); return c_apx; 
    }
    
    VORELL::Voronoi_circle_exact<ET>& circle_exact() { 
        init_exact(); return c_exact; 
    }
    
public:    
    Voronoi_circle() { }

    Voronoi_circle(const Ellipse_triplet<ET>& triplet_): 
        triplet(triplet_), apx_init(false), exact_init(false) { }
    
    Voronoi_circle(const Ellipse_2<ET>& e1, const Ellipse_2<ET>& e2, 
                   const Ellipse_2<ET>& e3): triplet(e1, e2, e3),
                       apx_init(false), exact_init(false)  { }
    
    Sign order_on_boundary(Voronoi_circle<ET>& v) {
        // TODO: use always exact/ adaptive prec. etc.
        init_apx();
        Sign r = c_apx.order_on_boundary(v.circle_apx());
        //Sign r = ZERO;
        if (r == CGAL::ZERO) {
            init_exact();
            r = static_cast<Sign>(
                 triplet.get_bt12().CH_range().order(
                        c_exact(), v.circle_exact()()));
        }
        return r;
    }
};
    
} //namespace CGAL
#endif

