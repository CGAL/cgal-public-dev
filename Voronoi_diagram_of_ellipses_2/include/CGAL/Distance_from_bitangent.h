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

#ifndef CGAL_DISTANCE_FROM_BITANGENT_H
#define CGAL_DISTANCE_FROM_BITANGENT_H

#include<iostream>
#include<CGAL/basic.h>
#include<CGAL/vorell.h>

#include<CGAL/Ellipse_2.h>
#include<CGAL/Bitangent.h>
#include<CGAL/Voronoi_diagram_of_ellipses_2/Generated_code.h>

namespace CGAL {


template<class EllipseTraits>
class Distance_from_bitangent {
    typedef EllipseTraits ET;
    typedef typename ET::QQ QQ;
    typedef typename ET::AK AK;
    typedef typename ET::upolz_t upolz_t;
    typedef typename ET::upolq_t upolq_t;
    typedef typename AK::Algebraic_real_1 Root;
    typedef typename AK::Sign_at_1 SignAt;
    typedef VORELL::Bitangent<ET> Bitangent;
    typedef CGAL::Ellipse_2<ET> Ellipse_2;
    typedef CGAL::VORELL::Generated_code<ET> gencode;
    
    Bitangent bt12, bt13;

public:
    
    Distance_from_bitangent() { }
        
    Distance_from_bitangent(const Bitangent& bt12_, const Bitangent& bt13_):
            bt12(bt12_), bt13(bt13_) { }

    Sign evaluate() const {
        CGAL_assertion(bt12.relative_position() >= VORELL::PSEUDO_CIRCLES);
        
        Bounded_side relpos3 = bt12.relative_position_of_ellipse(bt13);
        if (relpos3 == CGAL::ON_BOUNDED_SIDE) return CGAL::NEGATIVE;
        
        if (bt13.relative_position() == VORELL::HIDDEN ||
            bt13.relative_position() == VORELL::HIDING) return CGAL::POSITIVE;
        
        Root tanpoint = bt12.external();
        
        CGAL_assertion_code( upolz_t tan11 = \
                gencode().tan_poly_xy(ELL_PARAM_COEFFS(bt12.get_e1()), \
                bt12.get_e1().x_center(), bt12.get_e1().y_center()); )
        CGAL_assertion( SignAt()(tan11, tanpoint) == CGAL::POSITIVE );
        
	// TODO: check if this is faster, or we can re-use some info from relpos3
	// (ordering of tan. points)
        upolz_t tan13 = 
                gencode().tan_poly_xy(ELL_PARAM_COEFFS(bt13.get_e1()), 
                bt13.get_e2().x_center(), bt13.get_e2().y_center());
        Sign side = SignAt()(tan13, tanpoint);

        // now side can't be 0 ;-)    

        if (side == CGAL::POSITIVE) {
            // 1 or 0 = -conflict;
            if (relpos3 == CGAL::ON_UNBOUNDED_SIDE) return CGAL::POSITIVE;
            return CGAL::ZERO;
        }
        return CGAL::NEGATIVE;
    }
    
    Sign evaluate(const Ellipse_2& e3) {
        bt13 = Bitangent(bt12.get_e1(), e3);
        return evaluate();
    }
    
    // ccw 2-1-inf vs 3
    Sign operator()(const Ellipse_2& e1, const Ellipse_2& e2, 
                    const Ellipse_2& e3) {
        bt12 = Bitangent(e1, e2);
        bt13 = Bitangent(e1, e3);
        return evaluate();
    }
};

} //namespace CGAL
#endif
