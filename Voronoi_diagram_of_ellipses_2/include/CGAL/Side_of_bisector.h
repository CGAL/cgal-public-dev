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

#ifndef CGAL_SIDE_OF_BISECTOR_H
#define CGAL_SIDE_OF_BISECTOR_H

#include<iostream>
#include<CGAL/basic.h>
#include<CGAL/vorell.h>

#include<CGAL/Ellipse_2.h>
#include<CGAL/Voronoi_diagram_of_ellipses_2/Generated_code.h>

//#include<CGAL/algebraic_kernel_1_tools.h>

namespace CGAL {


template<class EllipseTraits>
class Side_of_bisector {
    typedef EllipseTraits ET;
    typedef typename ET::QQ QQ;
    typedef typename ET::AK AK;
    typedef typename ET::upolz_t upolz_t;
    typedef typename CGAL::VORELL::Generated_code<ET> gencode;
    typedef typename ET::Root Root;

public:
    
    // return -1 = closer to e1
    //        +1 = closer to e2
    //         0 = equidistant
    Comparison_result operator()(const Ellipse_2<ET>& e1, const Ellipse_2<ET>& e2, 
                                    const QQ& v1, const QQ& v2) {
        bool positive_distance = true;
        Bounded_side b1, b2;
       
        b1 = e1.bounded_side(v1, v2);
        b2 = e2.bounded_side(v1, v2);
        
        if (b1 == ON_UNBOUNDED_SIDE) {
            if (b2 != b1) return POSITIVE;
        } else if (b1 == ON_BOUNDED_SIDE) {
            if (b2 != b1) return NEGATIVE;
            positive_distance = false;
        } else { // b1 == ON_BOUNDARY
            if (b2 == ON_BOUNDARY) return ZERO;
            else if (b2 == ON_UNBOUNDED_SIDE) return NEGATIVE;
            else return POSITIVE;
        } 

        QQ qd1(-1), qd2(-1);
        upolz_t dpoly1, dpoly2;

        if (e1.is_on_axes(v1,v2))
            qd1 = e1.squared_distance_from_edges(v1,v2);
        else {
            dpoly1 = gencode().distance_poly(ELL_PARAM_COEFFS(e1), v1, v2);
        //        std::cerr << dpoly1 << std::endl;
            CGAL_assertion( typename ET::PTZ::Degree()(dpoly1) % 2 == 0 );
        }

        if (e2.is_on_axes(v1,v2))
            qd2 = e2.squared_distance_from_edges(v1,v2);
        else {
            dpoly2 = gencode().distance_poly(ELL_PARAM_COEFFS(e2), v1, v2);
        //        std::cerr << dpoly2 << std::endl;
            CGAL_assertion( typename ET::PTZ::Degree()(dpoly1) % 2 == 0 );
        }

        std::vector<Root> sol1, sol2;
        // actually, the smallest root is nonnegative!
        if (qd1 == -1) sol1 = ET::Real_roots(dpoly1);
        if (qd2 == -1) sol2 = ET::Real_roots(dpoly2);

//        CGAL_assertion( sol1.size() > 0 );
//        CGAL_assertion( sol2.size() > 0 );
//        CGAL_assertion( is_positive(sol1[0]) );
//        CGAL_assertion( is_positive(sol2[0]) );

        Comparison_result r;
        if (qd1 > -1) {
            if (qd2 > -1) r = compare(qd1, qd2);
            else r = compare(qd1, sol2[0]);
        } else {
            if (qd2 > -1) r = compare(sol1[0], qd2);
            else r = compare(sol1[0], sol2[0]);
        }

#if VERBOSE > 2
        std::cerr << "Side_of_bisector(" << e1.get_id() << ',' << e2.get_id() << ")\n";
        typename std::vector<Root>::iterator it;
        std::cerr << "d1-sols: ";
        for (it = sol1.begin(); it != sol1.end(); it++) {
            std::cerr << CGAL::to_double(*it);
            std::cerr << ' ';
        }
        std::cerr << std::endl;
        std::cerr << "d2-sols: ";
        for (it = sol2.begin(); it != sol2.end(); it++) {
            std::cerr << CGAL::to_double(*it);
            std::cerr << ' ';
        }
        std::cerr << std::endl;
#endif

//        Comparison_result r =
//                compare_smallest_nonnegative_roots(AK(), dpoly1, dpoly2);
        if (positive_distance) return r; else return -r;
    }

};

} //namespace CGAL
#endif
