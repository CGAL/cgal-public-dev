// ============================================================================
//
// Copyright (c) 2001-2007 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : QdX
//
// File          : include/QdX/Surface_intersections_with_moebius_strip_traits_2.h
// QdX_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//
// ============================================================================

/*! \file include/QdX/Surface_intersections_with moebius_strip_traits_2
 * \brief Provides traits class to compute arrangement on a moebius strip
 * induced by other intersection curves with other surfaces
 */

#ifndef QdX_SURFACE_INTERSECTIONS_WITH_MOEBIUS_TRAITS_2
#define QdX_SURFACE_INTERSECTIONS_WITH_MOEBIUS_TRAITS_2

#include <QdX/basic.h>

// TODO remove GAPS_2 with CKvA_2
#include <SoX/GAPS/GAPS_2.h>
#include <SoX/GAPS/CGAL_Arrangement_2_for_GAPS_traits.h>

#include <AcX/Algebraic_curve_2.h>
#include <AcX/Algebraic_curve_pair_2.h>

#include <QdX/SfX/Algebraic_surface_3.h>

#include <CGAL/Arrangement_2.h>
#include <CGAL/Sweep_line_2/Sweep_line_2_utils.h>

// TODO replace with usual gaps #include <QdX/P_kernel_3.h>

namespace QdX {

/*!\brief
 * Derived ArrangmentTraits class for points and segments embedded on 
 * a given moebius strip
 */
template < class ArithmeticTraits >
class Surface_intersections_with_moebius_strip_traits_2 : 
        public SoX::CGAL_Arrangement_2_for_GAPS_traits< 
SoX::GAPS_2<  
AcX::Algebraic_curve_pair_2< 
AcX::Algebraic_curve_2< ArithmeticTraits > 
>
>,  
SoX::GAPS_2<  
AcX::Algebraic_curve_pair_2< 
AcX::Algebraic_curve_2< ArithmeticTraits > 
> 
>
> {
    
public:
    
    //! this instance's template parameter
    typedef ArithmeticTraits Arithmetic_traits;
    
    //! the class itself
    typedef 
    Surface_intersections_with_moebius_strip_traits_2< Arithmetic_traits > 
    Self;
    
private:
    typedef Arithmetic_traits AT;
    
    typedef AcX::Algebraic_curve_2< AT > Algebraic_curve_2;
    typedef AcX::Algebraic_curve_pair_2< Algebraic_curve_2 > 
    Algebraic_curve_pair_2;

public:
    
    //! type of kernel
    typedef SoX::GAPS_2< Algebraic_curve_pair_2 > Kernel;

    //! type of Base
    typedef SoX::CGAL_Arrangement_2_for_GAPS_traits< Kernel, Kernel > Base;
    
    // TODO class for moebius strip
    //! type of moebius strip
    typedef QdX::Algebraic_surface_3< AT > Moebius_strip_3;

    //! type of Surface
    typedef QdX::Algebraic_surface_3< AT > Surface_3;

    //! typedef of Curve_2
    typedef Surface_3 Curve_2;
    
    //! type of x-monotone curve
    typedef typename Base::X_monotone_curve_2 X_monotone_curve_2;
    
    //! type of point
    typedef typename Base::Point_2 Point_2;
  

public:
    //!\name Constructors
    //!@{

    // default constructor should be not supported!
    /*\brief 
     * Default constructor 
     */
    Surface_intersections_with_moebius_strip_traits_2() {
    };
    
    /*\brief 
     * Standard constructor that stores \c base as base moebius strip
     */
    Surface_intersections_with_moebius_strip_traits_2(
            const Moebius_strip_3& base) :
        _m_base(base) {
    };
    
    //@}

    //!\name Accessors
    //!@{
    
    /*!\brief
     * returns the stored base moebius strip
     */
    Moebius_strip_3 base() const {
        return this->_m_base;
    }
    
    //!@}

#if 0 // can be removed?!
    
private:
    /*! \brief Compares the x-coordinates of two events
     */
    template < class BaseTraits >
    class Moebius_strip_compare_x_2 {
    public:
        typedef BaseTraits Base_traits;
        
        Moebius_strip_compare_x_2(const Base_traits& base_traits) :
            _m_base_traits(base_traits) {
        }
        
        /*!
         * Compare the x-coordinates of two points.
         * \param p1 The first point.
         * \param p2 The second point.
         * \return LARGER if x(p1) > x(p2);
         *         SMALLER if x(p1) < x(p2);
         *         EQUAL if x(p1) = x(p2).
         */
        CGAL::Comparison_result operator() (
                const Point_2& p1, 
                const Point_2& p2) const {
            
            if (p1.is_at_infinity() && p2.is_at_infinity()) {
                // used for comparison in top traits
                return _compare_on_rim(p1,p2);
            }
            return p1.compare_x(p2);
        }
        
        /*!
         * Compare the relative positions of a vertical curve and another given
         * curves at y = +/- oo.
         * \param p A reference point; we refer to a vertical line 
         * incident to p.
         * \param cv The compared curve.
         * \param end MIN_END if we refer to cv's minimal end,
         *            MIN_END if we refer to its maximal end.
         * \pre cv's relevant end is defined at y = +/- oo.
         * \return SMALLER if p lies to the left of cv;
         *         LARGER if p lies to the right cv;
         *         EQUAL in case of an overlap.
         */
        CGAL::Comparison_result operator() (
                const Point_2& p,
                const X_monotone_curve_2& cv, 
                CGAL::Curve_end end) const {
            
            Point_2 p2 = 
                (end == CGAL::MIN_END ?
                 cv.min_endpoint() : cv.max_endpoint());
            
            // used in event-compare of sweep
            return _compare_on_rim(p, p2);
        }
        
        /*!
         * Compare the relative positions of two curves at y = +/- oo.
         * \param cv1 The first curve.
         * \param end1 MIN_END if we refer to cv1's minimal end,
         *             MIN_END if we refer to its maximal end.
         * \param cv2 The second curve.
         * \param end2 MIN_END if we refer to cv2's minimal end,
         *             MIN_END if we refer to its maximal end.
         * \pre The curves are defined at y = +/- oo.
         * \return SMALLER if cv1 lies to the left of cv2;
         *         LARGER if cv1 lies to the right cv2;
         *         EQUAL in case of an overlap.
         */
        CGAL::Comparison_result
        operator() (
                const X_monotone_curve_2& cv1, 
                CGAL::Curve_end end1,
                const X_monotone_curve_2& cv2, 
                CGAL::Curve_end end2)
            const {
            
            Point_2 p1 = 
                (end1 == CGAL::MIN_END ?
                 cv1.min_endpoint() : cv1.max_endpoint());
            
            Point_2 p2 = 
                (end2 == CGAL::MIN_END ?
                 cv2.min_endpoint() : cv2.max_endpoint());
            
            // used in event-compare of sweep
            return _compare_on_rim(p1, p2);
        }
        
    private:
        
        //! compares the actual x-coordinates
        CGAL::Comparison_result _compare_on_rim(
                const Point_2& p1,
                const Point_2& p2) const {
            return p1.x().number().finite().compare(p2.x().number().finite());
        }
        
        //! the store base_traits instance
        const Base_traits& _m_base_traits;
    };
    
public:
    //! type Compare_x_2 functor
    typedef Moebius_strip_compare_x_2< Self > Compare_x_2;
    /*! \brief Get a Compare_x_2 functor object. */
    Compare_x_2 compare_x_2_object () const {
        return Compare_x_2(*_m_base_traits);
    }
    
#endif
    

private:
    //! Compares two points lexicographically
    template < class BaseTraits >
    class Moebius_strip_compare_xy_2 {
    public:
        typedef BaseTraits Base_traits;
        
        //! standard constructor
        Moebius_strip_compare_xy_2(const Base_traits& base_traits) :
            _m_base_traits(base_traits) {
        }
        
        //! this instance itself
        typedef Moebius_strip_compare_xy_2< Base_traits > Self;
        
        //! function call operator
        CGAL::Comparison_result operator()(
                Point_2 p1,
                Point_2 p2) const {
            
            typedef typename AT::Integer Integer;
            typedef typename AT::Rational Rational;
            typedef typename AT::Poly_int1 Polynomial_1;
            typedef typename Algebraic_curve_2::Event1_info Event1_info;

            if (p1.is_at_infinity() && !p2.is_at_infinity() ||
                !p1.is_at_infinity() && p2.is_at_infinity()) {
                
                Moebius_strip_compare_x_2< Base_traits > 
                    compare_x(_m_base_traits);
                
                CGAL::Comparison_result res =
                    compare_x(p1, p2);
                
                if (res != CGAL::EQUAL) {
                    return res;
                }
                // else

                // X_coordinate is Algebraic real
                typedef typename Algebraic_curve_2::X_coordinate X_coordinate;
                X_coordinate y1, y2;
                
                if (p1.x().number().infty() == NiX::PLUS_INFTY) {
                    y1 = p1.curve().
                        horizontal_asymptote_for_arc_to_plus_infinity(
                                p1.arcno()
                        ).finite();
                } else {
                    CGAL_assertion(p1.x().number().is_finite());
                    Polynomial_1 poly1 = 
                        CGAL::CGALi::make_square_free(
                                NiX::substitute_x(p1.curve().f(),
                                                  Integer(0))
                        );
                    poly1.scale(Integer(-1),Integer(1));
                    Event1_info ei = 
                        p1.curve().event_info_at_x(X_coordinate(0));
                    y1 = X_coordinate(poly1, 
                                      -ei.upper_boundary(p1.arcno()),
                                      -ei.lower_boundary(p1.arcno()));
                }
                
                if (p2.x().number().infty() == NiX::PLUS_INFTY) {
                    y2 = p2.curve().
                        horizontal_asymptote_for_arc_to_plus_infinity(
                                p2.arcno()
                        ).finite();
                } else {
                    CGAL_assertion(p2.x().number().is_finite());
                    Polynomial_1 poly2 = 
                        CGAL::CGALi::make_square_free(
                                NiX::substitute_x(p2.curve().f(),
                                                  Integer(0))
                        );
                    poly2.scale(Integer(-1),Integer(1));
                    Event1_info ei = 
                        p2.curve().event_info_at_x(X_coordinate(0));
                    y2 = X_coordinate(poly2, 
                                      -ei.upper_boundary(p2.arcno()),
                                      -ei.lower_boundary(p2.arcno()));
                }
                
                // else
                // reverse here, as we computed -root for 0-side
                // used to compare on identification
                return NiX::compare(y2, y1);
            }

            // else 
            if (p1.is_at_infinity()) {
                CGAL_assertion(p2.is_at_infinity());
                return -p1.compare_xy(p2);
            }
            
            // else
            // geometric comparison in interior of parameter space
            return p1.compare_xy(p2);
        }
    private:
        const Base_traits& _m_base_traits;
        
    };

public:
    //! type Compare_xy_2 functor
    typedef Moebius_strip_compare_xy_2< Self > Compare_xy_2;
    /*! \brief Get a Compare_xy_2 functor object. */
    Compare_xy_2 compare_xy_2_object () const {
        return Compare_xy_2(*_m_base_traits);
    }
    

private:
    //! splits intersection curve into x-monotone curves or isolated points
    template < class BaseTraits >
    class Moebius_strip_make_x_monotone_2 {
    public:
        typedef BaseTraits Base_traits;

        //! standard constructor
        Moebius_strip_make_x_monotone_2(const Base_traits& base_traits) :
            _m_base_traits(base_traits) {
        }
        
        template<class OutputIterator>
        OutputIterator operator() (
                const Curve_2& cv, 
                OutputIterator oi) 
        {
            typedef typename Point_2::Compactified_x Cp_x;
            typedef typename Algebraic_curve_2::X_coordinate X_coordinate;
            typedef typename AT::Poly_int1 Polynomial_1;
            typedef typename AT::Poly_int2 Polynomial_2;
            typedef typename AT::Integer Integer;
            Polynomial_2 p;
            
            // TODO parametrization of cv into UV-space

            p = CGAL::CGALi::canonicalize_polynomial(p);

            // TODO what to do if not square-free?
            CGAL_assertion(CGAL::CGALi::is_square_free(p));
#if !NDEBUG
            CGAL::set_ascii_mode(std::cout);
            std::cout << "P(u,v)= " << p << std::endl;
            CGAL::set_pretty_mode(std::cout);
            std::cout << "P(u,v)= " << p << std::endl;
#endif
            
            // curve cache
            Algebraic_curve_2 curr_curve = 
                Algebraic_curve_2::get_curve_cache()(p);
            
            // create segments
            std::list< X_monotone_curve_2 > tmp;
            SoX::curve_to_segments< X_monotone_curve_2 >(
                    curr_curve,
                    std::back_inserter(tmp)
            );
            std::list< X_monotone_curve_2 > xcurves;
            std::list< Point_2 > xpoints;
            
            for (typename std::list< X_monotone_curve_2 >::iterator 
                     it = tmp.begin(); it != tmp.end(); ++it) {
                if (it->is_degenerate()) {
                    xpoints.push_back(it->source());
                } else {
                    xcurves.push_back(*it);
                }
            }
            
            // create intersection points with line y=+-1
            Polynomial_2 hlines(
                    Polynomial_1(Integer(-1)), //y^0
                    Polynomial_1(Integer(0)),  //y^1
                    Polynomial_1(Integer(1))   //y^2
            );
            
            // y=+-1 is not included in p 
            // TODO allow this input in the future
            CGAL_assertion(NiX::gcd(hlines, p).degree() == 0);
            
            Algebraic_curve_pair_2 pair(p, hlines);
            
            for (int i = 0; i < pair.num_events(); i++) {
                int fg = pair.event_indices(i).fg;
                if (fg >= 0) {
                    typename Algebraic_curve_pair_2::Event2_slice 
                        slice = pair.slice_at_event(i);
                    int arc;
                    arc = slice.arc_of_other_curve(hlines, 0);
                    if (arc != -1) {
                        Point_2 pt(pair.event_x(i), hlines, 0);
                        xpoints.push_back(pt);
                    }
                    arc = slice.arc_of_other_curve(hlines, 1);
                    if (arc != -1) {
                        Point_2 pt(pair.event_x(i), hlines, 1);
                        xpoints.push_back(pt);
                    }
                }
            }
            
            // compute arrangement
            typedef CGAL::Arrangement_2< Base > Base_arr;
            Base_arr arr;
            CGAL::insert_empty(arr, 
                               xcurves.begin(), xcurves.end(),
                               xpoints.begin(), xpoints.end());
            
            Point_2 hsrc0(X_coordinate(0),hlines,0);
            Point_2 htgt0(Cp_x(NiX::PLUS_INFTY), hlines,0);
            X_monotone_curve_2 hline0(hsrc0, htgt0, hlines,0,0,0);
            Point_2 hsrc1(X_coordinate(0),hlines,1);
            Point_2 htgt1(Cp_x(NiX::PLUS_INFTY), hlines,1);
            X_monotone_curve_2 hline1(hsrc1, htgt1, hlines,1,1,1);

            // TODO add hline0 and hline1 to arrangement instead of fictitious
            
            // and select points in stripe [-1,1]x[0,oo[
            for (typename Base_arr::Vertex_const_iterator 
                     vit = arr.vertices_begin();
                 vit != arr.vertices_end(); ++vit) {
                if (!vit->is_at_infinity() && vit->is_isolated()) {
                    Point_2 pt = vit->point();
                    if (hline0.compare_y_at_x(pt) != CGAL::SMALLER &&
                        hline1.compare_y_at_x(pt) != CGAL::LARGER) {
                        // TODO assign boundary conditions
                        *oi++ = CGAL::make_object(pt);
                    }
                }
            }
            
            // and select curves in stripe [-1,1]x[0,oo[
            for (typename Base_arr::Edge_const_iterator 
                     eit = arr.edges_begin();
                 eit != arr.edges_end(); ++eit) {
                if (hline0.compare_y_at_x(
                            eit->curve().min_endpoint()
                    ) != CGAL::SMALLER &&
                    hline1.compare_y_at_x(
                            eit->curve().min_endpoint()
                    ) != CGAL::LARGER &&
                    hline0.compare_y_at_x(
                            eit->curve().max_endpoint()
                    ) != CGAL::SMALLER &&
                    hline1.compare_y_at_x(
                            eit->curve().max_endpoint()
                    ) != CGAL::LARGER) {
                    if (eit->curve().min_endpoint().x().
                        number().finite().compare(X_coordinate(0)) ==
                        CGAL::EQUAL) {
                        eit->curve().set_boundary_in_x(
                                CGAL::MIN_END,
                                CGAL::AFTER_DISCONTINUITY
                        );
                    } else if (eit->curve().min_endpoint().x().
                               number().infty() == NiX::PLUS_INFTY) {
                        eit->curve().set_boundary_in_x(
                                CGAL::MIN_END,
                                CGAL::BEFORE_DISCONTINUITY
                        );
                    } else if (hline0.compare_y_at_x(
                                       eit->curve().min_endpoint()
                               ) == CGAL::EQUAL) {
                        eit->curve().set_boundary_in_y(
                                CGAL::MIN_END,
                                // TODO introduce FINITE boundary
                                CGAL::AFTER_DISCONTINUITY
                        );
                    } else if (hline1.compare_y_at_x(
                                       eit->curve().min_endpoint()
                               ) == CGAL::EQUAL) {
                        eit->curve().set_boundary_in_y(
                                CGAL::MIN_END,
                                // TODO introduce FINITE boundary
                                CGAL::BEFORE_DISCONTINUITY
                        );
                    }
                    if (eit->curve().max_endpoint().x().
                        number().finite().compare(X_coordinate(0)) ==
                        CGAL::EQUAL) {
                        eit->curve().set_boundary_in_x(
                                CGAL::MAX_END,
                                CGAL::AFTER_DISCONTINUITY
                        );
                    } else if (eit->curve().max_endpoint().x().
                               number().infty() == NiX::PLUS_INFTY) {
                        eit->curve().set_boundary_in_x(
                                CGAL::MAX_END,
                                CGAL::BEFORE_DISCONTINUITY
                        );
                    } else if (hline0.compare_y_at_x(
                                       eit->curve().max_endpoint()
                               ) == CGAL::EQUAL) {
                        eit->curve().set_boundary_in_y(
                                CGAL::MAX_END,
                                // TODO introduce FINITE boundary
                                CGAL::AFTER_DISCONTINUITY
                        );
                    } else if (hline1.compare_y_at_x(
                                       eit->curve().max_endpoint()
                               ) == CGAL::EQUAL) {
                        eit->curve().set_boundary_in_y(
                                CGAL::MAX_END,
                                // TODO introduce FINITE boundary
                                CGAL::BEFORE_DISCONTINUITY
                        );
                    }

                    *oi++ = CGAL::make_object(eit->curve());
                }
            }
            
            return oi;
        }
        
    private:
        // the stored base_traits instance
        const Base_traits& _m_base_traits;
        
    };

public:
    //! type of Make_x_monotone_2 functor
    typedef Moebius_strip_make_x_monotone_2< Self > Make_x_monotone_2;
    /*! \brief Get a Make_x_monotone_2 functor object. */
    Make_x_monotone_2 make_x_monotone_2_object () {
        return Make_x_monotone_2(*_m_base_traits);
    }

private:
    //! the stored base moebius_strip
    Surface_3 _m_base;
    
    //! instance of the base_traits class
    Self *_m_base_traits;
};

} // namespace QdX

#endif // QdX_SURFACE_INTERSECTIONS_WITH_TORUS_TRAITS_2
// EOF
