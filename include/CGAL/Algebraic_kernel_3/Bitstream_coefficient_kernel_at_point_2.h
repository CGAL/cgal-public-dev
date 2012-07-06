// Copyright (c) 2006-2012 Max-Planck-Institute Saarbruecken (Germany), 
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//                 Michael Kerber <mkerber@mpi-inf.mpg.de>

#ifndef CGAL_ALGEBRAIC_KERNEL_3_BITSTREAM_COEFFICIENT_KERNEL_AT_POINT_2_H
#define CGAL_ALGEBRAIC_KERNEL_3_BITSTREAM_COEFFICIENT_KERNEL_AT_POINT_2_H 1

/*!\file include/CGAL/Algebraic_kernel_3/Bitstream_coefficient_kernel_at_point_2.h
 * \brief Kernel class for CGAL::internal::Bitstream_descartes
 */

#include <CGAL/config.h>

#include <boost/optional.hpp>

#include <CGAL/Algebraic_kernel_d/Real_embeddable_extension.h>

#include <CGAL/Arrangement_2l/macros.h>

namespace CGAL {

// predeclaration
template < class SurfaceZAtXyIsolatorTraits >
class Bitstream_coefficient_kernel_at_point_2;

namespace internal { 

/*!\brief
 * Representation class for traits of bitstream tree
 */
template < class SurfaceZAtXyIsolatorTraits >
class Bitstream_coefficient_kernel_at_point_2_rep {
public:

    // types
    //! this instance's first template parameter
    typedef SurfaceZAtXyIsolatorTraits Surface_z_at_xy_isolator_traits;

    CGAL_SURFACE_Z_AT_XY_ISOLATOR_TRAITS_SNAP_TYPEDEFS(
            Surface_z_at_xy_isolator_traits
    );

#if DOXYGEN_RUNNING
    //! type of point
    typedef typename Surface_z_at_xy_isolator_traits::Point_2 Point_2;

    //! type of bivariate polynomial
    typedef typename P_curve_2::Polynomial_2 Polynomial_2;
#endif

    //! type of x-coordinate
    typedef typename Curve_kernel_2::Coordinate_1 Coordinate_1;
    
    //! type of rational numbers
    typedef typename Coordinate_1::Rational Rational;

    //! bigfloat interval
    typedef typename CGAL::Get_arithmetic_kernel< Rational >::
      Arithmetic_kernel::Bigfloat_interval Bigfloat_interval;

    //! type of curve analysis
    typedef typename Curve_analysis_2::Status_line_1 Status_line_1;

    //! type for intervals
    typedef boost::numeric::interval< Rational > Interval;

    //! type of coefficient for Bitstream_descartes_rndl_tree
    typedef Polynomial_2 Coefficient;
    
    //! type of integer for Bitstream_descartes_rndl_tree
    typedef typename CGAL::Fraction_traits< Rational >::Numerator_type Integer;
    
    //!\name Constructors
    //!@{
    
 private:
    /*!\brief
     * default constructor
     */
    Bitstream_coefficient_kernel_at_point_2_rep() {
    }
   
public:
    
    /*!\brief
     * standard constructor from a given refineable point
     */
    Bitstream_coefficient_kernel_at_point_2_rep(const Point_2& pt,
                                                bool use_artificial_x_interval) :
        _m_point(pt),
        _m_use_artificial_x_interval(use_artificial_x_interval),
        _m_prec_y(4) {
    };
    
    //!@}

private:
    // members

    //! the stored projected point
    mutable Point_2 _m_point;
    
    //! stores whether artificial x intervals is used
    mutable bool _m_use_artificial_x_interval;
    
    //! if x is rational this is its value
    mutable Rational _m_x_rational;
    
    //! current interval radius
    mutable Rational _m_eps;

    //! artificial x-interval
    mutable boost::optional< Interval > _m_x_interval;

    mutable int _m_prec_y;
    
    // friends
    friend class CGAL::Bitstream_coefficient_kernel_at_point_2< Surface_z_at_xy_isolator_traits >;
    
}; // Bitstream_coefficient_kernel_at_point_2_rep

} // namespace internal

/** *******************************************************************/

/*!\brief
 * Model of BitstreamCoefficientKernel concept to isolate the
 * real roots of a square-free polynomial defined by an algebraic
 * surface in 3d and a line parallel to the z-axis running through a
 * given point, whose coordinates are usally algebraic.
 */
template < class SurfaceZAtXyIsolatorTraits >
class Bitstream_coefficient_kernel_at_point_2 : public 
::CGAL::Handle_with_policy< CGAL::internal::Bitstream_coefficient_kernel_at_point_2_rep< 
SurfaceZAtXyIsolatorTraits > > {
    
public:
    
    // types
    //! this instance's first template parameter
    typedef SurfaceZAtXyIsolatorTraits Surface_z_at_xy_isolator_traits;

    CGAL_SURFACE_Z_AT_XY_ISOLATOR_TRAITS_SNAP_TYPEDEFS(
            Surface_z_at_xy_isolator_traits
    );

#if DOXYGEN_RUNNING
    //! type of point
    typedef typename Surface_z_at_xy_isolator_traits::Point_2 Point_2;
#endif
    
    //! type of instantiated class
    typedef Bitstream_coefficient_kernel_at_point_2< Surface_z_at_xy_isolator_traits >  Self;

    //! type of representation
    typedef internal::Bitstream_coefficient_kernel_at_point_2_rep< Surface_z_at_xy_isolator_traits > Rep;

    //! type of Base
    typedef CGAL::Handle_with_policy< Rep > Base;
    
private:

    //! type of x-coordinate
    typedef typename Rep::Coordinate_1 Coordinate_1;
    
    //! type of status line
    typedef typename Rep::Status_line_1 Status_line_1;

public:
    
    //! type for rational numbers
    typedef typename Rep::Rational Rational;

    //! type for bifgloat intervals
    typedef typename Rep::Bigfloat_interval Bigfloat_interval;

    //! type for intervals
    typedef typename Rep::Interval Interval;

public:
    // BS-types
    //!\name Types for Bitstream Traits
    //!@{
    
    //! type of coefficient for Bitstream_descartes_rndl_tree
    typedef Polynomial_2 Coefficient;
    
    //! type of integer for Bitstream_descartes_rndl_tree
    typedef typename CGAL::Fraction_traits< Rational >::Numerator_type Integer;
    
    //! type of boundary for Bitstream_descartes_rndl_tree
    typedef Rational Bound;

    //! type for the box
    typedef std::pair<Interval, Interval> Box;
    
    //!@}

    //!\name Constructors
    //!@{
    
    /*!\ brief
     * Default constructor
     */
    // TODO check whether default constructor is useful
    Bitstream_coefficient_kernel_at_point_2() :
      Base(Rep()) {
    }

 public:

    /*!\brief
     * Standard constructor to define xy with the help of a point
     */
    Bitstream_coefficient_kernel_at_point_2(
            const Point_2& pt,
            bool use_artificial_x_interval = false) : 
        Base(Rep(pt, use_artificial_x_interval)) {
        CGAL_precondition(pt.arcno() >= 0);
        // checks if rational and if so, sets artifical interval
        _x_iv();
    }

    //!@}
    
    //!\name Access members
    //!@{
    // functions to collaborate with functors

    //! returns stored point
    inline 
    Point_2 point() const {
        return this->ptr()->_m_point;
    }
    
    /*!\brief
     * returns current approximation of p as interval
     *
     * Note that is_zero(p) is not public as it uses internally point_on_curve
     * that uses this class itself -> cyclic dependency
     */
    inline
    Interval approximation(const Coefficient& p) const {
        
        // initialize x_iv
        Interval x_iv(_x_iv());
        
        // initialize y_iv
        Interval y_iv(_y_iv());
        
        return _evaluate_iv_2(p, x_iv, y_iv);
    }

private:
    
    //! returns finite x-coordinate of stored point
    inline
    static
    Coordinate_1 _x(const Point_2& pt) {
        return pt.x();
    }

    //! returns arcno of stored point
    inline 
    static
    int _arcno(const Point_2& pt) {
        return pt.arcno();
    }

    /*!\brief
     * returns Event1_info of supporting curve of stored 
     * point at its x-coordinate
     */
    inline
    static 
    Status_line_1 _sl(const Point_2& pt) {
        return pt.curve().status_line_at_exact_x(_x(pt));
    }
    
    //!@}
    
    //!\name Intervals and their evaluations
    //!@{

public:

    //! returns current interval approximation for x
    inline
    Interval x_interval() const {
        return _x_iv();
    }

    //! returns current interval approximation for y
    inline
    Interval y_interval() const {
        return _y_iv();
    }
    
private:

    //! returns current interval for x
    inline
    Interval _x_iv(bool try_to_use_exact_rational = false) const {
        bool use_axiv = this->ptr()->_m_use_artificial_x_interval;
        if (use_axiv && this->ptr()->_m_x_interval) {
            if(try_to_use_exact_rational) {
                return Interval(this->ptr()->_m_x_rational,
                                this->ptr()->_m_x_rational);
            }
            return *this->ptr()->_m_x_interval;
        } else {
            if (use_axiv && _x(point()).is_rational()) {
                Interval y_iv = _y_iv();
                Rational mid = _x(point()).low();
                this->ptr()->_m_x_rational = mid;
                Rational eps = 
                    (y_iv.upper() - y_iv.lower()) / Rational(2);
                this->ptr()->_m_eps = eps;
                this->ptr()->_m_x_interval  =
                    Interval(mid - eps, mid + eps);
                if(try_to_use_exact_rational) {
                    return Interval(this->ptr()->_m_x_rational,
                                    this->ptr()->_m_x_rational);
                }
                return *this->ptr()->_m_x_interval;
            }
        }
        // else 
        return Interval(_x(point()).low(), _x(point()).high());
    }
    
    //! return current interval for y
    inline
    Interval _y_iv() const {
        
        typename Arrangement_traits_2::Curve_kernel_2::Approximate_relative_y_2
            approx_y = 
            Arrangement_traits_2::instance().kernel().
            approximate_relative_y_2_object();

	std::pair<Bound,Bound> bound_pair = approx_y(point().xy(),this->ptr()->_m_prec_y);
	
        return Interval(bound_pair.first,bound_pair.second);
	
    }
    
    //! interval evaluation of univariate polynomial 
    inline
    static
    Interval _evaluate_iv_1(const typename Coefficient::NT& p, 
                            const Interval& x_iv) {
        // TASK cache for iv(p)
        int n = p.degree();
        Interval ret(p[n],p[n]);
        for(int i = n - 1; i >= 0; i--) {
            ret *= x_iv;;
            ret += Interval(p[i],p[i]);
        }
        return ret;
    }
    
    //! interval evaluation of bivariate polynomial 
    inline
    static 
    Interval _evaluate_iv_2(const Coefficient& p, 
                            const Interval& x_iv, 
                            const Interval& y_iv) {
        // TASK cache for iv(p)
        int i = p.degree();
        Interval ret = _evaluate_iv_1(p[i--],x_iv);
        while (i >= 0) { 
            ret *= y_iv; 
            ret += _evaluate_iv_1(p[i--],x_iv); 
        }
        return ret;
    }
    
    //!@}
    
    //!\name Representation
    //!@{

public:    
    //! refines stored representation
    inline
    void refine() const {
        
        // Goal: keep xy-approx as quadratic as possible
        
        bool do_refine_y = true;

        if (this->ptr()->_m_use_artificial_x_interval || !is_x_rational()) {
            
            Interval x_iv = _x_iv();
            Interval y_iv = _y_iv();
            
            Rational x_len = x_iv.upper() - x_iv.lower();
            Rational y_len = y_iv.upper() - y_iv.lower();

//             if(x_len <= 0) {
//                 std::cerr << "x_len: " << x_len << "; " << CGAL::to_double(x_len) << "\n";
//             }
//             if(y_len <= 0) {
//                 std::cerr << "y_len: " << y_len << "; " << CGAL::to_double(y_len) << "\n";
//             }
            
            CGAL_assertion(x_len >= 0);
            CGAL_assertion(y_len >= 0);

            if (x_len > y_len) {
                do_refine_y = false;
            }
        }
        
        // Remark: it is not good to make the x- and y- length similar
        //         long, as this requires an initial lengthy refinement of
        //         y since x is usally very small.
        //         The point is different: One refinement step of x
        //         implies a large number of refinement steps for y
        //         until |x| ~ |y|.
        
        if (do_refine_y) {
            refine_y();
        } else {
            // refine x
            refine_x();
        }
    }

    //! refines x of stored representation
    inline
    void refine_x() const {
        if (this->ptr()->_m_use_artificial_x_interval && 
            this->ptr()->_m_x_interval) {
            this->ptr()->_m_eps /= Rational(2);
            const Rational& mid = this->ptr()->_m_x_rational;
            const Rational& eps = this->ptr()->_m_eps;
            this->ptr()->_m_x_interval = Interval(mid - eps, mid + eps);
        } else {
            _x(point()).refine();
        }
    }

    //! refines y of stored representation
    inline
    void refine_y() const {
      this->ptr()->_m_prec_y*=2;
    }
    
    //! refines approximation until \c pt does not lie in approximation
    //! \pre pt != stored point
    inline
    void refine_against(const Self& bck) {
        CGAL_precondition(bck.point() != this->point());
        
        Interval this_x_iv = x_interval();
        Interval this_y_iv = y_interval();
        
        Interval that_x_iv = bck.x_interval();
        Interval that_y_iv = bck.y_interval();
        
        while (boost::numeric::overlap(this_x_iv, that_x_iv) && 
               boost::numeric::overlap(this_y_iv, that_y_iv)) {
            this->refine();
            bck.refine();
        }
    }

    //! returns approximation as a box around the point
    Box approximation_box(long prec) {
        
        typename CGAL::Fraction_traits<Rational>::Compose compose;

        Rational bound = (prec < 0) 
            ? compose(CGAL::ipower(Integer(2),-prec),1)
            : compose(1,CGAL::ipower(Integer(2),prec));

        Rational bound_times_2 = (prec < 0) 
            ? compose(CGAL::ipower(Integer(2),-(prec+1)),1)
            : compose(1,CGAL::ipower(Integer(2),prec+1));

        

        // Refine until the approximation is inside the box
        while(y_interval().upper()-y_interval().lower() > bound_times_2 ||
              x_interval().upper()-x_interval().lower() > bound_times_2) {
            refine();
        }
        
        // now, scale to the correct size
        Interval x_iv = _rescale_to(x_interval(),prec);
        CGAL_assertion(x_iv.upper()-x_iv.lower() == bound);
        Interval y_iv = _rescale_to(y_interval(),prec);
        CGAL_assertion(y_iv.upper()-y_iv.lower() == bound);

        return std::make_pair(x_iv,y_iv);

    }

private:

    Interval _rescale_to(Interval iv, long prec) {
        
        // TODO Implement for non-integers

        typedef CGAL::Fraction_traits<Rational> Fraction_traits;
        typename Fraction_traits::Decompose decompose;
        Integer num, denom;
        Integer pow = CGAL::ipower(Integer(2),CGAL::abs(prec)+1);
        Integer low_result, high_result;
        Integer remainder;

        decompose(iv.lower(), num, denom);

        //std::cout << "Lower: " << num << ", " << denom << std::endl;

        if(prec<0) {
            denom = denom * pow;
        } else {
            num = num * pow;
        }
        CGAL::div_mod(num, denom, low_result,remainder);

        if(remainder<0) {
            low_result = low_result - 1;
        }

        decompose(iv.upper(), num, denom);

        //std::cout << "Upper: " << num << ", " << denom << std::endl;

        if(prec<0) {
            denom = denom * pow;
        } else {
            num = num * pow;
        }
        
        CGAL::div_mod(num, denom,high_result,remainder);
        if(remainder>0) {
            high_result = high_result + 1;
        }

        CGAL_assertion(high_result - low_result <=2);
        if(high_result==low_result) {
            
            low_result = low_result - 1;
            high_result = high_result + 1;
            

        } else if(high_result - low_result == 1) {
            if(CGAL::mod(low_result,Integer(2)) == 0) {
                low_result = CGAL::div(low_result,Integer(2));
                high_result = CGAL::div(high_result+1,Integer(2));
            } else {
                low_result = CGAL::div(low_result-1,Integer(2));
                high_result = CGAL::div(high_result,Integer(2));
            }
            pow = pow / 2;
        }

        typename Fraction_traits::Compose compose;

        if(prec < 0) {
            return Interval(compose(pow*low_result,Integer(1)),
                            compose(pow*high_result,1));
        }
        
        return Interval(compose(low_result,pow),
                        compose(high_result,pow));
        
    }
    

    // TODO add refine_y

public:

    //!\name Rational neighbors
    //!@{
    
    //! returns \c true iff x is rational
    bool is_x_rational() const {
        return _x(point()).is_rational();
    }
    
    //! return the rational x, if x is rational
    Rational rational_x() const {
        CGAL_precondition(is_x_rational());
        return _x(point()).low();
    }

    // TODO make use of Curve_kernel_2::Bound_between

    //! returns a rational to the left of x within the current x-range approx
    Rational left_x() const {
        CGAL_precondition(!is_x_rational() || 
                          this->ptr()->_m_use_artificial_x_interval);
        if (is_x_rational()) {
            CGAL_assertion(this->ptr()->_m_use_artificial_x_interval);
            this->_x_iv();
            return 
                (this->ptr()->_m_x_interval->lower() + 
                 this->ptr()->_m_x_rational) / Rational(2);
        }
        return _x(point()).rational_between(Coordinate_1(_x(point()).low()));
    }

    //! returns a rational to the right of x within the current x-range approx
    Rational right_x() const {
        CGAL_precondition(!is_x_rational() || 
                          this->ptr()->_m_use_artificial_x_interval);
        if (is_x_rational()) {
            CGAL_assertion(this->ptr()->_m_use_artificial_x_interval);
            this->_x_iv();
            return 
                (this->ptr()->_m_x_rational + 
                 this->ptr()->_m_x_interval->upper()) / Rational(2);
        }
        return _x(point()).rational_between(
                Coordinate_1(_x(point()).high())
        );
    }

    //! returns whether a given \c y is in current y-approximation
    bool is_in_x_interval_interior(Rational x) const {
        Interval x_iv = _x_iv();
        return (x_iv.lower() < x) && (x < x_iv.upper());
    }
    
    
    //! returns mean value of y-approximation
    Rational mean_y() const {
        Interval y_iv = _y_iv();
        return (y_iv.upper() - y_iv.lower()) / Rational(2);
    }
    
    //! returns whether a given \c y is in current y-approximation
    bool is_in_y_interval_interior(Rational y) const {
        Interval y_iv = _y_iv();
        return (y_iv.lower() < y) && (y < y_iv.upper());
    }
    
    //!@}


public:

    // TASK document functors
    /** *********************************************************************/

    struct Is_zero : public std::unary_function<Coefficient,bool> {

        Is_zero(const Self& bck) :
            _m_bck(bck) {
        }
        
        inline
        bool operator() (Coefficient f) const {
          typename Surface_z_at_xy_isolator_traits::Point_on_curve_2 
            point_on_curve
            // TODO (point_on_curve_object())
            ;
          bool on_curve = point_on_curve(_m_bck.point(), f);
          
          return on_curve;
        }

    private:
        // members
        mutable Self _m_bck;

    };

    Is_zero is_zero_object() const {
        return Is_zero(*this);
    }

// TASK convert_to_bfi should only use computations over bigfloats
// to avoid double conversion from rational to bigfloat and again
// to rational
    struct Convert_to_bfi : public std::unary_function< Coefficient, Bigfloat_interval > {

        Convert_to_bfi(const Self& bck) :
            _m_bck(bck) {
        }

        Bigfloat_interval operator() (Coefficient f) const {

          // TASK cleanup implementation

            long p = CGAL::get_precision(Bigfloat_interval());
//             long prec = 16, wbit = 0;

            //! can \c p be less than 0 ?? nope..
            Rational bound = (Integer)1 << std::abs(p);
            if(p > 0) {
                bound = Rational(1) / bound;
            }    
            
//             typedef typename CGAL::Polynomial_traits_d<Coefficient>
//                 ::template Rebind<Bigfloat_interval ,2>::Other P_traits_bfi;
//             typedef typename P_traits_bfi::Type Poly_bfi_2;
// 
//             typename P_traits_bfi::Substitute subs;
//             Bigfloat_interval pt_bfi;

            Interval f_eval_iv;
            while(true) {
//                 CGAL::set_precision(Bigfloat_interval(), prec);

                Interval x_iv(this->_m_bck._x_iv());
                Interval y_iv(this->_m_bck._y_iv());

#if 0
                typename CGAL::Coercion_traits<
                  Coefficient /* = Polynomial_2 */, Poly_bfi_2 >::Cast cvt;

                Poly_bfi_2 f_bfi = cvt(f);

                Bigfloat_interval s[] =
                    { CGAL::hull(CGAL::convert_to_bfi(x_iv.lower()),
                                CGAL::convert_to_bfi(x_iv.upper())),
                     CGAL::hull(CGAL::convert_to_bfi(y_iv.lower()),
                                CGAL::convert_to_bfi(y_iv.upper()))};

                pt_bfi = subs(f_bfi, s, s + 2);
//                  std::cerr << "prec: " << prec <<
//                     "; x-width: " << CGAL::convert_to_bfi(
//                         x_bounds.second - x_bounds.first) <<
//                      "; y-width: " << CGAL::convert_to_bfi(y_iv.upper() - y_iv.lower()) << "; pt_bfi: " << pt_bfi << "\n";
#else
                f_eval_iv = this->_m_bck._evaluate_iv_2(f, x_iv, y_iv);
#endif

//                if (CGAL::singleton(pt_bfi)) {
//                  break;
//                }
// 
   // TASK precision might be too high if already rational!!!
                if(f_eval_iv.upper() - f_eval_iv.lower() < bound) {
                    break;
                }
/*               long ceil = CGAL::internal::ceil_log2_abs(pt_bfi);
               long signi = CGAL::get_significant_bits(pt_bfi);
               wbit = ceil - signi + p;
            
               if (wbit < -5) {
                 break;
               }
                prec *= 2;*/
                this->_m_bck.refine();
            }
            return CGAL::hull(CGAL::convert_to_bfi(f_eval_iv.lower()),
                                CGAL::convert_to_bfi(f_eval_iv.upper()));
        }

    private:
        mutable Self _m_bck;
    };

    Convert_to_bfi convert_to_bfi_object() const {
      return Convert_to_bfi(*this);
    }

}; // class Bitstream_coefficient_kernel_at_point_2

} //namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_3_BITSTREAM_COEFFICIENT_KERNEL_AT_POINT_2_H
// EOF
