// Copyright (c) 2009, 2010, 2011, 2012 Max-Planck-Institut fuer Informatik (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>

/*!\file Range_analysis_2
 * \brief a collection of range analysis methods (1D and 2D) with the support
 * for complex arithmetic
 *  
 * 2D range methods:
 * eval_range_RT_2_aff - recursive Taylor bottom-up with affine muls
 * eval_range_RT_2 - recursive Taylor bottom-up with IA muls
 * eval_range_AF1_RD_2 - the first affine form with derivative information
 * eval_range_AF1_2/aff - the first affine form (a building block for
 *                                                   the advanced methods)
 * 1D range methods:
 * eval_range_AF1_RT_1_aff - recursive Taylor with affine muls
 * eval_range_AF1_RD_1 - the first affine form with derivative information
                            (specialized for complex case)
 * eval_range_AF1_1 - the first affine form (a building block for the
 *                                                  advanced methods)
 * 2D point methods:
 * eval_range_AF1_2 - evaluates a 2D polynomial with interval coefficients
 *                          at a point (a building block for advanced methods)
 */

#ifndef CGAL_ALGEBRAIC_KERNEL_2_RANGE_ANALYSIS_2
#define CGAL_ALGEBRAIC_KERNEL_2_RANGE_ANALYSIS_2

#include <CGAL/config.h>

#include <complex>
#include <boost/mpl/if.hpp>
#include <CGAL/Algebraic_kernel_2/Type_conversion_internal_traits.h>
#include <CGAL/Algebraic_kernel_2/Complex_sector.h>

// maximal default order of taylor expansion (bivariate case)
#define CGAL_IA_2_MIXED_DER_ORDER 2U
// maximal default order of derivatives f^(i) to consider (univariate case)
#define CGAL_IA_1_BINOM_SIZE 6U

// extern CGAL::Bbox_2 cplx_box[cplx_samples];
// 
// extern std::vector< std::complex< double > > cplx_pts[cplx_samples];

namespace CGAL {

namespace internal {

template < class NT >
struct Quad_ {
    Quad_() {  }
    Quad_(const NT& l_, const NT& r_, const NT& b_, const NT& t_) :
            l(l_), r(r_), b(b_), t(t_) { }
    NT l, r, b, t; // left right bottom top
};

template < class NT >
std::ostream& operator <<(std::ostream& os, const Quad_< NT >& q) {

    os << "[" << q.l << "; " << q.r << "] x [" << q.b << "; " << q.t << "]";
    return os;
}

} // internal


template < class PolyInt_2, class FloatCoeff, bool UseComplexArithmetic >
class Range_analysis_2
{
public: 
    //! \name public typedefs 
    //!@{ 
    //! this instance's template argument (concise explanation ??)
    typedef PolyInt_2 Poly_int_2; // integer polynomial 
    //! interval of doubles or BigFloat
    typedef FloatCoeff Float_coeff;
    //! univariate integer polynomial
    typedef typename Poly_int_2::NT Poly_int_1; 

    //! integer number type
    typedef typename CGAL::Polynomial_traits_d< Poly_int_2 >::
            Innermost_coefficient_type Integer;

    //! arithmetic kernel
    typedef typename CGAL::Get_arithmetic_kernel< Integer >::
         Arithmetic_kernel Arithmetic_kernel;
    //! rational NT
    typedef typename Arithmetic_kernel::Rational Rational;

    //! polynomials with floating-point coeffs for internal usage
    typedef CGAL::Polynomial< Float_coeff > Poly_1;
    typedef CGAL::Polynomial< Poly_1 > Poly_2;  

    //! transform traits
    typedef internal::Type_conversion_internal_traits< Float_coeff, Rational >
             TC_traits;

    //! floating-point NT for internal use
    typedef typename TC_traits::Float NT_;

    //! decide whether to use complex arithmetic internally
    typedef typename boost::mpl::if_c< UseComplexArithmetic, 
        std::complex< NT_ >, NT_ >::type NT;

    typedef internal::Complex_IA_traits< NT > IA_traits;
    //! real number type (not to confuse with NT that can be complex)
    typedef typename IA_traits::Real Real;

    //! represents a box !?
    typedef internal::Quad_< NT > Quad;
    //! sector in C
    typedef internal::Complex_sector< Real > Complex_sector;

    //! computes vector magnitude (as you might guess)
    typedef typename TC_traits::Vector_magnitude Vector_magnitude;
    //! approximates absolute value: looks like a versatile construction..
    //! no better solution found yet
    typedef typename IA_traits::template Approximate_abs< Vector_magnitude >
            Approximate_abs;
    //!@}
protected:
    //!\name protected typedefs
    //!@{

    //! to store binomial coefficients
    typedef std::vector< Real > Vector_real_1;
    typedef std::vector< Vector_real_1 > Vector_real_2;

    //! rational polynomials
    typedef CGAL::Polynomial< Rational > Poly_rat_1;
    typedef CGAL::Polynomial< Poly_rat_1 > Poly_rat_2;

    struct _Cache_entry_2 { //! caching bivariate data
        Poly_rat_2 poly_rat;
        Poly_2 poly_y;
        Poly_2 poly_x;
        std::vector< Poly_2 > mixed_dfs;
        unsigned mixed_dfs_depth;
    };

    typedef std::map< Poly_int_2, _Cache_entry_2,
            CGAL::Handle_id_less_than< Poly_int_2 > > Cache_2;

    struct _Cache_entry_1 { //! caching univariate data
        Poly_rat_1 poly_rat;
        Poly_1 poly_1;
        Vector_real_2 binoms;
        unsigned dfx_size;
    };

    typedef std::map< Poly_int_1, _Cache_entry_1, 
            CGAL::Handle_id_less_than< Poly_int_1 > > Cache_1;

    //! approximates the absolute value of \c NT
    Approximate_abs nt_abs;
    //! returns "spread" (disc around a center point)
    typename IA_traits::Disc disc;
    //! multiplies an interval by [-1;1]: not the best solution but we need
    //! to distinguish complex arithmetic from the real one here..
    typename IA_traits::template Mul_IA_exp< Approximate_abs, 1 > mul_IA_1;
    //! multiplies an interval by [-1;1]^2
    typename IA_traits::template Mul_IA_exp< Approximate_abs, 2 > mul_IA_2;
     //! converts \c Float_coeff to an affine form [center; spread]
    typename TC_traits::Coeff_to_affine fp2aff;
    //! computes a convex hull of \c NT interval (not used now)
//     typename IA_traits::template Hull< typename TC_traits::Bigfloat_hull >
//              nt_hull;

    //! affine "triple": represents: \hat(x) = x0 + x1\eps_1 + e\eps_{n+1}
    //! (one correlated unknown)
    struct Affine_triple {
        NT x0, x1; // center, first noise symbol 
        Real e; // non-affine error term (always real)
    };

    //! affine "quadruple": represents: 
    //! \hat(x) = x0 + x1\eps_1 + x2\eps_2 + e\eps_{n+1}
    //! (two correlated unknowns)
    struct Affine_quadruple {
        NT x0, x1, x2; // center, two noise symbols
        Real e; // non-affine error term (always real)
    };

    //!@}
public:
    //! \name Constructors
    //!@{

    //! default constructor
    Range_analysis_2() : poly_2_set(false), poly_1_set(false),
             mixed_dfs(0),binoms(0) {
        max_mixed_der_order = CGAL_IA_2_MIXED_DER_ORDER;
        max_dfx_order = CGAL_IA_1_BINOM_SIZE;
    }
 
    //!@}
    //! \name set & get
    //!@{ 

    //! sets up maximal order of mixed derivatives to be used for interval
    //! analysis: 0 - no mixed derivaties will be used
    void set_mixed_der_order(unsigned order) {
        bool need_precompute = (order > max_mixed_der_order);
        max_mixed_der_order = order;
        if(need_precompute && poly_2_set) {
            precompute_2(input_poly_2);
        }
    }

    //! sets up maximal order of derivatives (univariate case)
    void set_dfx_order(unsigned order) {
        bool need_precompute = (order > max_dfx_order);
        max_dfx_order = order;
        if(need_precompute && poly_1_set) {
            precompute_1(input_poly_1);
        }
    }
    
    //! sets up the polynomial (bivariate case)
    //! \c precision_change - indicates that internal polynomial needs to be
    //! recreated because FP precision has been changed
    //! \c normalize_coeffs - whether coefficients of the polynomial should be
    //! divided by geometric mean prior to conversion to floating-point
    void set_polynomial(const Poly_int_2& poly, bool precision_change = true,
                bool normalize_coeffs = true) {
        poly_2_set = true;
        precompute_2(poly, precision_change, normalize_coeffs);
    }

    //! sets up the polynomial (univariate case)
    //! \c precision_change - indicates that internal polynomial needs to be
    //! recreated because FP precision has been changed
    //! \c normalize_coeffs - whether coefficients of the polynomial should be
    //! divided by geometric mean prior to conversion to floating-point
    void set_polynomial(const Poly_int_1& poly, bool precision_change = true,
                    bool normalize_coeffs = true) {
        poly_1_set = true;
        precompute_1(poly, precision_change, normalize_coeffs);
    }

    //! univariate case: alternative version for a floating-point polynomial
    void set_polynomial(const Poly_1& poly) {
        poly_1_set = true;
        binoms = 0;
        precompute_1(poly);
    }
    
    //! returns bivariate polynomial (if set)
    Poly_int_2 polynomial_2() const {
        return input_poly_2;
    }

    //! bivariate polynomial with coefficients given by interval of doubles
    Poly_2 internal_poly_2() const {
        return poly_y;
    }

    //! returns univariate polynomial (if set)
    Poly_int_1 polynomial_1() const {
        return input_poly_1;
    }

    //! bivariate polynomial with coefficients given by interval of doubles
    Poly_1 internal_poly_1() const {
        return poly_1;
    }

    void clear_caches() const {
        cache_2.clear();
        cache_1.clear();
    }

    //!@}
    //! \name range evaluation methods
    //!@{ 

    //! Recursive Taylor bottom-up approach (unwinded recursion tree)
    //! using affine multiplications (works better in complex arithmetic)
    void eval_range_RT_2_aff(const Quad& q, NT& l, NT& h) {

        if(q.l.real() > q.r.real() || q.l.imag() > q.r.imag() ||
                q.b.real() > q.t.real() || q.b.imag() > q.t.imag()) {
            std::cerr << "FATAL: wrong quad: " << q << "\n";
            throw -1;
        }

        if(!poly_2_set) // top-level methods require precomputed data
            return;
        // NOTE: use const_reverse_iterator rbegin() ??
        typename std::vector< Poly_2 >::const_iterator di =
            mixed_dfs->end() - 1;
        int i, d = mixed_dfs_depth;
        if(d & 1) { // we start with an even tree level
                    // since we need only second-order derivaties
            di -= (d + 1), d--;
        }
        // bounds for 1st- and 2nd-order derivatives
        std::vector< Affine_triple > _1(d + 1);
        std::vector< Affine_quadruple > _2(d + 1);
        // evaluate at the highest level using AF1_2
        for(i = d; i >= 0; i--, di--) {
            eval_range_AF1_2_aff(*di, q, _2[i]);
        }

        if(d != 0) {

        d -= 1; // descend one level

        NT x0, x1, y0, y1;
        IA_to_AA(q.l, q.r, x0, x1);
        IA_to_AA(q.b, q.t, y0, y1);

        while(1) {
            // NOTE: *di now points to 'odd' derivative level
            // evaluate 1st-order derivatives at a point (x0, y0)
            for(i = d; i >= 0; i--, di--) {
                eval_point_AF1_2_aff(*di, x0, y0, _1[i]);
            }
            d -= 1; // descend one level

            // now compute bounds for level d using precomputed information
            for(i = d; i >= 0; i--, di--) {
                
                Affine_triple f;
                eval_point_AF1_2_aff(*di, x0, y0, f); // evaluate f(x0, y0)
                // *di uses first derivatives info from _1[i] & _1[i+1]

                Affine_triple fx = _1[i], fy = _1[i + 1];
    //! NOTE: important: observe that x1 & y1 are correlated in different
    //! variables !!!
                mul_AA_11(fx, x1);
                mul_AA_11(fy, y1);

                // second order derivatives are stored in _2 as affine forms
                Affine_quadruple fxx = _2[i],
                    fxy = _2[i + 1], fyy = _2[i + 2];

                mul_AA_01(fxx, x1 * Real(std::sqrt(0.5)));
                mul_AA_01(fyy, y1 * Real(std::sqrt(0.5)));
                mul_AA_11x11(fxy, x1, y1);
    
                Affine_quadruple r;
                r.x0 = f.x0 + fx.x0 + fy.x0 + fxx.x0 + fxy.x0 + fyy.x0;  
                r.e = f.e + fx.e + fy.e + fxx.e + fxy.e + fyy.e;  

                r.x1 = fx.x1 + fxx.x1 + fxy.x1 + fyy.x1; // correlations in
                r.x2 = fy.x1 + fxx.x2 + fxy.x2 + fyy.x2; // different vars !!
                
                 // i+2 to protect from overwriting
                _2[i + 2] = r;
            }
            // TODO: better solution without shifting ??
            for(i = 0; i < d + 1; i++) {
                _2[i] = _2[i+2]; // shift everything back
            }
            if(d == 0)
                break; 
            d -= 1;
        } // while(1)
        } // d != 0
        const Affine_quadruple& a = _2[0];
        Real spr = nt_abs(a.x1) + nt_abs(a.x2) + a.e;
        l = a.x0 - disc(spr), h = a.x0 + disc(spr);
    }

    //! Recursive Taylor bottom-up approach (unwinded recursion tree)
    //! using IA multplications (not sure which is better)
    void eval_range_RT_2(const Quad& q, NT& l, NT& h) {

        // NOTE: use const_reverse_iterator rbegin() ??
        typename std::vector< Poly_2 >::const_iterator di =
            mixed_dfs->end() - 1;
        int i, d = mixed_dfs_depth;
        if(d & 1) { // we start with an even tree level
                    // since we need only second-order derivaties
            di -= (d + 1), d--;
        }
        // bounds for 1st- and 2nd-order derivatives
        std::vector< std::pair< NT, NT > > _1(d + 1), _2(d + 1);
        // evaluate at the highest level using AF1_2
        for(i = d; i >= 0; i--, di--) {
            eval_range_AF1_2(*di, q, _2[i].first, _2[i].second);
        }
        if(d == 0) {
            l = _2[0].first, h = _2[0].second; // hurra!! done..!
            return;
        }
        d -= 1; // descend one level

        NT x0, x1, y0, y1;
        IA_to_AA(q.l, q.r, x0, x1);
        IA_to_AA(q.b, q.t, y0, y1);

        while(1) {
            // NOTE: *di now points to 'odd' derivative level
            // evaluate 1st-order derivatives at a point (x0, y0)
            for(i = d; i >= 0; i--, di--) {
                eval_point_AF1_2(*di, x0, y0, _1[i].first, _1[i].second);
            }
            d -= 1; // descend one level

            // now compute bounds for level d using precomputed information
            for(i = d; i >= 0; i--, di--) {
                NT fl, fh, fxl, fyl, fxh, fyh;
                eval_point_AF1_2(*di, x0, y0, fl, fh); // evaluate f(x0, y0)
                // *di uses first derivatives info from _1[i] & _1[i+1]
                fxl = _1[i].first, fxh = _1[i].second;
                fyl = _1[i+1].first, fyh = _1[i+1].second;

                mul_IA_1(fxl, fxh, x1);
                mul_IA_1(fyl, fyh, y1);
                fl += fxl + fyl, fh += fxh + fyh;

                // second order derivatives are stored in _2
                NT fxxl, fxxh, fxyl, fxyh, fyyl, fyyh;
                fxxl = _2[i].first, fxxh = _2[i].second;
                fxyl = _2[i+1].first, fxyh = _2[i+1].second;
                fyyl = _2[i+2].first, fyyh = _2[i+2].second;
    // in each iteration: fyy <- fxy, fxy <- fxx, fxx = new

                mul_IA_2(fxxl, fxxh, x1* Real(std::sqrt(0.5)));
                mul_IA_2(fyyl, fyyh, y1* Real(std::sqrt(0.5)));
//TODO TODO: probably lot more easier is to map the polynomial to
// symmetric range in real & imag parts: optimizations possible

// NOTE NOTE NOTE: this can be a potential problem because
    // of multipliying complex intervals..: [-1;1] * [-1;1] != [-1;1] !!!
                mul_IA_1(fxyl, fxyh, x1*y1);
                fxxl += fyyl + fxyl;
                fxxh += fyyh + fxyh;
                l = fl + fxxl, h = fh + fxxh;
                 // i+2 to protect from overwriting
                _2[i+2].first = l, _2[i+2].second = h;
            }
            // TODO: better solution without shifting ??
            for(i = 0; i < d + 1; i++) {
                _2[i] = _2[i+2]; // shift everything back
            }

            if(d == 0)
                break; 
            d -= 1;
        } // while(1)
        // not necessary to fetch from the array
        l = _2[0].first, h = _2[0].second; // hurra!! done..!
    }

    //! the same as AF1_2 but uses mixed derivatives to improve bounds
    //! NOTE: does not work for complex case
    void eval_range_AF1_RD_2(const Quad& q, NT& l, NT& h) {

        typename std::vector< Poly_2 >::const_iterator di =
            mixed_dfs->end() - 1;
        int i, d = mixed_dfs_depth;

        // bounds for 1st-order derivatives
        std::vector< std::pair< NT, NT > > _1(d + 1);
        // evaluate at the highest level using AF1_2
        for(i = d; i >= 0; i--, di--) {
            eval_range_AF1_2(*di, q, _1[i].first, _1[i].second);
        }
        if(d == 0) {
            l = _1[0].first, h = _1[0].second; // hurra!! done..!
            return;
        }
        d -= 1; // descend one level

        while(1) {
            // compute bounds for level d using precomputed information
            for(i = d; i >= 0; i--, di--) {
                // get precomputed bounds for derivatives in x and y
                NT fxl = _1[i].first, fxh = _1[i].second,
                   fyl = _1[i+1].first, fyh = _1[i+1].second;

                // if sign of partial derivaties agrees on q, we can use them
                // to improve bounds for f
        // TODO: use sign function instead of multiply (better for BigFloats)
                if(fxl * fxh > 0 && fyl * fyh > 0) {
                    NT xl(q.l), yl(q.b), xh(q.r), yh(q.t);
                    if(fxl < 0)
                        std::swap(xl, xh);
                    if(fyl < 0) 
                        std::swap(yl, yh);

                    NT _; // merge two intervals
                    eval_point_AF1_2(*di, xl, yl, l, _);
                    eval_point_AF1_2(*di, xh, yh, _, h);

                } else { // otherwise proceed in a usual way
                    eval_range_AF1_2(*di, q, l, h);
                }
                _1[i+1].first = l, _1[i+1].second = h;
            }
            for(i = 0; i < d + 1; i++) {
                _1[i] = _1[i + 1]; // shift everything back
            }
            if(d == 0)
                break; 
            d -= 1;
        } // while(1)
        l = _1[0].first, h = _1[0].second;
    }
     
    //! the same as below but converts AF to explicit interval
    inline void eval_range_AF1_2(const Poly_2& poly, const Quad& q, 
            NT& l, NT& h) {
            
        Affine_quadruple aq;
        eval_range_AF1_2_aff(poly, q, aq);
        //   NT spr = CGAL_ABS(z1) + CGAL_ABS(z2) + e; // spread
        Real s = nt_abs(aq.x1) + nt_abs(aq.x2) + aq.e;
        l = aq.x0 - disc(s), h = aq.x0 + disc(s);
    }

    //! evaluates the 2D range of a bivariate polynomial using 1st affine form
    //! returns the final affine form (with two correlated symbols)
    //! \c aq.x1 correlates with x-range, \c aq.x2 correlates with y-range
    void eval_range_AF1_2_aff(const Poly_2& poly, const Quad& q, 
            Affine_quadruple& aq) {
        //TODO: use CGAL::to_double_interval() method here ??

        NT y0, y2;
        IA_to_AA(q.b, q.t, y0, y2);
        Real a_y0 = nt_abs(y0), a_y2 = nt_abs(y2), s = a_y0 + a_y2, e;

        NT z0, z1, z2;
        // iterate over polynomials in x-variable
        typename Poly_2::const_iterator fxi = poly.end() - 1;
        Affine_triple a;
        eval_range_AF1_1_aff< false >(*fxi, q.l, q.r, a);

        z0 = a.x0, z1 = a.x1, z2 = NT(0), e = a.e;
        // TODO: certainly one can obtain better bounds when all powers of
        // x are artificially precomputed
        while((fxi--) != poly.begin()) {
//             e = a_y0 * e + y2 * (CGAL_ABS(z1) + CGAL_ABS(z2) + e);
            e = e * s + a_y2 * (nt_abs(z1) + nt_abs(z2));
            z1 = z1 * y0; z2 = z2 * y0 + z0 * y2;
            z0 = z0 * y0;
            // add the next x
            eval_range_AF1_1_aff< false >(*fxi, q.l, q.r, a);
            z0 += a.x0, z1 += a.x1, e += a.e;
        }
        aq.x0 = z0, aq.x1 = z1, aq.x2 = z2, aq.e = e;
    }

    // evaluates polynomial over the spherical boundary with radius
    // rad_x & rad_y
    Complex_sector eval_boundary_sector_2(const Poly_2& poly,
            const Complex_sector& x, const Complex_sector& y) const {
        
        Complex_sector r;

        typename Poly_2::const_iterator fxi = poly.end() - 1;
        r = eval_range_sector_1(*fxi, x);

        while((fxi--) != poly.begin()) {

            r = r.mul(y);
            Complex_sector z = eval_range_sector_1(*fxi, x);
            r = r.add(z);
        }
        return r;
    }

    Complex_sector eval_range_sector_1(const Poly_1& poly,
        const Complex_sector& x) const {

        typename Poly_1::const_iterator pi = poly.end() - 1;

        Real y0, y1;
        fp2aff(*pi, y0, y1);
        
        Complex_sector r(NT(y0, 0));
        while((pi--) != poly.begin()) {
            r = r.mul(x);
            fp2aff(*pi, y0, y1);
// 
//             r = r.add_interval(y0 - y1, y0 + y1);
            r = r.add(Complex_sector(y0 - y1, y0 + y1));
        }
        return r;
    }

//     Complex_sector eval_range_sector_1_alt(const Poly_1& poly,
//         const Complex_sector& x) const {
// 
//         typename Poly_1::const_iterator pi = poly.begin();
// 
//         Real y0, y1;
//         fp2aff(*pi, y0, y1);
//         
//         // TODO: special add for the case when the sector is just
//         // a real number
//         Complex_sector r(NT(y0, 0)), exp(x);
//         while((++pi) != poly.end()) {
//             fp2aff(*pi, y0, y1);
//             Complex_sector z = exp.mul(Complex_sector(NT(y0, 0)));
//             r = r.add(z);
//             exp = exp.mul(x);
//         }
//         return r;
//     }

    //! the same as below but converts AF to an explicit interval
    inline void eval_range_AF1_1(const Poly_1& poly, const NT& xl,
            const NT& xh, NT& l, NT& h) const {
            
        Affine_triple a;
        eval_range_AF1_1_aff< false >(poly, xl, xh, a);
        triple_to_IA(a, l, h);
    }

    template < bool UseBinomials >
    void eval_range_AF1_1_aff(const Poly_1& poly, const NT& xl, const NT& xh,
             Affine_triple& a) const {

//TODO: does it make sense to use AF2 or QF instead ??
        typename Poly_1::const_iterator pi = poly.end() - 1,
            begin = poly.begin();
        
        typename Vector_real_1::const_iterator bi;
        if(UseBinomials) {
            bi = binom_i->end() - 1;
        }

        Real _, e;  // error term (always real)
        fp2aff(*pi, _, e);
        if(UseBinomials) {
            _ *= (*bi), e *= (*bi), pi--;
        }
        NT y0(_), y1 = NT(0), x0, x1;

        if(poly.degree() != 0) {
            IA_to_AA(xl, xh, x0, x1);
            Real a_x0 = nt_abs(x0), a_x1 = nt_abs(x1), s = a_x0 + a_x1;

            while((UseBinomials && (bi--) != binom_i->begin()) ||
                  (!UseBinomials && (pi--) != begin)) {
                e = e * s + nt_abs(y1) * a_x1; // mul y by x0 + x1*e1
                y1 = y1 * x0 + x1 * y0;
                y0 = y0 * x0;
                Real p0, p1;
                fp2aff(*pi, p0, p1);
                if(UseBinomials) {
                    p0 *= (*bi), p1 *= (*bi), pi--;
                }
                // add p0 + p1*en to y
                y0 += p0, e += p1; // p1 must be positive
            }
        }
        a.x0 = y0, a.x1 = y1, a.e = e;
    }

    //! Recursive Taylor 1D case with affine multiplications & 
    //! variable Taylor \c ExpansionLength (2 or 3)
    //! \c rec_lenght controls the number of recursion steps
    template < unsigned ExpansionLength >
    void eval_range_RT_1_aff(const NT& xl, const NT& xh, NT& l, NT& h, 
                int rec_lenght = -1) const {

        if(!poly_1_set)
            return;

        if(ExpansionLength != 2 && ExpansionLength != 3)
            throw "NYI";

        NT x0, x1, y0, y1;
        IA_to_AA(xl, xh, x0, x1);

        if(poly_1.degree() == 0) {
            Real x0, x1;
            fp2aff(poly_1.lcoeff(), x0, x1);
            l = NT(x0 - x1), h = NT(x0 + x1);
            return;
        }

//      typename std::vector< Poly_1 >::const_iterator di = dfx_set.end() - 1;
    // NOTE: if rec_lenght > max_dfx_order => error
        int d = rec_lenght;
        if((unsigned)d >= dfx_size)
            d = (int)dfx_size - 1;

        //! we start with the level divisible by \c ExpansionLength
        d -= (d % ExpansionLength);
        binom_i = binoms->begin() + d;

        Affine_triple top_fx, f2, f1, f;
        eval_range_AF1_1_aff< true >(poly_1, xl, xh, top_fx);

    // NOTE: do we need d explicitly ??
        if(d != 0) {
        binom_i--, d--;

        while(1) { // here di points to derivative f2 level

            if(ExpansionLength == 3) { // compile-time decision
                eval_point_AF1_1< true >(poly_1, x0, f2);
                binom_i--, d--;
            }
            eval_point_AF1_1< true >(poly_1, x0, f1);
            binom_i--, d--;
            eval_point_AF1_1< true >(poly_1, x0, f);
    
            //! now compute: f + f1*[-1;1]*x1 + f2*[-1,1]^2*x1^2/2 +
            //!       + f3*[-1,1]^3*x1^3/6
            mul_AA_11(f1, x1);

//         mul_AA< false, true >(f1, NT(0), x1); // correlated
//         mul_AA< true, false >(f1, NT(0), x1); // zeronoise, uncorrelated
            typename IA_traits::template Symm_IA_exp< 2 >()
                (x1 * Real(std::sqrt(0.5)), y0, y1);
            
            if(ExpansionLength == 2) 
                mul_AA< false, false >(top_fx, y0, y1);
            else // shouldn't f2 be correlated ??
                mul_AA< true, false >(f2, y0, y1);

            if(ExpansionLength == 3) {
                typename IA_traits::template Symm_IA_exp< 3 >()
                    (x1 * Real(0.550321208149), y0, y1); // 6^(-1/3)
    // template < bool ZeroNoise, bool Correlated > inline void mul_AA
                // noise is non-zero for f3
                mul_AA< false, false >(top_fx, y0, y1); 
            }
            
            top_fx.x0 += f1.x0 + f.x0;
            top_fx.x1 += f1.x1 + f.x1;
            top_fx.e += f1.e + f.e;

            if(ExpansionLength == 3) {
                top_fx.x0 += f2.x0;
                top_fx.x1 += f2.x1;
                top_fx.e += f2.e;
            }
        
            if(d == 0)
                break;
            binom_i--, d--; // descend one more level & repeat
        } // while(1)
        }
        triple_to_IA(top_fx, l, h);
    }

    NT eval_point_2(const Poly_2& poly, const NT& x, const NT& y) {
        typename Poly_2::const_iterator fxi = poly.end() - 1;
        NT r = eval_point_1(*fxi, x);

        while((fxi--) != poly.begin()) {
            r = r * y + eval_point_1(*fxi, x);
        }
        return r;
    }

    NT eval_point_1(const Poly_1& poly, const NT& x) {
        typename Poly_1::const_iterator pi = poly.end() - 1;

        Real _, __;
        fp2aff(*pi, _, __);
        NT r(_);
        while((pi--) != poly.begin()) {
            fp2aff(*pi, _, __);
            r = r * x + NT(_);
        }
        return r;
    }

    //! the same as below but returns an explicit interval
    inline void eval_point_AF1_2(const Poly_2& poly, const NT& x, const NT& y,
            NT& l, NT& h) {

        Affine_triple a;
        eval_point_AF1_2_aff(poly, x, y, a);
        
        NT spr = disc(a.e);
        l = a.x0 - spr, h = a.x0 + spr;
    }

    //! evaluates a polynomial at a point \c (x,y)
    //! polynomial is assumed to have interval coefficients
    void eval_point_AF1_2_aff(const Poly_2& poly, const NT& x, const NT& y,
            Affine_triple& a) {

        typename Poly_2::const_iterator fxi = poly.end() - 1;
        // eval point returns x0 and e only (which is the same as interval)
        eval_point_AF1_1< false >(*fxi, x, a);
        NT z0 = a.x0;
        Real e = a.e, a_y = nt_abs(y);
        while((fxi--) != poly.begin()) {
            e = a_y * e; // mul by y
            z0 = z0 * y;
            // add the next x
            eval_point_AF1_1< false >(*fxi, x, a);
            z0 += a.x0, e += a.e;
        }
        a.x0 = z0, a.x1 = NT(0), a.e = e;
    }

    //! evaluates a polynomial at \c x
    //! polynomial is assumed to have interval coefficients
    template < bool UseBinomials >
    void eval_point_AF1_1(const Poly_1& poly, const NT& x, Affine_triple& a)
                 const {

// NOTE NOTE: this is not so useful unless polynomial coefficients
// are real intervals..
        typename Poly_1::const_iterator pi = poly.end() - 1,
            begin = poly.begin();
        
        typename Vector_real_1::const_iterator bi;
        if(UseBinomials) {
            bi = binom_i->end() - 1;
        }
        Real _, e;
        fp2aff(*pi, _, e);
        if(UseBinomials) {
            _ *= (*bi), e *= (*bi), pi--;
        }
        NT y0(_);

        if(poly.degree() != 0) {
            Real a_x = nt_abs(x);

             while((UseBinomials && (bi--) != binom_i->begin()) ||
                  (!UseBinomials && (pi--) != begin)) {
                e = a_x * e; // mul y by x, e is always positive
                y0 = y0 * x;
                Real p0, p1;
                fp2aff(*pi, p0, p1);
                 if(UseBinomials) {
                    p0 *= (*bi), p1 *= (*bi), pi--;
                }
                // add p0 + p1*e to y
                y0 += p0, e += p1; // p1 must be positive
            }
        }
        a.x0 = y0, a.x1 = NT(0), a.e = e;
    }
 
    void exact_range_2_complex(const Poly_2& poly,
        const Quad& q, NT& l, NT& h) {

// static unsigned ncalls = 0;
// static unsigned second_run = 0;

        l = NT(1e10,1e10), h = NT(-1e10,-1e10);

//         std::cout << "checking quad: " << q << "\n";

        unsigned nr = 5, ni = 5, i, j, k, s;
        NT dy = q.t - q.b, dx = q.r - q.l;
        Real syr = dy.real() / nr, syi = dy.imag() / ni;
        Real sxr = dx.real() / nr, sxi = dx.imag() / ni;

        Real me(1e10);

        Real yr = q.b.real();
        for(i = 0; i <= nr; i++, yr += syr) {
            Real yi = q.b.imag();
            for(j = 0; j <= ni; j++, yi += syi) {
                Real xr = q.l.real();
                for(k = 0; k <= nr; k++, xr += sxr) {
                    Real xi = q.l.imag();
                    for(s = 0; s <= ni; s++, xi += sxi) {
                        NT r = eval_point_2(poly, NT(xr,xi), NT(yr,yi));

//                 std::cout << "point: " << NT(xr,xi) << "; " << NT(yr,yi) <<
//                     "; r: " << r << "\n";
                me=std::min(me, std::abs(r));
//                         cplx_pts[ncalls].push_back(r);
                        l = NT(std::min(l.real(), r.real()),
                            std::min(l.imag(), r.imag()));
                        h = NT(std::max(h.real(), r.real()),
                            std::max(h.imag(), r.imag()));
                    }
                }       
            }
        }

// cplx_box[ncalls] = CGAL::Bbox_2(l.real(), l.imag(), h.real(), h.imag());
std::cerr << "magnitude min: " << me << " ********************\n\n";
//         
//         if(++ncalls == cplx_samples) {
//             if(second_run < 20) {
//                 second_run++, ncalls = 0;
// for(int ii=0;ii<cplx_samples;ii++)
//                 cplx_pts[ii].clear();
//             } else 
//             throw "done";
//         }
    }

    void exact_range_1_complex(const Poly_1& poly,
        const NT& xl, const NT& xh, NT& l, NT& h) {

        l = NT(1e10,1e10), h = NT(-1e10,-1e10);
// static unsigned ncalls = 0;

        unsigned nx = 70, ny = 70, i, j;
        NT diff = xh - xl;
        Real sx = diff.real() / Real((int)nx), sy = diff.imag() / Real((int)ny);

        Real x = xl.real();
        for(i = 0; i <= nx; i++, x += sx) {
            Real y = xl.imag();
            for(j = 0; j <= ny; j++, y += sy) {
// if(!(i==0||j==0||i==nx||i==ny)) continue;

                NT r = eval_point_1(poly, NT(x,y));
//                 cplx_pts[ncalls].push_back(r);

                l = NT(std::min(l.real(), r.real()),
                            std::min(l.imag(), r.imag()));
                h = NT(std::max(h.real(), r.real()),
                            std::max(h.imag(), r.imag()));
            }
        }
//  cplx_box[ncalls] = CGAL::Bbox_2(l.real(), l.imag(), h.real(), h.imag());

//         std::cout << "exact 1D range complex: [" << l << "; " << h << "]\n";
    }

    //!@}
protected:
    //!\name protected methods
    //!@{

    //! multiplies \c a = (x0, 0, e) by \c [-s,s]
    //! returns an affine triple (\c s is correlated with \c a.x1 )
    void mul_AA_11(Affine_triple& a, const NT& s) const {

        NT y0(0), y1 = s;
        NT s0 = a.x0 * y0, s1 = a.x0 * y1;

        // equivalent to mul_AA< true, true >(a, y0, y1)
        a.e = a.e * (nt_abs(y1) + nt_abs(y0));
        a.x0 = s0, a.x1 = s1;
    }

    //! multiplies \c aq by \c [-s,s]^2 ; no correlations assumed
    void mul_AA_01(Affine_quadruple& aq, const NT& s) const {

        NT y0, y1;
        typename IA_traits::template Symm_IA_exp< 2 > symm_IA_exp;
        symm_IA_exp(s, y0, y1);
        Real ey = nt_abs(y1);
        
        NT x0 = aq.x0, x1 = aq.x1, x2 = aq.x2;
        aq.e = aq.e * (nt_abs(y0) + ey) + ey * (nt_abs(x0) + nt_abs(x1) +
                 nt_abs(x2));
        aq.x0 = x0 * y0, aq.x1 = x1*y0, aq.x2 = x2*y0;
    }

    //! multiples \c aq by [-s;s] and [-t;t]. It is assumed that
    //! \c aq.x1 correlates with \c s while \c aq.x2 correlates with \c t
    void mul_AA_11x11(Affine_quadruple& aq, const NT& s, const NT& t) const {

//!NOTE: as an alternative we can use IA multiplication between
//! s and t , then multiply the resulting interval with \a aq affinnely
        NT x1 = aq.x1, x2 = aq.x2, y1 = s, y2 = t;
        aq.e = nt_abs(y1) * (nt_abs(x1) + nt_abs(x2) + aq.e);
        x1 = aq.x0 * y1, x2 = 0;

        aq.e = nt_abs(y2) * (nt_abs(x1) + aq.e);
        aq.x0 = 0, aq.x1 = 0, aq.x2 = 0;
    }

    //! multiplies two affine forms (no correlations assumed)
    //! \c y0 - center, \c y1 - spread
    template < bool ZeroNoise, bool Correlated >
    inline void mul_AA(Affine_triple& a, const NT& y0, const NT& y1) const {
//! a = x0 + x1\eps_1 + e 
//! y = y0 + y1\eps_2 (no correlated variables)
//! x0' = x0*y0, x1' = x1*y0
//! e' = |y0|*e + |y1|*(|x1| + e) + |x0*y1|
//! NOTE: the last term |x0*y1| is due to the fact that we do not introduce
//! new variable (just put it to the error term)
        Real ey = nt_abs(y1);
        NT x1 = (ZeroNoise ? NT(0) : a.x1); // and pray for compiler constant
                                            // propagation..
        a.e = a.e * (nt_abs(y0) + ey) + ey * nt_abs(x1);
        a.x1 = x1 * y0;
        if(Correlated) {
            a.x1 = a.x1 + a.x0 * y1;
        } else { // if no correlated vars: add to the error term
            a.e = a.e + ey * nt_abs(a.x0);
        }
        a.x0 = a.x0 * y0;
    }

    //! converst an interval to affine form
    inline void IA_to_AA(const NT& l, const NT& h, NT& c, NT& spread) const {
        c = (l + h) * Real(0.5), spread = (h - l) * Real(0.5);
    }

    inline void triple_to_IA(const Affine_triple& a, NT& l, NT& h) const {
        Real s = nt_abs(a.x1) + a.e;
        l = a.x0 - disc(s), h = a.x0 + disc(s);
    }

    //!@}
public:
    //!\name precompute methods
    //!@{

    //! initializes with univariate polynomial
    void precompute_1(const Poly_int_1& in, bool precision_change = true,
            bool normalize_coeffs = true) {

        input_poly_1 = in;
        typename Cache_1::iterator ci = cache_1.find(in);
       
        bool found = (ci != cache_1.end());
        if(!found) {
            typedef internal::Convert_and_normalize1< Poly_int_1,
                 Arithmetic_kernel > C_and_n;
            _Cache_entry_1 ce;

            Rational factor(1);
            ce.poly_rat = C_and_n()(in, factor, normalize_coeffs);
            ci = cache_1.insert(std::make_pair(in, ce)).first;
        } else 
             ;//std::cout << "prefetching from cache..\n";

        if(!found || precision_change) {
            typedef typename CGAL::Coercion_traits< Rational, 
                Float_coeff >::Cast Cast;        
            internal::Transform1< Poly_1, Poly_rat_1, Cast > convert_poly;
            // polynomial with outermost var y
            ci->second.poly_1 = convert_poly(ci->second.poly_rat);
        } 
        poly_1 = ci->second.poly_1;
        binoms = &ci->second.binoms;
        dfx_size = ci->second.dfx_size;

        if(!found) { // do not need to recompute binomials even for 
                     // precision change
            precompute_1(poly_1);
            ci->second.dfx_size = dfx_size;
        }
    }

    //! precomputes first binomial coefficients for polynomial derivatives 
    void precompute_1(const Poly_1& poly) {
//TODO: find out what's wrong with Gmpfis: 

        poly_1 = poly;
        if(binoms == 0) { // indicates that no int poly has been set
            static Poly_int_1 ps(Integer(0)); // reserve one entry for the case
                                        // when no integer poly is given
            _Cache_entry_1 ce;
            ce.poly_1 = poly;
            typename Cache_1::iterator ci = cache_1.insert(
                    std::make_pair(ps, ce)).first;
            binoms = &ci->second.binoms;
        }

        int degree = poly.degree();

        dfx_size = std::min((unsigned)degree, max_dfx_order);
        //! the i-th row is binomial coefficient C^(j)_(i+1)
        //! where j = [degree..i]

        binoms->clear();
        binoms->push_back(Vector_real_1(degree + 1, Real(1)));
        Vector_real_1 v(degree);

        for(int i = 0; i < (int)dfx_size; i++) {
            for(int j = 1; j <= degree - i; j++) {
                v[j - 1] = (i == 0 ? Real(j) : v[j] * Real(j));
            }
            binoms->push_back(v);
            v.pop_back(); // decrease length
        }

//         std::cout << "------------------ binoms:\n";
//         typename Vector_real_2::const_iterator di = binoms->begin();
//         for(; di != binoms->end(); di++) {
//             typename Vector_real_1::const_iterator ddi = di->begin();
//             for(; ddi != di->end(); ddi++) {
//                 std::cout << *ddi << " ";
//             }
//             std::cout <<  "\n";
//         }
    }

    void precompute_2(const Poly_int_2& in, bool precision_change = true,
            bool normalize_coeffs = true) {

        input_poly_2 = in;
        typename Cache_2::iterator ci = cache_2.find(in);
        
        bool found = (ci != cache_2.end());
        if(!found) {
            _Cache_entry_2 ce;
            typedef internal::Convert_and_normalize1< Poly_int_2,
                 Arithmetic_kernel > C_and_n;

            Rational factor(1);
            ce.poly_rat = C_and_n()(in, factor, normalize_coeffs);
            ci = cache_2.insert(std::make_pair(in, ce)).first;
        }
        
        if(!found || precision_change) {

            typedef typename CGAL::Coercion_traits< Rational, 
                    Float_coeff >::Cast Cast;        
            internal::Transform1< Poly_2, Poly_rat_2, Cast > convert_poly;
            typedef CGAL::Polynomial_traits_d< Poly_2 > P2_traits;

            poly_y = convert_poly(ci->second.poly_rat); // y is outer var
            poly_x = typename P2_traits::Swap()(poly_y, 0, 1); // x is outer
            
            unsigned i, j, idx;
            typename P2_traits::Degree degree;
            // maximum of degree in x- and degree in y-variable
            mixed_dfs_depth = std::max(degree(poly_y, 0), degree(poly_y, 1));
            if(mixed_dfs_depth > max_mixed_der_order)
                mixed_dfs_depth = max_mixed_der_order;

            mixed_dfs = &ci->second.mixed_dfs;
            // f, fx, fy, fxx, fxy, fyy, fxxx, fxxy, fxyy, fyyy, ...
            mixed_dfs->clear();
            mixed_dfs->push_back(poly_y);

            typename P2_traits::Differentiate diff;
            // level i has (i + 1) nodes
            for(i = 1, idx = 0; i < mixed_dfs_depth + 1; i++) {
//             NT inv = NT(1) / NT(j); // need make_exact in case of bigfloats
                Poly_2 py = diff((*mixed_dfs)[idx], 0); // f_x
                //! multiply by inverse ??  reduce_coeffs(py, NT(1));
                mixed_dfs->push_back(py);
                for(j = idx; j < idx + i; j++) {
                    Poly_2 py = diff((*mixed_dfs)[j], 1); // f_y
                //  reduce_coeffs(py, NT(1));
                    mixed_dfs->push_back(py);
                }
                idx += i;
            }
            ci->second.poly_y = poly_y;
            ci->second.poly_x = poly_x;
            ci->second.mixed_dfs_depth = mixed_dfs_depth;

        } else {
            poly_y = ci->second.poly_y;
            poly_x = ci->second.poly_x;
            mixed_dfs = &ci->second.mixed_dfs;
            mixed_dfs_depth = ci->second.mixed_dfs_depth;
        }
#if 0
        typedef CGAL::Polynomial< double > Poly_double_1;
        typedef CGAL::Polynomial< Poly_double_1 > Poly_double_2;
        typename std::vector< Poly_2 >::const_reverse_iterator dit =
             mixed_dfs.rbegin();
        for(i = mixed_dfs_depth; (int)i >= 0; i--) {

            std::cout << i << "th mixed derivatives: " << std::endl;
            for(j = 0; j < i+1; j++, dit++) {
                Poly_double_2 pd = CGAL::to_double(*dit);
                std::cout << pd << std::endl;
            }
        }
#endif
    }

    //!@}   
protected:
    //! \name protected properties 
    //!@{ 

    bool poly_2_set, poly_1_set;

    Poly_int_2 input_poly_2;      //! input polynomial 
    Poly_2 poly_y, poly_x;      //! f(y(x)) & f(x(y))

    unsigned max_mixed_der_order;
    //! depth of mixed derivaties tree (starting with 0)
    unsigned mixed_dfs_depth; //! accordingly, level i has (i + 1) nodes
    std::vector< Poly_2 > *mixed_dfs;
    mutable Cache_2 cache_2; // caching bivariate data

    Poly_int_1 input_poly_1;      //! input polynomial (univariate) 
    Poly_1 poly_1;                //! internal poly_1
    unsigned max_dfx_order; //! maximal order of derivatives to consider
    unsigned dfx_size;
    //std::vector< Poly_1 > dfx_set;  //! set of derivatives f^(i)
    Vector_real_2 *binoms;

    // binomial iterator    
    mutable typename Vector_real_2::const_iterator binom_i; 
    mutable Cache_1 cache_1; // caching univariate data

    //!@}
}; // class Range_analysis_2


} // namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_2_RANGE_ANALYSIS_2
// EOF

