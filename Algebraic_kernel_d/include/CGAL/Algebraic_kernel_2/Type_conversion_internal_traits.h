// Copyright (c) 2010, 2011, 2012 Max-Planck-Institut fuer Informatik (Germany).
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

#ifndef CGAL_ALGEBRAIC_KERNEL_2_TYPE_CONVERSION_INTERNAL_TRAITS_H
#define CGAL_ALGEBRAIC_KERNEL_2_TYPE_CONVERSION_INTERNAL_TRAITS_H

#include <CGAL/config.h>
#include <CGAL/function_objects.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/Arithmetic_kernel.h>
#include <complex>

/*! \file Type_conversion_internal_traits.h
 *  \brief various type conversions to support interval arithmetic over
 *  different number types
 */
 
namespace CGAL {

// transformation routine
namespace internal {

  // @Pavel: Replace with CGAL::abs and CGAL::sign
#ifndef CGAL_ABS
#define CGAL_ABS(x) ((x) < 0 ? -(x): (x))
#endif

#ifndef CGAL_SGN
#define CGAL_SGN(x) ((x) > 0 ? 1 : ((x) < 0 ? -1 : 0))
#endif

template <class NT>
struct Get_minmax1 {

    template <class X>
    void operator()(const CGAL::Polynomial<X>& p, NT& pmin, NT& pmax) const
    {
        typename CGAL::Polynomial<X>::const_iterator it = p.begin();
        Get_minmax1< NT > minmax;
        minmax(*it, pmin, pmax);
        while(++it != p.end()) {
            NT smin, smax;
            minmax(*it, smin, smax);
            if(pmax < smax)
                pmax = smax;
            if(pmin == NT(0) || (smin != NT(0) && smin < pmin))
                pmin = smin;
        }
    }
    void operator()(const NT& x, NT& xmin, NT& xmax) const
    {  xmin = CGAL_ABS(x), xmax = xmin; }
};

/*!\brief
 * divides input value by a contant value
 *
 * provided that there is a coercion between \c Input and \c Result types
 */
template <class Result, class Input>
struct Reduce_by1 {

    typedef Input argument_type;
    typedef Result result_type;

    Reduce_by1(const Input& denom_) :
        denom(cast(denom_)) {
    }
    
    Result operator()(const Input& x) {
        return (cast(x)/denom);
    }

    typename CGAL::Coercion_traits<Input, Result>::Cast cast;
    Result denom;
};

/*!\brief
 * transforms bivaritate polynomial of type \c InputPoly_2 to
 * \c OutputPoly_2 by recursively applying operation \c Op to all of its
 * coefficients
 *
 * <tt>Op: InputPoly_2::Inntermost_coefficient_type ->
 *             OutputPoly_2::Inntermost_coefficient_type</tt>
 */
template <class OutputPoly_2, class InputPoly_2, class Op>
struct Transform1 {

    typedef InputPoly_2  first_argument_type;
    typedef Op           second_argument_type;
    typedef OutputPoly_2 result_type;

    template <class X>
    OutputPoly_2 operator()(const CGAL::Polynomial<X>& p, Op op = Op()) const {

        Transform1<typename OutputPoly_2::NT, typename InputPoly_2::NT, Op> tr;
        return OutputPoly_2(
            ::boost::make_transform_iterator(p.begin(), std::bind2nd(tr, op)),
            ::boost::make_transform_iterator(p.end(), std::bind2nd(tr, op)));
    }

    OutputPoly_2 operator()(
        const typename Innermost_coefficient_type<InputPoly_2>::Type& x,
                Op op) const {
        return static_cast<OutputPoly_2>(op(x));
    }
};

//! converts an integer polynomial to rational one;
//! coefficients are normalized by dividing them by geometric mean of min and
//! max coefficient
template < class PolyInt_d, class ArithmeticKernel >
struct Convert_and_normalize1 {

    //! this first template argument
    // NOTE: apparently PolyInt_d coefficient type must match the one provided
    // by arithmetic kernel
    typedef PolyInt_d Poly_int_d;
    //! this second template argument
    typedef ArithmeticKernel Arithmetic_kernel;

    //! extract number types
    typedef typename Arithmetic_kernel::Integer Integer;
    typedef typename Arithmetic_kernel::Rational Rational;
    typedef typename Arithmetic_kernel::Bigfloat_interval BFI;

    //! dimensionality
    static const int d = CGAL::internal::Dimension< PolyInt_d >::value;
    //! polynomial over rationals
    typedef typename CGAL::Polynomial_type_generator< Rational, d >::Type
        Poly_rat_d;

    typedef Poly_rat_d result_type;

    //! returns normalization \c factor if \c normalize is true
    Poly_rat_d operator()(const Poly_int_d& in,
            Rational& factor, bool normalize = true) const {

        typedef typename CGAL::Bigfloat_interval_traits< BFI >::Bound BigFloat;
        typedef CGAL::Polynomial_traits_d< PolyInt_d > PT;
        
        factor = Rational(1);
        if(normalize) {
            Integer pmin(0), pmax(0); // pmin is guaranteed to be > 0
            Get_minmax1< Integer >()(in, pmin, pmax);
    
            // approximate square root with lower precision
            long old_prec = CGAL::set_precision(BFI(), 40);
            BigFloat bf(pmax*pmin);
            // divide by geometric mean: sqrt(pmin*pmax)
            bf = CGAL::sqrt(bf);
            factor = typename CGAL::Coercion_traits< BigFloat, Rational >::
                 Cast()(bf);
            CGAL::set_precision(BFI(), old_prec);
        }

        typedef Reduce_by1< Rational, Rational > Reduce_op;
        Transform1< Poly_rat_d, Poly_int_d, Reduce_op > transform;
        Reduce_op op(factor);

        Poly_rat_d res = transform(in, op); // divides by ratio inside
        return res;
    }
};

//! support for complex NTs in interval analysis
template < class NT_ >
struct Complex_IA_traits {

    //! this instance's template argument
    typedef NT_ NT;
    //! real NT: the same as NT by default
    typedef NT_ Real;

    //! raises symmetric interval \c [-s;s] to power \c e
    template < unsigned e >
    struct Symm_IA_exp {
        inline void operator()(const NT& s, NT& y0, NT& y1) const {

            y0 = NT(0);
            if(e == 1) {            
                y1 = s;
            } else if(e == 2) {
                y1 = s*s;
            } else if(e == 3) {
                y1 = s*s*s;
            } else
                throw "NYI";
        }
    };

    //! multiply an interval \c [a;b] by \c [-s;s]^e (s > 0)
    //! certainly not the best solution but partial templates are not
    //! supported for class methods..
    //! and we need to tell real NT from complex ones
    template < class ApproximateAbs, unsigned e >
    struct Mul_IA_exp {
        inline void operator()(NT& a, NT& b, const NT& s) const {
            if(e & 1) { // [-1;1]^(2n+1) = [-1;1]
                // min(a,-a,b,-b) = min(a, -b) since a < b
                // max(a,-a,b,-b) = max(-a, b) since a < b
                NT na(a), se(s);
                if(e == 3)
                    se = s*s*s;
                a = std::min(a, -b) * se;
                b = std::max(-na, b) * se;

            } else { // [-1;1]^(2n) = [0;1]
                NT sq = s * s;
                a = std::min(a, NT(0)) * sq;
                b = std::max(b, NT(0)) * sq;
            }
        }
    };
 
    //! returns disc radius (spread) around a center point
    struct Disc {
        NT operator()(const Real& x) const {
            return x;
        }
    };

    //! absolute value (approximate)
    template < class MagnitudeApproximator = void >
    struct Approximate_abs {
        Real operator()(const NT& x) const {
            return CGAL_ABS(x);
        }
    };

     //! hull of an interval (more detailed explanations will follow)
    template < class BigfloatHull >
    struct Hull {
        void operator()(NT& a, NT& b) const {
            BigfloatHull hull;
            hull(a, b);
        }
    };
};

template < class NT_ >
struct Complex_IA_traits< std::complex< NT_ > > {

    //! template argument
    typedef std::complex< NT_ > NT;
    //! real/imag number type
    typedef typename NT::value_type Real;

    struct Disc {
        NT operator()(const Real& x) const {
            return NT(x, x);
        }
    };

    //! raises a symmetric interval \c [-s;s] (s > 0) to power \c e
    //! returns a middle point \c y0 and a radius \c y1
    template < unsigned e >
    struct Symm_IA_exp {
        inline void operator()(const NT& s, NT& y0, NT& y1) const {

            if(e == 1) { // compile-time decision
                y0 = NT(0), y1 = s;

            } else if(e == 2) {
            // [-c-di; c+di]^2 = [-d^2-2cdi; c^2+2cdi]
                Real c = s.real(), d = s.imag();
                Real lr = -d*d, hr = c*c, irad = Real(2) * c * d;
                y0 = NT((lr + hr) * Real(0.5), Real(0));
                y1 = NT((hr - lr) * Real(0.5), irad);

            } else if(e == 3) {
            // [-c-di; c+di]^3 = [min(-c^3,c^3-3cd^2)+i*min(-d^3,d^3-3c^2d);..]
                Real c = s.real(), d = s.imag(), c2 = c*c, c3 = c2*c,
                        d2 = d*d, d3 = d2*d;
                Real p = Real(3.0)*c*d2 - c3, q = Real(3.0)*c2*d - d3;
    // Real: min(-c3,-p); max(c3, p); Imag: min(-d3,-q); max(d3, q)
                // final interval: [-p-qi; p+qi] (p > 0, q > 0)
                if(p < c3)
                    p = c3;
                if(q < d3)
                    q = d3;
                y0 = 0, y1 = NT(p, q);
            }
        }
    };
    
    //! multiply an interval \c [a;b] by \c [-s;s]^e (s > 0)
    template < class ApproximateAbs, unsigned e >
    struct Mul_IA_exp {
        inline void operator()(NT& a, NT& b, const NT& s) const {
            // multiples [a;b] by [-c-di;c+di]
            NT x0 = (a + b) * Real(0.5), x1 = (b - a) * Real(0.5);
            NT y0, y1;        

            Symm_IA_exp< e >(s, y0, y1);
            NT s0 = x0*y0, s1 = x0 * y1, s2 = x1 * y0, e = x1 * y1;

            ApproximateAbs nt_abs;
            Real re = nt_abs(s1) + nt_abs(s2) + nt_abs(e);
            NT spr(re, re);
            a = s0 - spr, b = s0 + spr;
        }
    };

    //! absolute value (approximate)
    template < class MagnitudeApproximator >
    struct Approximate_abs {

        Real operator()(const NT& x) const {
            return MagnitudeApproximator()(x.real(), x.imag());
        }
    };
    
    //! hull of an interval (more detailed explanations will follow)
    template < class BigfloatHull >
    struct Hull {
        void operator()(NT& a, NT& b) const {
            BigfloatHull hull;
            // separately for complex and reals
            hull(a.real(), b.real());
            hull(a.imag(), b.imag());
        }
    };
};

/*!\brief
 * class template \c Type_conversion_traits_base
 */
template <class Coeff_, class Integer_, class Rational_, class Float_ >
struct Type_conversion_traits_base
{ 
    //! type of innermost polynomial coefficients
    typedef Coeff_ Coeff;
    
    //! an integer number type
    typedef Integer_ Integer; 

    //! rational
    typedef Rational_ Rational;

    //! float
    typedef Float_ Float;

    typedef typename CGAL::Coercion_traits< Rational, Float >::Cast
        Rat_to_float;

    typedef typename CGAL::Coercion_traits< Float, Rational >::Cast
        Float_to_rat;

    //! provided for convenience when there exists an implicit coercion
    //! between number types
    template <class To>
    struct Implicit_coercion {
        typedef To result_type;
        
        template <class From> 
        inline To operator()(const From& x) const {
            return static_cast<To>(x);
        }
    };

    struct Float_to_int {
        typedef int result_type;
        
        template <class Float>
        inline int operator()(const Float& x) const
        { return static_cast<int>(std::floor(CGAL::to_double(x))); }
            //return static_cast<int>(std::floor(x)); }
    };

    //! converts \c Coeff to an affine form (center and radius)
    struct Coeff_to_affine {
        inline void operator()(const Coeff& c, Float& x0, Float& x1) const {
            x0 = static_cast< Float >(c); x1 = Float(0);
        }
    };

    //! vector magnitude: exact and approximate ones (defaults to exact)
    struct Exact_magnitude {
        inline Float operator()(const Float& a, const Float& b) const {
     // TODO: change this to
            return std::abs(std::complex< Float >(a, b));
//            return std::sqrt(a*a + b*b);
        }
    };

    struct Approximate_magnitude {

        Float operator()(const Float& _a, const Float& _b) const {

            static const double threshold = 1e-24;

            double ad = CGAL::to_double(_a), bd = CGAL::to_double(_b);
            ad *= ad, bd *= bd;
            // try with double precision first
            if(ad > threshold && bd > threshold) {
                return std::sqrt(ad + bd);
            }
            // rough approximation of a square root
            Float a = CGAL_ABS(_a), b = CGAL_ABS(_b);
            return a + b;
/*           return CGAL::sqrt(_a*_a + _b*_b);
            if(a < b)
                std::swap(a, b);
            Float z = Float(15.0 / 16.0)*a + Float(15.0 / 32.0)*b;
//             Float d = a*a + b*b;
//             if(y < 1e-10)
//                 return y;
//             Float z = y - (y*y - d) / (2.0 * y);
//             return z;*/
        }
    };

    typedef Exact_magnitude Vector_magnitude;

    //! computes a hull of interval [a; b]: relevant only for \c Floats which
    //! accumulate error internally, e.g., \c CORE::BigFloat
    struct Bigfloat_hull {
        inline void operator()(Float& a, Float& b) const {
            // bla bla.. nothing interesting happening here.. discouraged ?
        }
    };
};

template <class Float1, class Rational>
struct Type_conversion_internal_traits
{  };


#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL

//! Specialization for \c CGAL::Interval_nt<true>
template <>
struct Type_conversion_internal_traits<CGAL::Interval_nt<true>, CORE::BigRat >:
        public Type_conversion_traits_base< CGAL::Interval_nt<true>, int,
            CORE::BigRat, double > {

    //! converts \c Coeff to an affine form (center and radius)
    struct Coeff_to_affine {
        void operator()(const Coeff& c, Float& x0, Float& x1) const {
            x0 = (c.inf() + c.sup()) * Float(0.5);
            x1 = (c.sup() - c.inf()) * Float(0.5);
        }
    };
 
    struct Rat_to_float {
        typedef Float result_type;
                
        template <class Extended>
        Float operator()(const Extended& x) const {
            return CGAL::to_double(x); 
        }
    };

    typedef Implicit_coercion<double> Float_to_rat;

    struct Rat_to_integer {
        typedef Rational argument_type;
        typedef Integer result_type;

        Integer operator()(const Rational& x) const {
            return static_cast<int>(std::floor(CGAL::to_double(x)));
        }
    };
};

//! Specialization for \c CORE::BigFloat
template <>
struct Type_conversion_internal_traits<CORE::BigFloat, class CORE::BigRat>
         : public Type_conversion_traits_base< CORE::BigFloat, CORE::BigInt,
                CORE::BigRat, CORE::BigFloat> {

    struct Rat_to_integer {
        typedef Rational argument_type;
        typedef Integer result_type;
        
        Integer operator()(const Rational& x) const { 
            return x.BigIntValue(); 
        }
    };

    //! replace by approximate magnitude because square roots are expensive
    //! for bigfloats
    typedef Approximate_magnitude Vector_magnitude;

    //! \c c has some error introduced during conversion from rationals
    //! we can incorporate it explicitly as an affine noise term: this gives
    //! us better control over the algorithm accuracy
    struct Coeff_to_affine {
        void operator()(const Coeff& c, Float& x0, Float& x1) const {
            x0 = CGAL::median(c);
            x1 = Float(Integer(c.err()), 0, c.exp());
        }
    };

    //! expands an interval accounting for accumulated errors: specialization
    //! provided for \c CORE::BigFloat because it is represented by an interval
    struct Bigfloat_hull {
        inline void operator()(Float& a, Float& b) const {
            a = CGAL::lower(a), b = CGAL::upper(b);
        }
    };
};

//! Specialization for \c CORE::BigRat
template <>
struct Type_conversion_internal_traits<CORE::BigRat, CORE::BigRat> :
    public Type_conversion_traits_base<CORE::BigRat, CORE::BigInt,
        CORE::BigRat, CORE::BigRat> {

//     typedef Implicit_coercion<Float> Float_to_rat;

    struct Rat_to_integer {
        typedef Rational argument_type;
        typedef Integer result_type;
        
        Integer operator()(const Rational& x) const { 
            return x.BigIntValue(); 
        }
    };

};
#endif // CGAL_HAS_CORE_ARITHMETIC_KERNEL

#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL

// namespace {
// 
// typedef CGAL::GMP_arithmetic_kernel::Integer BigInt;
// typedef CGAL::GMP_arithmetic_kernel::Rational BigRat;
// typedef CGAL::GMP_arithmetic_kernel::BigFloat_interval BigFloat;
// 
// }

template <>
struct Type_conversion_internal_traits< CGAL::Interval_nt<true>, CGAL::Gmpq > :
        public Type_conversion_traits_base< CGAL::Interval_nt<true>, int,
            CGAL::Gmpq, double > {

    struct Rat_to_float {
        typedef Float result_type;
                
        template <class Extended>
        Float operator()(const Extended& x) const {
            return CGAL::to_double(x); 
        }
    };

    struct Rat_to_integer {
        typedef Rational argument_type;
        typedef Integer result_type;

        Integer operator()(const Rational& x) const {
            return static_cast<int>(std::floor(CGAL::to_double(x)));
        }
    };

    //! converts \c Coeff to an affine form (center and radius)
    struct Coeff_to_affine {
        void operator()(const Coeff& c, Float& x0, Float& x1) const {
            x0 = (c.inf() + c.sup()) * Float(0.5);
            x1 = (c.sup() - c.inf()) * Float(0.5);
        }
    };
};

template <>
struct Type_conversion_internal_traits< CGAL::Gmpfi, CGAL::Gmpq >
         : public Type_conversion_traits_base< CGAL::Gmpfi, CGAL::Gmpz,
                CGAL::Gmpq, CGAL::Gmpfr > {

    struct Rat_to_integer {
        typedef Rational argument_type;
        typedef Integer result_type;
        
        Integer operator()(const Rational& x) const { 
            return x.numerator(); 
        }
    };

    //! replace by approximate magnitude because square roots are expensive
    //! for bigfloats
    typedef Approximate_magnitude Vector_magnitude;

    struct Coeff_to_affine {
        void operator()(const Coeff& c, Float& x0, Float& x1) const {
    //NOTE: use CGAL::lower / upper ??
            x0 = (c.inf() + c.sup()) * Float(0.5);
            x1 = (c.sup() - c.inf()) * Float(0.5);
        }
    };

//     typedef Rat_to_float<Float> Rat_to_float;
    
//     struct Float_to_rat {
//         typedef Float argument_type;
//         typedef Rational result_type;
//
//         Rational operator()(const Float& x) const
//         {   std::cout << "NYI\n"; throw -1; }
//     };
};

template <>
struct Type_conversion_internal_traits< CGAL::Gmpq, CGAL::Gmpq > :
    public Type_conversion_traits_base< CGAL::Gmpq, CGAL::Gmpz,
        CGAL::Gmpq, CGAL::Gmpq > {

    struct Rat_to_integer {
        typedef Rational argument_type;
        typedef Integer result_type;
        
        Integer operator()(const Rational& x) const { 
            return x.numerator();  // assumes denominator is zero
        }
    };

};
#endif // CGAL_HAS_GMP_ARITHMETIC_KERNEL

#ifdef CGAL_HAS_LEDA_ARITHMETIC_KERNEL

template <>
struct Type_conversion_internal_traits<CGAL::Interval_nt<true>, leda::rational > :
        public Type_conversion_traits_base<CGAL::Interval_nt<true>, int,
            leda::rational, double > {
 
    struct Rat_to_float {
        typedef Float result_type;
                
        template <class Extended>
        Float operator()(const Extended& x) const {
            return CGAL::to_double(x); 
        }
    };

    struct Rat_to_integer {
        typedef Rational argument_type;
        typedef Integer result_type;

        Integer operator()(const Rational& x) const {
            return static_cast<int>(std::floor(CGAL::to_double(x)));
            //return static_cast<int>(x.to_double());
        }
    };

    //! converts \c Coeff to an affine form (center and radius)
    struct Coeff_to_affine {
        void operator()(const Coeff& c, Float& x0, Float& x1) const {
            x0 = (c.inf() + c.sup()) * Float(0.5);
            x1 = (c.sup() - c.inf()) * Float(0.5);
        }
    };
};

//! Specialization for \c leda::bigfloat
template <>
struct Type_conversion_internal_traits< leda_bigfloat_interval,
    class leda::rational>
    : public Type_conversion_traits_base<leda_bigfloat_interval,
         leda::integer, leda::rational, leda::bigfloat> {

    struct Rat_to_integer {
        typedef Rational argument_type;
        typedef Integer result_type;
        
        Integer operator()(const Rational& x) const { 
            return leda::floor(x); 
        }
    };

    struct Rat_to_float {
        typedef Float result_type;

        // no implicit coercion between leda rational and floats, therefore
        // decompose rational to compute the result
        Float operator()(const Rational& x) const {
            return (static_cast<Float>(x.numerator()) /
                    static_cast<Float>(x.denominator()));
        }
    };

     struct Coeff_to_affine {
        void operator()(const Coeff& c, Float& x0, Float& x1) const {
    //NOTE: use CGAL::lower / upper ??
            x0 = (c.lower() + c.upper()) * Float(0.5);
            x1 = (c.upper() - c.lower()) * Float(0.5);
        }
    };

    //! replace by approximate magnitude because square roots are expensive
    //! for bigfloats
    typedef Approximate_magnitude Vector_magnitude;
};

//! Specialization for \c leda::rational
template <>
struct Type_conversion_internal_traits<leda::rational, leda::rational> :
    public Type_conversion_traits_base<leda::rational, leda::integer,
        leda::rational, leda::rational> {

    struct Rat_to_integer {
        typedef Rational argument_type;
        typedef Integer result_type;
        
        Integer operator()(const Rational& x) const { 
            return leda::floor(x);  
        }
    };
};

#endif // CGAL_USE_LEDA

} // namespace internal

} // namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_2_TYPE_CONVERSION_INTERNAL_TRAITS_H
// EOF

