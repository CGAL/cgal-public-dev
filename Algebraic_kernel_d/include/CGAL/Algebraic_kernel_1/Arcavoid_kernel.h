#ifndef CGAL_ARCAVOID_KERNEL_H
#define CGAL_ARCAVOID_KERNEL_H

/*****************************************************************
 * Requirements for input to Arcavoid
 *
 * Coeff can be approximated to arbitrary relative precision p
 *   (if exact, say so)
 * TODO: How to deal with unknown zeroes (not for leading coefficient)?
 *       -> use absolute precision for other coeffs
 *****************************************************************/

#include <CGAL/basic.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/LEDA_arithmetic_kernel.h>
#endif // CGAL_USE_LEDA
#ifdef CGAL_USE_GMP
#include <CGAL/GMP_arithmetic_kernel.h>
#endif // CGAL_USE_GMP

#include <CGAL/Algebraic_kernel_1/Cartesian_complex.h>
#include <CGAL/Polynomial_type_generator.h>

#include <CGAL/Algebraic_kernel_1/Double_with_exponent.h>

namespace CGAL {

/*****************************************************************
 * provides:
 * bool operator() (const NT &IN, long P, Bigfloat &OUT);
 *   - set OUT to a bigfloat representation of IN up to relative
 *     precision >= P
 *     (i.e. >= P effective mantissa digits
 *      AND relative error <= 2^(-p))
 *   - return whether the conversion was exact
 *     AND the mantissa does not exceed P effective digits
 *****************************************************************/

template< class NT, class Bigfloat >
class Convert_to_bigfloat_base;

template< class Coeff, class Bigfloat >
struct Convert_polynomial_to_bigfloat_base {
  const bool operator()
  (const typename Polynomial_type_generator< Coeff, 1 >::Type &in,
   long p,
   typename Polynomial_type_generator< Bigfloat, 1>::Type &out) const {
    typedef typename Polynomial_type_generator< Bigfloat, 1>::Type OutPoly;

    const int n = in.degree();

    if (n < 0) {
      out = OutPoly();
      CGAL_postcondition (out.degree() == in.degree());
      return true;
    }

    Convert_to_bigfloat_base< Coeff, Bigfloat > convert;
    std::vector< Bigfloat > cs (n+1, Bigfloat(0));
    bool exact = true;
    for (int i = 0; i <= n; ++i)
      exact &= convert (in[i], p, cs[i]);
    out = OutPoly (cs.begin(), cs.end());
    CGAL_postcondition (out.degree() == in.degree());
    return exact;
  }
};

#ifdef CGAL_USE_LEDA

template<>
struct Convert_to_bigfloat_base< leda_integer, leda_bigfloat > {
  const bool operator() (const leda_integer &n, long p,
                         leda_bigfloat &x) const {
    // leda_integer are reference counted, so just set x to n
    x = n;
    return false;
    // return (p >= n.length()
    //         || p >= x.get_significant_length());
  }
};

template<>
struct Convert_to_bigfloat_base< leda_rational, leda_bigfloat > {
  const bool operator() (const leda_rational &q, long p,
                         leda_bigfloat &x) const {
    bool exact;
    x = leda::div (leda_bigfloat (q.numerator()),
                   leda_bigfloat (q.denominator()),
                   p,
                   leda::TO_NEAREST,
                   exact);
    return exact;
  }
};

template<>
struct Convert_to_bigfloat_base< leda_bigfloat, leda_bigfloat > {
  const bool operator() (const leda_bigfloat &in, long p,
                         leda_bigfloat &x) const {
    x = in;
    return (p >= in.get_significant_length());
  }
};

template<>
struct Convert_to_bigfloat_base< leda_real, leda_bigfloat > {
  const bool operator() (const leda_real &in, long p,
                         leda_bigfloat &x) const {
    in.guarantee_relative_error (p);
    x = in.to_bigfloat();
    return (p >= x.get_significant_length()
            // TODO: is anything cheaper than get_bigfloat_error()?
            && is_zero (in.get_bigfloat_error()));
  }
};

// template< class LEDA_number_type >
// struct Convert_to_bigfloat_base< LEDA_number_type, Double_with_exponent > {
//   const bool operator() (const LEDA_number_type &in, long p,
//                          Double_with_exponent &x) const {
//     leda_bigfloat bf = 0;
//     bool exact = Convert_to_bigfloat_base< LEDA_number_type, leda_bigfloat >() (in, 53, bf);
//     x = Double_with_exponent (bf);
//     return exact;
//   }
// };

#endif // CGAL_USE_LEDA

#ifdef CGAL_USE_GMP
#ifdef CGAL_USE_MPFR

template<>
struct Convert_to_bigfloat_base< Gmpz, Gmpfr > {
  const bool operator() (const Gmpz &n, size_t p,
                         Gmpfr &x) const {
    // Gmpzf mantissae are reference counted, so just set x to n
    x = n;
    return (p >= n.bit_size());
  }
};

template<>
struct Convert_to_bigfloat_base< Gmpq, Gmpfr > {
  const bool operator() (const Gmpq &q, size_t p,
                         Gmpfr &x) const {
    x = Gmpfr::div (Gmpfr (q.numerator()),
                    Gmpfr (q.denominator()),
                    p,
                    std::round_to_nearest);
    return ! Gmpfr::inex_flag();
  }
};

#endif // CGAL_USE_MPFR
#endif // CGAL_USE_GMP

template< class NT >
struct Convert_to_bigfloat_base< NT, Double_with_exponent > {
  const bool operator() (const NT &in, long p,
                         Double_with_exponent &x) const {
    typedef typename CGAL::Get_arithmetic_kernel< NT >::Arithmetic_kernel::Bigfloat BF;
    BF bf = 0;
    bool exact = Convert_to_bigfloat_base< NT, BF >() (in, 53, bf);
    x = Double_with_exponent (bf);
    return exact;
  }
};


template< class Polynomial_,
          class RRFloat_ = typename CGAL::Get_arithmetic_kernel< typename Polynomial_::NT >::Arithmetic_kernel::Bigfloat >
class Arcavoid_kernel {
public:
  typedef Polynomial_ Polynomial;
  typedef typename Polynomial::NT Coeff;
private:
  typedef typename CGAL::Get_arithmetic_kernel< Coeff >::Arithmetic_kernel AK;
  typedef RRFloat_ RRFloat;
  // typedef typename AK::Bigfloat RRFloat;
  typedef Cartesian_complex< RRFloat > CCFloat;
  typedef typename Polynomial_type_generator< RRFloat, 1 >::Type RRPoly;
  typedef typename Polynomial_type_generator< CCFloat, 1 >::Type CCPoly;
public:
  typedef AK Arithmetic_kernel;
  typedef RRFloat Bigfloat;
  typedef CCFloat Complex_bigfloat;
  typedef RRPoly Bigfloat_polynomial;
  typedef CCPoly Complex_bigfloat_polynomial;

  typedef Convert_to_bigfloat_base< Coeff, RRFloat > Convert_to_bigfloat;
  typedef Convert_polynomial_to_bigfloat_base< Coeff, RRFloat > Convert_polynomial_to_bigfloat;

  const Convert_to_bigfloat convert_to_bigfloat_object () const {
    return Convert_to_bigfloat();
  }
  const Convert_polynomial_to_bigfloat convert_polynomial_to_bigfloat_object () const {
    return Convert_polynomial_to_bigfloat();
  }
};

} // namespace CGAL

#endif // CGAL_ARCAVOID_KERNEL_H
