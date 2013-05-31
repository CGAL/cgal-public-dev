#ifndef CGAL_DOUBLE_WITH_EXPONENT_H
#define CGAL_DOUBLE_WITH_EXPONENT_H

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>

#include <CGAL/basic.h>
#include <CGAL/Algebraic_kernel_1/Bigfloat_traits.h>

#ifdef CGAL_USE_GMP
#include <CGAL/GMP_arithmetic_kernel.h>
#endif // CGAL_USE_GMP

#ifdef CGAL_USE_LEDA
#include <CGAL/LEDA_arithmetic_kernel.h>
#endif // CGAL_USE_LEDA

//#include <boost/preprocessor/seq/adt.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/seq/for_each_product.hpp>
#include <boost/preprocessor/stringize.hpp>

#include <boost/numeric/conversion/bounds.hpp>
#include <boost/tuple/tuple.hpp>
//#include <boost/numeric_limits.hpp>

#ifndef CGAL_DOUBLE_WITH_EXPONENT_WRAP_DPE
#define CGAL_DOUBLE_WITH_EXPONENT_WRAP_DPE 1
#endif

#if !CGAL_DOUBLE_WITH_EXPONENT_WRAP_DPE
#error Currently no support for Double_with_exponent without dpe library
#else

#define DPE_USE_LONG true
#include <CGAL/Algebraic_kernel_1/dpe.h>

#define CGAL_DOUBLE_WITH_EXPONENT_BUILTIN_FLOAT_TYPES \
  (float) (double) (long double)
#define CGAL_DOUBLE_WITH_EXPONENT_BUILTIN_SIGNED_INTEGRAL_TYPES \
  (signed char) (signed short) (signed int) (signed long)
#define CGAL_DOUBLE_WITH_EXPONENT_BUILTIN_UNSIGNED_INTEGRAL_TYPES \
  (unsigned char) (unsigned short) (unsigned int) (unsigned long)
#define CGAL_DOUBLE_WITH_EXPONENT_BUILTIN_INTEGRAL_TYPES    \
  CGAL_DOUBLE_WITH_EXPONENT_BUILTIN_SIGNED_INTEGRAL_TYPES   \
  CGAL_DOUBLE_WITH_EXPONENT_BUILTIN_UNSIGNED_INTEGRAL_TYPES
#define CGAL_DOUBLE_WITH_EXPONENT_BUILTIN_TYPES     \
  CGAL_DOUBLE_WITH_EXPONENT_BUILTIN_FLOAT_TYPES     \
  CGAL_DOUBLE_WITH_EXPONENT_BUILTIN_INTEGRAL_TYPES

namespace CGAL {

template< class BF > class Bigfloat_traits;
template< class BF, class T > class Convert_to_bigfloat;

//template< typename MantissaType = DPE_DOUBLE, typename ExponentType = DPE_EXP_T >
class Double_with_exponent {
public:
  // typedef MantissaType Mantissa_type;
  // typedef ExponentType Exponent_type;
  typedef DPE_DOUBLE Mantissa_type;
  typedef DPE_EXP_T Exponent_type;
private:
  mutable dpe_t rep;
  typedef Double_with_exponent Self;

  typedef Mantissa_type MT;     // Mantissa type
  typedef Exponent_type ET;     // Exponent type

public:
  Double_with_exponent () {
    dpe_init (rep);
    dpe_set_d (rep, 0.);
  }
  Double_with_exponent (const Double_with_exponent &rhs) {
    dpe_init (rep);
    dpe_set (rep, rhs.rep);
  }
  Double_with_exponent (dpe_t rhs_rep) {
    dpe_set (rep, rhs_rep);
  }

#define CGAL_DOUBLE_WITH_EXPONENT_BUILTIN_TYPE_CONSTRUCTOR(_,__,MT) \
  Double_with_exponent (MT m) {                                     \
    dpe_init (rep);                                                 \
    DPE_MANT (rep) = static_cast< Mantissa_type > (m);              \
    DPE_EXP (rep) = static_cast< Exponent_type > (0);               \
    dpe_normalize (rep);                                            \
  }

  BOOST_PP_SEQ_FOR_EACH(CGAL_DOUBLE_WITH_EXPONENT_BUILTIN_TYPE_CONSTRUCTOR, ,
                        CGAL_DOUBLE_WITH_EXPONENT_BUILTIN_TYPES);
#undef CGAL_DOUBLE_WITH_EXPONENT_BUILTIN_TYPE_CONSTRUCTOR

#define CGAL_DOUBLE_WITH_BOTH_BUILTIN_TYPES_CONSTRUCTOR(_,MT_ET) \
  Double_with_exponent (BOOST_PP_SEQ_HEAD (MT_ET) m,             \
                        BOOST_PP_SEQ_ELEM (1,MT_ET) e) {         \
    dpe_init (rep);                                              \
    DPE_MANT (rep) = static_cast< Mantissa_type > (m);           \
    DPE_EXP (rep) = static_cast< Exponent_type > (e);            \
    dpe_normalize (rep);                                         \
  }

  BOOST_PP_SEQ_FOR_EACH_PRODUCT(CGAL_DOUBLE_WITH_BOTH_BUILTIN_TYPES_CONSTRUCTOR,
                                (CGAL_DOUBLE_WITH_EXPONENT_BUILTIN_TYPES)
                                (CGAL_DOUBLE_WITH_EXPONENT_BUILTIN_INTEGRAL_TYPES));
#undef CGAL_DOUBLE_WITH_BOTH_BUILTIN_TYPE_CONSTRUCTOR

  static const Self ipow2 (long e) {
    dpe_t x;
    DPE_MANT (x) = 0.5;
    DPE_EXP (x) = e+1;
    return Double_with_exponent (x);
  }

  ~Double_with_exponent () {
    dpe_clear (rep);
  }

  operator long () const {
    return dpe_get_si (rep);
  }
  operator unsigned long () const {
    return dpe_get_ui (rep);
  }
  operator double () const {
    return dpe_get_d (rep);
  }
  operator long double () const {
    return dpe_get_ld (rep);
  }

#ifdef CGAL_USE_LEDA
  Double_with_exponent (const leda_bigfloat &x) {
    dpe_init (rep);
    leda_integer le = x.get_exponent();

    if (DPE_UNLIKELY (! le.is_long()))
      *this = Double_with_exponent (1./0.);
    else {
      long e = le.to_long() + x.get_significant_length();
      DPE_MANT (rep) = (x * leda::ipow2 (-e)).to_double();
      DPE_EXP (rep) = e;
    }
  }

  operator leda_bigfloat () const {
    return to_leda_bigfloat();
  }

  const leda_bigfloat to_leda_bigfloat () const {
    return leda_bigfloat (DPE_MANT (rep))
      * leda::ipow2 (static_cast< long > (DPE_EXP (rep)));
  }
#endif // CGAL_USE_LEDA

#ifdef CGAL_USE_GMP
#ifdef CGAL_USE_MPFR

  Double_with_exponent (const Gmpfr &x) {
    dpe_init (rep);
    boost::tie (DPE_MANT (rep), DPE_EXP (rep)) = x.to_double_exp();

    // from dpe.h: dpe_normalize
    if (DPE_UNLIKELY (DPE_MANT (rep) == 0.0))
      DPE_EXP (rep) = DPE_EXPMIN;

    CGAL_postcondition_code (dpe_t n;               \
                             dpe_init (n);          \
                             dpe_set (n, rep);      \
                             dpe_normalize (n); );
    CGAL_postcondition_msg (DPE_MANT (rep) == DPE_MANT (n) && DPE_EXP (rep) == DPE_EXP (n), \
                            "Check Gmpfr::to_double_exp()");
  }

  operator Gmpfr () const {
    return to_gmpfr();
  }

  const Gmpfr to_gmpfr () const {
    Gmpfr x (DPE_MANT (rep));
    mpfr_set_exp (x.fr(), DPE_EXP (rep));
    return x;
    // return CGAL::Gmpfr (DPE_MANT (rep))
    //   * CGAL::ipow2< CGAL::Gmpfr > (static_cast< long > (DPE_EXP (rep)));
  }

#endif // CGAL_USE_MPFR
#endif // CGAL_USE_GMP

  const Self & operator+() const {
    return *this;
  }
  const Self operator-() const {
    Self result;
    dpe_neg (result.rep, rep);
    return result;
  }

  Self & operator+= (const Self &rhs) {
    dpe_add (rep, rep, rhs.rep);
    return *this;
  }
  Self & operator-= (const Self &rhs) {
    dpe_sub (rep, rep, rhs.rep);
    return *this;
  }
  Self & operator*= (const Self &rhs) {
    dpe_mul (rep, rep, rhs.rep);
    return *this;
  }
  Self & operator/= (const Self &rhs) {
    dpe_div (rep, rep, rhs.rep);
    return *this;
  }

  const Self operator+ (Self rhs) const {
    return rhs += *this;
  }
  const Self operator- (const Self &rhs) const {
    Self result = *this;
    return result -= rhs;
  }
  const Self operator* (Self rhs) const {
    return rhs *= *this;
  }
  const Self operator/ (const Self &rhs) const {
    Self result = *this;
    return result /= rhs;
  }

  Self & operator<<= (unsigned long e) {
    dpe_mul_2exp (rep, rep, e);
    return *this;
  }
  Self & operator>>= (unsigned long e) {
    dpe_div_2exp (rep, rep, e);
    return *this;
  }

  const Self operator<< (unsigned long e) const {
    Self result;
    result <<= e;
    return result;
  }
  const Self operator>> (unsigned long e) const {
    Self result;
    result >>= e;
    return result;
  }

  const bool operator<= (const Self &rhs) const {
    return dpe_cmp (rep, rhs.rep) <= 0;
  }
  const bool operator< (const Self &rhs) const {
    return dpe_cmp (rep, rhs.rep) == -1;
  }
  const bool operator== (const Self &rhs) const {
    return dpe_cmp (rep, rhs.rep) == 0;
  }
  const bool operator!= (const Self &rhs) const {
    return dpe_cmp (rep, rhs.rep) != 0;
  }
  const bool operator>= (const Self &rhs) const {
    return dpe_cmp (rep, rhs.rep) >= 0;
  }
  const bool operator> (const Self &rhs) const {
    return dpe_cmp (rep, rhs.rep) == 1;
  }
  const int sign () const {
    return DPE_SIGN (rep);
  }

  dpe_t & dpe() {
    return rep;
  }
  const dpe_t & dpe() const {
    return rep;
  }

  const Self abs() const {
    Self result;
    dpe_abs (result.rep, rep);
    return result;
  }
  const Self sqrt() const {
    Self result;
    dpe_sqrt (result.rep, rep);
    return result;
  }

  /****************************************************************
   * New functionality
   ****************************************************************/
  const Exponent_type ilog2() const {
    return DPE_EXP (rep);
  }
  const Self root_d (int d) const {
    assert (d > 0);
    std::ldiv_t qr = std::ldiv (DPE_EXP (rep), d);
    if (DPE_UNLIKELY (qr.rem == 0))
      return Self (std::pow (DPE_MANT (rep), 1. / static_cast< MT > (d)),
                   qr.quot);
    else
      return Self (std::pow (std::pow (2., static_cast< double > (qr.rem)) * DPE_MANT (rep),
                             1. / static_cast< MT > (d)),
                   qr.quot);
  }
};

inline
const Double_with_exponent & operator>> (std::istream &in, Double_with_exponent & dwe) {
  Double_with_exponent::Mantissa_type d;
  in >> d;
  dwe = d;
  return dwe;
}

inline
std::ostream & operator<< (std::ostream &out, const Double_with_exponent & dwe) {
  return out << static_cast< Double_with_exponent::Mantissa_type > (dwe);
  // return out << DPE_MANT (dwe.dpe()) << "b" << DPE_EXP (dwe.dpe());
}

template <>
class Algebraic_structure_traits< Double_with_exponent >
  : public Algebraic_structure_traits_base< Double_with_exponent,
                                            Field_with_kth_root_tag >  {
public:
  typedef Tag_false            Is_exact;
  typedef Tag_true             Is_numerical_sensitive;

  class Sqrt
    : public std::unary_function< Type, Type > {
  public:
    Type operator()( const Type& x ) const {
      return x.sqrt();
    }
  };

  class Kth_root
    : public std::binary_function<int, Type, Type> {
  public:
    Type operator()( int k,
                     const Type& x) const {
          CGAL_precondition_msg( k > 0, "'k' must be positive for k-th roots");
          return x.root_d (k);
    }
  };
};

template <>
class Real_embeddable_traits< Double_with_exponent >
  : public INTERN_RET::Real_embeddable_traits_base< Double_with_exponent , CGAL::Tag_true> {
public:
  class Abs
    : public std::unary_function< Type, Type > {
  public:
    Type operator()( const Type& x ) const {
      return x.abs();
    }
  };
};

template<>
class Bigfloat_traits< Double_with_exponent > {
public:
  typedef Double_with_exponent Type;

  struct Integral_pow2
    : public std::unary_function< long, Double_with_exponent > {
    const Double_with_exponent operator() (long e) const {
      return Double_with_exponent::ipow2 (e);
    }
  };

  // TODO
  struct Set_precision;
  //    : public std::binary_function< leda_bigfloat, long, void >;

  struct Set_default_precision
    : public std::unary_function< long, void > {
    void operator() (long) const {}
  };

  struct Multiply_by_pow2
    : public std::binary_function< Double_with_exponent, long, void > {
    void operator() (Double_with_exponent &x, long e) const {
      x *= Integral_pow2() (e);
    }
  };

  struct Abs_integral_log2
    : public std::unary_function< Double_with_exponent, long > {
    const long operator() (const Double_with_exponent &x) const {
      return x.ilog2();
    }
  };

  struct Abs_round_up_to_pow2
    : public std::unary_function< Double_with_exponent, Double_with_exponent > {
    const Double_with_exponent operator() (const Double_with_exponent &x) const {
      return Integral_pow2() (Abs_integral_log2() (x));
    }
  };

  struct Root_d
    : public std::binary_function< Double_with_exponent, int, Double_with_exponent > {
    const Double_with_exponent operator() (const Double_with_exponent &x, int d, long) const {
      return x.root_d (d);
    }
    const Double_with_exponent operator() (const Double_with_exponent &x, int d) const {
      return x.root_d (d);
    }
  };
};

#ifdef CGAL_USE_LEDA
template<>
struct Convert_to_bigfloat< Double_with_exponent, leda_bigfloat >
  : public std::binary_function< leda_bigfloat, long, Double_with_exponent > {
  const Double_with_exponent operator() (const leda_bigfloat &x, long) {
    return Double_with_exponent (x);
  }
};

template<>
struct Convert_to_bigfloat< Double_with_exponent, leda_integer >
  : public std::binary_function< leda_integer, long, Double_with_exponent > {
  const Double_with_exponent operator() (const leda_integer &n, long p) {
    return Double_with_exponent (leda_bigfloat (n));
  }
};

template<>
struct Convert_to_bigfloat< Double_with_exponent, leda_rational >
  : public std::binary_function< leda_rational, long, Double_with_exponent > {
  const Double_with_exponent operator() (const leda_rational &q, long) {
    return Double_with_exponent (leda_bigfloat (q.numerator()) / leda_bigfloat (q.denominator()));
  }
};
#endif // CGAL_USE_LEDA

} // namespace CGAL

#endif // CGAL_DOUBLE_WITH_EXPONENT_WRAP_DPE

#endif // CGAL_DOUBLE_WITH_EXPONENT_H
