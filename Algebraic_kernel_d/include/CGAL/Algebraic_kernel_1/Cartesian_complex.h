/*
** Complex.h
**
** Made by Alexander Kobel
** Login   <perpeduumimmobile@lidinoid>
**
** Started on  Sun Jan  3 13:36:50 2010 Alexander Kobel
** Last update Sun Jan  3 13:36:50 2010 Alexander Kobel
*/

/****************************************************************
 * TODO:
 * - istream & operator>> (...) does not work as expected for leda_bigfloat
 * - implement missing functors in AST (see below)
 * - implement multiply_by_I
 * - implement abs to user-given precision.
 *   -> Interface?
 * - implement arg() (Which return type? Template? To user-given precision?)
 * - constructor from (radius, phase angle) pair?
 *   -> how to distinguish from Cartesian constructor?
 *   -> as a global/static function (polar() in std::complex)?
 * - should Cartesian_complex have additional template params
 *   for Modulus_type and Azimuth_type?
 * - check if inheritation of Algebraic_category from T is correct
 *   in each situation
 * - should global operators and AST/CET be friends and use x and y
 *   instead of real() and imag()?
 *
 * - should we offer a Homogenous_cartesian_complex type?
 ****************************************************************/

#ifndef   	CGAL_CARTESIAN_COMPLEX_H_
# define   	CGAL_CARTESIAN_COMPLEX_H_

#include <CGAL/basic.h>
#include <complex>
#include <CGAL/Algebraic_kernel_1/Complex_embeddable_traits.h>
#include <CGAL/Interval_traits.h>

namespace CGAL {

//#define CGAL_CARTESIAN_COMPLEX_USE_GAUSS_MULTIPLICATION

template< typename T >
class Cartesian_complex {
public:
  typedef Represents_cartesian_complex_tag      Representation;
  typedef T                                     Cartesian_type;
  // for interchangeability with std::complex< T >
  typedef T                                     value_type;

private:
  typedef Cartesian_complex< T >                Self;
  T x,y;

#define CST_TP typename CGAL::Coercion_traits< T, T2 >::Type
#define CSTR typename CGAL::Coercion_traits< T, T2 >::Cast()

public:
  Cartesian_complex () : x(0), y(0) {}

  Cartesian_complex (const T &x)
    : x(x), y(0) {}
  Cartesian_complex (const T &x, const T& y)
    : x(x), y(y) {}
  Cartesian_complex (const Self &rhs)
    : x(rhs.real()), y(rhs.imag()) {}
  Cartesian_complex (const ::std::complex< T > &rhs)
    : x(rhs.real()), y(rhs.imag()) {}
  Cartesian_complex (const ::std::pair< T, T > &rhs)
    : x(rhs.first), y(rhs.second) {}

  template< typename T2 >
  Cartesian_complex (const T2 &x)
    : x(CSTR (x)), y(0) {
    // BOOST_STATIC_ASSERT ((boost::is_same< CST_TP, T >::value));
  }
  template< typename T2 >
  Cartesian_complex (const T2 &x, const T2& y)
    : x(CSTR (x)), y(CSTR (y)) {
    // BOOST_STATIC_ASSERT ((boost::is_same< CST_TP, T >::value));
  }
  template< typename T2 >
  Cartesian_complex (const Cartesian_complex< T2 > &rhs)
    : x(CSTR (rhs.real())), y(CSTR (rhs.imag())) {
    // BOOST_STATIC_ASSERT ((boost::is_same< CST_TP, T >::value));
  }
  template< typename T2 >
  Cartesian_complex (const ::std::complex< T2 > &rhs)
    : x(CSTR (rhs.real())), y(CSTR (rhs.imag())) {
    // BOOST_STATIC_ASSERT ((boost::is_same< CST_TP, T >::value));
  }
  template< typename T2 >
  Cartesian_complex (const ::std::pair< T2, T2 > &rhs)
    : x(CSTR (rhs.first)), y(CSTR (rhs.second)) {
    // BOOST_STATIC_ASSERT ((boost::is_same< CST_TP, T >::value));
  }

  template< typename T2 >
  operator const Cartesian_complex< T2 > () {
    // BOOST_STATIC_ASSERT ((boost::is_same< CST_TP, T2 >::value));
    return Cartesian_complex< T2 > (CSTR (x), CSTR (y));
  }

  operator const ::std::complex< T > () {
    return ::std::complex< T > (x, y);
  }

  const T & real() const { return x; }
  T & real() { return x; }

  const T & imag() const { return y; }
  T & imag() { return y; }

  Self & operator= (const T &rhs) {
    x = rhs;
    y = 0;
    return *this;
  }
  Self & operator= (const Self& rhs) {
    x = rhs.real();
    y = rhs.imag();
    return *this;
  }
  Self & operator= (const std::complex< T >& rhs) {
    x = rhs.real();
    y = rhs.imag();
    return *this;
  }

  template< typename T2 >
  Self & operator= (const T2 &rhs) {
    BOOST_STATIC_ASSERT ((boost::is_same< CST_TP, T >::value));
    x = CSTR (rhs);
    y = 0;
    return *this;
  }
  template< typename T2 >
  Self & operator= (const Cartesian_complex< T2 >& rhs) {
    BOOST_STATIC_ASSERT ((boost::is_same< CST_TP, T >::value));
    x = CSTR (rhs.real());
    y = CSTR (rhs.imag());
    return *this;
  }
  template< typename T2 >
  Self & operator= (const std::complex< T2 >& rhs) {
    BOOST_STATIC_ASSERT ((boost::is_same< CST_TP, T >::value));
    x = CSTR (rhs.real());
    y = CSTR (rhs.imag());
    return *this;
  }

  Self & operator+= (const T &rhs) {
    x += rhs;
    return *this;
  }
  Self & operator+= (const Self &rhs) {
    x += rhs.real();
    y += rhs.imag();
    return *this;
  }
  template< typename T2 >
  Self & operator+= (const T2 &rhs) {
    BOOST_STATIC_ASSERT ((boost::is_same< CST_TP, T >::value));
    x += CSTR (rhs);
    return *this;
  }
  template< typename T2 >
  Self & operator+= (const Cartesian_complex< T2 > &rhs) {
    BOOST_STATIC_ASSERT ((boost::is_same< CST_TP, T >::value));
    x += CSTR (rhs.real());
    y += CSTR (rhs.imag());
    return *this;
  }

  Self & operator-= (const T &rhs) {
    x -= rhs;
    return *this;
  }
  Self & operator-= (const Self &rhs) {
    x -= rhs.real();
    y -= rhs.imag();
    return *this;
  }
  template< typename T2 >
  Self & operator-= (const T2 &rhs) {
    BOOST_STATIC_ASSERT ((boost::is_same< CST_TP, T >::value));
    x -= CSTR (rhs);
    return *this;
  }
  template< typename T2 >
  Self & operator-= (const Cartesian_complex< T2 > &rhs) {
    BOOST_STATIC_ASSERT ((boost::is_same< CST_TP, T >::value));
    x -= CSTR (rhs.real());
    y -= CSTR (rhs.imag());
    return *this;
  }

  Self & operator*= (const T &rhs) {
    x *= rhs;
    y *= rhs;
    return *this;
  }
  template< typename T2 >
  Self & operator*= (const T2 &rhs) {
    BOOST_STATIC_ASSERT ((boost::is_same< CST_TP, T >::value));
    x *= CSTR (rhs);
    y *= CSTR (rhs);
    return *this;
  }
  Self & operator*= (const Self &rhs) {
    /* Write
     *   lhs = a + b*I
     *   rhs = c + d*I.
     * Then,
     *   Re (lhs * rhs) = ac - bd
     *   Im (lhs * rhs) = ad + bc = (a + b)*(c + d) - ac - bd
     */

#ifdef CGAL_CARTESIAN_COMPLEX_USE_GAUSS_MULTIPLICATION
    const T ac = x * rhs.real();
    const T bd = y * rhs.imag();
    const T c_plus_d = rhs.real() + rhs.imag();

    // compute Im
    y += x;
    y *= c_plus_d;
    y -= ac;
    y -= bd;

    // compute Re
    x = ac - bd;
#else // CGAL_CARTESIAN_COMPLEX_USE_GAUSS_MULTIPLICATION
    const T ac = x * rhs.real();
    const T bd = y * rhs.imag();
    const T ad = x * rhs.imag();
    const T bc = y * rhs.real();

    x = ac - bd;
    y = ad + bc;
#endif // CGAL_CARTESIAN_COMPLEX_USE_GAUSS_MULTIPLICATION

    return *this;
  }
  template< typename T2 >
  Self & operator*= (const Cartesian_complex< T2 > &rhs) {
    return operator*= (Self (rhs));
  }

  const Self square() const {
    const T aa = CGAL::square (x);
    const T bb = CGAL::square (y);

    return Self (aa - bb, 2*x*y);
  }

  const Self reciprocal() const {
    // returns 1/z = 1/(x + iy) = (x - iy) / (x^2 + y^2)
    const T sqr_norm = squared_norm_2();
    return Self (x / sqr_norm, -y / sqr_norm);
  }

  Self & operator/= (const T &rhs) {
    x /= rhs;
    y /= rhs;
    return *this;
  }
  template< typename T2 >
  Self & operator/= (const T2 &rhs) {
    BOOST_STATIC_ASSERT ((boost::is_same< CST_TP, T >::value));
    x /= CSTR (rhs);
    y /= CSTR (rhs);
    return *this;
  }
  Self & operator/= (const Self &rhs) {
    /* Write
     *   lhs = a + b*I
     *   rhs = c + d*I.
     * Then,
     *   Re (lhs / rhs) * (c2+d2) = ac + bd
     *   Im (lhs / rhs) * (c2+d2) = bc - ad = (a + b)*(c - d) - ac + bd
     */

    const T ac = x * rhs.real();
    const T bd = y * rhs.imag();
    const T c_minus_d = rhs.real() - rhs.imag();

    // compute Im * (c2+d2)
    y += x;
    y *= c_minus_d;
    y -= ac;
    y += bd;

    // compute Re * (c2+d2)
    x = ac + bd;

    // compute (c2+d2)
    const T den = CGAL::square (rhs.real()) + CGAL::square (rhs.imag());

    x /= den;
    y /= den;

    return *this;
  }
  template< typename T2 >
  Self & operator/= (const Cartesian_complex< T2 > &rhs) {
    return operator/= (Self (rhs));
  }

  const Self operator+() const { return Self ( x,  y); }
  const Self operator~() const { return Self ( x, -y); }
  const Self operator-() const { return Self (-x, -y); }

  // for interchangeability with std::complex
  const Self conj() const { return operator~(); }
  const T norm() const { return squared_norm_2(); }
  const T abs() const { return CGAL::sqrt (norm()); }

  const T norm_1() const { return CGAL::abs (x) + CGAL::abs (y); }
  const T squared_norm_2() const { return CGAL::square (x) + CGAL::square (y); }
  const T norm_inf() const { return CGAL::max (CGAL::abs (x), CGAL::abs (y)); }

  void multiply_by_I() {
    std::swap (x, y);
    x = -x;
  }
  void multiply_by_minus_I() {
    std::swap (x,y);
    y = -y;
  }

  static const Self ZERO()      { return Self ( T(0),  T(0)); }
  static const Self ONE()       { return Self ( T(1),  T(0)); }
  static const Self I()         { return Self ( T(0),  T(1)); }
  static const Self MINUS_ONE() { return Self (-T(1),  T(0)); }
  static const Self MINUS_I()   { return Self ( T(0), -T(1)); }
};

template< typename T >
const T abs (const Cartesian_complex< T > &z) {
  return z.abs();
}

#define CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_RETURN_TYPE          \
  template< typename T1, typename T2 >                                \
  inline                                                              \
  const Cartesian_complex< typename Coercion_traits< T1, T2 >::Type >
#define CGAL_CARTESIAN_COMPLEX_OPERATOR_BOOL_RETURN_TYPE \
  template< typename T1, typename T2 >                   \
  inline                                                 \
  const bool
#define CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_TYPEDEFS \
  typedef Coercion_traits< T1, T2 >    CT;                \
  typedef typename CT::Cast            CT_cast;           \
  typedef typename CT::Type            CT_type;           \
  typedef Cartesian_complex< CT_type > Cartesian
#define CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_AB \
  const CT_type a = CT_cast() (lhs.real());         \
  const CT_type b = CT_cast() (lhs.imag())
#define CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_CD \
  const CT_type c = CT_cast() (rhs.real());         \
  const CT_type d = CT_cast() (rhs.imag())
#define CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_ABCD \
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_AB;        \
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_CD

/****************************************************************
 * OPERATOR (Complex< T1 >, Complex< T2 >)
 ****************************************************************/

#define CGAL_CARTESIAN_COMPLEX_OPERATOR_ARGUMENTS \
  const Cartesian_complex< T1 > &lhs,             \
  const Cartesian_complex< T2 > &rhs

CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_RETURN_TYPE
operator+ (CGAL_CARTESIAN_COMPLEX_OPERATOR_ARGUMENTS) {
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_TYPEDEFS;
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_ABCD;
  return Cartesian (a + c, b + d);
}

CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_RETURN_TYPE
operator- (CGAL_CARTESIAN_COMPLEX_OPERATOR_ARGUMENTS) {
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_TYPEDEFS;
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_ABCD;
  return Cartesian (a - c, b - d);
}

CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_RETURN_TYPE
operator* (CGAL_CARTESIAN_COMPLEX_OPERATOR_ARGUMENTS) {
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_TYPEDEFS;
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_ABCD;

#ifdef CGAL_CARTESIAN_COMPLEX_USE_GAUSS_MULTIPLICATION
  /* Write
   *   lhs = a + b*I
   *   rhs = c + d*I.
   * Then,
   *   Re (lhs * rhs) = ac - bd
   *   Im (lhs * rhs) = ad + bc = (a + b)*(c + d) - ac - bd
   */

  const CT_type ac = a * c;
  const CT_type bd = b * d;
  const CT_type c_plus_d = c + d;

  // compute Im
  CT_type y = a + b;
  y *= c_plus_d;
  y -= ac;

  // compute Re, finalize computation of Im, and return
  return Cartesian (ac - bd, y - bd);
#else // CGAL_CARTESIAN_COMPLEX_USE_GAUSS_MULTIPLICATION
  const CT_type ac = a * c;
  const CT_type bd = b * d;
  const CT_type ad = a * d;
  const CT_type bc = b * c;

  return Cartesian (ac - bd, ad + bc);
#endif // CGAL_CARTESIAN_COMPLEX_USE_GAUSS_MULTIPLICATION
}

CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_RETURN_TYPE
operator/ (CGAL_CARTESIAN_COMPLEX_OPERATOR_ARGUMENTS) {
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_TYPEDEFS;
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_ABCD;

  /* Write
   *   lhs = a + b*I
   *   rhs = c + d*I.
   * Then,
   *   Re (lhs / rhs) * (c2+d2) = ac + bd
   *   Im (lhs / rhs) * (c2+d2) = bc - ad = (a + b)*(c - d) - ac + bd
   */

  const CT_type ac = a * c;
  const CT_type bd = b * d;
  const CT_type c_minus_d = c - d;

  // compute Re * (c2+d2)
  const CT_type x = ac + bd;

  // compute Im * (c2+d2)
  CT_type y = a + b;
  y *= c_minus_d;
  y -= ac;
  y += bd;

  // compute (c2+d2)
  const CT_type den = CGAL::square (c) + CGAL::square (d);

  return Cartesian (x / den, y / den);
}

CGAL_CARTESIAN_COMPLEX_OPERATOR_BOOL_RETURN_TYPE
operator== (CGAL_CARTESIAN_COMPLEX_OPERATOR_ARGUMENTS) {
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_TYPEDEFS;
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_ABCD;
  return CGAL::certainly (a == c) && CGAL::certainly (b == d);
}

CGAL_CARTESIAN_COMPLEX_OPERATOR_BOOL_RETURN_TYPE
operator!= (CGAL_CARTESIAN_COMPLEX_OPERATOR_ARGUMENTS) {
  return ! (lhs == rhs);
  // CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_TYPEDEFS;
  // CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_ABCD;
  // return CGAL::possibly (a != c) || CGAL::possibly (b != d);
}

/****************************************************************
 * OPERATOR (T1, Complex< T2 >)
 ****************************************************************/

#undef CGAL_CARTESIAN_COMPLEX_OPERATOR_ARGUMENTS
#define CGAL_CARTESIAN_COMPLEX_OPERATOR_ARGUMENTS \
  const T1 &lhs, \
  const Cartesian_complex< T2 > &rhs

CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_RETURN_TYPE
operator+ (CGAL_CARTESIAN_COMPLEX_OPERATOR_ARGUMENTS) {
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_TYPEDEFS;
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_CD;
  return Cartesian (CT_cast() (lhs) + c, d);
}

CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_RETURN_TYPE
operator- (CGAL_CARTESIAN_COMPLEX_OPERATOR_ARGUMENTS) {
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_TYPEDEFS;
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_CD;
  return Cartesian (CT_cast() (lhs) - c, -d);
}

CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_RETURN_TYPE
operator* (CGAL_CARTESIAN_COMPLEX_OPERATOR_ARGUMENTS) {
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_TYPEDEFS;
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_CD;
  const CT_type a = CT_cast() (lhs);
  return Cartesian (a * c, a * d);
}

CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_RETURN_TYPE
operator/ (CGAL_CARTESIAN_COMPLEX_OPERATOR_ARGUMENTS) {
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_TYPEDEFS;
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_CD;
  const CT_type a = CT_cast() (lhs);

  const CT_type ac = a * c;
  const CT_type minus_ad = - (a * d);
  const CT_type den = CGAL::square (c) + CGAL::square (d);

  return Cartesian (ac / den, minus_ad / den);
}

CGAL_CARTESIAN_COMPLEX_OPERATOR_BOOL_RETURN_TYPE
operator== (CGAL_CARTESIAN_COMPLEX_OPERATOR_ARGUMENTS) {
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_TYPEDEFS;
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_CD;
  return CGAL::is_zero (d) && (CT_cast() (lhs) == c);
}

CGAL_CARTESIAN_COMPLEX_OPERATOR_BOOL_RETURN_TYPE
operator!= (CGAL_CARTESIAN_COMPLEX_OPERATOR_ARGUMENTS) {
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_TYPEDEFS;
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_CD;
  return (! CGAL::is_zero (d)) || (CT_cast() (lhs) != c);
}

/****************************************************************
 * OPERATOR (Complex< T1 >, T2)
 ****************************************************************/

#undef CGAL_CARTESIAN_COMPLEX_OPERATOR_ARGUMENTS
#define CGAL_CARTESIAN_COMPLEX_OPERATOR_ARGUMENTS \
  const Cartesian_complex< T1 > &lhs,             \
  const T2 &rhs

CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_RETURN_TYPE
operator+ (CGAL_CARTESIAN_COMPLEX_OPERATOR_ARGUMENTS) {
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_TYPEDEFS;
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_AB;
  return Cartesian (a + CT_cast() (rhs), b);
}

CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_RETURN_TYPE
operator- (CGAL_CARTESIAN_COMPLEX_OPERATOR_ARGUMENTS) {
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_TYPEDEFS;
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_AB;
  return Cartesian (a - CT_cast() (rhs), b);
}

CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_RETURN_TYPE
operator* (CGAL_CARTESIAN_COMPLEX_OPERATOR_ARGUMENTS) {
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_TYPEDEFS;
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_AB;
  const CT_type c = CT_cast() (rhs);
  return Cartesian (a * c, b * c);
}

CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_RETURN_TYPE
operator/ (CGAL_CARTESIAN_COMPLEX_OPERATOR_ARGUMENTS) {
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_TYPEDEFS;
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_AB;
  const CT_type c = CT_cast() (rhs);
  return Cartesian (a / c, b / c);
}

CGAL_CARTESIAN_COMPLEX_OPERATOR_BOOL_RETURN_TYPE
operator== (CGAL_CARTESIAN_COMPLEX_OPERATOR_ARGUMENTS) {
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_TYPEDEFS;
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_AB;
  return CGAL::is_zero (b) && (a == CT_cast() (rhs));
}

CGAL_CARTESIAN_COMPLEX_OPERATOR_BOOL_RETURN_TYPE
operator!= (CGAL_CARTESIAN_COMPLEX_OPERATOR_ARGUMENTS) {
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_TYPEDEFS;
  CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_AB;
  return (! CGAL::is_zero (b)) || (a != CT_cast() (rhs));
}

#undef CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_RETURN_TYPE
#undef CGAL_CARTESIAN_COMPLEX_OPERATOR_BOOL_RETURN_TYPE
#undef CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_TYPEDEFS
#undef CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_AB
#undef CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_CD
#undef CGAL_CARTESIAN_COMPLEX_OPERATOR_COERCION_ABCD
#undef CGAL_CARTESIAN_COMPLEX_OPERATOR_ARGUMENTS

template< class T >
inline
std::ostream & operator<< (std::ostream &os,
                           const Cartesian_complex< T > &z) {
  return os << '(' << z.real() << ',' << z.imag() << ')';
}

template< class T >
inline
std::istream & operator>> (std::istream &is,
                           Cartesian_complex< T > &z) {
  std::complex< T > intermediate;
  is >> intermediate;

  std::cerr << intermediate << std::endl;

  if (! is.fail())
    z = intermediate;

  return is;
}

CGAL_DEFINE_COERCION_TRAITS_FOR_SELF_TEM(Cartesian_complex< T >, typename T)

template< typename T1, typename T2 >
struct Coercion_traits< Cartesian_complex< T1 >, T2 > {
private:
  typedef Coercion_traits< T1, T2 > Base_CT;
  typedef typename Base_CT::Type    T;

  typedef Cartesian_complex< T1 > Complex1;
  typedef Cartesian_complex< T2 > Complex2;

public:
  typedef typename Base_CT::Are_implicit_interoperable Are_implicit_interoperable;
  typedef typename Base_CT::Are_explicit_interoperable Are_explicit_interoperable;
  typedef Cartesian_complex< T > Type;
  struct Cast {
    typedef Type result_type;
    const Type operator() (const Complex1 &z) const {
      return static_cast< const Type >(z);
    }
    const Type operator() (const T2 &z) const {
      return static_cast< const Type >(Complex2 (z));
    }
  };
};

template< typename T1, typename T2 >
struct Coercion_traits< T1, Cartesian_complex< T2 > > {
private:
  typedef Coercion_traits< T1, T2 > Base_CT;
  typedef typename Base_CT::Type    T;

  typedef Cartesian_complex< T1 > Complex1;
  typedef Cartesian_complex< T2 > Complex2;

public:
  typedef typename Base_CT::Are_implicit_interoperable Are_implicit_interoperable;
  typedef typename Base_CT::Are_explicit_interoperable Are_explicit_interoperable;
  typedef Cartesian_complex< T > Type;
  struct Cast {
    typedef Type result_type;
    const Type operator() (const T1 &z) const {
      return static_cast< const Type >(Complex1 (z));
    }
    const Type operator() (const Complex2 &z) const {
      return static_cast< const Type >(z);
    }
  };
};

template< typename T1, typename T2 >
struct Coercion_traits< Cartesian_complex< T1 >, Cartesian_complex< T2 > > {
private:
  typedef Coercion_traits< T1, T2 > Base_CT;
  typedef typename Base_CT::Type    T;

  typedef Cartesian_complex< T1 > Complex1;
  typedef Cartesian_complex< T2 > Complex2;

public:
  typedef typename Base_CT::Are_implicit_interoperable Are_implicit_interoperable;
  typedef typename Base_CT::Are_explicit_interoperable Are_explicit_interoperable;
  typedef Cartesian_complex< T > Type;
  struct Cast {
    typedef Type result_type;
    const Type operator() (const Complex1 &z) const {
      return static_cast< const Type >(z);
    }
    const Type operator() (const Complex2 &z) const {
      return static_cast< const Type >(z);
    }
  };
};

template< typename T >
class Complex_embeddable_traits< Cartesian_complex< T > >
  : public INTERN_CET::Complex_embeddable_traits_base< Cartesian_complex< T >,
                                                       Tag_true,
                                                       Represents_cartesian_complex_tag > {};

template< typename T >
class Algebraic_structure_traits< Cartesian_complex< T > > {
public:
  typedef Cartesian_complex< T >                 Type;

private:
  typedef Algebraic_structure_traits< T >        AST_T;
  typedef Complex_embeddable_traits< Type >      CET;

public:
  typedef typename AST_T::Algebraic_category     Algebraic_category;
  typedef typename AST_T::Is_exact               Is_exact;
  typedef typename AST_T::Is_numerical_sensitive Is_numerical_sensitive;

  typedef typename CET::Is_zero                  Is_zero;
  struct Is_one : public std::unary_function< Type, bool > {
    const bool operator() (const Type &z) const
    { return typename AST_T::Is_one() (z.real()) && typename AST_T::Is_zero() (z.imag()); }
  };
  struct Square : public std::unary_function< Type, Type > {
    const Type operator() (const Type &z) const
    { return z.square(); }
  };
  struct Simplify : public std::unary_function< Type, void > {
    void operator() (Type &z) const {
      typename AST_T::Simplify() (z.real());
      typename AST_T::Simplify() (z.imag());
    }
  };

  // struct Unit_part
  //   : public std::unary_function< Type,
  //                                 typename Coercion_traits< Type,
  //                                                           Cartesian_complex< typename Complex_embeddable_traits< Type >::Modulus_type > >
  //                                 ::Type > {
  //   const typename Coercion_traits< Type,
  //                                   Cartesian_complex< typename Complex_embeddable_traits< Type >::Modulus_type > >::Type
  //   operator() (const Type &z) {
  //     return z / z.abs();
  //   }
  // };
  typedef Null_functor Unit_part;

  struct Integral_division
    : public std::binary_function< Type, Type, Type > {
    const Type operator() (const Type &lhs, const Type &rhs) const
    { return lhs / rhs; }
  };

  typedef Null_functor Divides;
  typedef Null_functor Is_square;
  typedef Null_functor Gcd;
  typedef Null_functor Mod;
  typedef Null_functor Div;
  typedef Null_functor Div_mod;

  struct Sqrt
    : public std::unary_function< Type, Type > {
    const Type operator() (const Type &z) const {
      const T abs = z.abs();

      T sqr_x = abs + z.real();
      T sqr_y_den = sqr_x;
      sqr_x /= T(2);
      sqr_y_den *= T(2);
      const T y_den = CGAL::sqrt (sqr_y_den);

      return Type (CGAL::sqrt (sqr_x), z.imag() / y_den);
    }
  };

  typedef Null_functor Kth_root;
  typedef Null_functor Root_of;
};

template< class Real_embeddable >
struct To_double< Cartesian_complex< Real_embeddable > >
  : public std::unary_function< Cartesian_complex< Real_embeddable >,
                                typename Complex_embeddable_traits< Cartesian_complex< Real_embeddable > >::To_double::result_type > {
  typename Complex_embeddable_traits< Cartesian_complex< Real_embeddable > >::To_double::result_type
  operator() (const Cartesian_complex< Real_embeddable > &z) const {
    typename Complex_embeddable_traits< Cartesian_complex< Real_embeddable > >::To_double to_double;
    return to_double (z);
  }
};

template< class Real_embeddable >
struct To_interval< Cartesian_complex< Real_embeddable > >
  : public std::unary_function< Cartesian_complex< Real_embeddable >,
                                typename Complex_embeddable_traits< Cartesian_complex< Real_embeddable > >::To_interval::result_type > {
  typename Complex_embeddable_traits< Cartesian_complex< Real_embeddable > >::To_interval::result_type
  operator() (const Cartesian_complex< Real_embeddable > &z) const {
    typename Complex_embeddable_traits< Cartesian_complex< Real_embeddable > >::To_interval to_interval;
    return to_interval (z);
  }
};

template< class Real_embeddable >
inline
typename Complex_embeddable_traits< Cartesian_complex< Real_embeddable > >::To_double::result_type
// std::complex< double >
to_double (const Cartesian_complex< Real_embeddable > &z) {
  typename Complex_embeddable_traits< Cartesian_complex< Real_embeddable > >::To_double to_double;
  return to_double (z);
}

template< class Real_embeddable >
inline
typename Complex_embeddable_traits< Cartesian_complex< Real_embeddable > >::To_interval::result_type
//std::complex< Interval_nt< true > >
to_interval (const Cartesian_complex< Real_embeddable > &z) {
  typename Complex_embeddable_traits< Cartesian_complex< Real_embeddable > >::To_interval to_interval;
  return to_interval (z);
}

namespace INTERN_CET {

template< class T, class Is_interval_, class Representation_ >
class Interval_traits_base {
public:
  typedef Cartesian_complex< T > Type;
  typedef CGAL::Null_tag Is_interval;
  typedef CGAL::Null_tag With_empty_interval;

  typedef CGAL::Null_functor Lower;
  typedef CGAL::Null_functor Upper;
  typedef CGAL::Null_functor Width;
  typedef CGAL::Null_functor Median;
  typedef CGAL::Null_functor Norm;
  typedef CGAL::Null_functor Empty;
  typedef CGAL::Null_functor Singleton;
  typedef CGAL::Null_functor In;
  typedef CGAL::Null_functor Zero_in;
  typedef CGAL::Null_functor Equal;
  typedef CGAL::Null_functor Overlap;
  typedef CGAL::Null_functor Subset;
  typedef CGAL::Null_functor Proper_Subset;
  typedef CGAL::Null_functor Intersection;
  typedef CGAL::Null_functor Hull;
};

template< class T >
class Interval_traits_base< Cartesian_complex< T >, CGAL::Tag_true, Represents_cartesian_complex_tag > {
public:
  typedef Cartesian_complex< T >                Type;

private:
  typedef Interval_traits< T >                  IT_T;
  typedef typename Interval_traits< T >::Bound  T_Bound;
  typedef Complex_embeddable_traits< Type >     CET;

  typedef typename IT_T::Median                 Median_T;
  typedef typename IT_T::Empty                  Empty_T;
  typedef typename IT_T::Singleton              Singleton_T;
  typedef typename IT_T::Zero_in                Zero_in_T;

public:
  typedef typename IT_T::Is_interval            Is_interval;
  typedef typename IT_T::With_empty_interval    With_empty_interval;
  typedef Cartesian_complex< T_Bound >          Bound;

  typedef CGAL::Null_functor Lower;
  typedef CGAL::Null_functor Upper;
  typedef CGAL::Null_functor Width;

  struct Median : public std::unary_function< Type, Bound > {
    const Bound operator() (const Type &I) const
    { return Bound (Median_T() (I.real()), Median_T() (I.imag())); }
  };

  typedef CGAL::Null_functor Norm;

  struct Empty : public std::unary_function< Type, bool > {
    const bool operator() (const Type &I) const
    { return (Empty_T() (I.real()) || Empty_T() (I.imag())); }
  };

  struct Singleton : public std::unary_function< Type, bool > {
    const bool operator() (const Type &I) const
    { return (Singleton_T() (I.real()) && Singleton_T() (I.imag())); }
  };

  typedef CGAL::Null_functor In; // TODO

  struct Zero_in : public std::unary_function< Type, bool > {
    const bool operator() (const Type &I) const
    { return (Zero_in_T() (I.real()) && Zero_in_T() (I.imag())); }
  };

  // TODO
  typedef CGAL::Null_functor Equal;
  typedef CGAL::Null_functor Overlap;
  typedef CGAL::Null_functor Subset;
  typedef CGAL::Null_functor Proper_Subset;
  typedef CGAL::Null_functor Intersection;
  typedef CGAL::Null_functor Hull;
};

} // namespace INTERN_CET

template< typename T >
class Interval_traits< Cartesian_complex< T > >
  : public INTERN_CET::Interval_traits_base< Cartesian_complex< T >,
                                             typename Interval_traits< T >::Is_interval,
                                             Represents_cartesian_complex_tag > {};

} // namespace CGAL

#endif 	    /* !CGAL_CARTESIAN_COMPLEX_H_ */
