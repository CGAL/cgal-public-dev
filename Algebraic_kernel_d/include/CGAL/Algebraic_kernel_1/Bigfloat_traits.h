#ifndef CGAL_BIGFLOAT_TRAITS_H
#define CGAL_BIGFLOAT_TRAITS_H

#include <CGAL/Algebraic_kernel_1/Cartesian_complex.h>
#include <CGAL/Algebraic_kernel_1/Double_with_exponent.h>
#include <CGAL/Algebraic_kernel_d/Float_traits.h>

namespace CGAL {

// Don't define default, results in more convinient compiler messages
template< class Bigfloat >
class Bigfloat_traits;

template< class Bigfloat, class Type >
struct Convert_to_bigfloat;

// Functor adapting functions
template< class Bigfloat >
const Bigfloat ipow2 (long e) {
  return typename Bigfloat_traits< Bigfloat >::Integral_pow2() (e);
}

template< class Bigfloat >
void set_precision (Bigfloat &x, long prec) {
  return typename Bigfloat_traits< Bigfloat >::Set_precision() (x, prec);
}

template< class Bigfloat >
void set_default_precision (long prec) {
  return typename Bigfloat_traits< Bigfloat >::Set_default_precision() (prec);
}

template< class Bigfloat >
const long get_precision (Bigfloat &x) {
  return typename Bigfloat_traits< Bigfloat >::Get_precision() (x);
}

template< class Bigfloat >
const long get_default_precision () {
  return typename Bigfloat_traits< Bigfloat >::Get_default_precision() ();
}

template< class Bigfloat >
void mult_by_pow2 (Bigfloat &x, long e) {
  return typename Bigfloat_traits< Bigfloat >::Multiply_by_pow2() (x, e);
}

template< class Bigfloat >
const long abs_ilog2 (const Bigfloat& x) {
  return typename Bigfloat_traits< Bigfloat >::Abs_integral_log2() (x);
}

template< class Bigfloat >
const Bigfloat abs_round_up_to_pow2 (const Bigfloat& x) {
  return typename Bigfloat_traits< Bigfloat >::Abs_round_up_to_pow2() (x);
}

template< class Bigfloat >
const Bigfloat root_d (const Bigfloat &x, int d) {
  return typename Bigfloat_traits< Bigfloat >::Root_d() (x, d);
}
template< class Bigfloat >
const Bigfloat root_d (const Bigfloat &x, int d, long p) {
  return typename Bigfloat_traits< Bigfloat >::Root_d() (x, d, p);
}

template< class Bigfloat >
const std::pair< double, long > to_double_exponent (const Bigfloat &x) {
  return typename Bigfloat_traits< Bigfloat >::To_double_exponent() (x);
}

// Specializations

#ifdef CGAL_USE_LEDA
template<>
class Bigfloat_traits< leda_bigfloat > {
public:
  typedef leda_bigfloat Type;

  struct Integral_pow2
    : public std::unary_function< long, leda_bigfloat > {
    const leda_bigfloat operator() (long e) const {
      return leda::ipow2 (e);
    }
  };

  // TODO
  struct Set_precision;
  //    : public std::binary_function< leda_bigfloat, long, void >;

  struct Set_default_precision
    : public std::unary_function< long, void > {
    void operator() (long prec) const {
      leda::bigfloat::set_precision (prec);
    }
  };

  // TODO
  struct Get_precision;

  struct Get_default_precision
    : public std::unary_function< void, long > {
    const long operator() () const {
      return leda::bigfloat::get_precision();
    }
  };

  struct Multiply_by_pow2
    : public std::binary_function< leda_bigfloat, long, void > {
    void operator() (leda_bigfloat &x, long e) const {
      x *= Integral_pow2() (e);
    }
  };

  struct Abs_integral_log2
    : public std::unary_function< leda_bigfloat, long > {
    const long operator() (const leda_bigfloat &x) const {
      CGAL_precondition (leda::ilog2 (x).is_long());
      return leda::ilog2 (x).to_long();
    }
  };

  struct Abs_round_up_to_pow2
    : public std::unary_function< leda_bigfloat, leda_bigfloat > {
    const leda_bigfloat operator() (const leda_bigfloat &x) const {
      return Integral_pow2() (Abs_integral_log2() (x));
    }
  };

  struct Root_d
    : public std::binary_function< leda_bigfloat, int, leda_bigfloat > {
    const leda_bigfloat operator() (const leda_bigfloat &x, int d, long p) const {
      return leda::sqrt_d (x, p, d);
    }
    const leda_bigfloat operator() (const leda_bigfloat &x, int d) const {
      return operator() (x, d, leda::bigfloat::get_precision());
    }
  };

  struct To_double_exponent
    : public std::unary_function< leda_bigfloat, std::pair< double, long > > {
    const std::pair< double, long > operator() (const leda_bigfloat &x) const {
      std::pair< double, long > res;
      res.second = Abs_integral_log2() (x);
      res.first = (x * CGAL::ipow2< leda_bigfloat > (-res.second)).to_double();
      return res;
    }
  };
};

template<>
struct Convert_to_bigfloat< leda_bigfloat, leda_bigfloat >
  : public std::binary_function< leda_bigfloat, long, leda_bigfloat > {
  const leda_bigfloat & operator() (const leda_bigfloat &x, long) {
    return x;
  }
};

template<>
struct Convert_to_bigfloat< leda_bigfloat, leda_integer >
  : public std::binary_function< leda_integer, long, leda_bigfloat > {
  const leda_bigfloat operator() (const leda_integer &n, long p) {
    return leda_bigfloat (n);
    // TODO: Early return?
    if (n.length() <= p || p <= 0)
      return leda_bigfloat (n);
    else
      return leda_bigfloat (n >> (n.length() - p), n.length() - p);
  }
};

template<>
struct Convert_to_bigfloat< leda_bigfloat, leda_rational >
  : public std::binary_function< leda_rational, long, leda_bigfloat > {
  const leda_bigfloat operator() (const leda_rational &q, long) {
    return leda_bigfloat (q.numerator()) / leda_bigfloat (q.denominator());
  }
};

template<>
class Bigfloat_traits< CGAL::LEDA_arithmetic_kernel::Bigfloat_interval > {
public:
  typedef CGAL::LEDA_arithmetic_kernel::Bigfloat_interval Type;

private:
  typedef Bigfloat_traits< leda_bigfloat > Base_BFT;
  typedef ::CGAL::internal::Float_traits< Type > FT;

public:
  struct Abs_integral_log2
    : public std::unary_function< Type, long > {
    const long operator() (const Type &x) const {
      return max (Base_BFT::Abs_integral_log2() (CGAL::lower (x)),
                  Base_BFT::Abs_integral_log2() (CGAL::upper (x)));
    }
  };

  struct Integral_pow2
    : public std::unary_function< long, Type > {
    const Type operator() (long e) const {
      return Type (Base_BFT::Integral_pow2() (e));
    }
  };

  struct Set_default_precision
    : public std::unary_function< long, void > {
    void operator() (long prec) const {
      leda::bigfloat::set_precision (prec);
    }
  };

  // TODO
  struct Get_precision;

  struct Get_default_precision
    : public std::unary_function< void, long > {
    const long operator() () const {
      leda::bigfloat::get_precision();
    }
  };

  struct Multiply_by_pow2
    : public std::binary_function< Type, long, void > {
    void operator() (Type &I, long e) const {
      I *= CGAL::ipow2< leda_bigfloat > (e);
    }
  };

  struct Root_d;

  struct To_double_exponent;
};
#endif // CGAL_USE_LEDA

#ifdef CGAL_USE_GMP
#ifdef CGAL_USE_MPFR
template<>
class Bigfloat_traits< Gmpfr > {
public:
  typedef Gmpfr Type;

private:
  typedef ::CGAL::internal::Float_traits< Gmpfr > FT;

public:
  struct Abs_integral_log2
    : public std::unary_function< Gmpfr, long > {
    const long operator() (const Gmpfr &x) const {
      return mpfr_get_exp (x.fr());
    }
  };

  struct Integral_pow2
    : public std::unary_function< long, Gmpfr > {
    const Gmpfr operator() (long e) const {
      return FT::Mul_by_pow_of_2() (Gmpfr (1), e);
    }
  };

  struct Set_default_precision
    : public std::unary_function< long, void > {
    void operator() (long prec) const {
      Gmpfr::set_default_precision (prec);
    }
  };

  // TODO
  struct Get_precision;

  struct Get_default_precision
    : public std::unary_function< void, long > {
    const long operator() () const {
      return Gmpfr::get_default_precision();
    }
  };

  struct Multiply_by_pow2
    : public std::binary_function< Gmpfr, long, void > {
    void operator() (Gmpfr &x, long e) const {
      x = FT::Mul_by_pow_of_2() (x, e);
    }
  };

  struct Root_d
    : public std::binary_function< Gmpfr, int, Gmpfr > {
    const Gmpfr operator() (const Gmpfr &x, int d, long p) const {
      Gmpfr xp = x;
      mpfr_prec_round (xp.fr(), p, MPFR_RNDN);
      Gmpfr r (0, p);
      mpfr_root (r.fr(), x.fr(), d, MPFR_RNDN);
      return r;
    }
    const Gmpfr operator() (const Gmpfr &x, int d) const {
      Gmpfr r (0);
      mpfr_root (r.fr(), x.fr(), d, MPFR_RNDN);
      return r;
    }
  };

  struct To_double_exponent
    : public std::unary_function< Gmpfr, std::pair< double, long > > {
    const std::pair< double, long > operator() (const Gmpfr &x) const {
      return x.to_double_exp();
    }
  };
};
#endif // CGAL_USE_MPFR

#ifdef CGAL_USE_MPFI
template<>
class Bigfloat_traits< Gmpfi > {
public:
  typedef Gmpfi Type;

private:
  typedef Bigfloat_traits< Gmpfr > Base_BFT;
  typedef ::CGAL::internal::Float_traits< Gmpfi > FT;

public:
  struct Abs_integral_log2
    : public std::unary_function< Gmpfi, long > {
    const long operator() (const Gmpfi &x) const {
      return max (Base_BFT::Abs_integral_log2() (x.inf()),
                  Base_BFT::Abs_integral_log2() (x.sup()));
    }
  };

  struct Integral_pow2
    : public std::unary_function< long, Gmpfi > {
    const Gmpfi operator() (long e) const {
      return Gmpfi (Base_BFT::Integral_pow2() (e));
    }
  };

  struct Set_default_precision
    : public std::unary_function< long, void > {
    void operator() (long prec) const {
      Gmpfi::set_default_precision (prec);
    }
  };

  // TODO
  struct Get_precision;

  struct Get_default_precision
    : public std::unary_function< void, long > {
    const long operator() () const {
      return Gmpfi::get_default_precision();
    }
  };

  struct Multiply_by_pow2
    : public std::binary_function< Gmpfi, long, void > {
    void operator() (Gmpfi &I, long e) const {
      I *= CGAL::ipow2< Gmpfr > (e);
    }
  };

  struct Root_d;

  struct To_double_exponent;
};
#endif // CGAL_USE_MPFI

// TODO
// template<>
// class Bigfloat_traits< Gmpzf >;
#endif // CGAL_USE_GMP

template<>
class Bigfloat_traits< double > {
public:
  typedef double Type;

  struct Abs_integral_log2
    : public std::unary_function< double, long > {
    const long operator() (const double &x) const {
      int e;
      std::frexp (x, &e);
      return e;
    }
  };

  struct Integral_pow2
    : public std::unary_function< long, double > {
    const double operator() (long e) const {
      return std::ldexp (1., e);
    }
  };

  struct Set_default_precision
    : public std::unary_function< long, void > {
    void operator() (long) const {}
  };

  struct Get_precision
    : public std::unary_function< double, long > {
    const long operator() (const double &x) { return 53; }
  };

  struct Get_default_precision
    : public std::unary_function< void, long > {
    const long operator() () const { return 53; }
  };

  struct Multiply_by_pow2
    : public std::binary_function< double, long, void > {
    void operator() (double &x, long e) const {
      x = std::ldexp (x, e);
    }
  };

  struct Root_d
    : public std::binary_function< double, int, double > {
    const double operator() (const double &x, int d, long) const {
      return std::pow (x, 1. / d);
    }
    const double operator() (const double &x, int d) const {
      return std::pow (x, 1. / d);
    }
  };

  struct To_double_exponent
    : public std::unary_function< double, std::pair< double, long > > {
    const std::pair< double, long > operator() (const double &x) const {
      std::pair< double, int > p;
      p.first = std::frexp (x, &(p.second));
      return p;
    }
  };
};

template< typename T >
class Bigfloat_traits< Cartesian_complex< T > > {
  typedef Bigfloat_traits< T > Base_BFT;
public:
  typedef Cartesian_complex< T > Type;

  struct Abs_integral_log2
    : public std::unary_function< Type, long > {
    const long operator() (const Type &z) const {
      return typename Base_BFT::Abs_integral_log2() (z.squared_norm_2()) - 1;
    }
  };

  struct Integral_pow2
    : public std::unary_function< long, Type > {
    const Type operator() (long e) const {
      return Type (typename Base_BFT::Integral_pow2() (e));
    }
  };

  struct Set_default_precision
    : public std::unary_function< long, void > {
    void operator() (long e) const {
      typename Base_BFT::Set_default_precision() (e);
    }
  };

  struct Get_default_precision
    : public std::unary_function< void, long > {
    const long operator() () const {
      return typename Base_BFT::Get_default_precision()();
    }
  };

  struct Get_precision
    : public std::unary_function< Type, long > {
    const long operator() (const Type &z) const {
      return max (typename Base_BFT::Get_precision() (z.real()),
                  typename Base_BFT::Get_precision() (z.imag()));
    }
  };
};

} // namespace CGAL

#endif // CGAL_BIGFLOAT_TRAITS_H
