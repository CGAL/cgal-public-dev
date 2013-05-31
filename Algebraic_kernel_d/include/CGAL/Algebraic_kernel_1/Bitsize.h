/*
** Bitsize.h
**
** Made by Alexander Kobel
** Login   <perpeduumimmobile@lidinoid>
**
** Started on  Thu Apr 22 12:59:47 2010 Alexander Kobel
** Last update Thu Apr 22 12:59:47 2010 Alexander Kobel
*/

#ifndef   	BITSIZE_H_
# define   	BITSIZE_H_

#include <CGAL/basic.h>
#include <CGAL/Polynomial.h>
#include <CGAL/LEDA_arithmetic_kernel.h>
#include <CGAL/GMP_arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_1/Cartesian_complex.h>
#include <gmpxx.h>

namespace CGAL {

template< typename T >
const size_t bitsize (const T &);

template< typename InputIterator >
const size_t bitsize (InputIterator, InputIterator);

template< class NT >
struct Bitsize {};

template< typename NT >
struct Bitsize< Polynomial< NT > > {
  const size_t operator() (const Polynomial< NT > &f) const {
    return bitsize (f.begin(), f.end());
  }
};

template< typename T >
struct Bitsize< CGAL::Cartesian_complex< T > > {
  const size_t operator() (const CGAL::Cartesian_complex< T > &z) const {
    return max (bitsize (z.real()), bitsize (z.imag()));
  }
};

#ifdef CGAL_USE_LEDA
/* LEDA specializations */
template<>
struct Bitsize< leda_integer > {
  const size_t operator() (const leda_integer &n) const {
    return n.length();
  }
};
template<>
struct Bitsize< leda_rational > {
  const size_t operator() (const leda_rational &q) const {
    return q.numerator().length() + q.denominator().length();
  }
};
template<>
struct Bitsize< leda_bigfloat > {
  const size_t operator() (const leda_bigfloat &x) const {
    return x.get_significant_length();
  }
};
template<>
struct Bitsize< leda_real > {
  const size_t operator() (const leda_real &x) const {
    if (! x.is_rational())
      return std::numeric_limits< size_t >::max();
    else
      return bitsize (x.to_rational());
  }
};
#endif // CGAL_USE_LEDA

/* GMP specializations */
#ifdef CGAL_USE_GMP
template<>
struct Bitsize< mpz_t > {
  const size_t operator() (const mpz_t &n) const {
    return ::mpz_sizeinbase (n, 2);
  }
};
template<>
struct Bitsize< mpz_class > {
  const size_t operator() (const mpz_class &n) const {
    return ::mpz_sizeinbase (n.get_mpz_t(), 2);
  }
};
template<>
struct Bitsize< Gmpz > {
  const size_t operator() (const Gmpz &n) const {
    return ::mpz_sizeinbase (n.mpz(), 2);
  }
};

template<>
struct Bitsize< mpq_t > {
  const size_t operator() (const mpq_t &q) const {
    return ::mpz_sizeinbase ( mpq_numref (q), 2)
      + ::mpz_sizeinbase ( mpq_denref (q), 2);
  }
};
template<>
struct Bitsize< mpq_class > {
  const size_t operator() (const mpq_class &q) const {
    return ::mpz_sizeinbase (mpq_numref (q.get_mpq_t()), 2)
      + ::mpz_sizeinbase (mpq_denref (q.get_mpq_t()), 2);
  }
};
template<>
struct Bitsize< Gmpq > {
  const size_t operator() (const Gmpq &q) const {
    return ::mpz_sizeinbase (mpq_numref (q.mpq()), 2)
      + ::mpz_sizeinbase (mpq_denref (q.mpq()), 2);
  }
};
#endif // CGAL_USE_GMP

#ifdef CGAL_USE_MPFR
template<>
struct Bitsize< Gmpfr > {
  const size_t operator() (const Gmpfr &x) const {
    return x.get_precision();
  }
};
#endif // CGAL_USE_MPFR

template< typename T >
inline
const size_t bitsize (const T &t) {
  return Bitsize< T >() (t);
}

template< typename InputIterator >
inline
const size_t bitsize (InputIterator begin, InputIterator end) {
  size_t L = 0;
  while (begin != end) {
    L = max (L, bitsize (*begin));
    ++begin;
  }
  return L;
}

} // namespace CGAL

#endif 	    /* !BITSIZE_H_ */
