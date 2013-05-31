#ifndef CGAL_ARCAVOID_FORMATTERS_H
#define CGAL_ARCAVOID_FORMATTERS_H

#include <iostream>
#include <string>
#include <boost/lexical_cast.hpp>

#include <CGAL/Cartesian_complex.h>
#include <CGAL/GMP_arithmetic_kernel.h>

#ifdef __GNUC__
#include <cxxabi.h>
template< class T >
std::string type_name (const T &t) {
  int status;
  return abi::__cxa_demangle (typeid (t).name(), 0, 0, &status);
}
#else
template< class T >
std::string type_name (const T &t) {
  return typeid (t).name();
}
#endif

namespace CGAL {

class Question_mark_interval_formatter {
  int ed;
  bool sc;
  char ec;
  int b;
public:
  Question_mark_interval_formatter (int error_digits = 0,
                                    bool scientific = false,
                                    char exp_char = 'e',
                                    int base = 10)
    : ed (error_digits), sc (scientific), ec (exp_char), b (base) {}
  
  const int error_digits () const { return ed; }
  const bool scientific () const { return sc; }
  const char exponent_character () const { return ec; }
  const int base () const { return b; }
};

#if defined(CGAL_USE_GMP) && defined(CGAL_USE_MPFR)

template<>
struct Output_rep< Gmpfr, Parens_as_product_tag > {
  const Gmpfr &x;
  
  Output_rep (const Gmpfr &x) : x(x) {}

  std::ostream & operator() (std::ostream &out) const {
    const int width = 40;
    char buf [width+20];
    std::string s;
    
    const double abs_x = CGAL::abs (CGAL::to_double (x));
    
    if (false && // always use scientific format
        1e-3 < abs_x && abs_x < 10) {
      int n = mpfr_sprintf (buf, "%.*RNf", width-3, x.fr());
      s = buf;
    } else {
      int n = mpfr_sprintf (buf, "%.*RNe", width-8, x.fr());
      s = buf;
      size_t ind_e = s.find_last_of ('e');
      if (ind_e != std::string::npos) {
        std::string exp_str = s.substr (++ind_e);
        long exp = boost::lexical_cast< long > (exp_str);
        n = mpfr_sprintf (buf, "%+04i", exp);
        s = s.substr (0, ind_e);
        s += std::string (buf);
      }
    }
    
    if (s.length() < width)
      s.insert (s.begin(), width - s.length(), ' ');

    return out << s;
  }
};

template<>
struct Output_rep< Gmpfr, Null_tag > {
  const Gmpfr &x;
  
  Output_rep (const Gmpfr &x) : x(x) {}

  std::ostream & operator() (std::ostream &out) const {
    return Output_rep< Gmpfr, Parens_as_product_tag > (x) (out);
  }
};

template< class F >
struct Output_rep< Gmpfr, F > {
  const Gmpfr &x;
  
  Output_rep (const Gmpfr &x) : x(x) {}

  std::ostream & operator() (std::ostream &out) const {
    return Output_rep< Gmpfr, Parens_as_product_tag > (x) (out);
  }
};

template<>
struct Output_rep< Gmpfi, Question_mark_interval_formatter > {
  const Gmpfi &x;
  const Question_mark_interval_formatter formatter;

  Output_rep (const Gmpfi& x) : x(x) {}
  Output_rep (const Gmpfi& x,
              const Question_mark_interval_formatter &formatter)
    : x (x), formatter (formatter) {}

  std::ostream & operator() (std::ostream &out) const {
    const int base = formatter.base();
    const char exp_char = formatter.exponent_character();
    int error_digits = formatter.error_digits();
    bool scientific = formatter.scientific();

    const Gmpfr a = CGAL::lower (x);
    const Gmpfr b = CGAL::upper (x);

    const bool exact = x.is_point();

    return out << "[" << CGAL::oformat (a) << " .. " << CGAL::oformat (b) << "]";

    if (b < a)
      return out << "[..]";

    if (mpfr_number_p (a.fr()) == 0
        || mpfr_number_p (b.fr()) == 0)
      return out << "[" << a << " .. " << b << "]";

    char *a_cstr, *b_cstr;
    mp_exp_t a_exp, b_exp;
    mpz_t a_mpz, b_mpz;

    a_cstr = mpfr_get_str (NULL, &a_exp, base, 0, a.fr(), MPFR_RNDD);
    CGAL_assertion (a_cstr != NULL);
    b_cstr = mpfr_get_str (NULL, &b_exp, base, 0, b.fr(), MPFR_RNDU);
    CGAL_assertion (b_cstr != NULL);

    int digits = strlen (a_cstr);
    if (a_cstr[0] == '-')
      --digits;
    a_exp -= digits;

    digits = strlen (b_cstr);
    if (b_cstr[0] == '-')
      --digits;
    b_exp -= digits;

    mpz_init_set_str (a_mpz, a_cstr, base);
    mpz_init_set_str (b_mpz, b_cstr, base);
    mpfr_free_str (a_cstr);
    mpfr_free_str (b_cstr);

    mpz_t tmp;
    mpz_init (tmp);

    if (mpfr_zero_p (a.fr()))
      a_exp = b_exp;
    if (mpfr_zero_p (b.fr()))
      b_exp = a_exp;

    int exp_delta = CGAL::abs (a_exp - b_exp);
    if (a_exp < b_exp) {
      if (mpz_sizeinbase (a_mpz, base) < exp_delta) {
        if (mpz_sgn (a_mpz) < 0)
          mpz_set_si (a_mpz, -1);
        else
          mpz_set_ui (a_mpz, 0);
      } else {
        mpz_ui_pow_ui (tmp, base, exp_delta);
        mpz_fdiv_q (a_mpz, a_mpz, tmp);
      }
      a_exp = b_exp;
    } else if (b_exp < a_exp) {
      if (mpz_sizeinbase (b_mpz, base) < exp_delta) {
        if (mpz_sgn (a_mpz) > 0)
          mpz_set_si (a_mpz, 1);
        else
          mpz_set_ui (a_mpz, 0);
      } else {
        mpz_ui_pow_ui (tmp, base, exp_delta);
        mpz_fdiv_q (b_mpz, b_mpz, tmp);
      }
      b_exp = a_exp;
    }

    int exp = a_exp;

    mpz_t cur_error, max_error;
    mpz_init (cur_error);
    mpz_sub (cur_error, b_mpz, a_mpz);
    mpz_init (max_error);

    if (error_digits == 0)
      mpz_set_ui (max_error, 2);
    else {
      mpz_ui_pow_ui (max_error, 10, error_digits);
      mpz_sub_ui (max_error, max_error, 1);
      mpz_mul_2exp (max_error, max_error, 1);
    }
    
    int cur_error_digits = mpz_sizeinbase (cur_error, base);
    int max_error_digits = mpz_sizeinbase (max_error, base);
    int k = cur_error_digits - 1 - max_error_digits;

    if (k > 0) {
      mpz_ui_pow_ui (tmp, base, k);
      mpz_fdiv_q (a_mpz, a_mpz, tmp);
      mpz_cdiv_q (b_mpz, b_mpz, tmp);
      exp += k;
      mpz_sub (cur_error, b_mpz, a_mpz);
    }

    while (mpz_cmp (cur_error, max_error) > 0) {
      mpz_fdiv_q_ui (a_mpz, a_mpz, base);
      mpz_cdiv_q_ui (b_mpz, b_mpz, base);
      ++exp;
      mpz_sub (cur_error, b_mpz, a_mpz);
    }

    mpz_add (a_mpz, a_mpz, b_mpz);

    if (mpz_sgn (a_mpz) >= 0)
      mpz_cdiv_q_2exp (a_mpz, a_mpz, 1);
    else
      mpz_fdiv_q_2exp (a_mpz, a_mpz, 1);

    mpz_cdiv_q_2exp (cur_error, cur_error, 1);

    char *tmp_cstr;

    tmp_cstr = mpz_get_str (NULL, base, a_mpz);
    digits = strlen (tmp_cstr);
    std::string mant_str = tmp_cstr;
    std::string sign_str = "";
    std::string error_str = "";
    if (mant_str[0] == '-') {
      mant_str = mant_str.substr (1);
      sign_str = "-";
    }

    free (tmp_cstr);

    if (error_digits != 0) {
      tmp_cstr = mpz_get_str (NULL, base, cur_error);
      error_str = tmp_cstr;
      free (tmp_cstr);
    }

    mpz_clear (a_mpz);
    mpz_clear (b_mpz);
    mpz_clear (tmp);
    mpz_clear (cur_error);
    mpz_clear (max_error);

    if (exp > 0)
      scientific = true;
    int sci_exp = exp + digits - 1;
    if (CGAL::abs (sci_exp) >= 6)
      scientific = true;

    const std::string mark = exact ? "?" : "!";
    
    if (scientific)
      return out << sign_str << mant_str[0] << "." << mant_str.substr (1) << mark
                 << error_str << exp_char << sci_exp;

    if (exp + digits <= 0)
      return out << sign_str << "0." << std::string (-(exp+digits), '0')
                 << mant_str << mark << error_str;

    return out << sign_str << mant_str.substr (0, exp+digits) << "." << mant_str.substr (exp+digits)
               << mark << error_str;
  }
};

template< class F >
struct Output_rep< Gmpfi, F > {
  const Gmpfi &x;

  Output_rep (const Gmpfi &x) : x(x) {}
  std::ostream & operator() (std::ostream &out) const {
    return Output_rep< Gmpfi, Question_mark_interval_formatter > (x) (out);
  }
};

template<>
struct Output_rep< Gmpfi, Parens_as_product_tag > {
  const Gmpfi &x;

  Output_rep (const Gmpfi &x) : x(x) {}
  std::ostream & operator() (std::ostream &out) const {
    return Output_rep< Gmpfi, Question_mark_interval_formatter > (x) (out);
  }
};

#endif // defined(CGAL_USE_GMP) && defined(CGAL_USE_MPFR)

template< class T, class F >
struct Output_rep< CGAL::Cartesian_complex< T >, F > {
  const CGAL::Cartesian_complex< T > &z;

  Output_rep (const CGAL::Cartesian_complex< T > &z) : z(z) {}
  std::ostream & operator() (std::ostream &out) const {
    return out << Output_rep< T, F > (z.real()) << " + I*" << Output_rep< T, F > (z.imag());
  }
};

template< class T >
struct Output_rep< CGAL::Cartesian_complex< T >, Parens_as_product_tag > {
  const CGAL::Cartesian_complex< T > &z;

  Output_rep (const CGAL::Cartesian_complex< T > &z) : z(z) {}
  std::ostream & operator() (std::ostream &out) const {
    return out << "(" << Output_rep< T, Parens_as_product_tag > (z.real()) << ", "
               << Output_rep< T, Parens_as_product_tag > (z.imag()) << ")";
  }
};

}

#endif // CGAL_ARCAVOID_FORMATTERS_H
