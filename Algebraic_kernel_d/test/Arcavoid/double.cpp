#include <cassert>
#include <fenv.h>
#include <iostream>
#include <cmath>
#include <cfloat>
#include <string>
#include <bitset>
#include <limits>
#include <boost/lexical_cast.hpp>

template< typename T, int N >
struct To_binary {
  std::string operator() (const T &value) {
    typedef std::bitset< N > Bitset;
    return Bitset (reinterpret_cast< const long & > (value)).to_string();
  }
};

template< typename T >
struct To_binary< T, 64 > {
  std::string operator() (const T &value) {
    typedef std::bitset< 64 > Bitset;
    std::string s = Bitset (reinterpret_cast< const long & > (value)).to_string();
    s.insert (12, " ");
    s.insert (1, " ");
    return s;
  }
};

template< typename T >
const std::string to_binary (const T & value) {
  return To_binary< T, sizeof (T) << 3 >() (value);
}

const double trunc_prec_RND_TO_INF (const double &x, int prec_loss) {
  long x_lng = reinterpret_cast< const long & > (x) | (~(-1L << prec_loss));
  assert (std::abs (reinterpret_cast< const double & > (x_lng)) > std::abs (x));
  return reinterpret_cast< const double & > (x_lng);

  if (x != 0.) {
    long x_lng = reinterpret_cast< const long & > (x);
    x_lng &= (-1L << prec_loss);
    unsigned long x_ulng = reinterpret_cast< const unsigned long & > (x_lng);
    x_ulng += (1 << prec_loss);
    const double ret_val = reinterpret_cast< const double & > (x_ulng);
    assert (abs(x) < abs(ret_val));
    return ret_val;
  } else {
    return DBL_MIN;
  }
}
const double trunc_prec_RND_TO_ZERO (const double &x, int prec_loss) {
  if (x != 0.) {
    long x_lng = reinterpret_cast< const long & > (x);
    x_lng &= (-1L << prec_loss);
    unsigned long x_ulng = reinterpret_cast< const unsigned long & > (x_lng);
    x_ulng += (1 << prec_loss);
    const double ret_val = reinterpret_cast< const double & > (x_ulng);
    assert (abs(x) < abs(ret_val));
    return ret_val;
  } else {
    return DBL_MIN;
  }
}

int main () {
  using namespace std;

  double inf = -1. / 0.;
  double min_dbl = DBL_MIN;
  double max_dbl = DBL_MAX;
  double nan = +0. / 0.;
  unsigned long lnan = reinterpret_cast< unsigned long & > (nan);
  lnan ^= (1UL << 63);
  nan = reinterpret_cast< double & > (lnan);
  double a = -sqrt(2.) + 0.25;
  a = -0.999999999999999;
  a = 0;
  float tmp = a;
  double f_a = tmp;

  int a_exp;
  double a_mant = frexp (a, &a_exp);
  a_exp += (1 << 10);

  cerr.precision (20);

  cerr << "a: " << a << endl;
  cerr << "a_mant: " << a_mant << endl;

  //cerr << "(long)a: " << reinterpret_cast< long & > (a) << endl;
  cerr << "                 " << 'S' << ' ' << string (11, 'e') << ' ' << string (52, 'M') << endl;
  cerr << "(binary)inf:     " << to_binary (inf) << endl;
  cerr << "inf:             " << inf << endl;
  cerr << "(binary)nan:     " << to_binary (nan) << endl;
  cerr << "nan:             " << nan << endl;
  cerr << "(binary)min_dbl: " << to_binary (min_dbl) << endl;
  cerr << "(binary)max_dbl: " << to_binary (max_dbl) << endl;
  cerr << "(binary)f_a:     " << to_binary (f_a) << endl;
  cerr << "(binary)a:       " << to_binary (a)     << "\t" << a << endl;
  double rnd_a = trunc_prec_RND_TO_INF (a, 5);
  cerr << "(binary)rnd_a:   " << to_binary (rnd_a) << "\t" << rnd_a << "\t" << "| | > |a|? " << (abs(rnd_a) > abs(a) ? "true  " : "false") << "\trel_error: " << (abs (a-rnd_a) / abs (a)) << endl;
  rnd_a = ldexp (rnd_a, 48);
  cerr << "(binary)rnd_a:   " << to_binary (rnd_a) << "\t" << rnd_a << "\t" << "| | > |a|? " << (abs(rnd_a) > abs(a) ? "true  " : "false") << "\trel_error: " << (abs (a-rnd_a) / abs (a)) << endl;
  cerr << "(binary)a_mant:  " << to_binary (a_mant) << endl;
  cerr << "(binary)a_exp:   " << (reinterpret_cast< long & > (a) >> 52) - 1023 << endl;

  double huge = 1.23456789e10;
  double smaller = std::pow (huge, 0.1);
  double small = 1.23456789;
  double larger = std::pow (small, 10.);

  fesetround (FE_UPWARD);
  double smaller_rnd_up = std::pow (huge, 0.1);
  double larger_rnd_up = std::pow (small, 10.);

  cerr << "huge:           " << huge << endl
       << "smaller:        " << smaller << endl
       << "smaller_rnd_up: " << smaller_rnd_up << endl;

  cerr << "small:         " << small << endl
       << "larger:        " << larger << endl
       << "larger_rnd_up: " << larger_rnd_up << endl;

  return 0;
}
