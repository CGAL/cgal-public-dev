#ifndef CGAL_ARCAVOID_READ_POLY_H
#define CGAL_ARCAVOID_READ_POLY_H

#include <iostream>
#include <locale>
#include <string>
#include <vector>

#include <boost/lexical_cast.hpp>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d/Polynomial_parser_d.h>

inline
void ltrim (std::string &str, const std::locale &loc) {
  std::string::size_type pos = 0;
  while (pos < str.size() && std::isspace (str[pos], loc))
    ++pos;
  str.erase (0, pos);
}
 
inline
void rtrim (std::string &str, const std::locale &loc) {
  std::string::size_type pos = str.size();
  while (pos > 0 && isspace (str[pos - 1], loc))
    --pos;
  str.erase (pos);
}
 
inline
void btrim (std::string &str, const std::locale &loc = std::locale()) {
  ltrim (str, loc);
  rtrim (str, loc);
}

inline
const std::string & next (std::istream &in, std::string &line) {
  if (! in.good())
    throw std::ios_base::failure ("End of input");
  do {
    std::getline (in, line);
    btrim (line);
  } while (line.empty()
           || line[0] == '!' || line[0] == '#' || line[0] == '%');

  return line;
}

template< class Poly >
bool read_poly (std::istream &in, Poly &f) {
  using namespace std;

  typedef typename Poly::NT CT;

  try {
    string line;
    next (in, line);

    if (line == "dri") {
      // MPSOLVE dense real integer
      clog << "Read MPSOLVE dense real integer polynomial" << endl;
      
      long p = boost::lexical_cast< long > (next (in, line));
      // clog << "p = " << p << endl;
      int n = boost::lexical_cast< int > (next (in, line));
      // clog << "n = " << n << endl;

      typedef typename CGAL::Get_arithmetic_kernel< CT >::Arithmetic_kernel AK;
      typedef typename AK::Integer IT;
      vector< IT > coeffs;
      // clog << "p = " << p << endl;
      while (coeffs.size() <= n) {
#if USE_LEDA_AK
        coeffs.push_back (IT (next (in, line).c_str()));
#else
        coeffs.push_back (boost::lexical_cast< IT > (next (in, line)));
#endif
      }
      // clog << "n = " << n << endl;

      // for (int i = 0; i <= n; ++i)
      //   std::cerr << "f[" << i << "] = " << coeffs[i] << std::endl;
      
      f = Poly (coeffs.begin(), coeffs.end());
      return true;
    }
    
    if (line == "drq") {
      // MPSOLVE dense real rational
      clog << "Read MPSOLVE dense real rational polynomial" << endl;
      
      long p = boost::lexical_cast< long > (next (in, line));
      int n = boost::lexical_cast< int > (next (in, line));
      
      typedef typename CGAL::Get_arithmetic_kernel< CT >::Arithmetic_kernel AK;
      typedef typename AK::Integer IT;
      vector< CT > coeffs;
      while (coeffs.size() <= n) {
        CT c = boost::lexical_cast< IT > (next (in, line));
        c /= boost::lexical_cast< IT > (next (in, line));
        coeffs.push_back (c);
      }
      
      f = Poly (coeffs.begin(), coeffs.end());
      return true;
    }

    if (line[0] == 'P') {
      in >> f;
      return true;
    }
    
#if USE_LEDA_AK
    return (CGAL::Polynomial_parser_d< Poly > () (line, f));
#else
    return (CGAL::Polynomial_parser_d< Poly > () (line, f)
            || CGAL::Polynomial_parser_d< Poly, CGAL::Mixed_rational_parser_policy< Poly > > () (line, f));
#endif
  } catch (...) {
    return false;
  }
}

#endif // CGAL_ARCAVOID_READ_POLY_H
