#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <boost/lexical_cast.hpp>

#include <CGAL/basic.h>
#include <CGAL/GMP_arithmetic_kernel.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial_type_generator.h>

#include <CGAL/Algebraic_kernel_d/Polynomial_parser_d.h>
#include <CGAL/Bitsize.h>

/*
namespace CGAL {
template<>
struct Coercion_traits< Gmpq, Gmpfr > {
  struct Cast : public std::unary_function< Gmpq, Gmpfr > {
    const Gmpfr operator() (const Gmpq &q) {
      return Gmpfr (q.numerator()) / Gmpfr (q.denominator());
    }
  };
};
}
*/

typedef CGAL::GMP_arithmetic_kernel AK;
typedef AK::Integer ZZ;
typedef CGAL::Polynomial_type_generator< ZZ, 1 >::Type IPoly_1;
typedef AK::Rational QQ;
typedef CGAL::Polynomial_type_generator< QQ, 1 >::Type QPoly_1;
typedef AK::Bigfloat BF;
typedef CGAL::Polynomial_type_generator< BF, 1 >::Type BFPoly_1;

// Original from Polynomial_type.h
template <class NT> inline
void org_print_maple_monomial(std::ostream& os, const NT& coeff,
    const char *var, int expn)
{
  using namespace CGAL;

  if (expn == 0 || coeff != NT(1)) {
    os << CGAL::oformat(coeff, Parens_as_product_tag());
    if (expn >= 1) os << "*";
  }
  if (expn >= 1) {
    os << var;
    if (expn > 1) os << "^" << CGAL::oformat(expn);
  }
}

template <class NT>
void org_output_maple(std::ostream& os, const CGAL::Polynomial< NT > &p) {
  using namespace CGAL;

  const char *varname;
  char vnbuf[42];
    
  // use variable names x, y, z, w1, w2, w3, ...
  if (Polynomial_traits_d<NT>::d < 3) {
    static const char *varnames[] = { "x", "y", "z" };
    varname = varnames[Polynomial_traits_d<NT>::d];
  } else {
    sprintf(vnbuf, "w%d", Polynomial_traits_d<NT>::d - 2);
    varname = vnbuf;
  }
    
  int i = p.degree();
  org_print_maple_monomial(os, p[i], varname, i);
  while (--i >= 0) {
    if (! CGAL::certainly (CGAL::is_zero (p[i]))) { // TODO: repair
      os << " + ";
      org_print_maple_monomial(os, p[i], varname, i);
    }
  }
}


void print_mpsolve (std::ostream &out,
                    const IPoly_1 &f,
                    const std::vector< std::string > &comments) {
  for (int i = 0; i < comments.size(); ++i)
    out << "!" << comments[i] << '\n';

  out << "dri\n"
      << "0\n"
      << f.degree() << '\n'
      << '\n';

  for (int i = 0; i <= f.degree(); ++i)
    out << f[i] << '\n';
}

void print_mpsolve (std::ostream &out,
                    const QPoly_1 &f,
                    const std::vector< std::string > &comments) {
  for (int i = 0; i < comments.size(); ++i)
    out << "!" << comments[i] << '\n';

  out << "drq\n"
      << "0\n"
      << f.degree() << '\n'
      << '\n';

  for (int i = 0; i <= f.degree(); ++i)
    out << f[i].numerator() << '\n'
        << f[i].denominator() << '\n';
}

void print_mpsolve (std::ostream &out,
                    const BFPoly_1 &f,
                    const std::vector< std::string > &comments) {
  for (int i = 0; i < comments.size(); ++i)
    out << "!" << comments[i] << '\n';

  out << "drf\n"
      << "100000\n" // TODO: get precision
      << f.degree() << '\n'
      << '\n';

  for (int i = 0; i <= f.degree(); ++i)
    out << f[i] << '\n';
}

template< class Poly >
void print_cgal (std::ostream &out,
                 const Poly &f,
                 const std::vector< std::string > &comments) {
  for (int i = 0; i < comments.size(); ++i)
    out << "#" << comments[i] << '\n';
  CGAL::set_ascii_mode (out);
  out << f << '\n';
}

template< class Poly >
void print_maple (std::ostream &out,
                  const Poly &f,
                  const std::vector< std::string > &comments) {
  for (int i = 0; i < comments.size(); ++i)
    out << "#" << comments[i] << '\n';
  org_output_maple (out, f);
  out << '\n';
}

template< class Poly >
void print_gp_script (std::ostream &out,
                      const Poly &f,
                      const std::vector< std::string > &comments) {
  for (int i = 0; i < comments.size(); ++i)
    out << "\\" << comments[i] << '\n';

  out << "f = ";
  org_output_maple (out, f);
  out << ";\n"
      << "sols = polroots (f);\n"
      << "##\n"
      << "print (\"Solutions:\");\n"
      << "for (i = 1, poldegree (f), print (real (sols[i]), \"\\t\", imag (sols[i])));\n";
}

template< class Poly >
void print_sage_complex_script (std::ostream &out,
                                const Poly &f,
                                const std::vector< std::string > &comments) {
  for (int i = 0; i < comments.size(); ++i)
    out << "#" << comments[i] << '\n';

  out << "R.<x> = QQ[];\n"
      << "from sage.rings.polynomial.complex_roots import *;\n"
      << "f = ";
  org_output_maple (out, f);
  out << ";\n"
      << "print \"\";\n"
      << "%time sols = complex_roots (f);\n"
      << "print \"\";\n"
      << "for i in sols :\n"
      << "\tprint \"ROOT: [\", real(i[0]), \"\\t\", imag(i[0]), \"]\"\n\n"
    //<< "\tprint i\n\n"
      << "quit\n";
}

template< class Poly >
void print_sage_real_script (std::ostream &out,
                             const Poly &f,
                             const std::vector< std::string > &comments) {
  for (int i = 0; i < comments.size(); ++i)
    out << "#" << comments[i] << '\n';

  out << "R.<x> = QQ[];\n"
      << "from sage.rings.polynomial.real_roots import *;\n"
      << "f = ";
  org_output_maple (out, f);
  out << ";\n"
      << "print \"\";\n"
      << "%time sols = real_roots (f);\n"
      << "print \"\";\n"
      << "for i in sols :\n"
      << "\tprint i\n\n"
      << "quit\n";
}


template< class Poly >
void write_all (const std::string &basename,
                const Poly &f,
                const std::vector< std::string > &comments) {
  using namespace std;

  {
    ofstream out ((basename + ".cgal").c_str());
    print_cgal (out, f, comments);
  }
  {
    ofstream out ((basename + ".maple").c_str());
    print_maple (out, f, comments);
  }
  {
    ofstream out ((basename + ".mpsolve").c_str());
    print_mpsolve (out, f, comments);
  }
  {
    ofstream out ((basename + ".gp").c_str());
    print_gp_script (out, f, comments);
  }
  {
    ofstream out ((basename + ".sage_real").c_str());
    print_sage_real_script (out, f, comments);
  }
  {
    ofstream out ((basename + ".sage_complex").c_str());
    print_sage_complex_script (out, f, comments);
  }

  std::cout << "Degree: " << f.degree() << std::endl
            << "Bitsize: " << CGAL::bitsize (f) << std::endl;
}

int main (int argc, char **argv) {
  using namespace CGAL;
  using namespace std;

  vector< string > args (argv, argv+argc);
  istream *in_ptr = &cin;
  if (args.size() > 1 && args.back() != "--")
    in_ptr = new ifstream (args.back().c_str());

  string line;
  vector< string > comments;
  while (line.empty()
         || line[0] == '!'
         || line[0] == '#') {
    if (! line.empty())
      comments.push_back (line.substr (1));
    getline (*in_ptr, line);
  }

  string basename = "polynomial";
  if (args.size() > 1 && args.back() != "--") {
    const size_t ext_pos = args.back().find_last_of ('.');
    basename = args.back().substr (0, ext_pos);
  }

  CGAL::set_pretty_mode (cout);

  if (line == "dri") {

    vector< ZZ > coeffs;
    int n;
    long p;
    do { getline (*in_ptr, line); } while (line.empty() && in_ptr->good());
    p = boost::lexical_cast< long > (line);
    do { getline (*in_ptr, line); } while (line.empty() && in_ptr->good());
    n = boost::lexical_cast< int > (line);
    
    while (n-- >= 0) {
      do { getline (*in_ptr, line); } while (line.empty() && in_ptr->good());
      coeffs.push_back (boost::lexical_cast< ZZ > (line));
    }
    
    IPoly_1 f (coeffs.begin(), coeffs.end());
    write_all (basename, f, comments);

  } else if (line == "drq") {

    vector< QQ > coeffs;
    int n;
    long p;
    do { getline (*in_ptr, line); } while (line.empty() && in_ptr->good());
    p = boost::lexical_cast< long > (line);
    do { getline (*in_ptr, line); } while (line.empty() && in_ptr->good());
    n = boost::lexical_cast< int > (line);
    
    while (n-- >= 0) {
      do { getline (*in_ptr, line); } while (line.empty() && in_ptr->good());
      coeffs.push_back (boost::lexical_cast< ZZ > (line));
      do { getline (*in_ptr, line); } while (line.empty() && in_ptr->good());
      coeffs.back() /= boost::lexical_cast< ZZ > (line);
    }
    
    QPoly_1 f (coeffs.begin(), coeffs.end());
    write_all (basename, f, comments);

  } else if (line == "drf") {

    vector< BF > coeffs;
    int n;
    long p;
    do { getline (*in_ptr, line); } while (line.empty() && in_ptr->good());
    p = boost::lexical_cast< long > (line);
    do { getline (*in_ptr, line); } while (line.empty() && in_ptr->good());
    n = boost::lexical_cast< int > (line);
    
    while (n-- >= 0) {
      do { getline (*in_ptr, line); } while (line.empty() && in_ptr->good());
      coeffs.push_back (boost::lexical_cast< BF > (line));
    }
    
    BFPoly_1 f (coeffs.begin(), coeffs.end());
    write_all (basename, f, comments);

  } else {

    CGAL::Polynomial_parser_d< IPoly_1 > parser;
    IPoly_1 f;
    if (parser (line, f)) {
      write_all (basename, f, comments);
    } else {
      CGAL::Polynomial_parser_d< QPoly_1,
        CGAL::Mixed_rational_parser_policy< QPoly_1 > > parser;
      QPoly_1 f;
      if (parser (line, f)) {
        write_all (basename, f, comments);
      } else {
        CGAL::Polynomial_parser_d< BFPoly_1 > parser;
        BFPoly_1 f;
        if (parser (line, f)) {
          write_all (basename, f, comments);
        }
      }
    }
  }

  // CLEAN UP
  if (in_ptr && in_ptr != &cin)
    delete in_ptr;
  return 0;
}
