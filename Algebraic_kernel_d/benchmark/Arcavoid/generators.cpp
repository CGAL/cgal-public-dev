#include <iostream>
#include <ctime>
#include <boost/lexical_cast.hpp>

#include <CGAL/Gmpz.h>
#include <CGAL/GMP_arithmetic_kernel.h>

#include <CGAL/Bitsize.h>
#include <CGAL/Bounds.h>
#include "include/generators.h"

template< class Polynomial >
const std::string to_maple (const Polynomial &f) {
  std::stringstream ss;
  ss << "f := ";

  int i = f.degree();
  ss << "(" << f[i] << ")*z^" << i << "\n";
  for (--i; i >= 2; --i)
    ss << "   + (" << f[i] << ")*z^" << i << "\n";
  ss << "   + (" << f[1] << ")*z\n";
  ss << "   + (" << f[0] << ")";

  return ss.str();
}

template< class Polynomial >
const std::string to_mathematica (const Polynomial &f) {
  std::stringstream ss;
  ss << "f = ";

  int i = f.degree();
  ss << "(" << f[i] << ")*z^" << i;
  for (--i; i >= 2; --i)
    ss << " + (" << f[i] << ")*z^" << i;
  ss << " + (" << f[1] << ")*z";
  ss << " + (" << f[0] << ")";

  return ss.str();
}

template< class Polynomial >
const std::string to_matlab (const Polynomial &f) {
  std::stringstream ss;
  ss << "f = [ ";

  for (int i = f.degree(); i >= 0; --i)
    ss << f[i] << " ";

  ss << "]";
  return ss.str();
}

template< class Polynomial >
const std::string to_mpsolve (const Polynomial &f) {
  std::stringstream ss;
  ss << "dri\n"; // dense, real coefficients, integers
  ss << "0\n"; // infinite input precision
  ss << f.degree() << "\n";

  for (int i = 0; i <= f.degree(); ++i)
    ss << f[i] << "\n";

  return ss.str();
}

template< class Polynomial >
const std::string to_pari (const Polynomial &f) {
  std::stringstream ss;
  ss << "f = "; // dense, real coefficients, integers
  int i = f.degree();
  ss << "(" << f[i] << ")*z^" << i;
  for (--i; i >= 2; --i)
    ss << " + (" << f[i] << ")*z^" << i;
  ss << " + (" << f[1] << ")*z";
  ss << " + (" << f[0] << ")";

  return ss.str();
}


int main (int argc, char **argv) {
  using namespace std;

  typedef CGAL::Gmpz                              NT;
  typedef CGAL::Random_polynomial_generator< NT > RPG;
  typedef RPG::Poly                               Poly;
  typedef RPG::mp_bitcnt_t                        mp_bitcnt_t;

  std::string mode = "u";
  std::string of = "cgal";
  size_t n = 8;
  mp_bitcnt_t L = 16;
  unsigned long seed = time (NULL);

  if (argc > 1) mode = argv[1];
  if (argc > 2) of = argv[2];
  if (argc > 3) n = boost::lexical_cast< size_t > (std::string (argv[3]));
  if (argc > 4) L = boost::lexical_cast< mp_bitcnt_t > (std::string (argv[4]));
  if (argc > 5) seed = boost::lexical_cast< unsigned long > (std::string (argv[5]));

  map< string, string > desc;
  desc["1" ] = "cyclotomic";
  desc["u" ] = "random uniform (dense)";
  desc["mu"] = "monic uniform (dense)";
  desc["us"] = "random uniform (sparse)";
  desc["ms"] = "monic uniform (sparse)";
  desc["m" ] = "Mignotte";
  desc["m3"] = "Mignotte3";
  desc["m4"] = "Mignotte4";
  desc["m5"] = "Mignotte5";
  desc["m6"] = "Mignotte6";
  desc["m7"] = "Mignotte7";
  desc["bm"] = "\"bad\" Mignotte-like (product of two)";
  desc["pm"] = "product of n Mignotte-like of degree 16 and bitsize 8";
  desc["pM"] = "product of n Mignotte-like of degree 15 and bitsize 8";
  desc["r2"] = "resultant of two random uniform dense";
  desc["rs"] = "resultant of two random uniform sparse";
  desc["rm"] = "resultant of Mandelbrot polynomials of levels n and n-1";
  desc["rl"] = "resultant of two sparse lemniscate polynomials";
  desc["e" ] = "Taylor expansion of exp at 0 until degree n";
  desc["b" ] = "Bezier [TODO: document]";
  desc["f" ] = "Fahramand [TODO: document]";
  desc["w" ] = "Wilkinson polynomial of degree n";

  if (mode == "-h") {
    cout << "Usage: generators [mode] [of] [n] [L] [seed]" << endl
         << "       n : degree   L : target bitsize" << endl
         << "Modes:" << endl
         << "       1 : " << desc["1" ] << endl
         << "       u : " << desc["u" ] << endl
         << "       mu: " << desc["mu"] << endl
         << "       us: " << desc["us"] << endl
         << "       ms: " << desc["ms"] << endl
         << "       m : " << desc["m" ] << endl
         << "       m3: " << desc["m3"] << endl
         << "       m4: " << desc["m4"] << endl
         << "       m5: " << desc["m5"] << endl
         << "       m6: " << desc["m6"] << endl
         << "       m7: " << desc["m7"] << endl
         << "       bm: " << desc["bm"] << endl
         << "       pm: " << desc["pm"] << endl
         << "       pM: " << desc["pM"] << endl
         << "       r2: " << desc["r2"] << endl
         << "       rs: " << desc["rs"] << endl
         << "       rm: " << desc["rm"] << endl
         << "       rl: " << desc["rl"] << endl
         << "       e:  " << desc["e" ] << endl
         << "       b:  " << desc["b" ] << endl
         << "       f:  " << desc["f" ] << endl
         << "       w:  " << desc["w" ] << endl << endl
         << "       [mand,lemn,bez]" << endl << endl
         << "Output formats:" << endl
         << "       cgal       " << endl
         << "       maple      " << endl
         << "       mathematica" << endl
         << "       matlab     " << endl
         << "       maxima     " << endl
         << "       mpsolve    " << endl
         << "       pari       " << endl
         << "       scilab     " << endl;
    return 1;
  }

  RPG rpg (seed);

  Poly f;
  if      (mode == "1" ) f = rpg.unit (n);
  else if (mode == "u" ) f = rpg.uniform (n, L);
  else if (mode == "mu") f = rpg.monic_uniform (n, L);
  else if (mode == "us") f = rpg.uniform_sparse (n, L);
  else if (mode == "ms") f = rpg.monic_uniform_sparse (n, L);
  else if (mode == "m" ) f = rpg.mignotte (n, L);
  else if (mode == "m3") f = rpg.k_mignotte (n, L, 3);
  else if (mode == "m4") f = rpg.k_mignotte (n, L, 4);
  else if (mode == "m5") f = rpg.k_mignotte (n, L, 5);
  else if (mode == "m6") f = rpg.k_mignotte (n, L, 6);
  else if (mode == "m7") f = rpg.k_mignotte (n, L, 7);
  else if (mode == "bm") f = rpg.bad_mignotte (n, L);
  else if (mode == "pm") {
    f = rpg.bad_mignotte (16, 8);
    for (int i = 1; i <= n; ++i)
      f *= rpg.bad_mignotte (16, 8);
  }
  else if (mode == "pM") {
    f = rpg.bad_mignotte (15, 8);
    for (int i = 1; i <= n; ++i)
      f *= rpg.bad_mignotte (15, 8);
  }
  else if (mode == "r2") f = rpg.resultant_of_two (n, L);
  else if (mode == "rs") f = rpg.resultant_of_two_sparse (n, L);
  else if (mode == "rm") f = rpg.mandelbrot_resultant (n);
  else if (mode == "rl") f = rpg.lemniscate_resultant (n, L);
  else if (mode == "e" ) f = rpg.exp_taylor (n);
  else if (mode == "b" ) f = rpg.bezier_resultant (n, L);
  else if (mode == "f" ) f = rpg.fahramand (n, L);
  else if (mode == "w" ) f = rpg.wilkinson (n);

  else if (mode == "mand") {
    RPG::Poly_2 f = rpg.mandelbrot (n);
    CGAL::set_pretty_mode (cout);
    cerr << "# Mandelbrot lemniscate (iteration " << n << ")" << endl
         << "# y-degree:  " << f.degree() << endl
         << "# Bitsize: " << CGAL::bitsize (f) << endl;
    cout << f << endl;
    return 0;
  }

  else if (mode == "lemn") {
    RPG::Poly_2 f = rpg.lemniscate (n, L);
    CGAL::set_pretty_mode (cout);
    cerr << "# lemniscate (n = " << n << ", L = " << L << ")" << endl
         << "# y-degree:  " << f.degree() << endl
         << "# Bitsize: " << CGAL::bitsize (f) << endl;
    cout << f << endl;
    return 0;
  }

  else if (mode == "bez") {
    RPG::Poly_2 f = rpg.bezier_2 (n, L);
    CGAL::set_pretty_mode (cout);
    cerr << "# bezier_grid (n = " << n << ", L = " << L << ")" << endl
         << "# y-degree:  " << f.degree() << endl
         << "# Bitsize: " << CGAL::bitsize (f) << endl;
    cout << f << endl;
    return 0;
  }
  
  else {
    cerr << "Unknown generator (see " << argv[0] << " -h)" << endl;
    return 1;
  }
  
  const int fb = CGAL::internal::Fujiwara_upper_bound_log2() (f);

  if (of == "cgal") {

    CGAL::set_pretty_mode (cout);
    cout << "# " << desc[mode] << endl
         << "# Seed:    " << seed << endl
         << "# Degree:  " << f.degree() << "\t(target: " << n << ")" << endl
         << "# Bitsize: " << CGAL::bitsize (f) << "\t(target: " << L << ")" << endl
         << "# Fujiwara bound: \t2^" << fb << endl
         << "# " << f << endl;
    CGAL::set_ascii_mode (cout);
    cout << f << endl;

  } else if (of == "maple") {

    cout << "# " << desc[mode] << endl
         << "# Seed:    " << seed << endl
         << "# Degree:  " << f.degree() << "\t(target: " << n << ")" << endl
         << "# Bitsize: " << CGAL::bitsize (f) << "\t(target: " << L << ")" << endl
         << "# Fujiwara bound: \t2^" << fb << endl
         << endl
         << to_maple (f) << ":" << endl
         << "showtime();" << endl
         << "RootFinding[Isolate] (f);" << endl
         << "RootFinding[Analytic] (f, z, " 
         << "re=-2^(" << fb << ")..2^(" << fb << "),"
         << "im=-2^(" << fb << ")..2^(" << fb << "));" << endl;

  } else if (of == "mathematica") {

    cout << "(* " << desc[mode] << " *)" << endl
         << "(* Seed:    " << seed << " *)" << endl
         << "(* Degree:  " << f.degree() << "\t(target: " << n << ") *)" << endl
         << "(* Bitsize: " << CGAL::bitsize (f) << "\t(target: " << L << ") *)" << endl
         << "(* Fujiwara bound: \t2^" << fb << " *)" << endl
         << endl
         << to_mathematica (f) << ";" << endl
         << "Timing [NSolve [f, z]]" << endl;

  } else if (of == "matlab") {

    cout << "% " << desc[mode] << endl
         << "% Seed:    " << seed << endl
         << "% Degree:  " << f.degree() << "\t(target: " << n << ")" << endl
         << "% Bitsize: " << CGAL::bitsize (f) << "\t(target: " << L << ")" << endl
         << "% Fujiwara bound: \t2^" << fb << endl
         << endl
         << to_matlab (f) << endl
         << "tic; roots (f); toc;" << endl;

  } else if (of == "mpsolve") {

    cout << "! " << desc[mode] << endl
         << "! Seed:    " << seed << endl
         << "! Degree:  " << f.degree() << "\t(target: " << n << ")" << endl
         << "! Bitsize: " << CGAL::bitsize (f) << "\t(target: " << L << ")" << endl
         << "! Fujiwara bound: \t2^" << fb << endl
         << endl
         << to_mpsolve (f) << endl;

  } else if (of == "pari") {

    cout << "\\\\ " << desc[mode] << endl
         << "\\\\ Seed:    " << seed << endl
         << "\\\\ Degree:  " << f.degree() << "\t(target: " << n << ")" << endl
         << "\\\\ Bitsize: " << CGAL::bitsize (f) << "\t(target: " << L << ")" << endl
         << "\\\\ Fujiwara bound: \t2^" << fb << endl
         << endl
         << "#" << endl
         << to_pari (f) << ";" << endl
         << "polroots(f,1);" << endl
         << "polroots(f);" << endl;

  }

  return 0;
}
