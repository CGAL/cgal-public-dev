#define CGAL_PROFILE 1
//#define USE_LEDA_AK 1

#include "read_poly.h"
#include "formatters.h"
#include <CGAL/Absolute_void.h>
#include <CGAL/Arcavoid_root_isolator.h>

#include <vector>
#include <iostream>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include <boost/tuple/tuple.hpp>

#ifdef CGAL_USE_LEDA
#include <CGAL/LEDA_arithmetic_kernel.h>
#endif
#ifdef CGAL_USE_GMP
#include <CGAL/GMP_arithmetic_kernel.h>
#endif

#include <CGAL/Timer.h>

#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_quadratic_refinement_rep_bfi.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_rndl_tree_traits.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_coefficient_kernel.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_coefficient_kernel_at_alpha.h>

int main (int argc, char **argv) {
  using namespace std;
  using boost::format;
  
  vector< string > args (argv, argv+argc);
  istream *in_ptr = &cin;
  if (args.size() > 1 && args[1] != "--")
    in_ptr = new ifstream (args[1].c_str());
  
  int deg_gcd = -1;
  if (args.size() > 2)
    deg_gcd = boost::lexical_cast< int > (args[2]);

  int nr_real_roots = -1;
  if (args.size() > 3)
    nr_real_roots = boost::lexical_cast< int > (args[3]);

  bool complex = true;
  if (nr_real_roots > -1 && deg_gcd < 0)
    complex = false;
  bool real = !complex;

  if (args.size() > 4) {
    complex = real = false;
    if (args[4].find ('c') != string::npos)
      complex = true;
    if (args[4].find ('r') != string::npos)
      real = true;
  }

#if USE_LEDA_AK
  typedef CGAL::LEDA_arithmetic_kernel AK;
#else
  typedef CGAL::GMP_arithmetic_kernel AK;
#endif
  typedef AK::Integer I;
  typedef AK::Rational Q;
  typedef AK::Bigfloat_interval BFI;
  typedef AK::Bigfloat BF;
  typedef CGAL::Polynomial_type_generator< I, 1>::Type IPoly_1;
  typedef CGAL::Polynomial_traits_d< IPoly_1 > IPT_1;
  typedef CGAL::Polynomial_type_generator< Q, 1>::Type QPoly_1;
  typedef CGAL::Polynomial_traits_d< QPoly_1 > QPT_1;
  typedef CGAL::Polynomial_type_generator< Q, 2>::Type QPoly_2;
  typedef CGAL::Polynomial_type_generator< BFI, 1>::Type BFIPoly_1;
  typedef CGAL::Polynomial_type_generator< BFI, 2>::Type BFIPoly_2;
  typedef CGAL::Polynomial_type_generator< BF, 1>::Type BFPoly_1;
  typedef CGAL::Polynomial_type_generator< BF, 2>::Type BFPoly_2;
  
#if false && USE_LEDA_AK
  typedef I          CT;
  typedef IPoly_1    Poly;
  typedef IPT_1      PT;
#else
  typedef Q          CT;
  typedef QPoly_1    Poly;
  typedef QPT_1      PT;
#endif
  CGAL::Polynomial_parser_d< Poly > parser;

  const Poly x = PT::Shift() (Poly (1), 1);
  Poly f = x - I(1);
  read_poly (*in_ptr, f);

  if (f.degree() > 0) {
    CT max = CGAL::abs (f[0]);
    for (int i = 0; i <= f.degree(); ++i)
      if (CGAL::abs (f[i]) > max)
        max = CGAL::abs (f[i]);
    f /= max;
  }
    
  const int n = f.degree();
  
  typedef CGAL::internal::Bitstream_coefficient_kernel< CT > BCK;
  BCK bck;
  
  CGAL::set_pretty_mode (cerr);
  // cerr << "Solving f = " << f << endl;

  typedef CGAL::internal::Arcavoid_list< BCK > Active_interval_set;
  // Active_interval_set set (bck, f);
  typedef Active_interval_set::Cluster_iterator Cluster_iterator;
  typedef Active_interval_set::Cluster_range Cluster_range;
  
  typedef CGAL::Arcavoid< BCK, CGAL::Arcavoid_real_root_isolator_tag > Real_root_isolator;
  typedef CGAL::Arcavoid< BCK, CGAL::Arcavoid_complex_root_isolator_tag > Complex_root_isolator;

  if (complex) {
    CGAL::Timer t_isol;
    t_isol.start();
    Complex_root_isolator isolator
      = (deg_gcd <= 0
         ? Complex_root_isolator (CGAL::Square_free_arcavoid_tag(), f)
         : Complex_root_isolator (CGAL::K_arcavoid_tag(), f, deg_gcd));
    t_isol.stop();
      
    cerr << "================================================================" << endl
         << "Complex roots of f =" << endl
      // << " " << f << endl
         << "================================================================" << endl;

    for (int i = 0; i < isolator.number_of_complex_roots(); ++i) {
      if (isolator[i].multiplicity() == 1)
        cout << "complex root " << i << ": (Mult: 1) " << *(isolator[i].cit->begin()) << endl;
      else
        cout << "complex multiple root " << i << endl << isolator[i] << endl;
    }
    cerr << "================================================================" << endl;

    for (int i = 0; i < isolator.number_of_complex_roots(); ++i) {
      for (int j = i+1; j < isolator.number_of_complex_roots(); ++j) {
        if ((isolator[i].center() - isolator[j].center()).abs()
            <= isolator[i].radius() + isolator[j].radius())
          cerr << "[ERROR]: Root regions " << i << " and " << j << " overlap" << endl;
      }
    }

    for (int i = 0; i < isolator.number_of_complex_roots(); ++i) {
      // touch real axis?
      if (CGAL::abs (isolator[i].center().imag()) > isolator[i].radius())
        continue;
    
      // check if conjugate root exists
      bool conj_exists = false;
      for (int j = 0; j < isolator.number_of_complex_roots(); ++j) {
        if (i == j)
          continue;
      
        if ((isolator[i].center().conj() - isolator[j].center()).abs()
            <= isolator[i].radius() + isolator[j].radius()) {
          conj_exists = true;
          cerr << "Roots " << i << " and " << j << " are complex conjugates?!" << endl;
          break;
        }
      }

      if (conj_exists)
        continue;

      // check if sign change is in interval
#if USE_LEDA_AK
      const Q a = (isolator[i].center().real() - isolator[i].radius()).to_rational();
      const Q b = (isolator[i].center().real() + isolator[i].radius()).to_rational();
#else
      const Q a = isolator[i].center().real() - isolator[i].radius();
      const Q b = isolator[i].center().real() + isolator[i].radius();
#endif
      const Q fa = f.evaluate (a);
      const Q fb = f.evaluate (b);
      if (! (CGAL::sign (fa) != CGAL::sign (fb)
             || (a == b && CGAL::is_zero (fa))))
        cerr << "Sign-change test failed for (real) root interval " << i << endl;
    }

    cerr << "Validity check finished" << endl
         << "================================================================" << endl;

    cout << "TIME FOR ISOLATION: " << t_isol.time() << endl;
  }

  if (real) {
    CGAL::Timer t_isol;
    t_isol.start();
    Real_root_isolator isolator
      = (deg_gcd <= 0
         ? (nr_real_roots < 0
            ? Real_root_isolator (CGAL::Square_free_arcavoid_tag(), f)
            : Real_root_isolator (CGAL::M_arcavoid_tag(), f, nr_real_roots))
         : (nr_real_roots < 0
            ? Real_root_isolator (CGAL::K_arcavoid_tag(), f, deg_gcd)
            : Real_root_isolator (CGAL::M_k_arcavoid_tag(), f, nr_real_roots, deg_gcd)));
    t_isol.stop();

    cerr << "================================================================" << endl
         << "Real roots of f =" << endl
         << " " << f << endl
         << "================================================================" << endl;
    
    for (int i = 0; i < isolator.number_of_real_roots(); ++i)
      cout
        << "real root in interval [" << CGAL::to_double (isolator.left_bound(i))
        << "|" << CGAL::to_double (isolator.right_bound(i)) << "]" << endl;
    cerr << "================================================================" << endl;
    
    for (int i = 1; i < isolator.number_of_real_roots(); ++i) {
      if (! (isolator.right_bound(i-1) < isolator.left_bound(i))) {
        cerr << "root intervals " << i-1 << " and " << i << " intersect" << endl;
      }
    }
    for (int i = 0; i < isolator.number_of_real_roots(); ++i) {
      const Q a = isolator.left_bound (i);
      const Q b = isolator.right_bound (i);
      const Q fa = f.evaluate (a);
      const Q fb = f.evaluate (b);
      if (! (CGAL::sign (fa) != CGAL::sign (fb)
             || (a == b && CGAL::is_zero (fa))))
        cerr << "Sign-change test failed for root interval " << i << endl;
    }
    cerr << "Validity check finished" << endl
         << "================================================================" << endl;

    cout << "Real_isolation_time: " << t_isol.time() << endl;
  }
  
  // Clean up
  mpfr_free_cache();
  if (in_ptr && in_ptr != &cin)
    delete in_ptr;
  
  return 0;
}
