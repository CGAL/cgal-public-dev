#include "formatters.h"
#include <CGAL/Arcavoid.h>
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

  {
    typedef CGAL::GMP_arithmetic_kernel AK;
    typedef AK::Rational Q;
    typedef AK::Bigfloat_interval BFI;
    typedef AK::Bigfloat BF;
    typedef CGAL::Polynomial_type_generator< Q, 1>::Type QPoly_1;
    typedef CGAL::Polynomial_traits_d< QPoly_1 > QPT_1;
    typedef CGAL::Polynomial_type_generator< Q, 2>::Type QPoly_2;
    typedef CGAL::Polynomial_type_generator< BFI, 1>::Type BFIPoly_1;
    typedef CGAL::Polynomial_type_generator< BFI, 2>::Type BFIPoly_2;
    typedef CGAL::Polynomial_type_generator< BF, 1>::Type BFPoly_1;
    typedef CGAL::Polynomial_type_generator< BF, 2>::Type BFPoly_2;
    Q a (1,5);
    Q Qfcs[7] = {
      /*      Q(1),
      2*a,
      a*a,
      Q(0),
      Q(0),
      Q(0),
      Q(1), */
      2, 1, -201, -99, 97, -100, 300
    };
    QPoly_1 Qf (Qfcs, Qfcs+7);
    std::vector< QPoly_1 > Qfs;
    Qfs.push_back (QPoly_1 (Qfcs, Qfcs+3));
    Qfs.push_back (QPoly_1 (Qfcs, Qfcs+3));
    Qfs.push_back (Qf);
    
    QPoly_1 Qx = QPT_1::Shift() (QPoly_1 (1), 1);
    QPoly_1 Qx50 = QPT_1::Shift() (QPoly_1 (1), 50);
    QPoly_1 Qm = Qx50 + 11449 * Qx*Qx - 214 * Qx + 1;

    Q M = 1234567;
    M *= M;
    M *= M;
    M *= M;
    M += 1;

    Qm = CGAL::ipower (Qx, 31) + M*M*Qx*Qx + 2*M* Qx + 1;

    // Qm = QPoly_1 (1);
    // for (int i = -10; i < 10; ++i)
    //   Qm *= (M+i)*Qx-1;
    // Qm += Qx50*Qx50*Qx;

    // Qm = Qm * Qx*Qx;
    //Qm = Qx*Qx*Qx*Qx - Qx*Qx*Qx;

    QPoly_2 Qg (Qfs.begin(), Qfs.end());

    typedef CGAL::internal::Bitstream_coefficient_kernel< Q > QBCK;
    QBCK qbck;
    
    CGAL::set_precision (BFI(), 160);
    
    std::cerr << Qf << std::endl;
    std::vector< BFI > fcoeffs;
    for (int i = 0; i <= Qf.degree(); ++i)
      fcoeffs.push_back (qbck.convert_to_bfi_object() (Qf[i]));
    std::cerr << BFIPoly_1 (fcoeffs.begin(), fcoeffs.end()) << std::endl;

    {
      typedef CGAL::internal::Arcavoid_list< QBCK > Active_interval_set;
      Active_interval_set set (qbck, Qm);
      typedef Active_interval_set::Cluster_iterator Cluster_iterator;
      typedef Active_interval_set::Cluster_range Cluster_range;

      typedef CGAL::Arcavoid< QBCK, CGAL::Arcavoid_real_root_isolator_tag > Complex_root_isolator;
      Complex_root_isolator isolator (qbck, Qm);
      
      std::cerr << "================================================================" << std::endl;
      std::cerr << "Complex roots of " << Qm << std::endl;
      std::cerr << "================================================================" << std::endl;
      for (int i = 0; i < isolator.number_of_real_roots(); ++i) {
        std::cout //<< isolator[i]
          << "real root in interval [" << CGAL::to_double (isolator.left_bound(i))
          << "|" << CGAL::to_double (isolator.right_bound(i)) << "]" << std::endl;
        std::cerr << "================================================================" << std::endl;
      }
      for (int i = 1; i < isolator.number_of_real_roots(); ++i) {
        if (! (isolator.right_bound(i-1) < isolator.left_bound(i))) {
          std::cerr << "root intervals " << i-1 << " and " << i << " intersect" << std::endl;
        }
      }
      for (int i = 0; i < isolator.number_of_real_roots(); ++i) {
        Q a = isolator.left_bound (i);
        Q b = isolator.right_bound (i);
        Q fa = Qm.evaluate (a);
        Q fb = Qm.evaluate (b);
        if (! (CGAL::sign (fa) != CGAL::sign (fb)
               || (a == b && CGAL::is_zero (fa))))
          std::cerr << "Sign-change test failed for root interval " << i << std::endl;
        // TODO: sign-change test
      }
      std::cerr << "Validity check finished" << std::endl
                << "================================================================" << std::endl;

      /*
      Cluster_range main = set.cluster_range().first;

      Cluster_range range = set.subdivide_cluster (main);
      for (Cluster_iterator it = range.first; it != range.second; ++it)
        std::cerr << *it << std::endl;

      Cluster_iterator mign;

      for (int i = 0; i < 5; ++i) {
        mign = range.second;
        std::advance (mign, -1);
        range = set.subdivide_cluster (mign);
        for (Cluster_iterator it = range.first; it != range.second; ++it)
          std::cerr << *it << std::endl;
      }
      */
    }
    
    typedef CGAL::internal::Algebraic_real_quadratic_refinement_rep_bfi < Q, Q > Rep_class;
    typedef CGAL::internal::Bitstream_descartes<
    CGAL::internal::Bitstream_descartes_rndl_tree_traits<
    CGAL::internal::Bitstream_coefficient_kernel<Q> > >
      Isolator;
    typedef CGAL::Algebraic_kernel_d_1<Q,Q,Rep_class, Isolator> Algebraic_kernel_d_1;
    typedef Algebraic_kernel_d_1::Algebraic_real_1 Algebraic_real_1;
    
    typedef CGAL::internal::Bitstream_coefficient_kernel_at_alpha< Algebraic_kernel_d_1 > QBCKalpha;
    Algebraic_real_1 alpha = 3;
    Algebraic_kernel_d_1 kernel;
    QBCKalpha qbcka (&kernel, alpha);

    CGAL::set_pretty_mode (std::cerr);
    std::cerr << Qg << std::endl;
    std::vector< BF > gcoeffs;
    for (int i = 0; i <= Qg.degree(); ++i)
      gcoeffs.push_back (CGAL::median (qbcka.convert_to_bfi_object() (Qg[i])));
    std::cerr << BFPoly_1 (gcoeffs.begin(), gcoeffs.end()) << std::endl;

    {
      typedef CGAL::internal::Arcavoid_list< QBCKalpha > Solver;
      Solver solver (qbcka, Qg);
    }
    const int n = Qf.degree();
  }

#ifdef CGAL_ARCAVOID_USE_OLD_ARCAVOID
  typedef CGAL::GMP_arithmetic_kernel AK;
  typedef AK::Integer Integer;
  typedef AK::Rational Rational;
  typedef AK::Bigfloat Bigfloat;
  typedef CGAL::Polynomial_type_generator< Integer, 1 >::Type IPoly_1;

  typedef CGAL::Timer Timer;

  // READ POLYNOMIAL
  vector< string > args (argv, argv+argc);
  istream *in_ptr = &cin;
  if (args.size() > 1 && args.back() != "--")
    in_ptr = new ifstream (args.back().c_str());

  string line;
  while (line.empty() || line[0] == '#')
    getline (*in_ptr, line);

  IPoly_1 f_int = boost::lexical_cast< IPoly_1 > (line);

  typedef Bigfloat F;
  typedef CGAL::Polynomial_type_generator< F, 1>::Type FPoly;

  F A = 123456789;
  F Z = A - A;
  vector< F > Fcs;
  Fcs.push_back (F (1));
  Fcs.push_back (A * 2);
  Fcs.push_back (A * A);
  for (int i = 3; i < 30; ++i)
    Fcs.push_back (Z);
  Fcs.push_back (F (1));
  FPoly  Ff (Fcs.begin(), Fcs.end());

  typedef Rational Q;
  typedef CGAL::Polynomial_type_generator< Q, 1>::Type QPoly;
  Q a (1,10000001);
  Q Qfcs[7] = {
    Q(1),
    2*a*a*a,
    a*a*a*a*a*a,
    Q(0),
    Q(0),
    Q(0),
    Q(1),
  };
  QPoly Qf (Qfcs, Qfcs+7);

  // typedef leda_real R;
  // typedef CGAL::Polynomial_type_generator< R, 1>::Type RPoly;
  // R Rfcs1[3] = {
  //   leda::sqrt (leda_real (2)),
  //   - 2 * leda::root (leda_real (2), 4),
  //   R(1.),
  // };
  // R Rfcs2[3] = {
  //   leda::sqrt (leda_real (3)),
  //   - 2 * leda::root (leda_real (3), 4),
  //   R(1.),
  // };
  // RPoly Rf1 (Rfcs1, Rfcs1+3);
  // RPoly Rf2 (Rfcs2, Rfcs2+3);
  // RPoly Rf = Rf1 * Rf2;
  
#ifdef ARCAVOID
  // typedef CGAL::Arcavoid< AK > Solver;

  Timer t;
  t.start();

  typedef CGAL::Arcavoid< IPoly_1 > Solver;

  bool all = false;
  if (args.size() > 1)
    if (args[1] == "ISOLATE_ALL") all = true;

  Solver solver (f_int,
                 all ? Solver::ISOLATE_ALL : Solver::ISOLATE_REAL,
                 1, 1);
  const int n = f_int.degree();

  /*
  typedef CGAL::Arcavoid< FPoly > Solver;
  Solver solver (Ff);
  const int n = Ff.degree();
  */

  /*
  typedef CGAL::Arcavoid< RPoly > Solver;
  Solver solver (Rf, Solver::ISOLATE_REAL);
  const int n = Rf.degree();
  */

  // typedef CGAL::Arcavoid< QPoly > Solver;
  // Solver solver (Qf);
  // const int n = Qf.degree();
  
  // solver.solve();

  t.stop();

  boost::format def_form ("%# .20f %|28t|");

  for (int i = 0; i < n; ++i) {
    if (CGAL::is_zero (solver.rad[i]))
      cerr << "Rad[" << i << "] = 0 !!!" << endl;
    if (solver[i].real() != solver.cntr[i].real())
      cerr << "z[" << i << "].real() = cntr[" << i << "].real() !!!" << endl;
    if (solver[i].imag() != solver.cntr[i].imag())
      cerr << "z[" << i << "].imag() = cntr[" << i << "].imag() !!!" << endl;

    cout << def_form % CGAL::to_double (solver[i].real());
#ifndef USE_DOUBLE_WITH_EXPONENT
    if (solver.tch_real[i])
      cout << "        0                   ";
    else
#endif // USE_DOUBLE_WITH_EXPONENT
      cout << def_form % CGAL::to_double (solver[i].imag());

    cout << def_form % CGAL::to_double (solver.rad[i]);
    cout << "( ";
    cout << (def_form % CGAL::to_double (solver.cntr[i].real()));
    cout << " , ";
    cout << (def_form % CGAL::to_double (solver.cntr[i].imag()));
    cout << " )"
         << std::endl;

    // cout << format ("%|# 24.16| %|# 24.16| %|# 24.16|")
    //   % solver[i].real()
    //   % solver[i].imag()
    //   % solver.rad[i];

    // cout << solver[i].real() << "\t";
    // // if (solver.tch_real[i])
    // //   cout << "\t0\t\t";
    // // else
    //   cout << solver[i].imag() << "\t";
    // cout << solver.rad[i] << "\t"
    //           << solver.cntr[i] << endl;
  }

  /*
  // SVG
  for (int i = 0; i < n; ++i) {
    cerr << "<circle cx=\"" << solver[i].real() << "\""
         << "cy=\"" << solver[i].imag() << "\""
         << "r=\"" << solver.rad[i] << "\""
         << "stroke=\"black\" fill=\"red\"/>" << endl;
    cerr << "<circle cx=\"" << solver[i].real() << "\""
         << "cy=\"" << solver[i].imag() << "\""
         << "r=\"" << solver.rad[i] / (n*n) << "\""
         << "stroke=\"black\" fill=\"black\"/>" << endl;
  }
  */

  cerr << "Time: " << t.time() << " sec" << endl;

  // leda_integer Iasdf = 1234567;
  // Iasdf *= Iasdf * Iasdf;
  // Iasdf *= Iasdf * Iasdf;
  // leda_bigfloat Fasdf (Iasdf, 1234);
  // cout << Fasdf.get_exponent() << endl
  //           << Fasdf.get_significant() << endl
  //           << Fasdf.get_significant_length() << endl
  //           << Fasdf.get_effective_significant_length() << endl;
#endif // ARCAVOID

  // CLEAN UP
  if (in_ptr && in_ptr != &cin)
    delete in_ptr;

#endif // CGAL_ARCAVOID_USE_OLD_ARCAVOID
  return 0;
}
