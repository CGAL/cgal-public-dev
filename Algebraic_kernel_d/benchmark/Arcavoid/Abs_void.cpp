#include "read_poly.h"
#include "formatters.h"
#include <fenv.h>
#include <csignal>

#include <CGAL/Arcavoid.h>

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



const double trunc_prec_RND_TO_INF (const double &x, int prec_loss) {
  long x_lng = reinterpret_cast< const long & > (x);
  x_lng &= (-1L << prec_loss);
  unsigned long x_ulng = reinterpret_cast< const unsigned long & > (x_lng);
  x_ulng += (1 << prec_loss);
  return reinterpret_cast< const double & > (x_ulng);
}

class SomeException {};
void handle_fpe (int signal) {
  switch (signal) {
  case (FE_DIVBYZERO):
    std::clog << "FE_DIVBYZERO" << std::endl;
    break;
  case (FE_INVALID):
    std::clog << "INVALID" << std::endl;
    break;
  case (FE_OVERFLOW):
    std::clog << "OVERFLOW" << std::endl;
    break;
  default:
    std::clog << "UNKNOWN FP EXCEPTION" << std::endl;
  }
  throw SomeException();
}

int main (int argc, char **argv) {
  using namespace std;
  using boost::format;

  typedef CGAL::GMP_arithmetic_kernel AK;
  typedef AK::Integer  ZZ;
  typedef AK::Rational QQ;
  typedef AK::Bigfloat BF;
  typedef AK::Bigfloat_interval BFI;

  typedef BF  RR;
  typedef BFI RRI;
  typedef CGAL::Cartesian_complex< RR > CC;
  typedef CGAL::Cartesian_complex< RRI > CCI;

  typedef double      DD;
  typedef CGAL::Interval_nt< false > DDI;
  typedef DDI::Protector Protect_fpu_rounding;
  typedef CGAL::Cartesian_complex< DD > DDCC;
  typedef CGAL::Cartesian_complex< DDI > DDCCI;

  typedef CGAL::Polynomial_type_generator< DD,  1 >::Type DDPoly_1;
  typedef CGAL::Polynomial_type_generator< DDI, 1 >::Type DDIPoly_1;
  typedef CGAL::Polynomial_type_generator< DDCC,  1 >::Type DDCCPoly_1;
  typedef CGAL::Polynomial_type_generator< DDCCI, 1 >::Type DDCCIPoly_1;
  typedef CGAL::Polynomial_type_generator< ZZ,  1 >::Type ZZPoly_1;
  typedef CGAL::Polynomial_type_generator< QQ,  1 >::Type QQPoly_1;
  typedef CGAL::Polynomial_type_generator< RR,  1 >::Type RRPoly_1;
  typedef CGAL::Polynomial_type_generator< RRI, 1 >::Type RRIPoly_1;
  typedef CGAL::Polynomial_type_generator< CC,  1 >::Type CCPoly_1;
  typedef CGAL::Polynomial_type_generator< CCI, 1 >::Type CCIPoly_1;
  typedef CGAL::Polynomial_traits_d< DDPoly_1  > DDPT_1;
  typedef CGAL::Polynomial_traits_d< DDIPoly_1 > DDIPT_1;
  typedef CGAL::Polynomial_traits_d< DDCCPoly_1  > DDCCPT_1;
  typedef CGAL::Polynomial_traits_d< DDCCIPoly_1 > DDCCIPT_1;
  typedef CGAL::Polynomial_traits_d< ZZPoly_1  > ZZPT_1;
  typedef CGAL::Polynomial_traits_d< QQPoly_1  > QQPT_1;
  typedef CGAL::Polynomial_traits_d< RRPoly_1  > RRPT_1;
  typedef CGAL::Polynomial_traits_d< RRIPoly_1 > RRIPT_1;
  typedef CGAL::Polynomial_traits_d< CCPoly_1  > CCPT_1;
  typedef CGAL::Polynomial_traits_d< CCIPoly_1 > CCIPT_1;

  int max_count = 100;
  int ref_count = 0;

  vector< string > args (argv, argv+argc);
  istream *in_ptr = &cin;
  if (args.size() > 1 && args.back() != "--")
    in_ptr = new ifstream (args.back().c_str());

  if (args.size() > 2)
    max_count = boost::lexical_cast< int > (args[1]);
  if (args.size() > 3)
    ref_count = boost::lexical_cast< int > (args[2]);

  CGAL::set_pretty_mode (clog);

  // Rational Mignotte polynomial
  {
    typedef QQ          CT;
    typedef QQPoly_1    Poly;
    typedef QQPT_1      PT;
    CGAL::Polynomial_parser_d< Poly > parser;
    
    const Poly x = PT::Shift() (Poly (1), 1);

    CT a = 1;
    a *= 99999;
    
    Poly f_exact = CGAL::ipower (x, 25) + CGAL::ipower (a*x + 1, 3);
    read_poly (*in_ptr, f_exact);

    clog << "Solving f_exact = " << f_exact << endl
         << "  w/ coeff type " << type_name (CT()) << endl;

    const int n = f_exact.degree();
    long e_min, e_max;
    e_min = e_max = CGAL::abs_ilog2 (CGAL::to_double (f_exact[n]));
    for (int i = 0; i < n; ++i) {
      if (! CGAL::is_zero (f_exact[i])) {
        long e = CGAL::abs_ilog2 (CGAL::to_double (f_exact[i]));
        e_min = CGAL::min (e_min, e);
        e_max = CGAL::max (e_max, e);
      }
    }
    const long e_norm = (e_min + e_max) / 2;

    clog << "Renormalizing by 2^" << e_norm << endl;
    const CT normalizer = CGAL::ipower (CT(2), CGAL::abs (e_norm));
    if (e_norm < 0)
      f_exact *= normalizer;
    else
      f_exact /= normalizer;
    clog << "  to " << f_exact << endl;

    DDPoly_1 f = CGAL::to_double (f_exact);
    DDPoly_1 df = CGAL::differentiate (f);


    std::vector< DDCC > z;
    /* DEFAULT INITIALIZATION ON UNIT CIRCLE */

    /*
    for (int i = 0; i < n; ++i)
      z.push_back (DDCC (std::cos (M_PI * (2*i+1) / n), std::sin (M_PI * (2*i+1) / n)));
    */

    /* CONVEX HULL MPSOLVE INITIALIZATION */

    typedef CGAL::Simple_cartesian< double > K;
    typedef K::Point_2 Point_2;
    std::vector< Point_2 > points, ch;
    for (int i = 0; i <= n; ++i)
      if (! CGAL::is_zero (f[i]))
        points.push_back (Point_2 (i, CGAL::abs_ilog2 (f[i])));
    CGAL::upper_hull_points_2 (points.begin(), points.end(), std::back_inserter (ch));
    ch.push_back (points.front());
    std::reverse (ch.begin(), ch.end());

    const double offset = (double)(std::rand()) / RAND_MAX * 2 * M_PI;

    for (int i = 1; i < ch.size(); ++i) {
      int k = (int)(CGAL::to_double (ch[i-1].x()));
      int K = (int)(CGAL::to_double (ch[i].x()));
      DD u = std::abs (f[k] / f[K]);
      u = std::pow (u, 1./(K-k));

      for (int j = 0; j < K-k; ++j) {
        const double angle = 2. * M_PI * (((double)j)/(K-k)  + ((double)i)/n) + offset;
        z.push_back (u * DDCC (std::cos (angle), std::sin (angle)));
        // const DDCC z = u * DDCC (std::cos (angle), std::sin (angle));
        // clsts.front().discs.push_front (Approximation (z, *this));
      }
    }

    /* END OF INITIALIZATION */

    Protect_fpu_rounding protector;
    feenableexcept (FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    signal (SIGFPE, handle_fpe);
    DDCCIPoly_1 fi (boost::make_transform_iterator (f_exact.begin(), CGAL::To_interval< CT >()),
                    boost::make_transform_iterator (f_exact.end(), CGAL::To_interval< CT >()));

    clog << "Solving f = " << f << endl
         << "       fi = " << fi << endl
         << "  w/ coeff type " << type_name (DD()) << endl;

    // widen intervals
    std::vector< DDCCI > fics (fi.begin(), fi.end());
    /*
    unsigned long fmin_l = n;
    double fmin = reinterpret_cast< double & > (fmin_l);
    fmin = FLT_MIN;
    fmin = ldexp (1., min (min (e_min, -e_max), -128L));
    clog << "fmin: " << fmin << endl;
    for (int i = 0; i <= n; ++i) {
      // DD l = CGAL::lower (fics[i]);
      // DD r = CGAL::upper (fics[i]);
      // volatile float fl = l;
      // volatile float fr = r;
      // volatile double dl = CGAL::abs (l - fl);
      // volatile double dr = CGAL::abs (r - fr);
      // l = l - dl;
      // r = r + dr;
      // fics[i] = DDI (l, r);

      DD lr = CGAL::lower (fics[i].real());
      DD rr = CGAL::upper (fics[i].real());
      DD li = CGAL::lower (fics[i].imag());
      DD ri = CGAL::upper (fics[i].imag());

      volatile float flr = lr;
      volatile float frr = rr;
      volatile float fli = li;
      volatile float fri = ri;

      volatile double dlr = std::max (CGAL::abs (lr - flr), fmin);
      volatile double drr = std::max (CGAL::abs (rr - frr), fmin);
      volatile double dli = std::max (CGAL::abs (li - fli), fmin);
      volatile double dri = std::max (CGAL::abs (ri - fri), fmin);

      lr = lr - dlr;
      rr = rr + drr;
      li = li - dli;
      ri = ri + dri;

      fics[i] = DDCCI (DDI (lr, rr), DDI (li, ri));
    }
    */

    const DD delta = std::ldexp (.5, -49); // 2^-50;
    const DDCCI delta_iv = DDCCI (DDI (-delta, delta), DDI (-delta, delta));
    for (int i = 0; i <= n; ++i)
      fics[i] += f[i] * delta_iv;

    fi = DDCCIPoly_1 (fics.begin(), fics.end());
    clog << "   widened to: " << fi << endl;

    std::vector< bool > in_nbh (n, false);
    bool all_in_nbh = false;

    int count = 0;
    while (! all_in_nbh && count < max_count) {
      // clog << "Approximations after " << count++ << " iterations" << endl;
      // for (int i = 0; i < n; ++i)
      //   clog << "  " << (in_nbh[i] ? '+' : '-') << ' ' << z[i] << " (fi[z] : " << fi.evaluate (z[i]) << ")"
      //        << "\t" << "ratio: " << (CGAL::width (fi.evaluate (z[i]).real()) / CGAL::width (fi.evaluate (z[i]).imag()))
      //        << endl;

      // Aberth
      for (int i = 0; i < n; ++i) {
        //if (in_nbh[i]) continue;

        DDCC abcorr = 0;
        for (int j = 0; j < n; ++j) {
          if (i == j) continue;
          abcorr += (z[i] - z[j]).reciprocal();
        }
        const DDCC f_df = f.evaluate (z[i]) / df.evaluate (z[i]);
        const DDCC corr = f_df / (DDCC (1) - f_df * abcorr);
        z[i] -= corr;

        DDCCI ziv = z[i];
        DDCCI iv = fi.evaluate (z[i]);
        in_nbh[i] = CGAL::zero_in (iv.norm_1());
      }

      all_in_nbh = true;
      for (int i = 0; i < n; ++i)
        all_in_nbh &= in_nbh[i];

      ++count;

      {
        DDCCIPoly_1 tight = DDCCIPoly_1 (boost::make_transform_iterator (f_exact.begin(), CGAL::To_interval< CT >()),
                                         boost::make_transform_iterator (f_exact.end(), CGAL::To_interval< CT >()));
        std::vector< DD > r (n);
        for (int i = 0; i < n; ++i) {
          DDI deniv = 1.;
          for (int j = 0; j < n; ++j)
            if (i != j)
              deniv *= CGAL::abs (DDCCI(z[i]) - DDCCI(z[j]));
          DDI riv = DDI(n) * CGAL::abs (tight.evaluate (z[i])) / tight[n].real();
          riv /= deniv;
          r[i] = CGAL::upper (riv);
        }

        // std::ofstream plotfile ("tmp.plot", std::ios_base::out | std::ios_base::trunc);
        // // plotfile << "set terminal pdf" << endl;
        // // plotfile << "set output 'tmp.pdf'" << endl;
        // plotfile.precision (80);
        // for (int i = 0; i < n; ++i)
        //   plotfile << "set object " << (i+1) << " ellipse center " << z[i].real() << " , " << z[i].imag()
        //            << " size " << r[i] << "," << r[i] << " front fs empty border 1" << endl;
        // plotfile << "plot 'tmp.plot' using 6:8 with points lc 1;" << endl
        //          << "replot 'mps.out' using 1:2 with points;" << endl;
        // plotfile << "pause mouse key;" << endl;
      }
      // int ret = std::system ("gnuplot tmp.plot");

      // if (count == 1)
      //   ret = std::system ("mv tmp.pdf plot.pdf");
      // else {
      //   ret = std::system ("pdftk plot.pdf tmp.pdf cat output plot-new.pdf");
      //   ret = std::system ("mv plot-new.pdf plot.pdf");
      // }
    }

    clog << "count: " << count << endl;

    clog << "AFTER INITIAL ITERATIONS:" << endl;
    for (int i = 0; i < n; ++i)
      clog << "  " << (in_nbh[i] ? '+' : '-') << ' ' << z[i] << " (fi[z] : " << fi.evaluate (z[i]) << ")" << endl;

    for (count = 0; count < ref_count; ++count) {
      for (int i = 0; i < n; ++i) {
        DDCC abcorr = 0;
        for (int j = 0; j < n; ++j)
          if (i != j)
            abcorr += (z[i] - z[j]).reciprocal();
        const DDCC f_df = f.evaluate (z[i]) / df.evaluate (z[i]);
        const DDCC corr = f_df / (DDCC (1) - f_df * abcorr);
        z[i] -= corr;

        DDCCI iv = fi.evaluate (z[i]);
        in_nbh[i] = CGAL::zero_in (iv.real()) && CGAL::zero_in (iv.imag());
      }
    }

    fi = DDCCIPoly_1 (boost::make_transform_iterator (f_exact.begin(), CGAL::To_interval< CT >()),
                      boost::make_transform_iterator (f_exact.end(), CGAL::To_interval< CT >()));
    std::vector< DD > r (n);
    for (int i = 0; i < n; ++i) {
      DDI deniv = 1.;
      for (int j = 0; j < n; ++j)
        if (i != j)
          deniv *= CGAL::abs (DDCCI(z[i]) - DDCCI(z[j]));
      DDI riv = DDI(n) * CGAL::abs (fi.evaluate (z[i])) / fi[n].real();
      riv /= deniv;
      r[i] = CGAL::upper (riv);
    }

    clog << "refinements: " << count << endl;

    clog << "AFTER REFINEMENTS:" << endl;
    for (int i = 0; i < n; ++i)
      clog << "  " << (in_nbh[i] ? '+' : '-') << ' ' << z[i] << " (fi[z] : " << fi.evaluate (z[i]) << ")" << endl;

    clog << "Solutions:" << endl;
    cout.precision (80);
    for (int i = 0; i < n; ++i)
    //   clog << "  " << z[i].real() << "\t" << z[i].imag() << "\t" << r[i] << endl;
      cout << "set object " << (i+1) << " ellipse center " << z[i].real() << " , " << z[i].imag()
           << " size " << r[i] << "," << r[i] << " front fs empty border 1" << endl;
    cout << "plot 'abs.out' using 6:8 with points lc 1;" << endl
         << "replot 'mps.out' using 1:2 with points;" << endl;
  }

  // CLEAN UP
  if (in_ptr && in_ptr != &cin)
    delete in_ptr;

  // RR a = RR(1) / RR(9);
  // RRPoly_1 x = RRPT_1::Shift() (RRPoly_1 (1), 1);
  // RRPoly_1 f = CGAL::ipower (x,13) + x*x + 2*a*x + a*a;

  // RRPoly_1 Rf (f.begin(), f.end());
  // cout << "Rf: " << Rf << endl;
  // RRIPoly_1 RIf (f.begin(), f.end());
  // cout << "RIf: " << RIf << endl;

  // CCPoly_1 Cf (f.begin(), f.end());
  // cout << "Cf: " << Cf << endl;
  // CCIPoly_1 CIf (f.begin(), f.end());
  // cout << "CIf: " << CIf << endl;
  // cout << "RIf (sqrt(2)+2.5*%i) = "
  //      << CGAL::oformat (RIf.evaluate (CC (CGAL::sqrt(2.), 2.5))) << endl;
  // cout << "  type: "
  //      << typeid (RIf.evaluate (CC (CGAL::sqrt(2.), 2.5))).name() << endl;

  return 0;
}
