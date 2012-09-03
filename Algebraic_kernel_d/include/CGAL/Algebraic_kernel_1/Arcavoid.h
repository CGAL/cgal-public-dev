#ifndef CGAL_USE_ABSOLUTE_VOID
#define CGAL_USE_ABSOLUTE_VOID 1
#endif

#if CGAL_USE_ABSOLUTE_VOID
#include <CGAL/Absolute_void.h>
#else // CGAL_USE_ABSOLUTE_VOID

#ifndef CGAL_ARCAVOID_H
#define CGAL_ARCAVOID_H

#ifndef Bisolve_telemetry_code
#define Bisolve_telemetry_code(x) 
#endif

#ifndef CGAL_ARCAVOID_DISABLE_DOUBLE_WITH_EXPONENT
#define CGAL_ARCAVOID_DISABLE_DOUBLE_WITH_EXPONENT 1 // default: TODO
#endif

#define CGAL_ARCAVOID_DISABLE_DOUBLE 1 // default: TODO
#ifndef CGAL_ARCAVOID_DISABLE_DOUBLE
#define CGAL_ARCAVOID_DISABLE_DOUBLE CGAL_ARCAVOID_DISABLE_DOUBLE_WITH_EXPONENT
#endif

#ifndef CGAL_ARCAVOID_USE_NEWTON_REFINEMENT
#define CGAL_ARCAVOID_USE_NEWTON_REFINEMENT 0 // default: 1 [TODO: after bugfix in Newton precision handling]
#endif

#ifndef CGAL_DEBUG_ARCAVOID
#define CGAL_DEBUG_ARCAVOID 0
#endif

#if CGAL_DEBUG_ARCAVOID
#define DBG_ARCA(x) x
#else
#define DBG_ARCA(x) static_cast< void > (0);
#endif

#define CGAL_DEBUG_ARCAVOID_PRECISION 0
#ifndef CGAL_DEBUG_ARCAVOID_PRECISION
#define CGAL_DEBUG_ARCAVOID_PRECISION CGAL_DEBUG_ARCAVOID
#endif

#if CGAL_DEBUG_ARCAVOID_PRECISION
#define DBG_ARCA_PREC(x) x
#else
#define DBG_ARCA_PREC(x) static_cast< void > (0);
#endif

#include <CGAL/Arcavoid_kernel.h>
#include <CGAL/Bitsize.h>
#include <CGAL/Bigfloat_traits.h>

#if defined(__GNUC__) && (__GNUC__ >= 3)
#define CGAL_ARCAVOID_LIKELY(x) (__builtin_expect(!!(x),1))
#define CGAL_ARCAVOID_UNLIKELY(x) (__builtin_expect((x),0))
#else
#define CGAL_ARCAVOID_LIKELY(x) (x)
#define CGAL_ARCAVOID_UNLIKELY(x) (x)
#endif

#if CGAL_ARCAVOID_DISABLE_DOUBLE_WITH_EXPONENT
#warning not using Double_with_exponent
#else
#warning using Double_with_exponent
#include <CGAL/Double_with_exponent.h>
#endif

#include <CGAL/basic.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Cartesian_complex.h>

//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/convex_hull_2.h>

#include <list>
#include <vector>

#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/shared_ptr.hpp>

namespace CGAL {

#define CGAL_ARCAVOID_DOUBLE_OUTPUT 0
#if CGAL_ARCAVOID_DOUBLE_OUTPUT
template< class T >
const std::string uformat (const T &t) {
  return boost::lexical_cast< std::string > (CGAL::to_double (t));
}
#else // CGAL_ARCAVOID_DOUBLE_OUTPUT
template< class T >
const std::string uformat (const T &t) {
  return boost::lexical_cast< std::string > (t);
}

#define CGAL_ARCAVOID_NESTED_QUOTEME(x) #x
#define CGAL_ARCAVOID_QUOTEME(x) CGAL_ARCAVOID_NESTED_QUOTEME(x)
#define CGAL_ARCAVOID_FLOAT_OUTPUT_LENGTH 40
#define CGAL_ARCAVOID_FLOAT_OUTPUT_LENGTH_M_3 37
#define CGAL_ARCAVOID_FLOAT_OUTPUT_LENGTH_M_8 32
// -1.0089657199268891041299411881482222222 // %f
// |3|^--------------- 37 ----------------^
// -1.0089657199268891041299411881482222    // %f
// |3|^-------------- 34 --------------^|3|
// -1.00896571992688910412994118814822e+002 // %e
// |3|^------------- 32 -------------^|-5-|
template<>
const std::string uformat (const Gmpfr &x) {
  char buf [64];
  std::string s;
  
  const double abs_x = CGAL::abs (CGAL::to_double (x));

  if (1e-3 < abs_x && abs_x < 10) {
    int n = mpfr_sprintf (buf, "%." CGAL_ARCAVOID_QUOTEME(CGAL_ARCAVOID_FLOAT_OUTPUT_LENGTH_M_3) "RNf", x.fr());
    s = std::string (buf);
  } else {
    int n = mpfr_sprintf (buf, "%." CGAL_ARCAVOID_QUOTEME(CGAL_ARCAVOID_FLOAT_OUTPUT_LENGTH_M_8) "RNe", x.fr());
    s = buf;
    size_t ind_e = s.find_last_of ('e');
    if (ind_e != std::string::npos) {
      std::string exp_str = s.substr (++ind_e);
      long exp = boost::lexical_cast< long > (exp_str);
      n = mpfr_sprintf (buf, "%+04i", exp);             // hardcoded extent, too
      s = s.substr (0, ind_e);
      s += std::string (buf);
    }
  }

  if (s.length() < CGAL_ARCAVOID_FLOAT_OUTPUT_LENGTH)
    s.insert (s.begin(), CGAL_ARCAVOID_FLOAT_OUTPUT_LENGTH - s.length(), ' ');

  return s;
  
  /*
  Gmpfr r = x;
  mpfr_prec_round (r.fr(), 80, MPFR_RNDA);
  std::string s = boost::lexical_cast< std::string > (r);
  std::replace (s.begin(), s.end(), 'e', 'b');
  return s;
  */
}

template< class T >
const std::string uformat (const Cartesian_complex< T > &z) {
  std::string s = "(";
  s += uformat (z.real());
  s += ", ";
  s += uformat (z.imag());
  s += ")";
  return s;
}
#endif // CGAL_ARCAVOID_DOUBLE_OUTPUT

#ifndef CGAL_ARCAVOID_USE_OLD_ARCAVOID
class Arcavoid_real_root_isolator_tag {};
class Arcavoid_complex_root_isolator_tag {};

template< class BCK, class ArcavoidTag >
class Arcavoid;

#endif // CGAL_ARCAVOID_USE_OLD_ARCAVOID

namespace internal {

template< class BCK >
class Arcavoid_list {
#ifndef CGAL_ARCAVOID_USE_OLD_ARCAVOID
  friend class CGAL::Arcavoid< BCK, Arcavoid_complex_root_isolator_tag >;

  // template< class ArcavoidTag >
  // friend class ::CGAL::Arcavoid< BCK, ArcavoidTag >;
#endif // CGAL_ARCAVOID_USE_OLD_ARCAVOID

public:
  class Approximation;
  class Cluster;

  typedef BCK                               Bitstream_coefficient_kernel;
  typedef typename Bitstream_coefficient_kernel::Coefficient        Coefficient;
  typedef typename Bitstream_coefficient_kernel::Arithmetic_kernel  Arithmetic_kernel;
  typedef typename Arithmetic_kernel::Bigfloat          Bigfloat;
  typedef typename Bitstream_coefficient_kernel::Bigfloat_interval  Bigfloat_interval;

protected:
  typedef std::list< Approximation >                    Approximation_list;
  typedef std::list< Cluster >                          Cluster_list;

public:
  typedef typename Approximation_list::iterator         Approximation_iterator;
  typedef typename Approximation_list::const_iterator   Approximation_const_iterator;
  typedef typename Cluster_list::iterator               Cluster_iterator;
  typedef typename Cluster_list::const_iterator         Cluster_const_iterator;
  typedef std::pair< Cluster_iterator,
                     Cluster_iterator >                 Cluster_range;
  typedef std::pair< Cluster_const_iterator,
                     Cluster_const_iterator >           Cluster_const_range;

protected:
  typedef Arcavoid_list< Bitstream_coefficient_kernel >             List;
  typedef typename Bitstream_coefficient_kernel::Convert_to_bfi     Convert_to_bfi;
  typedef typename Bitstream_coefficient_kernel::Is_zero            Is_zero;
  typedef typename CGAL::Polynomial_type_generator< Coefficient, 1 >::Type Input_polynomial;
  typedef Bigfloat              BF;
  typedef Bigfloat_interval     BFI;
  typedef Bigfloat                      RR;
  typedef CGAL::Cartesian_complex< RR > CC;
  typedef typename CGAL::Polynomial_type_generator< RR, 1 >::Type RRPoly;
  typedef typename CGAL::Polynomial_type_generator< CC, 1 >::Type CCPoly;

  typedef CGAL::Polynomial_traits_d< Input_polynomial > Input_PT;
  typedef CGAL::Polynomial_traits_d< RRPoly >           RRPT;
  typedef CGAL::Polynomial_traits_d< CCPoly >           CCPT;

#if ! CGAL_ARCAVOID_DISABLE_DOUBLE_WITH_EXPONENT
  typedef CGAL::Double_with_exponent            RRDwE;
  typedef CGAL::Cartesian_complex< RRDwE >      CCDwE;
  typedef typename CGAL::Polynomial_type_generator< RRDwE, 1 >::Type RRDwEPoly;
  typedef typename CGAL::Polynomial_type_generator< CCDwE, 1 >::Type CCDwEPoly;

#if ! CGAL_ARCAVOID_DISABLE_DOUBLE
  typedef CGAL::Double_with_exponent            RRD;
  typedef CGAL::Cartesian_complex< RRD >        CCD;
  typedef typename CGAL::Polynomial_type_generator< RRD, 1 >::Type RRDPoly;
  typedef typename CGAL::Polynomial_type_generator< CCD, 1 >::Type CCDPoly;
#endif // ! CGAL_ARCAVOID_DISABLE_DOUBLE
#endif // ! CGAL_ARCAVOID_DISABLE_DOUBLE_WITH_EXPONENT


  /****************************************************************
   * Precision handling
   ****************************************************************/

  // RAII precision lock for RR Bigfloat types
  class Precision_guard {
    long old_prec;
  public:
    Precision_guard (long prec)
      : old_prec (CGAL::get_default_precision< RR >()) {
      DBG_ARCA_PREC (std::cerr << "Requested precision " << prec << " (max so far: " << "[TODO]" << ")" << std::endl);
      //DBG_ARCA_PREC (max_prec = max (prec, max_prec););
      CGAL::set_default_precision< RR > (prec);
    }
    ~Precision_guard () {
      CGAL::set_default_precision< RR > (old_prec);
      DBG_ARCA_PREC (std::cerr << "Restored precision to " << old_prec << std::endl);
    }
    
    // make non-copyable
  private:
    Precision_guard (const Precision_guard &);
    Precision_guard & operator= (const Precision_guard &);
  };

  // RRAI precision lock for Bigfloat intervals
  class Interval_precision_guard {
    long old_prec;
  public:
    Interval_precision_guard (long prec)
      : old_prec (CGAL::get_precision (Bigfloat_interval())) {
      DBG_ARCA_PREC (std::cerr << "Requested interval precision " << prec << " (max so far: " << "[TODO]" << ")" << std::endl);
      CGAL::set_precision (Bigfloat_interval(), prec);
    }
    ~Interval_precision_guard () {
      CGAL::set_precision (Bigfloat_interval(), old_prec);
      DBG_ARCA_PREC (std::cerr << "Restored interval precision to " << old_prec << std::endl);
    }
    
    // make non-copyable
  private:
    Interval_precision_guard (const Interval_precision_guard &);
    Interval_precision_guard & operator= (const Interval_precision_guard &);
  };

  class Approximator : public std::unary_function< Coefficient, RR > {
    Convert_to_bfi convert_to_bfi;
    Is_zero coefficient_is_zero;
    long prec;
  public:
    Approximator()
      : convert_to_bfi (Bitstream_coefficient_kernel().convert_to_bfi_object()),
        coefficient_is_zero (Bitstream_coefficient_kernel().is_zero_object()) {
    }
    Approximator (const Bitstream_coefficient_kernel &bck, long prec)
      : convert_to_bfi (bck.convert_to_bfi_object()),
        coefficient_is_zero (bck.is_zero_object()),
        prec (prec) {}
    const RR operator() (const Coefficient &c) const {
      Bisolve_telemetry_code(t_arca_transform.start();)
      bool is_zero = coefficient_is_zero (c);
      Bisolve_telemetry_code(t_arca_transform.stop();)
      if (is_zero)
        return RR (0);
      
      Precision_guard guard (prec);
      RR m;
      int p = 0;
      const RR mu = CGAL::ipow2< RR > (-2 - prec);
      while (true) {
        Interval_precision_guard intv_guard (prec << p);
        Precision_guard guard (prec << p);
        const Bigfloat_interval I = convert_to_bfi (c);

        ++p;
        if (CGAL::sign (I.inf()) != CGAL::sign (I.sup()))
          continue;

        if (I.sup() - I.inf()
            < mu * CGAL::min (CGAL::abs (I.inf()), CGAL::abs (I.sup()))) {
          m = CGAL::median (I);
          break;
        }
      }
      return m;
    }
    void set_precision (long p) {
      CGAL_precondition (prec < p);
      prec = p;
    }
    const long get_precision () { return prec; }
  };

public:
  class Cluster {
    friend class Arcavoid_list< Bitstream_coefficient_kernel >;

  protected:
    /* Data members */
    Approximation_list discs;
    CC scntr;
    RR srad;
    size_t k;
    long prec;

  public:
    Cluster (long prec)
      : prec (prec) {}

    const size_t multiplicity () const { return k; }
    const CC & center () const { return scntr; }
    const RR & radius () const { return srad; }

    const Approximation_iterator begin () { return discs.begin(); }
    const Approximation_const_iterator begin () const { return discs.begin(); }
    const Approximation_iterator end () { return discs.end(); }
    const Approximation_const_iterator end () const { return discs.end(); }
    
    const bool touch_real () const {
      return ! (radius() < CGAL::abs (center().imag()));
    }
    const bool touch_origin () const {
      return CGAL::abs (center()) <= radius();
    }

    struct Compare_real : public std::binary_function< bool, Cluster, Cluster > {
      const bool operator() (const Cluster &c1, const Cluster &c2) const {
        // TODO: check for overlaps etc.
        return c1.center().real() < c2.center().real();
      }
    };
    struct Touch_real : public std::unary_function< bool, Cluster > {
      const bool operator() (const Cluster &c) const {
        return c.touch_real();
      }
    };
    struct Touch_origin : public std::unary_function< bool, Cluster > {
      const bool operator() (const Cluster &c) const {
        return c.touch_origin();
      }
    };

    struct Cluster_multiplicity_add {
      int operator()(int result, const Cluster& c) const {
        return result + c.multiplicity();
      }
    };

    struct Cluster_critical_multiplicity {
      const bool operator()(const Cluster& c) const {
        return (c.multiplicity() > 1);
      }
    };


  protected:
    void compute_k_center_radius () {
      // Precondition: discs have center and radius computed
      k = discs.size();
      CGAL_precondition (k > 0);
      scntr = CC (0);
      for (Approximation_iterator it = discs.begin(); it != discs.end(); ++it)
        scntr += it->center();
      scntr /= RR (static_cast< long > (k));
      srad = RR (-1);
      for (Approximation_iterator it = discs.begin(); it != discs.end(); ++it)
        srad = CGAL::max (srad, CGAL::abs (it->z - scntr) + it->rad);

      CGAL_postcondition (! CGAL::is_negative (srad));
    }


    friend
    std::ostream & operator<< (std::ostream &out,  const Cluster &c) {
      out << "  Cluster" << std::endl
          << "    Multiplicity: " << c.multiplicity() << std::endl
          << "    Center: " << CGAL::uformat (c.center()) << std::endl
          << "    Radius: " << CGAL::uformat (c.radius())
          << " <= 2^" << CGAL::abs_ilog2 (c.radius()) << std::endl
          << "    Touch real: " << c.touch_real() << std::endl
          << "    Discs:" << std::endl;
      for (Approximation_const_iterator jt = c.discs.begin(); jt != c.discs.end(); ++jt)
        out << "    " << *jt << std::endl;
      return out;
    }
  };

  class Approximation {
    friend class Arcavoid_list< Bitstream_coefficient_kernel >;

  protected:
    Approximation ()
      : z (0),
        fz (0),
        abs_fz (0),
        dfz (0),
        abs_z (0),
        sz (0),
        rad (0) {}

    Approximation (const CC &z, const List &list)
      : z (z),
        fz (list.f.evaluate (z)),
        abs_fz (CGAL::abs (fz)),
        dfz (list.df.evaluate (z)),
        abs_z (CGAL::abs (z)),
        sz (list.s.evaluate (abs_z)) {}

    const CC & center () const { return z; }
    const RR & radius () const { return rad; }

    friend
    std::ostream & operator<< (std::ostream &out,  const Approximation &D) {
      return out << "  Disc (" << CGAL::uformat (D.center()) << "; "
                 << CGAL::uformat (D.rad) << " <= 2^" << CGAL::abs_ilog2 (D.rad) << ")";
    }

    /* Data members */
    CC z, fz, dfz;
    RR abs_z, abs_fz, sz;
    RR rad;
  };

  struct Overlapping_regions {
    template< class OutputIterator >
    void operator() (const Arcavoid_list &lhs, const Arcavoid_list &rhs,
                     OutputIterator out) const {
      for (Cluster_const_iterator it = lhs.clsts.begin(); it != lhs.clsts.end(); ++it) {
        for (Cluster_const_iterator jt = rhs.clsts.begin(); jt != rhs.clsts.end(); ++jt) {
          if (CGAL::abs (it->center() - jt->center()) <= it->radius() + jt->radius()) {
            *out++ = std::make_pair (it, jt);
          }
        }
      }
    }
  };

private:
  /* Input */
  Input_polynomial input;        // TODO: better to store a "raw" std::vector< Coefficient >?

  // TODO: are these three necessary for Approximator only? If so, remove from Arcavoid_list.
  Bitstream_coefficient_kernel bck;
  Convert_to_bfi convert_to_bfi;
  Is_zero coefficient_is_zero;

  /* Working variables */
  long prec;
  RR mu;
  Approximator approximator;

  RRPoly f, df, s;
  int n, N;
  int mult_zero;
  Cluster_list clsts;
  
  bool first_stage_finished;

  /* temporaries */
  Approximation_list eps_nbh;

public:
  /****************************************************************
   * Constructors
   ****************************************************************/
  Arcavoid_list ()
    : bck(),
      convert_to_bfi (bck.convert_to_bfi_object()),
      coefficient_is_zero (bck.is_zero_object()),
      input(),
      prec (1),
      approximator (bck, prec),
      n (input.degree()),
      N (n),
      mult_zero (0),
      first_stage_finished (true)
  {}            // TODO: make reasonable defaults

  template< class InputIterator >
  Arcavoid_list (const Bitstream_coefficient_kernel &_bck, InputIterator first, InputIterator beyond)
    : bck (_bck),
      convert_to_bfi (bck.convert_to_bfi_object()),
      coefficient_is_zero (bck.is_zero_object()),
      input (first, beyond),
      first_stage_finished (false) {
    init();
  }

  /* convenience constructor */
  Arcavoid_list (const Bitstream_coefficient_kernel &_bck, const Input_polynomial &poly)
    : bck (_bck),
      convert_to_bfi (bck.convert_to_bfi_object()),
      coefficient_is_zero (bck.is_zero_object()),
      input (poly),
      first_stage_finished (false) {
    init();
  }

  ~Arcavoid_list () {
    std::cerr << "~Arcavoid_list - max precision: " << prec << std::endl;
  }

private:
  // Define Horner-style evaluation of f and it's derivatives
  // (Knowing the evaluation scheme is essential for guarantees about the evaluation error)

  //! evaluate f at t
  template< typename T, typename U, class Polynomial >
  static void horner (const Polynomial &f, const T& t, U& f0) {
    int d = f.degree();

    if (CGAL_ARCAVOID_UNLIKELY (d < 0)) {
      f0 = 0;
      return;
    }

    f0 = f[d];
    while (--d >= 0) {
      f0 *= t;
      f0 += f[d];
    }
  }

  //! evaluate f and f' at t
  template< typename T, typename U, class Polynomial >
  static void horner_1 (const Polynomial &f, const T& t, U& f0, U& f1) {
    int d = f.degree();

    if (CGAL_ARCAVOID_UNLIKELY (d < 1)) {
      if (d < 0) {
        f0 = f1 = 0;
      } else { // d == 0
        f0 = f[0];
        f1 = 0;
      }
      return;
    }

    f0 = f[d];
    f1 = f[d];
    while (--d >= 1) {
      f0 *= t;
      f0 += f[d];
      f1 *= t;
      f1 += f0;
    }
    f0 *= t;
    f0 += f[0];
  }

  //! evaluate f, f', and f'' at t
  template< typename T, typename U, class Polynomial >
  static void horner_2 (const Polynomial &f, const T& t, U& f0, U& f1, U& f2) {
    int d = f.degree();

    if (CGAL_ARCAVOID_UNLIKELY (d < 2)) {
      if (d < 0) {
        f0 = f1 = f2 = 0;
      } else if (d == 0) {
        f0 = f[0];
        f1 = f2 = 0;
      } else { // d == 1
        f1 = f0 = f[1];
        f0 *= t;
        f0 += f[0];
        f2 = 0;
      }
      return;
    }

    f0 = f[d];
    f1 = f[d];
    f2 = f[d];
    while (--d >= 2) {
      f0 *= t;
      f0 += f[d];
      f1 *= t;
      f1 += f0;
      f2 *= t;
      f2 += f1;
    }

    f0 *= t;
    f0 += f[1];
    f1 *= t;
    f1 += f0;

    f0 *= t;
    f0 += f[0];
  }

  void init () {
    prec = 53;
    mu = CGAL::ipow2< RR > (1-prec);
    approximator = Approximator (bck, prec);

    Precision_guard guard (prec);
    f = RRPoly (boost::make_transform_iterator (input.begin(), approximator),
                boost::make_transform_iterator (input.end(), approximator));

    // TODO: check for degree loss

    N = f.degree();

    mult_zero = 0;
    while (mult_zero <= N && CGAL::is_zero (f[mult_zero]))
      ++mult_zero;
    n = N - mult_zero;

    if (mult_zero != 0)
      f = RRPoly (f.begin() + mult_zero, f.end());

    CGAL_postcondition (f.degree() < 0 || ! CGAL::is_zero (f[0]));
    CGAL_postcondition (n == f.degree());
    
    df = CGAL::differentiate (f);
    typename CGAL::Real_embeddable_traits< RR >::Abs abs;
    s = RRPoly (boost::make_transform_iterator (f.begin(), abs),
                boost::make_transform_iterator (f.end(), abs));
    DBG_ARCA(std::cerr << "Constructed Arcavoid_rep for f = " << f << std::endl
             << "                     MULT_ZERO = " << mult_zero << std::endl
             << "                            df = " << df << std::endl
             << "                             s = " << s << std::endl
             << "                         input = " << input << std::endl);
    
    initialize_roots();
  }

  const bool in_eps_nbh (long p, const Approximation &D) const {
    //return D.fz.norm_inf() <= D.sz * RR (4*n+1) * RR (n) * CGAL::ipow2< RR > (1-p);
    const RR eps = RR (8*n+2) * CGAL::ipow2< RR > (1-p);
    return D.abs_fz <= eps * D.sz;
  }

  void aberth_in_cluster (Cluster &c, int max_iter = 250) {
    CGAL::Timer t_aberth;
    t_aberth.start();

    // std::cerr << "ABERTH in cluster" << std::endl
    //           << c << std::endl;

    CGAL_precondition (eps_nbh.empty());
    int count = 0, icount = 0;

    do {
      // std::cerr << "WHILE loop" << std::endl
      //           << c << std::endl;
      Approximation_iterator it = c.discs.begin();
      while (it != c.discs.end()) {
        ++icount;
        CC abcorr = 0;
        for (Approximation_const_iterator jt = c.discs.begin(); jt != c.discs.end(); ++jt) {
          if (it == jt) continue;
          abcorr += (it->z - jt->z).reciprocal();
        }
        for (Approximation_const_iterator jt = eps_nbh.begin(); jt != eps_nbh.end(); ++jt)
          abcorr += (it->z - jt->z).reciprocal();
        CGAL_assertion (! CGAL::is_zero (it->dfz));
        const CC f_df = it->fz / it->dfz;
        CGAL_assertion (! CGAL::is_zero (CC (1) - f_df * abcorr));
        const CC corr = f_df / (CC (1) - f_df * abcorr);
        *it = Approximation (it->z - corr, *this);
        //compute_radius (it, c.prec);
        CGAL_postcondition (! CGAL::is_zero (it->z));

        if (in_eps_nbh (c.prec, *it)) {
          Approximation_iterator old_it = it++;
          eps_nbh.splice (eps_nbh.begin(), c.discs, old_it);
        } else {
          ++it;
        }
      }

      ++count;
    } while (! c.discs.empty());
    DBG_ARCA (std::cerr << "Aberth iterations: " << count << " / " << icount << std::endl;);

    t_aberth.stop();
    std::cerr << "Time for Aberth: " << t_aberth.time() << " (" << t_aberth.time() / icount << " per approx" << std::endl;

    for (Approximation_iterator it = eps_nbh.begin(); it != eps_nbh.end(); ++it)
      compute_radius (it, c.prec);

    c.discs.splice (c.discs.begin(), eps_nbh);
  }


  void compute_radius (Approximation_iterator dit, long p) {
    compute_radius_thm14 (dit, p);
  }

  const RR product_of_distances (Approximation_const_iterator dit) const {
    CGAL_postcondition_code (int found_self = 0;);
    CGAL_postcondition_code (int count = 0;);
    CGAL_postcondition_code (int skipped = 0;);

    RR prod = RR (1);
    for (Cluster_const_iterator it = clsts.begin(); it != clsts.end(); ++it) {
      for (Approximation_const_iterator jt = it->begin(); jt != it->end(); ++jt) {
        CGAL_postcondition_code (if (dit == jt) ++found_self;);
        if (dit == jt) continue;
        
        CGAL_postcondition_code (if (CGAL::is_zero (jt->z))
                                   ++skipped;);
        if (CGAL::is_zero (jt->z)) continue;
        
        CGAL_postcondition_code (++count;);
        prod *= CGAL::abs (dit->z - jt->z);
      }
    }
    
    for (Approximation_const_iterator jt = eps_nbh.begin(); jt != eps_nbh.end(); ++jt) {
      CGAL_postcondition_code (if (dit == jt) ++found_self;);
      if (dit == jt) continue;
      
      CGAL_postcondition_code (if (CGAL::is_zero (jt->z))
                                 ++skipped;);
      if (CGAL::is_zero (jt->z)) continue;
      
      CGAL_postcondition_code (++count;);
      prod *= CGAL::abs (dit->z - jt->z);
    }

    //DBG_ARCA (CGAL_postcondition_code (std::cerr << "found_self: " << found_self << " skipped: " << skipped << " count: " << count << std::endl;));
    CGAL_postcondition (found_self == 1);
    CGAL_postcondition (skipped == mult_zero || skipped == 0);
    CGAL_postcondition (count == n-1);

    return prod;
  }

  void compute_radius_gershgorin (Approximation_iterator dit, long) {
    // TODO: efficiently handle approxs in other clusters

    DBG_ARCA (std::cerr << "COMPUTE_RADIUS_GERSHGORIN" << std::endl;);
    
    if (CGAL::is_zero (dit->z))
      return;

    dit->rad = RR (n) * dit->abs_fz;
    dit->rad /= s.lcoeff() * product_of_distances (dit);
  }

  void compute_radius_thm14 (Approximation_iterator dit, long p) {
    // TODO: efficiently handle approxs in other clusters
    
    if (CGAL::is_zero (dit->z))
      return;

    RR eps = RR (12*n+3) * CGAL::ipow2< RR > (1-p);

    /*
    DBG_ARCA (std::cerr << "z: " << CGAL::uformat (dit->z) << std::endl;);
    DBG_ARCA (std::cerr << "fz: " << CGAL::uformat (dit->fz) << std::endl;);
    DBG_ARCA (std::cerr << "sz: " << CGAL::uformat (dit->sz) << std::endl;);
    DBG_ARCA (std::cerr << "eps before: " << eps << std::endl;);
    */

    if (eps_nbh.empty()) {
      for (Cluster_const_iterator it = clsts.begin(); it != clsts.end(); ++it) {
        for (Approximation_const_iterator jt = it->begin(); jt != it->end(); ++jt) {
          if (CGAL::is_zero (jt->z)) continue;
          CGAL_assertion (! CGAL::is_zero (jt->sz));
          eps = CGAL::max (eps, jt->abs_fz / jt->sz * RR(1.5));
        }
      }
    }
    // else {
    //   for (Approximation_const_iterator jt = eps_nbh.begin(); jt != eps_nbh.end(); ++jt) {
    //     if (CGAL::is_zero (jt->z)) continue;
    //     eps = CGAL::max (eps, jt->abs_fz / jt->sz * RR(1.5));
    //   }
    // }

    // eps = CGAL::max (eps, dit->abs_fz / dit->sz * RR(1.5));

    // DBG_ARCA (std::cerr << "eps after: " << eps << std::endl;);
    
    CGAL_postcondition_code
      (if (! eps_nbh.empty()) {
      //   for (Cluster_const_iterator it = clsts.begin(); it != clsts.end(); ++it)
      //     for (Approximation_const_iterator jt = it->begin(); jt != it->end(); ++jt)
      //       { CGAL_postcondition (jt->abs_fz <= RR(2./3.) * eps * jt->sz); }
      // } else {
        for (Approximation_const_iterator jt = eps_nbh.begin(); jt != eps_nbh.end(); ++jt)
          { CGAL_postcondition (jt->abs_fz <= RR(2./3.) * eps * jt->sz); }
      });
                           
    dit->rad = RR (2 * n) * eps * dit->sz;
    dit->rad /= s.lcoeff() * product_of_distances (dit);
  }

  void compute_radius_mpsolve (Approximation_iterator dit, long p) {
    if (CGAL::is_zero (dit->z)) {
      CGAL_assertion (CGAL::is_zero (dit->rad));
      return;
    }

    const RR eps = RR (4) * RR (n) * CGAL::ipow2< RR > (1-p);     // 4n mu = 4n 2^(1-p)
    CGAL_assertion (! CGAL::is_zero (dit->dfz));
    dit->rad = RR (n) * (dit->abs_fz + eps * dit->sz) / CGAL::abs (dit->dfz);
    // TODO: set prec
    // TODO: IMPORTANT check CGAL::max (eps, ...)
  }

  const bool discs_pairwise_newton_isolated (const Approximation &D1, const Approximation &D2) const {
    return RR (n << 1) /* RR (n << 1) */ * (D1.rad + D2.rad)
      < CGAL::abs (D1.center() - D2.center()); // TODO: check rounding
  }

  const Cluster_range split_cluster_into_newton_connected_components (Cluster_iterator cit) {
    DBG_ARCA (std::cerr << "SPLIT_CLUSTER_INTO_NEWTON_CONNECTED_COMPONENTS" << std::endl;);

    const size_t k = cit->multiplicity();
    std::vector< Approximation_iterator > dits;
    for (Approximation_iterator it = cit->discs.begin(); it != cit->discs.end(); ++it)
      dits.push_back (it);
    
    /* initialize adjacency graph */
    boost::adjacency_matrix< boost::undirectedS > G (k);
    for (size_t i = 0; i < k; ++i)
      for (size_t j = 0; j < i; ++j)
        if (! discs_pairwise_newton_isolated (*dits[i], *dits[j]))
          add_edge (i, j, G);

    /* compute connected components */
    std::vector< int > comp (k);
    const size_t nr_comp = boost::connected_components (G, &comp[0]);
    
    Cluster_iterator beyond = cit;
    ++beyond;

    if (nr_comp == 1) {
      CGAL_postcondition_code (size_t kold = cit->multiplicity(););
      cit->compute_k_center_radius();
      CGAL_assertion (kold == cit->multiplicity());
      return Cluster_range (cit, beyond);
    }

    Cluster_list subclusters;
    CGAL_postcondition_code (size_t ksum = 0;);
    for (int i = 0; i < nr_comp; ++i) {
      Cluster c (cit->prec);
      for (int j = 0; j < k; ++j)
        if (comp[j] == i)
          c.discs.splice (c.discs.begin(), cit->discs, dits[j]);
      c.compute_k_center_radius();
      CGAL_postcondition_code (ksum += c.multiplicity(););
      subclusters.push_front (c); // TODO: avoid copy
    }
    CGAL_assertion (ksum == k);

    cit = clsts.erase (cit);
    CGAL_assertion (cit == beyond);
    const Cluster_range range (subclusters.begin(), beyond);
    clsts.splice (beyond, subclusters);

    CGAL_postcondition (std::distance (range.first, range.second) == nr_comp);

    return range;
  }

#if ! CGAL_ARCAVOID_DISABLE_DOUBLE_WITH_EXPONENT
  const Cluster_range first_stage () {
    CGAL_precondition (! first_stage_finished);
    CGAL_precondition (clsts.size() == 1);

    DBG_ARCA (std::cerr << "DOUBLE STAGE" << std::endl;);

    
    std::vector< RRDwE > coeffs (n+1);

    for (int i = 0; i <= n; ++i)
      coeffs[i] = f[i];
    RRDwEPoly dwef (coeffs.begin(), coeffs.end());
    RRDwEPoly dwedf = CGAL::differentiate (dwef);

    for (int i = 0; i <= n; ++i)
      coeffs[i] = s[i];
    RRDwEPoly dwes (coeffs.begin(), coeffs.end());

    Approximation_iterator begin = clsts.front().begin();
    std::advance (begin, mult_zero);
    std::vector< CCDwE > dwez (n);
    Approximation_iterator it = begin;
    for (int i = 0; i < n; ++i) {
      dwez[i] = it->z;
      ++it;
    }
    CGAL_assertion (it == clsts.front().end());

    DBG_ARCA (std::cerr << "- CONVERSION DONE" << std::endl;);

    std::vector< bool > in_eps_nbh (n, false);
    bool all_in_eps_nbh = false;
    const RRDwE eps = RRDwE (8*n+2) * CGAL::ipow2< RRDwE > (-52);
    int count = 0;
    while (++count < 250 && ! all_in_eps_nbh) {
      for (int i = 0; i < n; ++i) {
        if (in_eps_nbh[i])
          continue;

        const CCDwE fz = dwef.evaluate (dwez[i]);
        const RRDwE sz = dwes.evaluate (CGAL::abs (dwez[i]));
        in_eps_nbh[i] = (CGAL::abs (fz) <= eps * sz);
        if (in_eps_nbh[i])
          continue;

        const CCDwE dfz = dwedf.evaluate (dwez[i]);
        CGAL_assertion (! CGAL::is_zero (dfz));
        const CCDwE f_df = fz / dfz;
        
        CCDwE abcorr = 0;
        for (int j = 0; j < n; ++j) {
          if (i == j) continue;
          abcorr += (dwez[i] - dwez[j]).reciprocal();
        }
        CGAL_assertion (! CGAL::is_zero (CCDwE (1) - f_df * abcorr));
        const CCDwE corr = f_df / (CCDwE (1) - f_df * abcorr);
        dwez[i] -= corr;
        CGAL_postcondition (! CGAL::is_zero (dwez[i]));
      }

      all_in_eps_nbh = true;
      for (int i = 0; all_in_eps_nbh && i < n; ++i)
        all_in_eps_nbh &= in_eps_nbh[i];
    }

    DBG_ARCA (std::cerr << "- ABERTH DONE" << std::endl;);

    it = begin;
    for (int i = 0; i < n; ++i) {
      const CC z (dwez[i].real().to_gmpfr(), dwez[i].imag().to_gmpfr());
      *it++ = Approximation (z, *this);
    }

    for (it = begin; it != clsts.front().end(); ++it)
      compute_radius (it, 53);

    const Cluster_range range = split_cluster_into_newton_connected_components (clsts.begin());

#if CGAL_DEBUG_ARCAVOID
    for (Cluster_const_iterator it = range.first; it != range.second; ++it) {
      std::cerr << "  C:";
      for (Approximation_const_iterator jt = it->discs.begin(); jt != it->discs.end(); ++jt)
        std::cerr << "\t" << *jt << std::endl;
    }
#endif

    DBG_ARCA (std::cerr << "- END OF DOUBLE STAGE" << std::endl;);

    first_stage_finished = true;
    return range;
  }
#endif

public:
  const Cluster_range subdivide_cluster (Cluster_iterator cit) {
    DBG_ARCA (std::cerr << "SUBDIVIDE_CLUSTER" << std::endl;);

#if ! CGAL_ARCAVOID_DISABLE_DOUBLE_WITH_EXPONENT
    if (! first_stage_finished) {
      CGAL_precondition (cit == clsts.begin());
      return first_stage();
    }
#endif

    if (CGAL::is_zero (cit->radius())) {
      // CGAL_precondition (CGAL::is_zero (cit->center())); // TODO: check this
      // TODO: problem with precondition if exact root found
      Cluster_iterator beyond = cit;
      ++beyond;
      return Cluster_range (cit, beyond);
    }

    // TODO: refine then, but check for zeros in the front.
    // if (cit->multiplicity() > 1) {
    //   Approximation_const_iterator second = cit->begin();
    //   Approximation_const_iterator first = second++;
    //   if (first->z == second->z) {
    //     refine (cit);
    //     Cluster_iterator beyond = cit;
    //     ++beyond;
    //     return Cluster_range (cit, beyond);
    //   }
    // }

    Approximation_list zeros;
    if (mult_zero > 0 && cit->touch_origin()) {
      Approximation_iterator it = cit->discs.begin();
      while (it != cit->discs.end()) {
        if (CGAL::is_zero (it->radius())) {
          Approximation_iterator old_it = it;
          ++it;
          zeros.splice (zeros.begin(), cit->discs, old_it);
        } else {
          ++it;
        }
      }
      CGAL_assertion (zeros.size() == mult_zero);
    }

    cit->prec <<= 1;
    refine_bitstream (cit->prec);
    Precision_guard guard (cit->prec);
    aberth_in_cluster (*cit);
    
    cit->discs.splice (cit->discs.begin(), zeros);

    const Cluster_range range = split_cluster_into_newton_connected_components (cit);

#if CGAL_DEBUG_ARCAVOID
    for (Cluster_const_iterator it = range.first; it != range.second; ++it) {
      std::cerr << "  C:";
      for (Approximation_const_iterator jt = it->discs.begin(); jt != it->discs.end(); ++jt)
        std::cerr << "\t" << *jt << std::endl;
    }
#endif

    return range;
  }

  void refine_real (Cluster_iterator cit) {
    DBG_ARCA (std::cerr << "REFINE REAL CLUSTER" << std::endl;);

    CGAL_precondition (cit->touch_real());

    if (CGAL::is_zero (cit->radius()))
      return;

    //const long p = - (CGAL::abs_ilog2 (cit->radius()) << 1);
    const long p = cit->prec << 1;

    CGAL_postcondition_code (const RR rad_before = cit->radius(););

    DBG_ARCA (std::cerr << "prec = " << p << std::endl;);
    cit->prec = p;
    refine_bitstream (p);
    Precision_guard guard (p);
    
    Approximation_iterator it = cit->begin();
    const RR corr = cit->multiplicity()
      * f.evaluate (cit->center().real()) / df.evaluate (cit->center().real());
    DBG_ARCA (std::cerr << "corr = " << uformat (corr) << std::endl;);
    CGAL_assertion (! CGAL::is_zero (df.evaluate (cit->center().real())));
    *it = Approximation (cit->center().real() - corr, *this);

    DBG_ARCA (int count = 1;);
    // while (! in_eps_nbh (cit->prec, *it)) {
    //   ++count;
    //   const CC corr = cit->multiplicity() * it->fz / it->dfz;
    //   DBG_ARCA (std::cerr << "corr = " << corr << std::endl;);
    //   CGAL_assertion (! CGAL::is_zero (it->dfz));
    //   *it = Approximation (it->center() - corr, *this);
    // }
    DBG_ARCA (std::cerr << count << " Newton steps performed" << std::endl;);

    it->rad = CGAL::square (cit->srad);
    //compute_radius_mpsolve (it, cit->prec);

    // long pr = CGAL::abs_ilog2 (it->rad);
    // long pz = CGAL::abs_ilog2 (it->z.real());
    // cit->prec = pz - pr;

    // rad = 0.00000000000000000000000000000001
    //                                        ^er
    // z =   0.0000000000000000000010001010101000101010100101001
    //                             ^ez
    // z~ =  0.00000000000000000000100010101010

    // rad <= 2^er
    // z = m * 2^ez

    cit->scntr = it->center();
    cit->srad = it->radius();

    CGAL_postcondition (cit->radius() <= rad_before);

    for (++it; it != cit->end(); ++it)
      *it = *(cit->begin());
  }

  void refine (Cluster_iterator cit) {
#if ! CGAL_ARCAVOID_USE_NEWTON_REFINEMENT
    // TODO: dummy implementation (uses subdivide as refinement step)
    // TODO: CHECK WHY DUMMY DOES NOT WORK
    Cluster_range range = subdivide_cluster (cit);
    CGAL_postcondition (cit == range.first);
    CGAL_postcondition (std::distance (range.first, range.second) == 1);
    return;
#else // CGAL_ARCAVOID_USE_NEWTON_REFINEMENT
    // TODO: CHECK PRECISION HANDLING IN "CLEVER" Newton REFINE
    // Try: (y^4 - (15*y-1)^2)^2 + (x^4 - (20*x-1)^2)^2

    const Input_polynomial x = typename Input_PT::Shift() (Input_polynomial (1), 1, 0);
    DBG_ARCA (std::cerr << "REFINE CLUSTER" << std::endl
              << "alpha: " << CGAL::uformat (approximator (x[0])) << std::endl
              << *cit;);

    if (cit->touch_real())
      return refine_real (cit); // TODO: use real part only

    if (CGAL::is_zero (cit->radius()))
      return;

    CGAL_postcondition_code (const RR rad_before = cit->radius(););

    //const long p = - (CGAL::abs_ilog2 (cit->radius()) << 1);
    const long p = cit->prec << 1;

    cit->prec = p;
    refine_bitstream (p);
    Precision_guard guard (p);
    
    //const RR rad = CGAL::ipow2< RR > (-p);
    const RR rad = CGAL::square (cit->radius());

    Approximation_iterator it = cit->begin();
    const CC corr = cit->multiplicity()
      * f.evaluate (cit->center()) / df.evaluate (cit->center());
    CGAL_assertion (! CGAL::is_zero (df.evaluate (cit->center())));
    *it = Approximation (cit->center() - corr, *this);

    int count = 1;
    while (! in_eps_nbh (cit->prec, *it)) {
      ++count;
      CGAL_assertion (! CGAL::is_zero (it->dfz));
      const CC corr = cit->multiplicity() * it->fz / it->dfz;
      *it = Approximation (it->center() - corr, *this);
    }
    DBG_ARCA (std::cerr << count << " Newton steps performed" << std::endl;);
    //it->rad = rad;

    it->rad = CGAL::square (cit->srad);
    // compute_radius_mpsolve (it, cit->prec);

    cit->scntr = it->center();
    cit->srad = it->radius();

    CGAL_postcondition (cit->radius() <= rad_before);

    for (++it; it != cit->end(); ++it)
      *it = *(cit->begin());

    DBG_ARCA (std::cerr << "REFINED TO: " << *(cit->begin()) << std::endl;);

    // TODO: dummy implementation (uses subdivide as refinement step)
    // TODO: CHECK WHY DUMMY DOES NOT WORK
    // Cluster_range range = subdivide_cluster (cit);
    // CGAL_postcondition (cit == range.first);
    // CGAL_postcondition (std::distance (range.first, range.second) == 1);
#endif
  }

  void sort_real_cluster_range (Cluster_range &real_range) {
    if (real_range.first == real_range.second) // empty range
      return;
    
    Cluster_list real_clsts;
    real_clsts.splice (real_clsts.begin(), clsts, real_range.first, real_range.second);
    real_clsts.sort (typename Cluster::Compare_real());
    real_range.first = real_clsts.begin();
    clsts.splice (real_range.second, real_clsts);
  }

  const Cluster_range real_subdivide_cluster (Cluster_iterator cit) {
    DBG_ARCA (if (! cit->touch_real())
                std::cerr << "real_subdivide_cluster() called for cluster" << std::endl
                          << *cit << "not touching the real line" << std::endl);
    
    const Cluster_range range = subdivide_cluster (cit);
    CGAL_postcondition_code (int total = std::distance (range.first, range.second););

    const Cluster_iterator begin_imag
      = std::partition (range.first, range.second, typename Cluster::Touch_real());
    // TODO: triple-check if stable_partition is necessary
    Cluster_range real_range (range.first, begin_imag);
    CGAL_postcondition_code (Cluster_range imag_range (begin_imag, range.second););

    sort_real_cluster_range (real_range);

#if CGAL_DEBUG_ARCAVOID
    // std::cerr << "# TOTAL " << total << std::endl;
    // std::cerr << "# REAL: " << real << std::endl;
    // std::cerr << "# IMAG: " << imag << std::endl;
    std::cerr << "Refined real root approximations to:" << std::endl;
    for (Cluster_const_iterator it = real_range.first; it != real_range.second; ++it) {
      std::cerr << "  C:";
      for (Approximation_const_iterator jt = it->discs.begin(); jt != it->discs.end(); ++jt)
        std::cerr << "\t" << *jt << std::endl;
    }
    std::cerr << "Refined complex root approximations to:" << std::endl;
    for (Cluster_const_iterator it = real_range.second; it != range.second; ++it) {
      std::cerr << "  C:";
      for (Approximation_const_iterator jt = it->discs.begin(); jt != it->discs.end(); ++jt)
        std::cerr << "\t" << *jt << std::endl;
    }
#endif

    CGAL_postcondition_code (int real = std::distance (real_range.first, real_range.second););
    CGAL_postcondition_code (int imag = std::distance (imag_range.first, imag_range.second););
    CGAL_postcondition (real + imag == total);
    CGAL_postcondition (real_range.second == imag_range.first);

    return real_range;
  }

  // TODO: write p_a_s_r (first, beyond) and use for real_subdivide_cluster (cit)
  const Cluster_iterator partition_and_sort_reals () {
    const Cluster_iterator begin_imag
      = std::partition (clsts.begin(), clsts.end(), typename Cluster::Touch_real());
    // TODO: triple-check if stable_partition is necessary
    Cluster_range real_range (clsts.begin(), begin_imag);
    sort_real_cluster_range (real_range);
    return begin_imag;
  }

  const Cluster_iterator begin () { return clsts.begin(); }
  const Cluster_const_iterator begin () const { return clsts.begin(); }
  const Cluster_iterator end () { return clsts.end(); }
  const Cluster_const_iterator end () const { return clsts.end(); }

  const int multiplicity_of_zero () const { return mult_zero; }
  const int degree () const { return N; }

  const Cluster_range cluster_range () {
    CGAL_assertion (std::distance (clsts.begin(), clsts.end()) >= 1);
    CGAL_assertion (std::distance (clsts.begin(), clsts.end()) <= n);
    return Cluster_range (clsts.begin(), clsts.end());
  }
  const Cluster_const_range cluster_range () const {
    CGAL_assertion (std::distance (clsts.begin(), clsts.end()) >= 1);
    CGAL_assertion (std::distance (clsts.begin(), clsts.end()) <= n);
    return Cluster_const_range (clsts.begin(), clsts.end());
  }

  const Input_polynomial & polynomial () const { return input; }
  const Bitstream_coefficient_kernel & kernel () const { return bck; }

protected:
  void refine_bitstream (long p) {
    while (p > prec)
      refine_bitstream();
  }
  
  void refine_bitstream () {
    prec <<= 1;
    mu = CGAL::ipow2< RR > (1-prec);
    approximator.set_precision (prec);
    Precision_guard guard (prec + n);
    f = RRPoly (boost::make_transform_iterator (input.begin() + mult_zero, approximator),
                boost::make_transform_iterator (input.end(), approximator));
    df = CGAL::differentiate (f);
    typename CGAL::Real_embeddable_traits< RR >::Abs abs;
    s = RRPoly (boost::make_transform_iterator (f.begin(), abs),
                boost::make_transform_iterator (f.end(), abs));
  }

  void initialize_roots_simple () {
    CGAL_precondition (clsts.empty());

    if (N <= 0)
      return;

    clsts.push_front (Cluster (prec));

    const CC w = CC (3, 4) / RR(5);
    CC z = CC (1, 0);
    for (int i = 0; i < n; ++i) {
      z *= w;
      clsts.front().discs.push_front (Approximation (z, *this));
    }
    for (int i = 0; i < mult_zero; ++i)
      clsts.front().discs.push_front (Approximation ());

    CGAL_postcondition (clsts.front().discs.size() == N);

    const Approximation_const_iterator beyond = clsts.front().end();
    for (Approximation_iterator it = clsts.front().begin();
         it != beyond; ++it)
      compute_radius (it, prec);

    clsts.front().compute_k_center_radius();

#if CGAL_DEBUG_ARCAVOID
    std::cerr << "Initialized root approximations to:" << std::endl;
    for (Cluster_const_iterator it = clsts.begin(); it != clsts.end(); ++it) {
      std::cerr << "  C:";
      for (Approximation_const_iterator jt = it->discs.begin(); jt != it->discs.end(); ++jt)
        std::cerr << "\t" << *jt << std::endl;
    }
#endif
  }

  void initialize_roots () {
    CGAL_precondition (clsts.empty());

    std::cerr << "Init roots for f = " << CGAL::oformat (f) << std::endl;

    if (N <= 0)
      return;

    clsts.push_front (Cluster (prec));

    typedef CGAL::Simple_cartesian< double > K;
    typedef K::Point_2 Point_2;
    std::vector< Point_2 > points, ch;
    for (int i = 0; i <= n; ++i)
      if (! is_zero (f[i]))
        points.push_back (Point_2 (i, abs_ilog2 (abs (f[i]))));
    CGAL::upper_hull_points_2 (points.begin(), points.end(), std::back_inserter (ch));
    ch.push_back (points.front());
    std::reverse (ch.begin(), ch.end());

    const double offset = (double)(std::rand()) / RAND_MAX * 2 * M_PI;

    for (int i = 1; i < ch.size(); ++i) {
      int k = (int)(to_double (ch[i-1].x()));
      int K = (int)(to_double (ch[i].x()));
      RR u = abs (f[k]) / abs (f[K]);
      if ((K-k) == 2)
        u = sqrt (u);
      else if ((K-k) > 2)
        u = root_d (u, K-k, 53);

      for (int j = 0; j < K-k; ++j) {
        const double angle = 2. * M_PI * (((double)j)/(K-k)  + ((double)i)/n) + offset;
        const CC z = u * CC (std::cos (angle), std::sin (angle));
        clsts.front().discs.push_front (Approximation (z, *this));
      }
    }

    // add zeros
    for (int i = 0; i < mult_zero; ++i)
      clsts.front().discs.push_front (Approximation ());

    CGAL_postcondition (clsts.front().discs.size() == N);

    const Approximation_const_iterator beyond = clsts.front().end();
    for (Approximation_iterator it = clsts.front().begin();
         it != beyond; ++it)
      compute_radius (it, prec);

    clsts.front().compute_k_center_radius();

#if CGAL_DEBUG_ARCAVOID
    std::cerr << "Initialized root approximations to:" << std::endl;
    for (Cluster_const_iterator it = clsts.begin(); it != clsts.end(); ++it) {
      std::cerr << "  C:";
      for (Approximation_const_iterator jt = it->discs.begin(); jt != it->discs.end(); ++jt)
        std::cerr << "\t" << *jt << std::endl;
    }
#endif
  }

  /****************************************************************
   * Statistics
   ****************************************************************/
  const long precision () const { return prec; }
};

} // namespace internal

} // namespace CGAL

#undef DBG_ARCA
#undef DBG_ARCA_PREC

#endif // CGAL_ARCAVOID_H

#endif // CGAL_USE_ABSOLUTE_VOID
