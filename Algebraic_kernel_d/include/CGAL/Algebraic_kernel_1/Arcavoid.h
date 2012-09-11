#ifndef CGAL_ARCAVOID_H
#define CGAL_ARCAVOID_H

#include <CGAL/Algebraic_kernel_d/flags.h>

////////////////////////////////////////////////////////////////
// PARAMETERS
////////////////////////////////////////////////////////////////

#ifndef CGAL_ARCAVOID_DISABLE_DOUBLE_WITH_EXPONENT
# define CGAL_ARCAVOID_DISABLE_DOUBLE_WITH_EXPONENT 1 // default: 0
#endif

#ifndef CGAL_ARCAVOID_DISABLE_DOUBLE
# define CGAL_ARCAVOID_DISABLE_DOUBLE 0 // default: 0
#endif

#ifndef CGAL_ARCAVOID_DISABLE_NEWTON_REFINEMENT
# define CGAL_ARCAVOID_DISABLE_NEWTON_REFINEMENT 1 // default: 0
#endif

#ifndef CGAL_DEBUG_ARCAVOID
# define CGAL_DEBUG_ARCAVOID 0
#endif

#if CGAL_AK_D_SHOW_COMPILE_OPTIONS_AS_WARNING
#warning Arcavoid: Using Absolute_void with parameters:
#if CGAL_ARCAVOID_DISABLE_DOUBLE_WITH_EXPONENT
# warning - NOT using Double with exponent
#else
# warning - using Double with exponent
#endif
#if CGAL_ARCAVOID_DISABLE_DOUBLE
# warning - NOT using double initial stage
#else
# warning - using double initial stage
#endif
#if CGAL_ARCAVOID_DISABLE_NEWTON_REFINEMENT
# warning - NOT using Newton refinement
#else
# warning - using Newton refinement
#endif
#endif

////////////////////////////////////////////////////////////////
// MAGIC NUMBERS
////////////////////////////////////////////////////////////////

#ifndef CGAL_ARCAVOID_DOUBLE_PRECISION
#if defined (DBL_MANT_DIG)
# define CGAL_ARCAVOID_DOUBLE_PRECISION DBL_MANT_DIG
#else
# define CGAL_ARCAVOID_DOUBLE_PRECISION 53
#endif
#endif

#if CGAL_AK_D_SHOW_COMPILE_OPTIONS_AS_WARNING && CGAL_ARCAVOID_DOUBLE_PRECISION != 53
# warning - CGAL_ARCAVOID_DOUBLE_PRECISION (default: DBL_MANT_DIG) is not 53.
# warning   You are sure I poor implementation am supposed to work properly on such a strange environment?
#endif

#ifndef CGAL_ARCAVOID_INITIAL_PRECISION
# define CGAL_ARCAVOID_INITIAL_PRECISION 64
#endif

#ifndef CGAL_ARCAVOID_DEFAULT_MAX_ITER
# define CGAL_ARCAVOID_DEFAULT_MAX_ITER 64
#endif

////////////////////////////////////////////////////////////////
// DEBUGGING
////////////////////////////////////////////////////////////////

#ifndef CGAL_DEBUG_ARCAVOID
# define CGAL_DEBUG_ARCAVOID 0
#endif
#ifndef CGAL_DEBUG_ARCAVOID_PRECISION
# define CGAL_DEBUG_ARCAVOID_PRECISION CGAL_DEBUG_ARCAVOID
#endif
#if CGAL_AK_D_SHOW_COMPILE_OPTIONS_AS_WARNING
#if CGAL_DEBUG_ARCAVOID
# warning - DEBUG mode
#endif
#if CGAL_DEBUG_ARCAVOID
# warning - DEBUG PRECISION mode
#endif
#endif

#if CGAL_DEBUG_ARCAVOID
# define DBG_ARCA(x) x
// Question_mark interval formatting is expensive, so don't use in release mode
// # include "formatters.h"
#else
# define DBG_ARCA(x) static_cast< void > (0)
#endif

#if CGAL_DEBUG_ARCAVOID_PRECISION
# define DBG_ARCA_PREC(x) x
#else
# define DBG_ARCA_PREC(x) static_cast< void > (0)
#endif

////////////////////////////////////////////////////////////////
// PROFILING
////////////////////////////////////////////////////////////////

#ifndef CGAL_DEBUG_ARCAVOID_TIME
# define CGAL_DEBUG_ARCAVOID_TIME CGAL_DEBUG_ARCAVOID
#endif

#if CGAL_DEBUG_ARCAVOID_TIME
# define DBG_ARCA_TIME(x) x
#else
# define DBG_ARCA_TIME(x) static_cast< void > (0)
#endif

#include <CGAL/Profile_timer.h>

#ifndef Bisolve_telemetry_code
# define Bisolve_telemetry_code(x)
#endif

////////////////////////////////////////////////////////////////
// WORKING CODE
////////////////////////////////////////////////////////////////

#if defined(__GNUC__) && (__GNUC__ >= 3)
# define CGAL_ARCAVOID_LIKELY(x) (__builtin_expect(!!(x),1))
# define CGAL_ARCAVOID_UNLIKELY(x) (__builtin_expect((x),0))
#else
# define CGAL_ARCAVOID_LIKELY(x) (x)
# define CGAL_ARCAVOID_UNLIKELY(x) (x)
#endif

#include <CGAL/Algebraic_kernel_1/Arcavoid_kernel.h>
#include <CGAL/Algebraic_kernel_1/Bitsize.h>
#include <CGAL/Algebraic_kernel_1/Bigfloat_traits.h>
#include <CGAL/Algebraic_kernel_1/Cartesian_complex.h>

#if ! CGAL_ARCAVOID_DISABLE_DOUBLE_WITH_EXPONENT
# include <CGAL/Algebraic_kernel_1/Double_with_exponent.h>
# include <boost/numeric/interval.hpp>
#endif

#include <CGAL/basic.h>
#include <CGAL/Random.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d/Real_embeddable_extension.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/Gmpfi.h>

// Convex hull calculation for initial approximations
#include <CGAL/Simple_cartesian.h>
#include <CGAL/convex_hull_2.h>

// Exception handling for machine precision floating point
#include <fenv.h>
#include <csignal>

#include <cmath>
#include <list>
#include <vector>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility.hpp>
#include <boost/detail/select_type.hpp>

// Computation of connected components (root clusters)
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/connected_components.hpp>

namespace CGAL {

class Arcavoid_real_root_isolator_tag;
class Arcavoid_complex_root_isolator_tag;
template< class BCK, class ArcavoidTag > class Arcavoid;

namespace internal {

template< class T > struct _Distance_rnd_u;

template<>
struct _Distance_rnd_u< Cartesian_complex< Gmpfr > > {
  typedef Cartesian_complex< Gmpfr > GmpfrCC;
  const Gmpfr operator() (const GmpfrCC &x, const GmpfrCC &y) {
    GmpfrCC d;
#if defined (MPFR_VERSION_MAJOR) && MPFR_VERSION_MAJOR >= 3
    mpfr_sub (d.real().fr(), x.real().fr(), y.real().fr(), MPFR_RNDA);
    mpfr_sqr (d.real().fr(), d.real().fr(), MPFR_RNDU);
    mpfr_sub (d.imag().fr(), x.imag().fr(), y.imag().fr(), MPFR_RNDA);
    mpfr_sqr (d.imag().fr(), d.imag().fr(), MPFR_RNDU);
#else
    // Rounding mode MPFR_RNDA not available in MPFR <= 2.x

    // mpfr_sub (d.real().fr(), x.real().fr(), y.real().fr(), MPFR_RNDA);
    if (x.real() >= y.real())
      mpfr_sub (d.real().fr(), x.real().fr(), y.real().fr(), MPFR_RNDU);
    else
      mpfr_sub (d.real().fr(), y.real().fr(), x.real().fr(), MPFR_RNDU);
    mpfr_sqr (d.real().fr(), d.real().fr(), MPFR_RNDU);

    // mpfr_sub (d.imag().fr(), x.imag().fr(), y.imag().fr(), MPFR_RNDA);
    if (x.imag() >= y.imag())
      mpfr_sub (d.imag().fr(), x.imag().fr(), y.imag().fr(), MPFR_RNDU);
    else
      mpfr_sub (d.imag().fr(), y.imag().fr(), x.imag().fr(), MPFR_RNDU);

    mpfr_sqr (d.imag().fr(), d.imag().fr(), MPFR_RNDU);
#endif
    mpfr_add (d.real().fr(), d.real().fr(), d.imag().fr(), MPFR_RNDU);
    mpfr_sqrt (d.real().fr(), d.real().fr(), MPFR_RNDU);
    return d.real();
  }
};

#if CGAL_USE_LEDA

template<>
struct _Distance_rnd_u< Cartesian_complex< leda_bigfloat > > {
  typedef Cartesian_complex< leda_bigfloat > ledaBFCC;
  const leda_bigfloat operator() (const ledaBFCC &x, const ledaBFCC &y) {
    ledaBFCC d;
    d.real() = leda::sub (x.real(), y.real(), leda::bigfloat::get_precision(), leda::TO_INF);
    d.real() = leda::mul (d.real(), d.real(), leda::bigfloat::get_precision(), leda::TO_INF);
    d.imag() = leda::sub (x.imag(), y.imag(), leda::bigfloat::get_precision(), leda::TO_INF);
    d.imag() = leda::mul (d.imag(), d.imag(), leda::bigfloat::get_precision(), leda::TO_INF);
    d.real() = leda::add (d.real(), d.imag(), leda::bigfloat::get_precision(), leda::TO_INF);
    d.real() = leda::sqrt (d.real(), leda::bigfloat::get_precision(), leda::TO_INF);
    return d.real();
  }
};

#endif // CGAL_USE_LEDA

template< class T > struct _Distance_rnd_d;

template<>
struct _Distance_rnd_d< Cartesian_complex< Gmpfr > > {
  typedef Cartesian_complex< Gmpfr > GmpfrCC;
  const Gmpfr operator() (const GmpfrCC &x, const GmpfrCC &y) {
    GmpfrCC d;
    mpfr_sub (d.real().fr(), x.real().fr(), y.real().fr(), MPFR_RNDZ);
    mpfr_sqr (d.real().fr(), d.real().fr(), MPFR_RNDD);
    mpfr_sub (d.imag().fr(), x.imag().fr(), y.imag().fr(), MPFR_RNDZ);
    mpfr_sqr (d.imag().fr(), d.imag().fr(), MPFR_RNDD);
    mpfr_add (d.real().fr(), d.real().fr(), d.imag().fr(), MPFR_RNDD);
    mpfr_sqrt (d.real().fr(), d.real().fr(), MPFR_RNDD);
    return d.real();
  }
};

#if CGAL_USE_LEDA

template<>
struct _Distance_rnd_d< Cartesian_complex< leda_bigfloat > > {
  typedef Cartesian_complex< leda_bigfloat > ledaBFCC;
  const leda_bigfloat operator() (const ledaBFCC &x, const ledaBFCC &y) {
    ledaBFCC d;
    d.real() = leda::sub (x.real(), y.real(), leda::bigfloat::get_precision(), leda::TO_ZERO);
    d.real() = leda::mul (d.real(), d.real(), leda::bigfloat::get_precision(), leda::TO_ZERO);
    d.imag() = leda::sub (x.imag(), y.imag(), leda::bigfloat::get_precision(), leda::TO_ZERO);
    d.imag() = leda::mul (d.imag(), d.imag(), leda::bigfloat::get_precision(), leda::TO_ZERO);
    d.real() = leda::add (d.real(), d.imag(), leda::bigfloat::get_precision(), leda::TO_ZERO);
    d.real() = leda::sqrt (d.real(), leda::bigfloat::get_precision(), leda::TO_ZERO);
    return d.real();
  }
};

#endif // CGAL_USE_LEDA

  //! evaluate at a point using the Ruffini-Horner scheme

#define CGAL_ARCAVOID_OPTIMIZED_HORNER 1
template< class Vector, class T >
struct _Horner {
  typedef CGAL::Coercion_traits< typename Vector::value_type, T > CT;
  typedef typename CT::Type Type;

  const long degree (const Vector &f) const {
    return f.size() - 1;
  }

  const Type horner (const Vector &f, const T &z) const {
    typename CT::Cast cast;

    int i = degree (f);
    if (CGAL_ARCAVOID_UNLIKELY (i < 0))
      return Type (0);

#if (! defined(CGAL_ARCAVOID_OPTIMIZED_HORNER)) || (! CGAL_ARCAVOID_OPTIMIZED_HORNER)
    Type res = cast (f[i]);
    for (--i; i >= 0; --i) {
      res *= z;
      res += f[i];
    }
#else
    const T zz = CGAL::square (z); // intentionally no cast (spares time for T real, Type complex)
    Type res = cast (f[i]);
    if (i % 2 == 1) {
      res *= z;
      res += f[--i];
    }
    while (i > 0) {
      res *= zz;
      res += f[--i] * z;
      res += f[--i];
    }
#endif

    return res;
  }

  //! evaluate f and derivative at a point using the parallel Ruffini-Horner scheme
  const boost::tuple< Type, Type > horner_2 (const Vector &f, const T &z) const {
    typename CT::Cast cast;

    int i = degree (f);

    if (CGAL_ARCAVOID_UNLIKELY (i < 1)) {
      if (i < 0)
        return boost::tuple< Type, Type > (0, 0);
      else
        return boost::tuple< Type, Type > (cast (f[0]), 0);
    }

    // TODO: optimization using CGAL::square()

    Type fz = cast (f[i]);
    Type dfz = cast (f[i]);
    for (--i; i > 0; --i) {
      fz *= z;
      fz += f[i];
      dfz *= z;
      dfz += fz;
    }
    fz *= z;
    fz += f[0];

    return boost::make_tuple (fz, dfz);
  }

  //! evaluate Newton correction at a point using the parallel Ruffini-Horner scheme
  const Type newton_correction (const Vector &f, const T &z) const {
    // if (boost::is_same< T, Cartesian_complex< double > >::value
    //     && z.squared_norm_2() > 1.)
    //   return newton_correction_via_inverse (f, z);

    const boost::tuple< Type, Type > fzs = horner_2 (f, z);
    return fzs.get<0>() / fzs.get<1>();
  }

  const Type newton_correction_via_inverse (const Vector &f, const T &z) const {
    const T z_inv = z.reciprocal();
    const Vector fR (f.rbegin(), f.rend());
    const boost::tuple< Type, Type > fzs_inv = horner_2 (fR, z_inv);

    const Type f_df_inv = fzs_inv.get<1>() / fzs_inv.get<0>();

    Type den = degree (f) * z_inv;
    den -= CGAL::square (z_inv) * f_df_inv;
    return den.reciprocal();
  }
};

template< class Interval >
struct _Rounder : public std::unary_function< Interval, void > {
  _Rounder (long) {}
  void operator() (const Interval &I) const {}
  // TODO: remove this; overload for LEDA
};

template<>
struct _Rounder< Gmpfi > : public std::unary_function< Gmpfi, void > {
  const long prec;
  _Rounder (long p) : prec (p) {}
  void operator() (Gmpfi &I) const {
    if (I.get_precision() > prec) {
      DBG_ARCA_PREC (std::cerr << "Rounding necessary (prec was: "
                     << I.get_precision() << ", should be " << prec << ")" << std::endl);
      I = I.round (prec);
    }
  }
};
template<>
struct _Rounder< CGAL::Cartesian_complex< Gmpfi > >
  : public std::unary_function< CGAL::Cartesian_complex< Gmpfi >, void > {
  typedef CGAL::Cartesian_complex< Gmpfi > GmpfiCC;
  const long prec;
  _Rounder (long p) : prec (p) {}
  void operator() (GmpfiCC &I) const {
    if (I.real().get_precision() > prec) {
      DBG_ARCA_PREC (std::cerr << "Rounding necessary (real.prec was: "
                     << I.real().get_precision() << ", should be " << prec << ")" << std::endl);
      I.real() = I.real().round (prec);
    }
    if (I.imag().get_precision() > prec) {
      DBG_ARCA_PREC (std::cerr << "Rounding necessary (imag.prec was: "
                     << I.imag().get_precision() << ", should be " << prec << ")" << std::endl);
      I.imag() = I.imag().round (prec);
    }
  }
};

template< class BCK >
class Arcavoid_list {
  //// public interface class
  friend class CGAL::Arcavoid< BCK, Arcavoid_complex_root_isolator_tag >;

  //// forward declarations
  template< class Real > class Approximation_template;
  template< class Real > class Cluster_template;

public:
  //// Self
  typedef Arcavoid_list< BCK >             List;

public:
  //// setup Bitstream_coefficient_kernel
  typedef BCK Bitstream_coefficient_kernel;
  typedef typename Bitstream_coefficient_kernel::Coefficient       Coefficient;
  typedef typename Bitstream_coefficient_kernel::Convert_to_bfi    Convert_to_bfi;
  typedef typename Bitstream_coefficient_kernel::Is_zero           Coefficient_is_zero;

  //// setup Arithmetic_kernel bigfloat types
  typedef typename Bitstream_coefficient_kernel::Arithmetic_kernel Arithmetic_kernel;
  typedef typename Arithmetic_kernel::Bigfloat                     Bigfloat;
  typedef typename Bitstream_coefficient_kernel::Bigfloat_interval Bigfloat_interval;

protected:
  //// Input types (possibly multivariate + alpha)
  typedef typename CGAL::Polynomial_type_generator< Coefficient, 1 >::Type Input_polynomial;
  typedef typename CGAL::Complex_embeddable_traits< Coefficient >::Is_complex_embeddable Input_has_complex_coefficients_tag;

#if ! CGAL_ARCAVOID_DISABLE_DOUBLE
  //// double precision FP types
  typedef double                         RR;
  typedef Interval_nt< false >           RRI;
  typedef CGAL::Cartesian_complex< RR >  CC;
  typedef CGAL::Cartesian_complex< RRI > CCI;
  typedef typename boost::detail::if_true< boost::is_same< Input_has_complex_coefficients_tag, CGAL::Tag_true >::value >::template then< CC, RR >::type   InCoeff;
  typedef typename boost::detail::if_true< boost::is_same< Input_has_complex_coefficients_tag, CGAL::Tag_true >::value >::template then< CCI, RRI >::type IInCoeff;
  typedef std::vector< RR >       RRVector;
  typedef std::vector< RRI >      RRIVector;
  typedef std::vector< CC >       CCVector;
  typedef std::vector< CCI >      CCIVector;
  typedef std::vector< InCoeff >  InVector;
  typedef std::vector< IInCoeff > IInVector;
#endif

#if ! CGAL_ARCAVOID_DISABLE_DOUBLE_WITH_EXPONENT
  //// double with exponent FP types
  typedef boost::numeric::interval< Double_with_exponent > Double_with_exponent_interval;
  typedef Double_with_exponent              DweRR;
  typedef Double_with_exponent_interval     DweRRI;
  typedef CGAL::Cartesian_complex< DweRR >  DweCC;
  typedef CGAL::Cartesian_complex< DweRRI > DweCCI;
  typedef typename boost::detail::if_true< boost::is_same< Input_has_complex_coefficients_tag, CGAL::Tag_true >::value >::template then< DweCC, DweRR >::type   DweInCoeff;
  typedef typename boost::detail::if_true< boost::is_same< Input_has_complex_coefficients_tag, CGAL::Tag_true >::value >::template then< DweCCI, DweRRI >::type DweIInCoeff;
  typedef std::vector< DweRR >       DweRRVector;
  typedef std::vector< DweRRI >      DweRRIVector;
  typedef std::vector< DweCC >       DweCCVector;
  typedef std::vector< DweCCI >      DweCCIVector;
  typedef std::vector< DweInCoeff >  DweInVector;
  typedef std::vector< DweIInCoeff > DweIInVector;
#endif

  //// multiprecision FP types
  typedef Bigfloat                         MpRR;
  typedef Bigfloat_interval                MpRRI;
  typedef CGAL::Cartesian_complex< MpRR >  MpCC;
  typedef CGAL::Cartesian_complex< MpRRI > MpCCI;
  typedef typename boost::detail::if_true< boost::is_same< Input_has_complex_coefficients_tag, CGAL::Tag_true >::value >::template then< MpCC, MpRR >::type   MpInCoeff;
  typedef typename boost::detail::if_true< boost::is_same< Input_has_complex_coefficients_tag, CGAL::Tag_true >::value >::template then< MpCCI, MpRRI >::type MpIInCoeff;
  typedef std::vector< MpRR >       MpRRVector;
  typedef std::vector< MpRRI >      MpRRIVector;
  typedef std::vector< MpCC >       MpCCVector;
  typedef std::vector< MpCCI >      MpCCIVector;
  typedef std::vector< MpInCoeff >  MpInVector;
  typedef std::vector< MpIInCoeff > MpIInVector;

  template< class T > struct MpIInCoeffToInterval
    : public std::unary_function< T, RRI > {
    const RRI operator() (const T &I) const {
      return RRI (CGAL::to_interval (I));
    }
  };
  template< class T > struct MpIInCoeffToInterval< CGAL::Cartesian_complex< T > >
    : public std::unary_function< CGAL::Cartesian_complex< T > , CCI > {
    const CCI operator() (const CGAL::Cartesian_complex< T > &I) const {
      return CCI (RRI (CGAL::to_interval (I.real())),
                  RRI (CGAL::to_interval (I.imag())));
    }
  };

protected:
  //// low precision clustering structures
  typedef Approximation_template< RR > LpApproximation;
  typedef Cluster_template< RR >       LpCluster;
  typedef std::list< LpApproximation > LpApproximation_list;
  typedef std::list< LpCluster >       LpCluster_list;
  typedef typename LpApproximation_list::iterator       LpApproximation_iterator;
  typedef typename LpApproximation_list::const_iterator LpApproximation_const_iterator;
  typedef typename LpCluster_list::iterator             LpCluster_iterator;
  typedef typename LpCluster_list::const_iterator       LpCluster_const_iterator;
  typedef std::pair< LpCluster_iterator, LpCluster_iterator >             LpCluster_range;
  typedef std::pair< LpCluster_const_iterator, LpCluster_const_iterator > LpCluster_const_range;

public:
  //// multiprecision clustering structures
  typedef Approximation_template< MpRR > Approximation;
  typedef Cluster_template< MpRR >       Cluster;
  typedef std::list< Approximation >     Approximation_list;
  typedef std::list< Cluster >           Cluster_list;
  typedef typename Approximation_list::iterator       Approximation_iterator;
  typedef typename Approximation_list::const_iterator Approximation_const_iterator;
  typedef typename Cluster_list::iterator             Cluster_iterator;
  typedef typename Cluster_list::const_iterator       Cluster_const_iterator;
  typedef std::pair< Cluster_iterator, Cluster_iterator >             Cluster_range;
  typedef std::pair< Cluster_const_iterator, Cluster_const_iterator > Cluster_const_range;

protected:
  //// RAII rounding mode protector for Interval_nt (low precision)
  typedef typename RRI::Protector Protector;

  //// RAII floating point environment lock
  class Floating_point_environment_handler : boost::noncopyable {
    enum { FE_EXCEPTIONS = (FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW) };
    fenv_t old_environment;
  public:
    Floating_point_environment_handler () {
      feholdexcept (&old_environment);
      DBG_ARCA (std::clog << "Installed floating point environment handler" << std::endl);
      CGAL_assertion (fetestexcept (FE_EXCEPTIONS) == 0);
    }
    ~Floating_point_environment_handler () {
      fesetenv (&old_environment);
      DBG_ARCA (std::clog << "Cleared floating point exception handling" << std::endl);
    }
    static const bool exception_occured () {
      DBG_ARCA (if (fetestexcept (FE_EXCEPTIONS) != 0)
                  std::clog << "Floating point exception " << fetestexcept (FE_ALL_EXCEPT) << " occured..." << std::endl);
      return fetestexcept (FE_EXCEPTIONS) != 0;
    }
  };

  //// RAII precision lock for bigfloats and BF intervals
  class Precision_guard : boost::noncopyable {
    long old_prec, old_iv_prec;
  public:
    Precision_guard (long prec)
      : old_prec (CGAL::get_default_precision< MpRR >()),
        old_iv_prec (CGAL::get_precision (Bigfloat_interval())) {
      DBG_ARCA_PREC (std::clog << "Requested precision " << prec << std::endl);
      CGAL::set_default_precision< MpRR > (prec);
      CGAL::set_precision (Bigfloat_interval(), prec);
    }
    ~Precision_guard () {
      CGAL::set_default_precision< MpRR > (old_prec);
      CGAL::set_precision (Bigfloat_interval(), old_iv_prec);
      DBG_ARCA_PREC (std::clog << "Restored precision to " << old_prec << " (BF) and " << old_iv_prec << " (BF interval)" << std::endl);
    }
  };

private:
  template< class Real > class Cluster_template {
    friend class Arcavoid_list< Bitstream_coefficient_kernel >;
    typedef Cluster_template< Real > Self;

  public:
    typedef Real                            Radius_type;
    typedef CGAL::Cartesian_complex< Real > Center_type;

  protected:
    typedef std::list< Approximation_template< Real > >  _Approximation_list;
    typedef typename _Approximation_list::iterator       _Approximation_iterator;
    typedef typename _Approximation_list::const_iterator _Approximation_const_iterator;

    //// Data members
    _Approximation_list discs;
    Center_type scntr;
    Radius_type srad;
    size_t k;
    long prec;

    /*
    unsigned char isolating : 1;
    unsigned char exact : 1;
    unsigned char needs_higher_precision : 1;
    unsigned char restarted : 1;
    unsigned char touch_real : 1;
    unsigned char touch_origin : 1;
    */

    Cluster_template (long prec)
      : srad (-1),
        k(0),
        prec (prec) {}

  public:
    const size_t multiplicity () const { return k; }
    const Center_type & center () const { return scntr; }
    const Radius_type & radius () const {
      CGAL_precondition (! CGAL::is_negative (srad));
      return srad;
    }
    const long precision () const { return prec; }

    const _Approximation_iterator begin () { return discs.begin(); }
    const _Approximation_const_iterator begin () const { return discs.begin(); }
    const _Approximation_iterator end () { return discs.end(); }
    const _Approximation_const_iterator end () const { return discs.end(); }

    const bool touch_real () const { // TODO: cache this
      return CGAL::abs (center().imag()) <= radius();
    }
    const bool touch_origin () const { // TODO: cache this
      return CGAL::abs (center()) <= radius(); // TODO: rounding
    }

    // TODO: rename to Less_real
    struct Compare_real : public std::binary_function< bool, Self, Self > {
      const bool operator() (const Self &c1, const Self &c2) const {
        // TODO: check for overlaps etc.
        return c1.center().real() < c2.center().real();
      }
    };
    struct Touch_real : public std::unary_function< bool, Self > {
      const bool operator() (const Self &c) const {
        return c.touch_real();
      }
    };
    struct Touch_origin : public std::unary_function< bool, Self > {
      const bool operator() (const Self &c) const {
        return c.touch_origin();
      }
    };
    struct Cluster_multiplicity_add : public std::binary_function< int, Self, int > {
      const int operator() (int result, const Self &c) const {
        return result + c.multiplicity();
      }
    };
    struct Cluster_critical_multiplicity : public std::unary_function< Self, bool > {
      const bool operator()(const Self &c) const {
        return (c.multiplicity() > 1);
      }
    };

  protected:
    void compute_k_center_radius () {
      // Precondition: discs have center and radius computed
      // TODO: rounding
      // TODO: weighted center?
      k = discs.size();
      CGAL_precondition (k > 0);
      scntr = Center_type::ZERO();
      for (_Approximation_const_iterator it = discs.begin(); it != discs.end(); ++it)
        scntr += it->center();
      scntr /= Radius_type (static_cast< long > (k));
      srad = Radius_type (-1);
      for (_Approximation_const_iterator it = discs.begin(); it != discs.end(); ++it)
        srad = CGAL::max (srad, CGAL::abs (it->z - scntr) + it->rad);

      CGAL_postcondition (! CGAL::is_negative (srad));
    }

    friend std::ostream & operator<< (std::ostream &out,  const Self &c) {
      out << "  Cluster" << std::endl
          << "    Multiplicity: " << c.multiplicity() << std::endl
          << "    Center: " << c.center() << std::endl
          << "    Radius: " << c.srad
          << " <= 2^" << CGAL::abs_ilog2 (c.srad) << std::endl
          << "    Touch real: " << c.touch_real() << std::endl
          << "    Discs:" << std::endl;
      for (_Approximation_const_iterator jt = c.discs.begin(); jt != c.discs.end(); ++jt)
        out << "    " << *jt << std::endl;
      return out;
    }
  };

  template< class Real > class Approximation_template {
    friend class Arcavoid_list< Bitstream_coefficient_kernel >;
    typedef Approximation_template< Real > Self;

  public:
    typedef Real                            Radius_type;
    typedef CGAL::Cartesian_complex< Real > Center_type;

  protected:
    /* Data members */
    Center_type z;
    Radius_type rad;
    bool exact;

    Approximation_template ()
      : z (Center_type::ZERO()),
        rad (-1),
        exact (false) {}

    Approximation_template (const Center_type &z)
      : z (z),
        rad (-1),
        exact (false) {}

    const Center_type & center () const { return z; }
    const Radius_type & radius () const {
      CGAL_precondition (! CGAL::is_negative (rad));
      return rad;
    }

    friend std::ostream & operator<< (std::ostream &out,  const Self &D) {
      return out << "  Disc (" << CGAL::oformat (D.center()) << "; "
                 << CGAL::oformat (D.rad) << " <= 2^" << CGAL::abs_ilog2 (D.rad) << ")";
    }
  };

public:
  struct Overlapping_regions {
    template< class OutputIterator >
    void operator() (const Arcavoid_list &lhs, const Arcavoid_list &rhs,
                     OutputIterator out) const {
      // TODO: check rounding
      for (Cluster_const_iterator it = lhs.clusters.begin(); it != lhs.clusters.end(); ++it)
        for (Cluster_const_iterator jt = rhs.clusters.begin(); jt != rhs.clusters.end(); ++jt)
          if (CGAL::abs (it->center() - jt->center()) <= it->radius() + jt->radius())
            *out++ = std::make_pair (it, jt);
    }
  };

private:
  MpIInVector Fiv;        // interval approximation of input & derivative
  MpCCIVector Fwiv;             // widened interval approximation of input
  MpInVector F;             // approximation of input & derivative

  IInVector fiv;
  CCIVector fwiv;
  InVector f;

  int N, n, mult_zero;          // N = deg(f)           n = deg(f) - mult_zero
                                // N = n = mult_zero = -1  iff  f = 0
  long prec;

  Cluster_list clusters;
  LpCluster_list lp_clusters;

  enum Stage {
    PRE_INITIALIZATION = 1,
    INITIALIZATION_COMPLETED = 2,
    LOW_PRECISION_STAGE_COMPLETED = 4,
    ISOLATION_COMPLETED = 8,
  };
  Stage stage;

  Random random;

  //// temporaries
  Approximation_list eps_nbh;

public:
  ////////////////////////////////////////////////////////////////
  // Constructors
  ////////////////////////////////////////////////////////////////
  Arcavoid_list ()
    : N (-1),
      n (-1),
      mult_zero (-1),
      prec (0),
      stage (ISOLATION_COMPLETED),
      approximation_cache() {
  }

  // TODO: make bck default to BCK()
  template< class InputIterator >
  Arcavoid_list (const Bitstream_coefficient_kernel &_bck, InputIterator first, InputIterator beyond)
    : prec (0),
      stage (PRE_INITIALIZATION),
      approximation_cache (Input_polynomial (first, beyond), _bck)
  {
      initialize();
  }

  // TODO: make bck default to BCK()
  // convenience constructor from polynomial
  template< class Poly >
  Arcavoid_list (const Bitstream_coefficient_kernel &_bck, const Poly &_f)
    : prec (0),
      stage (PRE_INITIALIZATION),
      approximation_cache (_f, _bck)
  {
    initialize();
  }

private:
  ////////////////////////////////////////////////////////////////
  // Polynomial_type rewrites for std::vector< Coeff >
  // (long-term) TODO:
  //   make CGAL::Polynomial< NT > support interval polynomials
  //   and replace these by the corresponding member functions
  ////////////////////////////////////////////////////////////////

  //! get degree
  template< class T > static const int degree (const std::vector< T > &f) {
    return f.size() - 1;
  }
  template< class T, int d > static const int degree (const typename CGAL::Polynomial_type_generator< T, d >::Type &f) {
    return f.degree();
  }

  //! get leading coefficient
  template< class T > static const T & lcoeff (const std::vector< T > &f) {
    CGAL_precondition (degree (f) >= 0);
    return f.back();
  }
  template< class T, int d > static const T & lcoeff (const typename CGAL::Polynomial_type_generator< T, d >::Type &f) {
    CGAL_precondition (degree (f) >= 0);
    return f.lcoeff();
  }

  //! differentiate
  template< class Vector > static const Vector differentiate (const Vector &f) {
    if (CGAL_ARCAVOID_UNLIKELY (degree (f) <= 0))
      return Vector();

    Vector df (CGAL::cpp0x::next(f.begin()), f.end());
    for (size_t i = 1; i < df.size(); ++i)
      df[i] *= typename Vector::value_type (i+1);
    return df;
  }

  //! evaluate at a point using the Ruffini-Horner scheme
  template< class Vector, class T >
  static const typename CGAL::Coercion_traits< typename Vector::value_type, T >::Type horner (const Vector &f, const T &z) {
    return _Horner< Vector, T >().horner (f, z);
  }

  //! evaluate f and derivative at a point using the parallel Ruffini-Horner scheme
  template< class Vector, class T >
  static const boost::tuple< typename CGAL::Coercion_traits< typename Vector::value_type, T >::Type,
                             typename CGAL::Coercion_traits< typename Vector::value_type, T >::Type > horner_2 (const Vector &f, const T &z) {
    return _Horner< Vector, T >().horner_2 (f, z);
  }

  //! evaluate Newton correction at a point using the parallel Ruffini-Horner scheme
  template< class Vector, class T >
  static const typename CGAL::Coercion_traits< typename Vector::value_type, T >::Type newton_correction (const Vector &f, const T &z) {
    return _Horner< Vector, T >().newton_correction (f, z);
  }

  //! compute truncated taylor shift to a point, using the Ruffini-Horner scheme
  template< class Vector, class T >
  static const std::vector< typename CGAL::Coercion_traits< typename Vector::value_type, T >::Type >
  truncated_taylor_shift (const Vector &f, const T &z, int k) {
    typedef typename CGAL::Coercion_traits< typename Vector::value_type, T > CT;
    typename CT::Cast cast;
    typedef typename CT::Type Type;

    CGAL_TIME_PROFILER ("truncated_taylor_shift()");

    // TODO: write optimized variant using CGAL::square()

    const int n = degree (f);
    std::vector< Type > coeffs (boost::make_transform_iterator (f.begin(), cast),
                                boost::make_transform_iterator (f.end(), cast));
    for (int j = n-1; j >= 0; --j)
      for (int i = j; i < j+k && i < n; ++i)
        coeffs[i] += z * coeffs[i+1]; // intentionally no cast (spares time for T real, Type complex)

    return std::vector< Type > (coeffs.begin(), coeffs.begin() + k);
  }

  template< class Vector >
  static void print_polynomial (std::ostream &out, const Vector &f) {
    int i = degree (f);

    if (CGAL_ARCAVOID_UNLIKELY (i < 1)) {
      if (i == 0)
        out << "(" << CGAL::oformat (f[0]) << ")";
      else
        out << "0";

      return;
    }

    out << "(" << CGAL::oformat (f[i]) << ")*x^" << i;
    for (--i; i > 1; --i) {
      if (! CGAL::certainly (CGAL::is_zero (f[i])))
        out << " + (" << CGAL::oformat (f[i]) << ")*x^" << i;
    }
    if (! CGAL::certainly (CGAL::is_zero (f[1])))
      out << " + (" << CGAL::oformat (f[1]) << ")*x";
    if (! CGAL::certainly (CGAL::is_zero (f[0])))
      out << " + (" << CGAL::oformat (f[0]) << ")";
  }

  ////////////////////////////////////////////////////////////////
  // Rounding-protected distance between complex values
  ////////////////////////////////////////////////////////////////

  static const MpRR distance_rnd_u (const MpCC &x, const MpCC &y) {
    return _Distance_rnd_u< MpCC >() (x, y);
  }
  static const MpRR distance_rnd_d (const MpCC &x, const MpCC &y) {
    return _Distance_rnd_d< MpCC >() (x, y);
  }

  ////////////////////////////////////////////////////////////////
  // Approximator
  ////////////////////////////////////////////////////////////////
private:
  struct Approximation_cache {
    Input_polynomial input;
    Bitstream_coefficient_kernel bck;
    Convert_to_bfi convert_to_bfi;
    Coefficient_is_zero coefficient_is_zero;

    int N, n, mult_zero, degree_loss;

    typename Input_polynomial::const_iterator input_begin;
    typename Input_polynomial::const_iterator input_end;

    struct Approximation_stage {
      MpIInVector Fiv;  // interval approximation of input & derivative
      MpCCIVector Fwiv; // widened interval approximation of input
      MpInVector F;     // bigfloat approximation of input & derivative
    };
    std::vector< Approximation_stage > stage_cache;

    long prec;

    Approximation_cache ()
      : input(),
        bck(),
        convert_to_bfi (bck.convert_to_bfi_object()),
        coefficient_is_zero (bck.is_zero_object()),
        input_begin (input.begin()),
        input_end (input.begin()),
        prec (0) {
      stage_cache.push_back (Approximation_stage());
    }

    Approximation_cache (const Input_polynomial &_f, const Bitstream_coefficient_kernel &_bck)
      : input (_f),
        bck (_bck),
        convert_to_bfi (bck.convert_to_bfi_object()),
        coefficient_is_zero (bck.is_zero_object()),
        input_begin (input.begin()),
        input_end (input.begin()),
        prec (0) {
      if (input.degree() > 0) {
        trim_zero_coefficients();
        approximate_input (CGAL_ARCAVOID_INITIAL_PRECISION);
      } else {
        N = n = mult_zero = input.degree();
        degree_loss = 0;
        stage_cache.push_back (Approximation_stage());
      }
    }

    static const size_t precision_to_index (long p) {
      CGAL_precondition (p % CGAL_ARCAVOID_INITIAL_PRECISION == 0);
      p /= CGAL_ARCAVOID_INITIAL_PRECISION;
      size_t log_p = 0;
      while (p > 1) {
        CGAL_precondition (p % 2 == 0);
        p >>= 1;
        ++log_p;
      }
      return log_p;
    }

    // TODO: promote ..._coeff_may_vanish to constructors
    void trim_zero_coefficients (bool leading_coeff_may_vanish = false,
                                 bool constant_coeff_may_vanish = true) {
      CGAL_precondition (leading_coeff_may_vanish
                         || input.degree() < 0
                         || ! coefficient_is_zero (input.lcoeff()));
      CGAL_precondition (constant_coeff_may_vanish
                         || input.degree() < 0
                         || ! coefficient_is_zero (input[0]));

      degree_loss = 0;
      N = input.degree();
      if (leading_coeff_may_vanish) {
        while (N >= 0 && coefficient_is_zero (input[N])) {
          --N;
          ++degree_loss;
        }
      }

      mult_zero = 0;
      if (constant_coeff_may_vanish) {
        // TODO: increase mult_zero adaptively during the algorithm
        while (mult_zero < N && coefficient_is_zero (input[mult_zero]))
          ++mult_zero;
      }

      input_begin = input.begin();
      std::advance (input_begin, mult_zero);
      input_end = input.end();
      std::advance (input_end, -degree_loss);

      n = N - mult_zero;
    }

    //! Set Fiv, F, to prec-approximations
    void approximate_input (long p) {
      CGAL_TIME_PROFILER ("approximate_input()");

      if (p <= prec)
        return;

      prec = p;
      Precision_guard guard (prec);
      // const size_t index = precision_to_index (p);
      // CGAL_assertion (index == stage_cache.size());
      CGAL_assertion (precision_to_index (p) == stage_cache.size());
      stage_cache.push_back (Approximation_stage());

      Approximation_stage &s = stage_cache.back();

      {
        CGAL_TIME_PROFILER ("approximate_input(): BCK::convert_to_bfi()");
        Bisolve_telemetry_code(t_arca_transform.start());
        s.Fiv = MpIInVector (boost::make_transform_iterator (input_begin, convert_to_bfi),
                             boost::make_transform_iterator (input_end, convert_to_bfi));
        _Rounder< MpIInCoeff > rounder (prec);
        std::for_each (s.Fiv.begin(), s.Fiv.end(), rounder);
        Bisolve_telemetry_code(t_arca_transform.stop());
      }

      const MpRR delta = CGAL::ipow2< MpRR > (-prec);
      const MpCCI delta_iv = MpCCI (MpRRI (-delta, delta), MpRRI (-delta, delta));
      DBG_ARCA (std::clog << "Using delta = " << delta << std::endl);

      // get widened interval approximation
      s.Fwiv = MpCCIVector (s.Fiv.begin(), s.Fiv.end());
      for (int i = 0; i <= n; ++i)
        s.Fwiv[i] += s.Fiv[i] * delta_iv;

      // get multiprecision approximations
      typename CGAL::Interval_traits< MpIInCoeff >::Median median;
      s.F = MpInVector (boost::make_transform_iterator (s.Fiv.begin(), median),
                        boost::make_transform_iterator (s.Fiv.end(), median));
    }
  };
  Approximation_cache approximation_cache;

  void _swap_approximations (typename Approximation_cache::Approximation_stage &stage) {
    Fiv.swap  (stage.Fiv);
    Fwiv.swap (stage.Fwiv);
    F.swap    (stage.F);
  }
  void swap_approximations (long p) {
    if (p == prec) return;

    approximation_cache.approximate_input (p);
    const size_t old_index = Approximation_cache::precision_to_index (prec);
    _swap_approximations (approximation_cache.stage_cache[old_index]);
    const size_t new_index = Approximation_cache::precision_to_index (p);
    _swap_approximations (approximation_cache.stage_cache[new_index]);

    prec = p;
  }

public:
  //! scale coefficients of multiprecision approximations s.t. exponent range is centered around 0
  void normalize_exponent_range () {
    long e_min = CGAL::abs_ilog2 (Fiv[0]);
    long e_max = e_min;
    for (int i = 1; i <= n; ++i) {
      if (! CGAL::certainly (CGAL::is_zero (Fiv[i]))) {
        const long e = CGAL::abs_ilog2 (Fiv[i]);
        e_min = min (e_min, e);
        e_max = max (e_max, e);
      }
    }

    // magic constant for max exponent after normalization
    // (allows underflow for near-zero coefficients)
    // const long e_norm = min (-(e_min + e_max) / 2, -e_max + 512);
    const long e_norm = -(e_min + e_max) / 2;
    if (e_norm == 0) return;

    DBG_ARCA (std::cerr << "Rescaling by 2^" << e_norm << std::endl);
    const MpRR scale = CGAL::ipow2< MpRR > (e_norm);

    for (int i = 0; i <= n; ++i) Fiv[i] *= scale;
    for (int i = 0; i <= n; ++i) Fwiv[i] *= scale;
    for (int i = 0; i <= n; ++i) F[i] *= scale;

    DBG_ARCA (std::cerr << "scaled F = ";
              print_polynomial (std::cerr, F);
              std::cerr << std::endl << "Fiv = ";
              print_polynomial (std::cerr, Fiv);
              std::cerr << std::endl
              << "Mult_zero = " << mult_zero << std::endl);
  }

  void initialize () {

    N = approximation_cache.N;
    n = approximation_cache.n;
    mult_zero = approximation_cache.mult_zero;

    if (N <= 0) {
      // constant polynomial
      stage = ISOLATION_COMPLETED;
      return;
    }

    if (N == mult_zero) {
      stage = ISOLATION_COMPLETED;
      clusters.push_front (Cluster (prec = 0));
      for (int i = 0; i < mult_zero; ++i) {
        clusters.front().discs.push_front (Approximation (MpCC::ZERO()));
        clusters.front().discs.front().rad = MpRR (0);
        clusters.front().discs.front().exact = true;
      }

      clusters.front().k = mult_zero;
      clusters.front().scntr = MpCC::ZERO();
      clusters.front().srad = MpRR (0);
      return;
    }

    _swap_approximations (approximation_cache.stage_cache.front());

    DBG_ARCA (std::cerr << "Initialize F = ";
              print_polynomial (std::cerr, F);
              std::cerr << std::endl << "Fiv = ";
              print_polynomial (std::cerr, Fiv);
              std::cerr << std::endl
              << "Mult_zero = " << mult_zero << std::endl);

    normalize_exponent_range();

    // before Protector is installed, since FE_UPWARD renders std::pow not working
    // on certain (my) architectures...
    initialize_roots();

    Protector protector;
    CGAL::To_double< MpInCoeff > incoeff_to_double;
    MpIInCoeffToInterval< MpIInCoeff > iincoeff_to_interval;
    MpIInCoeffToInterval< MpCCI > mpcci_to_interval;
    std::vector< CCI > CCIs;

    fiv =  IInVector (boost::make_transform_iterator (Fiv.begin(), iincoeff_to_interval),
                      boost::make_transform_iterator (Fiv.end(), iincoeff_to_interval));
    fwiv = CCIVector (boost::make_transform_iterator (Fwiv.begin(), mpcci_to_interval),
                      boost::make_transform_iterator (Fwiv.end(), mpcci_to_interval));
    f =  InVector (boost::make_transform_iterator (F.begin(), incoeff_to_double),
                   boost::make_transform_iterator (F.end(), incoeff_to_double));

    DBG_ARCA (std::cerr << "double approximation f of F = ";
              print_polynomial (std::cerr, f);
              std::cerr << std::endl);

    prec = 0;
    stage = INITIALIZATION_COMPLETED;

    print_state();
  }

  void initialize_roots () {
    CGAL_precondition (clusters.empty());

    if (n <= 0)
      return;

    clusters.push_front (Cluster (prec));

    get_initial_approximations (Fiv, std::front_inserter (clusters.front().discs));
    CGAL_assertion (clusters.front().discs.size() == static_cast< size_t >(n));

    const long fb_log2 = log2_fujiwara_bound (Fiv);
    const MpRR fb = CGAL::ipow2< MpRR > (fb_log2);
    const MpRR fb_2 = CGAL::ipow2< MpRR > (fb_log2 - 1);

    const Approximation_iterator beyond = clusters.front().end();
    for (Approximation_iterator it = clusters.front().begin(); it != beyond; ++it)
      it->rad = fb_2; // may not be correct, but is not publicly visible, so it doesn't matter

    // add zeros
    for (int i = 0; i < mult_zero; ++i) {
      clusters.front().discs.push_front (Approximation (MpCC::ZERO()));
      clusters.front().discs.front().rad = MpRR (0);
      clusters.front().discs.front().exact = true;
    }
    CGAL_postcondition (clusters.front().discs.size() == static_cast< size_t >(N));

    clusters.front().k = N;
    clusters.front().scntr = MpCC::ZERO();
    clusters.front().srad = fb;

    DBG_ARCA (std::cerr << "Initialized root approximations" << std::endl);
  }

  template< class Vector >
  static const long log2_fujiwara_bound (const Vector &f) {
    typename Real_embeddable_extension< MpRRI >::Ceil_log2_abs ceil_log2_abs;
    typename Real_embeddable_extension< MpRRI >::Floor_log2_abs floor_log2_abs;

    int n = degree (f);
    CGAL_precondition (n > 0);
    CGAL_precondition (! CGAL::possibly (CGAL::is_zero (f[n])));

    const long log2_lcoeff = floor_log2_abs (CGAL::abs (f[n]));

    long fb_log2 = std::numeric_limits< long >::min();

    for (int i = 1; i < n; ++i)
      if (! CGAL::certainly (CGAL::is_zero (CGAL::abs (f[n-i]))))
        fb_log2 = max (fb_log2, (ceil_log2_abs (CGAL::abs (f[n-i])) - log2_lcoeff) / i);
    if (! CGAL::certainly (CGAL::is_zero (CGAL::abs (f[0]))))
      fb_log2 = max (fb_log2, (ceil_log2_abs (CGAL::abs (f[0])) - log2_lcoeff) / n + 1);

    return fb_log2 + 2;
  }

  template< class Poly, class OutputIterator > static void get_initial_approximations (const Poly &f, OutputIterator out) {
    CGAL_TIME_PROFILER ("get_initial_approximations()");

    CGAL_precondition (degree (f) > 0);
    size_t n = degree (f);

    Precision_guard guard (CGAL_ARCAVOID_DOUBLE_PRECISION);

    // 1. generate point set (i, log_2 |f_i|)
    typedef CGAL::Simple_cartesian< double > K;
    typedef K::Point_2 Point_2;
    std::vector< Point_2 > points, ch;
    for (size_t i = 0; i <= n; ++i)
      if (! CGAL::certainly (CGAL::is_zero (f[i])))
        points.push_back (Point_2 (i, abs_ilog2 (f[i])));

    // 2. compute upper hull
    CGAL::upper_hull_points_2 (points.begin(), points.end(), std::back_inserter (ch));
    ch.push_back (points.front());
    std::reverse (ch.begin(), ch.end());

    // 3. choose approximations on concentrical circles
    const double offset = (double)(std::rand()) / RAND_MAX * 2 * M_PI;

    for (size_t i = 1; i < ch.size(); ++i) {
      int k = (int)(to_double (ch[i-1].x()));
      int K = (int)(to_double (ch[i].x()));

      // TODO: rewrite when using Double_with_exponent
      // detour for rounding to double precision
      MpRR fkK = CGAL::median (CGAL::abs (f[k]) / CGAL::upper (CGAL::abs (f[K])));
      if (CGAL_ARCAVOID_UNLIKELY (CGAL::is_zero (fkK)))
        fkK = CGAL::upper (CGAL::abs (f[k]) / CGAL::abs (f[K]));
      if (CGAL_ARCAVOID_UNLIKELY (fkK != fkK))
        fkK = CGAL::lower (CGAL::abs (f[k]) / CGAL::abs (f[K]));

      const std::pair< double, long > d_e = CGAL::to_double_exponent (fkK);
      MpRR u = d_e.first;
      CGAL::mult_by_pow2< MpRR > (u, d_e.second);

      if ((K-k) == 2)
        u = sqrt (u);
      else if ((K-k) > 2)
        u = root_d (u, K-k, CGAL_ARCAVOID_DOUBLE_PRECISION);

      for (int j = 0; j < K-k; ++j) {
        const double angle = 2. * M_PI * (((double)j)/(K-k)  + ((double)i)/n) + offset;
        const MpCC z = u * MpCC (std::cos (angle), std::sin (angle));
        *out++ = z;
      }
    }
  }

  template< class VectorIn, class VectorIIn, class VectorCCI >
  void aberth_in_cluster (Cluster &c, const VectorIn &F, const VectorIIn &Fiv, const VectorCCI &Fwiv,
                          bool skip_radius_computation, bool check_first, int max_iter = 1000) {
#ifdef CGAL_DEBUG_ARCAVOID_TIME
    CGAL::Timer t_aberth;
    t_aberth.start();
#endif // CGAL_DEBUG_ARCAVOID_TIME
    CGAL_TIME_PROFILER ("aberth_in_cluster ()");

    CGAL_precondition (eps_nbh.empty());
    int count = 0, icount = 0;

    DBG_ARCA_PREC (std::cerr << "aberth_in_cluster() of precision c.prec = " << c.prec << std::endl);
    Precision_guard guard (c.prec);

    if (check_first) {
      CGAL_TIME_PROFILER ("  aberth_in_cluster () -> pre1. check for eps nbh");
      Approximation_iterator it = c.discs.begin();
      while (it != c.discs.end()) {
        if (CGAL::possibly (CGAL::is_zero (horner (Fwiv, it->z)))) {
          Approximation_iterator old_it = it++;
          eps_nbh.splice (eps_nbh.begin(), c.discs, old_it);
        } else {
          ++it;
        }
      }
    }

    do {
      Approximation_iterator it = c.discs.begin();
      while (it != c.discs.end()) {
        ++icount;

        // 1. Compute Newton correction
        MpCC fz, dfz;
        {
          CGAL_TIME_PROFILER ("  aberth_in_cluster () -> 1. compute Newton correction");
          boost::tie (fz, dfz) = horner_2 (F, it->z);

          if (CGAL_ARCAVOID_UNLIKELY (CGAL::is_zero (dfz))) {
            MpCCI fziv, dfziv;
            boost::tie (fziv, dfziv) = horner_2 (Fiv, it->z);
            if (CGAL::certainly (CGAL::is_zero (fziv))) {
              it->exact = true;
              Approximation_iterator old_it = it++;
              eps_nbh.splice (eps_nbh.begin(), c.discs, old_it);
              continue;
            } else if (CGAL::certainly (CGAL::is_zero (dfziv))) {
              // TODO: perturb z
              // const bool perturbed = false;
              // CGAL_assertion (perturbed);
              CGAL_assertion (false);
            } else {
              const MpRR &dfzr1 = CGAL::lower (dfziv.real());
              const MpRR &dfzr2 = CGAL::upper (dfziv.real());
              const MpRR &dfzi1 = CGAL::lower (dfziv.imag());
              const MpRR &dfzi2 = CGAL::upper (dfziv.imag());
              dfz = MpCC (CGAL::abs (dfzr1) > CGAL::abs (dfzr2) ? dfzr1 : dfzr2,
                          CGAL::abs (dfzi1) > CGAL::abs (dfzi2) ? dfzi1 : dfzi2);
            }
          }
        }

        const MpCC f_df = fz / dfz;

        // 2. Compute Aberth correction
        MpCC abcorr = 0;
        {
          CGAL_TIME_PROFILER ("  aberth_in_cluster () -> 2. compute Aberth correction");
          for (Approximation_const_iterator jt = c.discs.begin(); jt != c.discs.end(); ++jt)
            if (it != jt)
              abcorr += (it->z - jt->z).reciprocal();
          for (Approximation_const_iterator jt = eps_nbh.begin(); jt != eps_nbh.end(); ++jt)
            abcorr += (it->z - jt->z).reciprocal();
        }

        // 3. Apply correction
        CGAL_assertion (! CGAL::is_zero (CC (1) - f_df * abcorr)); // TODO: gracefully handle this
        const MpCC corr = f_df / (CC (1) - f_df * abcorr);      // TODO: check if inf can occur even if assertion holds
        *it = Approximation (it->z - corr);

        // 4. Check for eps neighborhoood
        {
          CGAL_TIME_PROFILER ("  aberth_in_cluster () -> 4. check for eps nbh");
          if (count % 4 == 0 // TODO: check magic constant
              && CGAL::possibly (CGAL::is_zero (horner (Fwiv, it->z)))) {
            Approximation_iterator old_it = it++;
            eps_nbh.splice (eps_nbh.begin(), c.discs, old_it);
          } else {
            ++it;
          }
        }
      }

      ++count;
    } while (! c.discs.empty() && count < max_iter);
    DBG_ARCA_TIME (std::cerr << "# Aberth iterations: " << count << " / " << icount << std::endl);

#ifdef CGAL_DEBUG_ARCAVOID_TIME
    t_aberth.stop();
    DBG_ARCA_TIME (std::cerr << "Time for Aberth: " << t_aberth.time() << " (" << t_aberth.time() / icount << " per approx" << std::endl);
#endif // CGAL_DEBUG_ARCAVOID_TIME

    // 5. Compute inclusion radii
    if (! skip_radius_computation) {
      CGAL_TIME_PROFILER ("  aberth_in_cluster () -> 5. compute inclusion radii");
      eps_nbh.splice (eps_nbh.begin(), c.discs);
      for (Approximation_iterator it = eps_nbh.begin(); it != eps_nbh.end(); ++it)
        compute_radius (it, c.prec);
    }

    c.discs.splice (c.discs.begin(), eps_nbh);
  }

  void compute_radius (Approximation_iterator dit, long p) {
    // TODO: use Tk-test as filter
    // TODO: try Neumaier-Gershgorin
    compute_radius_gershgorin (dit);
  }

  void compute_radius_gershgorin (Approximation_iterator dit) {
    // TODO: efficiently handle approxs in other clusters
    // TODO: protect rounding (by caller?)

    CGAL_TIME_PROFILER ("compute_radius_gershgorin()");

    if (CGAL_ARCAVOID_UNLIKELY (dit->exact)) {
      dit->rad = 0;
      return;
    }

    MpRR riv = CGAL::upper (CGAL::abs (horner (Fiv, dit->z)) * MpRRI(n));
    riv /= CGAL::lower (CGAL::abs (lcoeff (Fiv)));

    CGAL_postcondition_code (int found_self = 0);
    CGAL_postcondition_code (int count = 0);
    CGAL_postcondition_code (int skipped = 0);

    // TODO: compute -prod for correct rounding?
    MpRR prod = 1.;
    for (Cluster_const_iterator it = clusters.begin(); it != clusters.end(); ++it) {
      for (Approximation_const_iterator jt = it->begin(); jt != it->end(); ++jt) {
        CGAL_postcondition_code (if (dit == jt) ++found_self);
        if (dit == jt) continue;

         CGAL_postcondition_code (if (CGAL::is_zero (jt->z)) ++skipped);
         if (CGAL::is_zero (jt->z)) continue;

        CGAL_postcondition_code (++count);
        ///prod *= CGAL::abs (dit->z - jt->z); // TODO: check rounding
        //prod *= CGAL::lower (CGAL::abs (MpCCI (dit->z) - MpCCI (jt->z)));
        prod *= distance_rnd_d (dit->z, jt->z);
      }
    }

    for (Approximation_const_iterator jt = eps_nbh.begin(); jt != eps_nbh.end(); ++jt) {
      CGAL_postcondition_code (if (dit == jt) ++found_self);
      if (dit == jt) continue;

      CGAL_postcondition_code (if (CGAL::is_zero (jt->z)) ++skipped);
      if (CGAL::is_zero (jt->z)) continue;

      CGAL_postcondition_code (++count);
      //prod *= CGAL::abs (dit->z - jt->z); // TODO: check rounding
      //prod *= CGAL::lower (CGAL::abs (MpCCI (dit->z) - MpCCI (jt->z)));
      prod *= distance_rnd_d (dit->z, jt->z);
    }

    // DBG_ARCA (CGAL_postcondition_code (std::cerr << "found_self: " << found_self << " skipped: " << skipped << " count: " << count << std::endl));
    CGAL_postcondition (found_self == 1);
    CGAL_postcondition (skipped == mult_zero || skipped == 0);
    CGAL_postcondition (count == n-1);

    riv /= prod;

    dit->rad = riv;
  }

  const bool discs_pairwise_newton_isolated (const Approximation &D1, const Approximation &D2) const {
    return MpRR (2*n) * (D1.radius() + D2.radius()) // TODO: check rounding
      < distance_rnd_d (D1.center(), D2.center());
  }

  const Cluster_range split_cluster_into_newton_connected_components (Cluster_iterator cit) {
    DBG_ARCA (std::cerr << "SPLIT_CLUSTER_INTO_NEWTON_CONNECTED_COMPONENTS" << std::endl);
    CGAL_TIME_PROFILER ("split_cluster_into_newton_connected_components()");

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
    std::vector< size_t > comp (k);
    const size_t nr_comp = boost::connected_components (G, &comp[0]);

    Cluster_iterator beyond = cit;
    ++beyond;

    if (nr_comp == 1) {
      CGAL_postcondition_code (size_t kold = cit->multiplicity());
      cit->compute_k_center_radius();
      CGAL_assertion (kold == cit->multiplicity());
      return Cluster_range (cit, beyond);
    }

    Cluster_list subclusters;
    CGAL_postcondition_code (size_t ksum = 0);
    for (size_t i = 0; i < nr_comp; ++i) {
      Cluster c (cit->prec);
      for (size_t j = 0; j < k; ++j)
        if (comp[j] == i)
          c.discs.splice (c.discs.begin(), cit->discs, dits[j]);
      c.compute_k_center_radius();
      CGAL_postcondition_code (ksum += c.multiplicity());
      subclusters.push_front (c); // TODO: avoid copy
    }
    CGAL_assertion (ksum == k);

    cit = clusters.erase (cit);
    CGAL_assertion (cit == beyond);
    const Cluster_range range (subclusters.begin(), beyond);
    clusters.splice (beyond, subclusters);

    CGAL_postcondition (static_cast< size_t >(std::distance (range.first, range.second)) == nr_comp);

    return range;
  }

#if ! CGAL_ARCAVOID_DISABLE_DOUBLE
  const bool low_precision_stage () {
    CGAL_precondition (stage == INITIALIZATION_COMPLETED);
    CGAL_precondition (clusters.size() == 1);

    CGAL_TIME_PROFILER ("low_precision_stage()");

    DBG_ARCA (std::cerr << "LOW_PRECISION_STAGE" << std::endl);
    DBG_ARCA_TIME (CGAL::Timer t_aberth);
    DBG_ARCA_TIME (t_aberth.start());

    Protector protector;
    Precision_guard guard (CGAL_ARCAVOID_DOUBLE_PRECISION);
    Floating_point_environment_handler fenv_handler;

    // TODO: no conversion Multiprecision -> double precision
    Approximation_iterator begin = clusters.front().begin();
    std::advance (begin, mult_zero);
    std::vector< CC > z (n);
    std::vector< RR > r (n);
    Approximation_iterator it = begin;
    for (int i = 0; i < n; ++i) {
      z[i] = CGAL::to_double (it->z);
      ++it;
    }
    CGAL_assertion (it == clusters.front().end());

    if (CGAL_ARCAVOID_UNLIKELY (fenv_handler.exception_occured()))
      return false;

    DBG_ARCA (std::cerr << "- CONVERSION DONE" << std::endl);

    std::vector< bool > in_nbh (n, false);
    bool all_in_nbh = false;

    const int max_count = CGAL::min (4*n, CGAL_ARCAVOID_DEFAULT_MAX_ITER);
    int count = 0, icount = 0;

    while (! all_in_nbh && ++count < max_count) {
      for (int i = 0; i < n; ++i) {
        if (in_nbh[i]) continue;
        ++icount;

        CC abcorr = 0;
        for (int j = 0; j < n; ++j)
          if (i != j)
            abcorr += (z[i] - z[j]).reciprocal();
        // CC fz, dfz;
        // boost::tie (fz, dfz) = horner_2 (f, z[i]);
        // const CC f_df = fz / dfz;
        const CC f_df = newton_correction (f, z[i]);
        const CC corr = f_df / (CC (1) - f_df * abcorr);
        z[i] -= corr;
      }

      if (count >= (CGAL_ARCAVOID_DEFAULT_MAX_ITER / 8)
          && count % 8 == 0) {
        all_in_nbh = true;
        for (int i = 0; i < n; ++i) {
          if (! in_nbh[i]) {
            const CCI iv = horner (fwiv, z[i]);
            in_nbh[i] = CGAL::possibly (CGAL::is_zero (iv));
            all_in_nbh &= in_nbh[i];
          }
        }
      }

      if (CGAL_ARCAVOID_UNLIKELY (fenv_handler.exception_occured()))
        return false;
    }

    DBG_ARCA (std::cerr << "- ABERTH DONE" << std::endl);

#if defined (CGAL_PROFILE) && CGAL_PROFILE
    if (count == max_count) {
      CGAL_PROFILER ("max count in double phase reached");
    }
#endif
    DBG_ARCA_TIME (std::cerr << "# Aberth iterations: " << count << " / " << icount << std::endl);
    DBG_ARCA_TIME (t_aberth.stop());
    DBG_ARCA_TIME (std::cerr << "Time for Aberth: " << t_aberth.time() << " (" << t_aberth.time() / icount << " per approx" << std::endl);

    for (int i = 0; i < n; ++i) {
      RRI deniv = 1.;
      for (int j = 0; j < n; ++j)
        if (i != j)
          deniv *= CGAL::abs (CCI(z[i]) - CCI(z[j])); // TODO: use proper rounding here instead of interval arithmetics
      RRI riv = RRI(n) * CGAL::abs (horner (fiv, z[i])) / CGAL::abs (fiv[n]);
      riv /= deniv;
      r[i] = CGAL::upper (riv);

      // r[i] = CGAL::upper (RRI(n) * CGAL::abs (horner (fiv, z[i])) / CGAL::abs (horner (dfiv, z[i])));
    }

    if (CGAL_ARCAVOID_UNLIKELY (fenv_handler.exception_occured())) {
      stage = LOW_PRECISION_STAGE_COMPLETED;
      return false;
    }

    it = begin;
    for (int i = 0; i < n; ++i) {
      const MpCC mpcc (MpRR (z[i].real()), MpRR (z[i].imag()));
      *it = Approximation (mpcc);
      it->rad = r[i];
      ++it;
    }

    const Cluster_range range = split_cluster_into_newton_connected_components (clusters.begin());

    stage = LOW_PRECISION_STAGE_COMPLETED;

    DBG_ARCA (std::cerr << "- END OF LOW_PRECISION_STAGE" << std::endl);
    print_state();

    return true;
  }
#endif // CGAL_ARCAVOID_DISABLE_DOUBLE

public:
  const Cluster_range subdivide_cluster (Cluster_iterator cit) {
    DBG_ARCA (std::cerr << "SUBDIVIDE_CLUSTER" << std::endl);

    CGAL_TIME_PROFILER ("subdivide_cluster (total)");
    CGAL_HISTOGRAM_PROFILER ("subdivide_cluster", cit->prec);

    CGAL_assertion (stage != PRE_INITIALIZATION);

#if ! CGAL_ARCAVOID_DISABLE_DOUBLE
    if (stage == INITIALIZATION_COMPLETED) {
      CGAL_precondition (cit == clusters.begin());

      CGAL_PROFILER ("low_precision_stage called");

      if (low_precision_stage()) {
        CGAL_PROFILER ("low_precision_stage successful");
        return Cluster_range (clusters.begin(), clusters.end());
      } else {
        stage = LOW_PRECISION_STAGE_COMPLETED;
        cit->prec = 0;
      }
    }
#endif

    if (CGAL::is_zero (cit->radius())) {
      // CGAL_precondition (CGAL::is_zero (cit->center())); // TODO: check this
      // TODO: problem with precondition if exact root found
      Cluster_iterator beyond = cit;
      ++beyond;
      return Cluster_range (cit, beyond);
    }

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
      CGAL_assertion (zeros.size() == static_cast< size_t >(mult_zero));
    }

    if (cit->prec == 0)
      cit->prec = CGAL_ARCAVOID_INITIAL_PRECISION;
    else
      cit->prec <<= 1;
    //approximation_cache.approximate_input (cit->prec);
    swap_approximations (cit->prec);
    Precision_guard guard (cit->prec);

    bool restarted = true;
    do { // will be done exactly once; but ease break
      if (cit->multiplicity() == static_cast< size_t >(n)
          || cit->multiplicity() <= 1) { // TODO: check this condition
        restarted = false;
        DBG_ARCA (std::cerr << "No restart: cluster multiplicity = " << cit->multiplicity() << std::endl);
        break;
      }

      // if (cit->touch_origin()) {
      //   restarted = false;
      //   DBG_ARCA (std::cerr << "No restart: cluster touching origin" << std::endl);
      //   break;
      // }

      // try restart
      const size_t k = cit->multiplicity();
      MpCC center = cit->center();

      bool sufficient_separation = true;

      // TODO: Probably a gross overestimation - intuition: should be ( n choose k )?
      // But better err on the conservative side...

      MpRR m = 20L;
      for (size_t i = 0; i < k; ++i)
	m *= n-i;

      // The following minimal separation factor is from Bini & Fiorentino: Design of MPSolve...
      // seems not sufficient:
      //   see experimental-packages/Algebraic_kernel_2/examples/Ak_d/
      //   and run gmp-rsisol-gres-ggcd-num-bca-sc-cana on data/easy_curves/swinnerton_dyer

      // const MpRR m = 20L * (n-k+1) / MpRR (sqrt (static_cast< double > (k)));

      for (Cluster_const_iterator it = clusters.begin(); it != clusters.end(); ++it) {
        if (it == cit) continue;

        if (CGAL::abs (it->center() - center)
            <= m * (cit->radius()) + it->radius()) {
          sufficient_separation = false;
          break;
        }
      }
      if (! sufficient_separation) {
        restarted = false;
        DBG_ARCA (std::cerr << "No restart: separation too low" << std::endl);
        break;
      }

      {
        CGAL_TIME_PROFILER ("Newton iteration on cluster center");

        MpInVector k_m_1_dF = F;
        {
          CGAL_TIME_PROFILER ("  Newton iteration on cluster center -> differentiation");
          for (size_t i = 1; i < k; ++i)
            k_m_1_dF = differentiate (k_m_1_dF);
        }
        {
          CGAL_TIME_PROFILER ("  Newton iteration on cluster center -> actual evaluation");
          MpCC km1dFz, kdFz;
          for (int i = 0; i < 2; ++i) { // TODO: check magic constant 1
            boost::tie (km1dFz, kdFz) = horner_2 (k_m_1_dF, center);
            center -= km1dFz / kdFz;
          }
        }
      }

      if (CGAL::abs (center - cit->center()) >= cit->radius()) {
        restarted = false;
        DBG_ARCA (std::cerr << "No restart: refined center not in cluster" << std::endl);
        break;
      }

      DBG_ARCA (std::cerr << "Shift to " << CGAL::oformat (center) << std::endl);
      const MpCCIVector shifted = truncated_taylor_shift (Fiv, center, k+1);
      DBG_ARCA (std::cerr << "local F: ");
      DBG_ARCA (print_polynomial (std::cerr, shifted));
      DBG_ARCA (std::cerr << std::endl);

      std::vector< MpCC > appr;
      get_initial_approximations (shifted, std::back_inserter (appr));

      CGAL_assertion (appr.size() == k);
      if (! appr.size() == k) {
        restarted = false;
        DBG_ARCA (std::cerr << "DISCREPANCY: expected " << k << " approximations, got " << appr.size() << std::endl);
        break;
      }

      // get widened interval approximation
      // const MpRR delta = CGAL::ipow2< MpRR > (-prec * k/n + CGAL::abs_ilog2 (static_cast< double > (n)));
      // const MpCCI delta_iv = MpCCI (MpRRI (-delta, delta), MpRRI (-delta, delta));
      // DBG_ARCA (std::clog << "Using delta = " << delta << std::endl);
      // MpCCIVector shifted_widened (shifted.begin(), shifted.end());
      // for (int i = 0; i <= k; ++i)
      //   shifted_widened[i] += shifted_widened[i] * delta_iv;
      MpCCIVector shifted_widened = truncated_taylor_shift (Fwiv, center, k+1);
      // get multiprecision approximation
      typename CGAL::Interval_traits< MpCCI >::Median median;
      const MpCCVector median_shifted (boost::make_transform_iterator (shifted.begin(), median),
                                       boost::make_transform_iterator (shifted.end(), median));

      Approximation_iterator dit = cit->discs.begin();
      for (size_t i = 0; i < appr.size(); ++i) {
        dit->z = appr[i];
        DBG_ARCA (std::cerr << " ===> new approximation " << i << " at " << CGAL::oformat (dit->z + center) << std::endl);
        ++dit;
      }

      aberth_in_cluster (*cit, median_shifted, shifted, shifted_widened, true, true);

      {
        // DBG_ARCA (std::cerr << "Before backshift: " << std::endl);
        // DBG_ARCA (print_state (cit, CGAL::cpp0x::next(cit)));

        CGAL_TIME_PROFILER ("  aberth_in_cluster () -> backshift and 5. (REDO) compute inclusion radii");
        for (Approximation_iterator it = cit->begin(); it != cit->end(); ++it)
          it->z += center;

        DBG_ARCA (std::cerr << "After backshift: " << std::endl);
        DBG_ARCA (print_state (cit, CGAL::cpp0x::next(cit)));

        eps_nbh.splice (eps_nbh.begin(), cit->discs);
        for (Approximation_iterator it = eps_nbh.begin(); it != eps_nbh.end(); ++it)
          compute_radius (it, cit->prec);
        cit->discs.splice (cit->discs.begin(), eps_nbh);
      }
    } while (false);

    if (! restarted) {
      aberth_in_cluster (*cit, F, Fiv, Fwiv, false, false);
    }

    cit->discs.splice (cit->discs.begin(), zeros);

    const Cluster_range range = split_cluster_into_newton_connected_components (cit);

    print_state (range.first, range.second);
    return range;
  }

  /*
  void refine (Cluster_iterator cit) {
    // TODO: dummy implementation (uses subdivide as refinement step)
    // TODO: CHECK WHY DUMMY DOES NOT WORK
    Cluster_range range = subdivide_cluster (cit);
    CGAL_postcondition (cit == range.first);
    CGAL_postcondition (std::distance (range.first, range.second) == 1);
  }
  */

  const Cluster_range real_subdivide_cluster (Cluster_iterator cit) {
    DBG_ARCA (std::cerr << "REAL_SUBDIVIDE_CLUSTER in" << std::endl);
    DBG_ARCA (print_state (cit, CGAL::cpp0x::next(cit)));
    DBG_ARCA (std::cerr << "STATE:" << std::endl);
    DBG_ARCA (print_state ());
    CGAL_precondition (cit->touch_real());

    const Cluster_range range = subdivide_cluster (cit);
    CGAL_postcondition_code (int total = std::distance (range.first, range.second));

    const Cluster_iterator begin_imag
      = std::partition (range.first, range.second, typename Cluster::Touch_real());
    // TODO: triple-check if stable_partition is necessary
    Cluster_range real_range (range.first, begin_imag);
    CGAL_postcondition_code (Cluster_range imag_range (begin_imag, range.second));

    sort_real_cluster_range (real_range);

    CGAL_postcondition_code (int real = std::distance (real_range.first, real_range.second));
    CGAL_postcondition_code (int imag = std::distance (imag_range.first, imag_range.second));
    CGAL_postcondition (real + imag == total);
    CGAL_postcondition (real_range.second == imag_range.first);

    return real_range;
  }

  //! sorts clusters in real_range by increasing real part of cluster centers
  void sort_real_cluster_range (Cluster_range &real_range) {
    if (real_range.first == real_range.second) // empty range
      return;

    Cluster_list real_clusters;
    real_clusters.splice (real_clusters.begin(), clusters, real_range.first, real_range.second);
    real_clusters.sort (typename Cluster::Compare_real());
    real_range.first = real_clusters.begin();
    clusters.splice (real_range.second, real_clusters);
  }

  // TODO: write p_a_s_r (first, beyond) and use for real_subdivide_cluster (cit)
  const Cluster_iterator partition_and_sort_reals () {
    const Cluster_iterator begin_imag
      = std::partition (clusters.begin(), clusters.end(), typename Cluster::Touch_real());
    // TODO: triple-check if stable_partition is necessary
    Cluster_range real_range (clusters.begin(), begin_imag);
    sort_real_cluster_range (real_range);
    return begin_imag;
  }

  ////////////////////////////////////////////////////////////////
  // Interface & Statistics
  ////////////////////////////////////////////////////////////////

  void print_state () const {
    return print_state (clusters.begin(), clusters.end());
  }

  void print_state (Cluster_const_iterator it, Cluster_const_iterator it_end) const {
#if CGAL_DEBUG_ARCAVOID
    std::cerr << "N: " << N << " n: " << n << " mult_zero: " << mult_zero << " prec: " << prec << "\tPhase: ";
    switch (stage) {
    case PRE_INITIALIZATION:            std::cerr << "PRE_INITIALIZATION" << std::endl;            break;
    case INITIALIZATION_COMPLETED:      std::cerr << "INITIALIZATION_COMPLETED" << std::endl;      break;
    case LOW_PRECISION_STAGE_COMPLETED: std::cerr << "LOW_PRECISION_STAGE_COMPLETED" << std::endl; break;
    case ISOLATION_COMPLETED:           std::cerr << "ISOLATION_COMPLETED" << std::endl;           break;
    default:                            std::cerr << "UNKNOWN: " << stage << std::endl;            break;
    }
    for (; it != it_end; ++it) {
      std::cerr << (it->touch_real() ? "  RC:" : "  IC:");
      if (it->multiplicity() > 1)
        std::cerr << "\tCLST (" << CGAL::oformat (it->scntr) << "; " << CGAL::oformat (it->srad) << std::endl;
      for (Approximation_const_iterator jt = it->discs.begin(); jt != it->discs.end(); ++jt)
        std::cerr << "\t" << *jt << std::endl;
    }
#endif
  }

  const Cluster_iterator begin () { return clusters.begin(); }
  const Cluster_const_iterator begin () const { return clusters.begin(); }
  const Cluster_iterator end () { return clusters.end(); }
  const Cluster_const_iterator end () const { return clusters.end(); }

  const int multiplicity_of_zero () const { return mult_zero; }
  const int degree () const { return N; }

  const Cluster_range cluster_range () {
    CGAL_assertion (std::distance (clusters.begin(), clusters.end()) >= 1);
    CGAL_assertion (std::distance (clusters.begin(), clusters.end()) <= n);
    return Cluster_range (clusters.begin(), clusters.end());
  }
  const Cluster_const_range cluster_range () const {
    CGAL_assertion (std::distance (clusters.begin(), clusters.end()) >= 1);
    CGAL_assertion (std::distance (clusters.begin(), clusters.end()) <= n);
    return Cluster_const_range (clusters.begin(), clusters.end());
  }

  const Input_polynomial & polynomial () const { return approximation_cache.input; }
  const Bitstream_coefficient_kernel & kernel () const { return approximation_cache.bck; }
  const long precision () const { return prec; }
};

} // namespace internal

} // namespace CGAL

#undef DBG_ARCA
#undef DBG_ARCA_PREC

#endif // CGAL_ARCAVOID_H
