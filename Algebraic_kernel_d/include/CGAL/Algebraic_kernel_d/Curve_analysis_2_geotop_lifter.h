// Copyright (c) 2009, 2010, 2011, 2012 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     :  Eric Berberich <eric.berberich@cgal.org>

/*! \file CGAL/Algebraic_kernel_d/Curve_analysis_2_geotop_lifter.h
  \brief Defines Geotop_lifter for Curve_analysis_2
*/

#ifndef CGAL_ALGEBRAIC_KERNEL_D_CURVE_ANALYSIS_2_GEOTOP_LIFTER_H
#define CGAL_ALGEBRAIC_KERNEL_D_CURVE_ANALYSIS_2_GEOTOP_LIFTER_H

#include <CGAL/config.h>

#include <CGAL/Algebraic_kernel_2/Bi_solve_2_flags.h>

#include <CGAL/Handle_with_policy.h>

#include <CGAL/Algebraic_kernel_d/Generic_isolator.h>
#if CGAL_USE_RS3
#include <CGAL/Algebraic_kernel_d/RS_ak_1_isolator.h>
#else
// Comment: Descartes_isolator_rep would be nice, but does not (yet) provide "refine_interval"
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_isolator_reps.h>
#endif

#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Algebraic_kernel_d/Interval_evaluate_2.h>
#include <CGAL/Algebraic_kernel_d/Bi_diff_vanish_2.h>

#if CGAL_ACK_CURVE_ANALYSES_BISOLVE_USE_TEISSIER
#include <CGAL/Algebraic_kernel_2/Certifier_cofactor_arcavoid_traits.h>
#include <CGAL/Algebraic_kernel_1/Arcavoid.h>
#else
#include <CGAL/Algebraic_kernel_2/Certifier_cofactor_bitstream_traits.h>
#endif

#ifndef Bisolve_telemetry_code
# define Bisolve_telemetry_code(x)
#endif

#ifdef CGAL_XTri_USE_TIMERS
#if CGAL_AK_D_SHOW_COMPILE_OPTIONS_AS_WARNING
#warning Curve_analysis_2_geotop_lifter: Using timers
#endif
extern CGAL::Timer tm_external_timer;
#endif

namespace CGAL {

namespace internal {

  // TODO document Geotop_lifter_bisolve_rep
  template< class BitstreamCoefficientKernel, typename HandlePolicy > // no default on policy, should be decided in higher level
class Geotop_lifter_bisolve_rep :
  public Generic_isolator_rep< typename BitstreamCoefficientKernel::Polynomial,
                               typename BitstreamCoefficientKernel::Bound,
                               HandlePolicy,
                               CGAL::Tag_true > {

  //!\name Tags
  //!@{

public:
  //! "refine_interval" function is implemented
  typedef CGAL::Tag_true Refine_interval_provided_tag;

  //!@} // Tags

public:

  //!\name Public typedefs
  //!@{

  //! this instance's template parameter
  typedef BitstreamCoefficientKernel Bitstream_coefficient_kernel;

   //! the policy type
   typedef HandlePolicy Handle_policy;

  //! type of Bound
  typedef typename Bitstream_coefficient_kernel::Polynomial Polynomial;

  //! type of Bound
  typedef typename Bitstream_coefficient_kernel::Bound Bound;

  //! the base class
  typedef Generic_isolator_rep< Polynomial, Bound, Handle_policy, Refine_interval_provided_tag > Base;

  //! the class itself
  typedef Geotop_lifter_bisolve_rep< Bitstream_coefficient_kernel, Handle_policy > Self;

  //! type of algebraic kernel
  typedef typename Bitstream_coefficient_kernel::Algebraic_kernel_d_1 Algebraic_kernel_d_1;

  //! multiplicity type
  typedef typename Algebraic_kernel_d_1::Multiplicity_type Multiplicity_type;

#if CGAL_ACK_CURVE_ANALYSES_BISOLVE_USE_TEISSIER
  //! type of certifier traits
  typedef CGAL::internal::Certifier_cofactor_arcavoid_traits< Algebraic_kernel_d_1 > Certifier_traits;
#else
  typedef CGAL::internal::Certifier_cofactor_bitstream_traits< Algebraic_kernel_d_1 > Certifier_traits;
#endif

  //! type of bi diff vanish
  typedef CGAL::Bi_diff_vanish_2< Certifier_traits > Bi_diff_vanish_2;

  //! type of bi algebraic real
  typedef typename Bi_diff_vanish_2::Algebraic_real_2 Algebraic_real_2;

  //! type of bi algebraic real
  typedef typename Bi_diff_vanish_2::Polynomial_2 Polynomial_2;

#if !CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS
  //! active interval set const iterator
  typedef typename Bi_diff_vanish_2::Certifier::Active_intervals_set Active_intervals_set;

  //! active interval set iterator
  typedef typename Bi_diff_vanish_2::Certifier::Active_intervals_set_iterator
    Active_intervals_set_iterator;

  //! active interval set const iterator
  typedef typename Bi_diff_vanish_2::Certifier::Active_intervals_set_const_iterator
    Active_intervals_set_const_iterator;
#else
#error Geotop_lifter needs active intervals
#endif

  //!@} // Public typedefs

public:

  //!\name Constructors
  //!@{

  //! Default constructor
#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_DEFAULT_CONSTRUCTIBLE
protected:
  Geotop_lifter_bisolve_rep(); // = disable
#else
public:
  Geotop_lifter_bisolve_rep() : Base(new Rep()) {}
#endif

  // needs no special assignable-implementation as no pointers
  // and this is a Rep of a handle class
#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_ASSIGNABLE
  protected:
    // Geotop_lifter_bisolve_rep(const Self&) // = default
    Geotop_lifter_bisolve_rep& operator=(const Self&); // = disable
#endif

 public:
  //! standard constructor
  Geotop_lifter_bisolve_rep(const Bitstream_coefficient_kernel& bck,
                            Multiplicity_type mult,
                            const Bi_diff_vanish_2& bdv,
                            bool local) :
    Base(Polynomial(bdv.f().begin(), bdv.f().end())),
    _m_bck(bck),
    _m_bi_diff_vanish(bdv),
    _m_mult(mult),
    _m_local(local)
  {
#if CGAL_ACK_DEBUG_FLAG
      CGAL_ACK_DEBUG_PRINT << "Use bdv-isolator for status line at " << CGAL::to_double(this->_m_bck.alpha()) << std::endl;
#endif
    // std::cout << "FFYisol poly: " << _m_bi_diff_vanish.f() << std::endl;
  }

  //!@} // Constructors

public:

  //!\name Access functions
  //!@{

  void isolate() {

    if (Base::_m_is_isolated) {
      return;
    }

    // std::cout << "FFYisol 1" << std::endl;

    _m_bi_diff_vanish._achieve(Bi_diff_vanish_2::CGAL_BISOLVE_ALL_BUT_LIFTING);

    // std::cout << "FFYisol 1a" << std::endl;

    bool has_fiber = _m_bi_diff_vanish.has_fiber(this->_m_bck.alpha());

    // std::cout << "FFYisol has_fiber: " << has_fiber << std::endl;

    // std::cout << "FFYisol 1b" << std::endl;

    // _m_diff is not needed if non-event line
    if (has_fiber) {
      // std::cout << "FFYisol 2" << std::endl;
      _m_bi_diff_vanish(this->_m_bck.alpha(), _m_mult, std::back_inserter(_m_diff), _m_local);
      // std::cout << "FFYisol mdiffsz: " << _m_diff.size() << std::endl;

      // isolate "roots" wrt to 'simple' solutions
      for (size_t i = 0; i < _m_diff.size(); i++) {

        Polynomial_2 fym;
        // std::cout << "rx: " << CGAL::to_double(_m_diff[i].first.x()) << std::endl;
        // std::cout << "ry: " << CGAL::to_double(_m_diff[i].first.y()) << std::endl;
        // std::cout << "multr: " << _m_diff[i].second << std::endl;
        if (_m_diff[i].second + 1< _m_bi_diff_vanish.max_diff()) {
          fym = _m_bi_diff_vanish.diff(_m_diff[i].second + 1);
        } else {
          fym = CGAL::canonicalize(CGAL::differentiate(
                _m_bi_diff_vanish.diff(_m_diff[i].second)));
        }
        //std::cout << "fym: " << fym << std::endl;

        int prec = 2;

        while (true) { // refine wrt f_{yˆm}
          typename internal::Interval_evaluate_2< Polynomial_2, Bound > ieval_2;

          std::pair< Bound, Bound >
            rx = this->_m_bck.kernel()->approximate_absolute_1_object()(_m_diff[i].first.x(), prec),
            ry = this->_m_bck.kernel()->approximate_absolute_1_object()(_m_diff[i].first.y(), prec);

          CGAL::cpp0x::array<Bound,4> box =
            CGAL::make_array(rx.first, rx.second, ry.first, ry.second);

          std::pair< Bound, Bound > fym_intv = ieval_2(fym, box);
          if (CGAL::sign(fym_intv.first) == CGAL::sign(fym_intv.second) &&
              CGAL::sign(fym_intv.first) != CGAL::ZERO) {
            _m_diff_y_precs.push_back(prec);
            break;
          }

          // else
          prec *= 2;

          // std::cout << "prec: " << prec << std::endl;
        }
      }

    }
    // std::cout << "FFYisol 3" << std::endl;

    Bisolve_telemetry_code(t_ai_init.start();)

    bool has_vert_line = has_fiber && _m_bi_diff_vanish.has_vert_line_f(this->_m_bck.alpha());
    // std::cout << "FFYisol has_vert_line: " << has_vert_line << std::endl;

    // "active_intervals_set_at is "large enough", i.e., initial bound is well chosen
    _m_ais =
    _m_bi_diff_vanish.certifier_traits()->active_intervals_set_at
    (this->_m_bck.alpha(),
     _m_bi_diff_vanish.f().begin(), // TODO use PT_2 for range
     has_fiber || has_vert_line ? // TODO && !has_ver_line???
       _m_bi_diff_vanish.fiber(this->_m_bck.alpha()).f_end : _m_bi_diff_vanish.f().end()
    );

    Bisolve_telemetry_code(t_ai_init.stop();)

    Active_intervals_set_iterator ait = _m_ais.begin();

    size_t r = 0;
    int a = 0;
    r = 0;

    // std::cout << "FFYisol 4" << std::endl;

    Bisolve_telemetry_code(t_ffy_refine.start();)
    while (ait != _m_ais.end()) { // _m_diff not at end and _m_ai not at end

#if CGAL_BISOLVE_ENABLE_ARCAVOID
      if (!ait->touch_real()) {
        ait++;
        a++; // increase id, too
        continue;
      }
#endif

      bool overlap_with_diff = false;
      bool mult1_found = false;

      // refine ait and _m_diff[r] until latter is contained in former or they are disjoint
      // if disjoint, refine ait until mult = 1
      while (true) {

        // std::cout << "FFYisol L1" << std::endl;

        Bound al = _m_bi_diff_vanish.certifier_traits()->lower(_m_ais, ait);
        Bound ah = _m_bi_diff_vanish.certifier_traits()->upper(_m_ais, ait);

        // std::cout << "al: " << (al) << std::endl;
        // std::cout << "ah: " << (ah) << std::endl;
        // std::cout << "ald: " << CGAL::to_double(al) << std::endl;
        // std::cout << "ahd: " << CGAL::to_double(ah) << std::endl;

        Multiplicity_type multa = _m_bi_diff_vanish.certifier_traits()->upper_bound_on_multiplicity(_m_ais, ait);
        // TODO ??? Multiplicity_type multa = (_m_ais.max_var(ait) == _m_ais.min_var(ait) ? _m_ais.max_var(ait) : 0);
        // std::cout << "multa: " << multa << std::endl;

        Multiplicity_type multr = -1;

        // TODO FUTURE sum of mults of aits overlapping with rit == multr -> stop

        bool no_overlap = true;
        if (has_fiber && r < _m_diff.size()) {

          multr = _m_diff[r].second;
          // std::cout << "multr: " << multr << std::endl;

          std::pair< Bound, Bound >
            ry = this->_m_bck.kernel()->approximate_absolute_1_object()(_m_diff[r].first.y(),
                                                                        _m_diff_y_precs[r]);
          Bound rl = ry.first;
          Bound rh = ry.second;

          // std::cout << "rl: " << (rl) << std::endl;
          // std::cout << "rh: " << (rh) << std::endl;
          // std::cout << "rld: " << CGAL::to_double(rl) << std::endl;
          // std::cout << "rhd: " << CGAL::to_double(rh) << std::endl;

          if (al <= rl && rh <= ah) {
            if (multa <= multr + 1) {
//                  std::cout << "######## FFY AI includes Diff's interval ";
                overlap_with_diff = true;
                break;
            }
          } else if (rl <= al && ah <= rh) {
              if (multa <= multr + 1) {
//                  std::cout << "######## FFY Diff's includes AI\n";
                overlap_with_diff = true;
                break;
              }
          }
          no_overlap = (ah <= rl);
        }

        //! assuming that intervals are sorted, the only possibility that
        //! \c I=[al;ah] does not overlap with any diff's interval [rl;rh]
        //! is that \c I lies to the left of diff's intervals
        if (no_overlap && multa == 1) {
            mult1_found = true;
            break;
        }

        // std::cout << "FFYisol L6 (subdiv)" << std::endl;

        Bisolve_telemetry_code(t_ffy_ais_subdivide.start();)
        Active_intervals_set_iterator cb, ce;
        std::pair< Active_intervals_set_iterator, Active_intervals_set_iterator > cbe =
          _m_bi_diff_vanish.certifier_traits()->refine(_m_ais, ait);
        cb = cbe.first;
        ce = cbe.second;
        Bisolve_telemetry_code(t_ffy_ais_subdivide.stop();)
        // std::cout << "FFYisol L7 dist: " << std::distance(cb, ce) << std::endl;

        ait = cb;
        //break;
        if (cb == ce) {
          // std::cout << "FFYisol L8 leaf interval EXIT" << std::endl;
          break;
        }
      }

      // std::cout << "FFYisol L9" << std::endl;

      // not both conditions at the same time
      CGAL_assertion(!(overlap_with_diff && mult1_found));

      if (overlap_with_diff) {
        // std::cout << "TRUE " << r << std::endl;
        _m_container_and_ids.push_back(std::make_pair(true, r));
        r++;
        ait++;
        a++;
      }
      if (mult1_found) {
        // std::cout << "FALSE " << a << std::endl;
        _m_container_and_ids.push_back(std::make_pair(false, a));
        ait++;
        a++;
      }
    }
    Bisolve_telemetry_code(t_ffy_refine.stop();)
    // std::cout << "FFYisol L12" << std::endl;

#if 0 // TODO not sure whether this past-processing is needed
    // increase precision of approximations until they are disjoint
    for (size_t i = 0; i < number_of_real_roots(); i++) {

      bool sufficientleft = false;
      bool sufficientright = false;

      //std::cout << "\n=== i : " << i << std::endl;

      do {

        if (!sufficientleft && (i > 0)) {
          //std::cout << "mr " << CGAL::to_double(right_bound(i-1)) << std::endl;
          //std::cout << "il " << CGAL::to_double(left_bound(i)) << std::endl;
          sufficientleft = (CGAL::compare(right_bound(i-1), left_bound(i)) == CGAL::SMALLER);
        }
        if (!sufficientright && (i + 1 < number_of_real_roots())) {
          //std::cout << "ir " << CGAL::to_double(right_bound(i)) << std::endl;
          //std::cout << "pl " << CGAL::to_double(left_bound(i+1)) << std::endl;
          sufficientright = (CGAL::compare(right_bound(i), left_bound(i+1)) == CGAL::SMALLER);
        }

        std::cout << std::endl;

        if (!sufficientleft && (i > 0)) {
          refine_interval(i-1);
        }
        if (!sufficientleft || !sufficientright) {
          refine_interval(i);
        }
        if (!sufficientright && (i + 1 < number_of_real_roots())) {
          refine_interval(i+1);
        }

      } while (!sufficientleft && !sufficientright);

    }
#endif
    Base::_m_is_isolated = true;

  }

  int number_of_real_roots() const {
    return _m_container_and_ids.size();
  }

  Bound left_bound(size_t i) const {
    //std::cout << "lbi: "<< i << std::endl;
    //std::cout << "lbs: "<< _m_container_and_ids.size() << std::endl;
    CGAL_assertion(i < _m_container_and_ids.size());
    const std::pair< bool, size_t >& cai = _m_container_and_ids[i];
    if (cai.first) {
      CGAL_assertion(cai.second < _m_diff.size());
      return this->_m_bck.kernel()->approximate_absolute_1_object()(_m_diff[cai.second].first.y(),
                                                                    _m_diff_y_precs[cai.second]).first;
    }
    // else
    Active_intervals_set_const_iterator ait = _m_ais.begin();
    //std::cout << "lbe: "<< cai.second << std::endl;
    //std::cout << "lbt: "<< std::distance(_m_ais.begin(), _m_ais.end()) << std::endl;
    std::advance(ait, cai.second);
    Bound lower = _m_bi_diff_vanish.certifier_traits()->lower(_m_ais, ait);
    return lower;
  }

  Bound right_bound(size_t i) const {
    CGAL_assertion(i < _m_container_and_ids.size());
    const std::pair< bool, size_t > cai = _m_container_and_ids[i];
    if (cai.first) {
      CGAL_assertion(cai.second < _m_diff.size());
      return this->_m_bck.kernel()->approximate_absolute_1_object()(_m_diff[cai.second].first.y(),
                                                                    _m_diff_y_precs[cai.second]).second;
    }
    // else
    Active_intervals_set_const_iterator ait = _m_ais.begin();
    std::advance(ait, cai.second);
    Bound upper = _m_bi_diff_vanish.certifier_traits()->upper(_m_ais, ait);
    return upper;
  }

  //! returns whether i-th root is exact
  bool is_exact_root(size_t i) const {
    const std::pair< bool, size_t > cai = _m_container_and_ids[i];
    return cai.first;
  };

  //! \brief Returns whether the \c i th root is definitely a simple root of the isolated polynomial
  bool is_certainly_simple_root(size_t i) const {
    return (multiplicity_of_root(i) == 1);
  }

  //!\brief Returns whether the \c i th root is definitely a multiple root of the isolated polynomial
  bool is_certainly_multiple_root(size_t i) const {
    return (multiplicity_of_root(i) > 1);
  }

  int upper_bound_for_multiplicity(size_t i) const {
    return multiplicity_of_root(i);
  }

  //! returns the \c i th root multiplicity
  int multiplicity_of_root(size_t i) const {
    CGAL_assertion(i < _m_container_and_ids.size());
    const std::pair< bool, size_t > cai = _m_container_and_ids[i];
    if (cai.first) {
      CGAL_assertion(cai.second < _m_diff.size());
      return _m_diff[cai.second].second + 1;
    }
    // else
    return 1;
  }

  void refine_interval(size_t i) {
    CGAL_assertion(i < _m_container_and_ids.size());
    const std::pair< bool, size_t > cai = _m_container_and_ids[i];
    if (cai.first) {
      CGAL_assertion(cai.second < _m_diff.size());
      _m_diff_y_precs[cai.second] *= 2; // TODO is doubling a good strategy
    } else {
      Active_intervals_set_iterator ait = _m_ais.begin();
      std::advance(ait, cai.second);

      _m_bi_diff_vanish.certifier_traits()->refine(_m_ais, ait);
    }
  }

  //!@} // Access functions

private:

  //! Needed for reference counting
  virtual CGAL::Reference_counted_hierarchy<>* clone() {
    return new Geotop_lifter_bisolve_rep(*this);
  }

protected:

  //!\name Protected members
  //!@{

  //! pointer to kernel
  Bitstream_coefficient_kernel _m_bck; // supposed to be a handle-rep class

  //! local instance of bi diff vanish
  Bi_diff_vanish_2 _m_bi_diff_vanish;

  //! multiplicity of root of res
  Multiplicity_type _m_mult;

  //! local computation
  bool _m_local;

  //! maintenaince of active intervals
  Active_intervals_set _m_ais;

  //! storage for derivatives
  std::vector< std::pair < Algebraic_real_2, Multiplicity_type > > _m_diff;

  //! storage for current precisions of y-coordinates
  std::vector< long > _m_diff_y_precs;

  //! read from one or another container (multiple roots from _m_bi_diff_vanish and ordinary roots from _m_ais)
  std::vector< std::pair< bool, size_t > > _m_container_and_ids;

  //!@}

 };


#if CGAL_ACK_CURVE_ANALYSES_BISOLVE_USE_TEISSIER

  // TODO document Geotop_lifter_teissier_rep
template< class BitstreamCoefficientKernel, typename HandlePolicy > // no default on policy, should be decided in higher level
class Geotop_lifter_teissier_rep :
  // Comment: Multiple inheritance!
  public Generic_isolator_rep_base< typename BitstreamCoefficientKernel::Bound,
                                    HandlePolicy,
                                    CGAL::Tag_true >,
  public Arcavoid_list< BitstreamCoefficientKernel > {

  //!\name Tags
  //!@{

public:
  //! "refine_interval" function is implemented
  typedef CGAL::Tag_true Refine_interval_provided_tag;

  //!@} // Tags

public:

  //!\name Public typedefs
  //!@{

  //! this instance's template parameter
  typedef BitstreamCoefficientKernel Bitstream_coefficient_kernel;

  //! the policy class
  typedef HandlePolicy Handle_policy;

  //! type of Bound
  typedef typename Bitstream_coefficient_kernel::Bound Bound;

  //! the first base class
  typedef Generic_isolator_rep_base< Bound, Handle_policy, Refine_interval_provided_tag > Base;

  //! the class itself
  typedef Geotop_lifter_teissier_rep< Bitstream_coefficient_kernel, Handle_policy > Self;

  //! the second base class
  typedef CGAL::internal::Arcavoid_list< Bitstream_coefficient_kernel > Arcavoid_list;

  //! multiplicity type
  typedef typename Bitstream_coefficient_kernel::Algebraic_kernel_d_1::Multiplicity_type Multiplicity_type;

  //! type of Cluster
  typedef typename Arcavoid_list::Cluster Cluster;

  //! type of const iterator
  typedef typename Arcavoid_list::Cluster_iterator Cluster_iterator;

  //! type of const iterator
  typedef typename Arcavoid_list::Cluster_const_iterator Cluster_const_iterator;

  //!@} // Public typedefs

  //!\name Constructors
  //!@{

#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_DEFAULT_CONSTRUCTIBLE
 protected:
   Geotop_lifter_teissier_rep(); // = disable
#else
 public:
  //! Default constructor
   Geotop_lifter_teissier_rep() : Base(new Rep()) {}
#endif

  // needs no special assignable-implementation as no pointers
  // and this is a Rep of a handle class
#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_ASSIGNABLE
  protected:
    //Geotop_lifter_teissier_rep(const Self&); // = default
    Self& operator=(const Self&); // = disable
#endif


 public:
  //! standard constructor
  template< class CoefficientInputIterator >
  Geotop_lifter_teissier_rep(const Bitstream_coefficient_kernel& bck,
                             CoefficientInputIterator first, CoefficientInputIterator beyond,
                             Multiplicity_type mult,
                             Multiplicity_type mult_r_fx_fy,
                             bool no_stop) :
    Arcavoid_list(bck, first, beyond),
    _m_mult(mult),
    _m_mult_r_fx_fy(mult_r_fx_fy),
    _m_no_stop(no_stop) {

#if CGAL_ACK_DEBUG_FLAG
      CGAL_ACK_DEBUG_PRINT << "Use Teissier-isolator for status line at " << CGAL::to_double(this->kernel().alpha()) << std::endl;
#endif
    }

  //!@} // Constructors

public:

  //!\name Access functions
  //!@{

  void isolate() {

    const int d = CGAL::degree(static_cast<Arcavoid_list*>(this)->polynomial());

    if (d < 1) {
      _m_real_end = this->begin();
      Base::_m_is_isolated = true;
      return;
    }

    // if not full degree
    if (this->kernel().kernel()->sign_at_1_object()(CGAL::leading_coefficient(static_cast<Arcavoid_list*>(this)->polynomial()),
                                                    this->kernel().alpha()) == CGAL::ZERO) {
#if CGAL_ACK_DEBUG_FLAG
      CGAL_ACK_DEBUG_PRINT << "failed as polynomial as not full degree" << std::endl;
#endif
      // call backup
      Base::_m_is_isolated = false;
      return;
    }

    // if res(f_x, f_y) vanishes
    if (_m_mult_r_fx_fy < 0) {

#if CGAL_ACK_DEBUG_FLAG
      CGAL_ACK_DEBUG_PRINT << "failed as res_fx_fy vanishes" << std::endl;
#endif
      // call backup
      Base::_m_is_isolated = false;
      return;
    }

    int num_clusters = std::distance(Arcavoid_list::begin(), Arcavoid_list::end());

    // EBEB for debugging
    //typename Cluster::Cluster_multiplicity_add cmadd;
    //int mult_sum;

    int i = 0;

#ifdef CGAL_XTri_USE_TIMERS
    tm_external_timer.start();
#endif
    do {

      typename Arcavoid_list::Cluster_iterator ait = Arcavoid_list::begin();
      while (ait != Arcavoid_list::end()) {
        if (ait->multiplicity() > 1) {
          ait = Arcavoid_list::subdivide_cluster(ait).second;
        } else {
          ait++;
        }
      }

      num_clusters = std::distance(Arcavoid_list::begin(), Arcavoid_list::end());
      // for debugging
      // mult_sum = std::accumulate(Arcavoid_list::begin(), Arcavoid_list::end(), 0, cmadd);

      // We ensure the conjugate-circle property by 2nr safety area of clusters

      const int stop = 8;
      if (!_m_no_stop && i == stop) {
        // after "stop" iterations of num-solve let Bisolve do the job

#if CGAL_ACK_DEBUG_FLAG
      CGAL_ACK_DEBUG_PRINT << "failed after " << stop << " iterations" << std::endl;
#endif
        Base::_m_is_isolated = false;
#ifdef CGAL_XTri_USE_TIMERS
    tm_external_timer.stop();
#endif
        return;
      }

#if CGAL_ACK_DEBUG_FLAG
      CGAL_ACK_DEBUG_PRINT << "d : " << d << std::endl;
      CGAL_ACK_DEBUG_PRINT << "mR: " << _m_mult << std::endl;
      CGAL_ACK_DEBUG_PRINT << "mS: " << _m_mult_r_fx_fy << std::endl;
      CGAL_ACK_DEBUG_PRINT << "nc: " << num_clusters << std::endl;
#endif
      // EBEB for debugging
      //std::cout << "ms: " << mult_sum << std::endl;

      i++;

    } while (d - _m_mult + _m_mult_r_fx_fy != num_clusters);
#ifdef CGAL_XTri_USE_TIMERS
    tm_external_timer.stop();
#endif

    // obtain real intervals and store past-the-end iterator
    _m_real_end = Arcavoid_list::partition_and_sort_reals();

    // std::cout << "DISTr: " << std::distance(this->begin(), this->end()) << std::endl;
    // std::cout << "DISTc: " << std::distance(Arcavoid_list::begin(), Arcavoid_list::end()) << std::endl;

    // done
    Base::_m_is_isolated = true;
    return;
  }

  //! returns the number of detected isolating intervals
  int number_of_real_roots() const {
    CGAL_precondition(_m_real_end);
    return std::distance(this->begin(), this->end());
  }

  //! the lower bound of the \c i th root
  Bound left_bound(size_t i) const {
    CGAL_precondition(_m_real_end);
    Cluster_const_iterator it = this->begin();
    std::advance(it, i);
    return it->center().real() - it->radius();
  }

  //! the upper bound of the \c i th root
  Bound right_bound(size_t i) const {
    CGAL_precondition(_m_real_end);
    Cluster_const_iterator it = this->begin();
    std::advance(it, i);
    return it->center().real() + it->radius();
  }

  //! returns whether i-th root is exact
  bool is_exact_root(size_t /* i */) const {
    // radius == 0 is not possible (as it gets a bitstream coefficient kernel)
    return false;
  };

  //! \brief Returns whether the \c i th root is definitely a simple root of the isolated polynomial
  bool is_certainly_simple_root(size_t i) const {
    return (multiplicity_of_root(i) == 1);
  }

  //!\brief Returns whether the \c i th root is definitely a multiple root of the isolated polynomial
  bool is_certainly_multiple_root(size_t i) const {
    return (multiplicity_of_root(i) > 1);
  }

  //! returns an upper bound of the \c i th root multiplicity
  int upper_bound_for_multiplicity(size_t i) const {
    CGAL_precondition(_m_real_end);
    Cluster_const_iterator it = this->begin();
    std::advance(it, i);
    return it->multiplicity();
  }

  //! returns the \c i th root multiplicity
  int multiplicity_of_root(size_t i) const {
    CGAL_precondition(_m_real_end);
    Cluster_const_iterator it = this->begin();
    std::advance(it, i);
    return it->multiplicity();
  }

  //! Refine the <tt>i</tt>th isolating interval
  void refine_interval(size_t i) {
    CGAL_precondition(_m_real_end);
    Cluster_iterator it = Arcavoid_list::begin();
    std::advance(it, i);
    Arcavoid_list::subdivide_cluster(it);
  }

  //!@} // Access functions

protected:

  Cluster_const_iterator begin() const {
    return Arcavoid_list::begin();
  }

  // DO NOT DELETE. We have to return this, as this isolator is only referring to REAL roots
  Cluster_const_iterator end() const {
    CGAL_precondition(_m_real_end);
    return *_m_real_end;
  }

  //! Needed for reference counting
  virtual CGAL::Reference_counted_hierarchy<>* clone() {
    return new Geotop_lifter_teissier_rep(*this);
  }

protected:

  //!\name Protected members
  //!@{

  //! multiplicity of x as root of res
  Multiplicity_type _m_mult;

  //! multiplicity of x as root of res_fx_fy
  Multiplicity_type _m_mult_r_fx_fy;

  //! keep iterations running and do not stop "early"
  bool _m_no_stop;

  //! real end in "Arcavoid_list"
  mutable boost::optional< Cluster_const_iterator > _m_real_end;

  //!@} // Protected members

};

#endif // CGAL_ACK_CURVE_ANALYSES_BISOLVE_USE_TEISSIER

/*\brief Lifts critical and non-critical x-coordinates in a curve analysis.
 *
 * Lifting means to compute the real roots (with multiplicities) for a
 * polynomial f(x0,y) with x0 possibly non-rational.
 * One out of four options are possible:
 * 1) lifting over a rational x0 -> RS
 * 2) lifting over a non-rational x0 (given by bck-instance),
 *    but f(x0,y) stays square-free -> Bitstream_descartes
 * 3) lifting over a non-rational x0 (given by bck-instance),
 *    but f(x0,y) is not square-free -> iterated Bisolve
 * 4) lifting real and complex roots -> Teissier
 */
template< class BitstreamCoefficientKernel, typename HandlePolicy = CGAL::Handle_policy_no_union >
class Geotop_lifter :
  public internal::Generic_isolator< typename BitstreamCoefficientKernel::Polynomial,
                                     typename BitstreamCoefficientKernel::Bound,
                                     HandlePolicy,
#if CGAL_USE_RS3
                                     RS_ak_1_isolator_rep< typename BitstreamCoefficientKernel::Polynomial::NT, typename BitstreamCoefficientKernel::Bound >
#else
                                     Square_free_bitstream_descartes_isolator_rep< BitstreamCoefficientKernel, HandlePolicy >
#endif
> {

public:

  //!\name Public typedefs
  //!@{

  //! this instance's first template parameter
  typedef BitstreamCoefficientKernel Bitstream_coefficient_kernel;

  //! second template paramter
  typedef HandlePolicy Handle_policy;

  //! The polynomial type
  typedef typename Bitstream_coefficient_kernel::Polynomial Polynomial;

  //! How the boundaries of the isolating intervals are represented
  typedef typename Bitstream_coefficient_kernel::Bound Bound;

#if CGAL_USE_RS3
  //! the base type
  typedef internal::Generic_isolator< Polynomial, Bound, Handle_policy, RS_ak_1_isolator_rep< typename Polynomial::NT, Bound > /* in fact, provide all of the used reps! */ > Base;
#else
typedef internal::Generic_isolator< Polynomial, Bound, Handle_policy, Square_free_bitstream_descartes_isolator_rep< Bitstream_coefficient_kernel, Handle_policy >  /* TODO in fact, provide all of the used reps! */ > Base;
#endif

  //! the class itself
  typedef Geotop_lifter< Bitstream_coefficient_kernel, Handle_policy > Self;

  // the next two are needed for two constructors

  //! multiplicity type
  typedef typename Bitstream_coefficient_kernel::Algebraic_kernel_d_1::Multiplicity_type Multiplicity_type;

  //! Bi_diff_vanish
  typedef typename Geotop_lifter_bisolve_rep< Bitstream_coefficient_kernel, Handle_policy >::Bi_diff_vanish_2 Bi_diff_vanish_2;

  //!@} // Public typedefs

  private:

  // TODO move to a helper file?
  struct Substitute_x {

    template< typename Polynomial_2, typename Bound >
    typename Polynomial_2::NT /* TODO use 'auto' in order to avoid fraction traits problems */
    operator()(const Polynomial_2& p, const Bound& b) {

      typedef CGAL::Fraction_traits< Bound > FT;
      typename FT::Decompose decompose;

      typename FT::Numerator_type num;
      typename FT::Denominator_type denom;
      decompose(b, num, denom);

      typename CGAL::Polynomial_traits_d<Polynomial_2>::
        Evaluate_homogeneous evh;
      typename CGAL::Polynomial_traits_d<Polynomial_2>::
        Move move;

      return evh(move(p,0,1),num,denom);
    }

  };

    public:
  //!\name Constructors
  //!@{

#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_DEFAULT_CONSTRUCTIBLE
 protected:
  // not recommended to have default constructor
  Geotop_lifter(); // = disable
#else
 public:
  //! Default constructor
  Geotop_lifter() : Base(new Rep()) {}
#endif

#if 0 // TODO 2012 assignable/shared_ptr
#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_ASSIGNABLE
  protected:
    Geotop_lifter(const Self&); // = disable
    Geotop_lifter& operator=(const Self&); // = disable
#endif
#endif

  public:

  // TODO 2012 replace all constructors by one ctor for event and one ctor for internal

  //! bound constructor (uses RS directly)
  template< class CoefficientInputIterator >
  Geotop_lifter(const Bitstream_coefficient_kernel& bck,
                CoefficientInputIterator first, CoefficientInputIterator beyond,
                const Bound& b /* b is kind of redundant here, as stored as bck.alpha().low() */) :
#if CGAL_USE_RS3
    Base(new RS_ak_1_isolator_rep< typename Polynomial::NT, Bound >((Substitute_x()(Polynomial(first, beyond),b)), false))
#else
Base(new Square_free_bitstream_descartes_isolator_rep< Bitstream_coefficient_kernel, Handle_policy >(Polynomial(first, beyond), bck)) // note: no subs takes place here!
#endif
  {
    CGAL_precondition(bck.alpha().low() == b);
  }

  //! bitstream constructor (gets a bivariate bck with x = alpha stored)
  Geotop_lifter(Polynomial f,
                const Bitstream_coefficient_kernel& bck) :
    // use square-free bitstream descartes
    Base(new Square_free_bitstream_descartes_isolator_rep< Bitstream_coefficient_kernel, Handle_policy >(f, bck))
  {
  }

  //! bisolve constructor
  Geotop_lifter(const Bi_diff_vanish_2& bdv,
                const Bitstream_coefficient_kernel& bck,
                Multiplicity_type mult,
                bool local = false) :
    Base(new Geotop_lifter_bisolve_rep< Bitstream_coefficient_kernel, Handle_policy >(bck, mult, bdv, local))
  {
  }

#if CGAL_ACK_CURVE_ANALYSES_BISOLVE_USE_TEISSIER
  //! arcavoid constructor
  template< class CoefficientInputIterator >
  Geotop_lifter(CoefficientInputIterator first, CoefficientInputIterator beyond,
                const Bitstream_coefficient_kernel& bck,
                Multiplicity_type mult,
                Multiplicity_type mult_r_fx_fy,
                bool no_stop = false) :
    Base(new Geotop_lifter_teissier_rep< Bitstream_coefficient_kernel, Handle_policy >(bck, first, beyond, mult, mult_r_fx_fy, no_stop))
  {
  }
#endif // CGAL_ACK_CURVE_ANALYSES_BISOLVE_USE_TEISSIER

  //!@} // Constructors

  //!\name Destructor
  //!@{

  virtual ~Geotop_lifter() {
  }

  //!@} // Destructor

  // all member functions are derived from Generic_isolator

};

} // namespace internal

} // namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_D_CURVE_ANALYSIS_2_GEOTOP_LIFTER_H
// EOF
