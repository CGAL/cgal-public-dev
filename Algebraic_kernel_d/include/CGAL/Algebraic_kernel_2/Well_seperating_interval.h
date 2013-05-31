// Copyright (c) 2010, 2012 Max-Planck-Institut fuer Informatik (Germany).
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
// $URL: svn+ssh://asm@scm.gforge.inria.fr/svn/cgal/branches/experimental-packages/Algebraic_kernel_2/include/CGAL/Bi_slice_certify.h $
// $Id: Bi_slice_certify.h 57108 2010-06-25 13:54:53Z eric $
//
//
// Author(s): Eric Berberich <eric.berberich@cgal.org>


#ifndef CGAL_ALGEBRAIC_KERNEL_2_WELL_SEPERATING_INTERVAL_H
#define CGAL_ALGEBRAIC_KERNEL_2_WELL_SEPERATING_INTERVAL_H

/*! \file Well_seperating_interval.h
 * \brief Compute a well seperating interval for a root.
 */


#include <CGAL/config.h>
#include <CGAL/Algebraic_kernel_2/Bi_solve_2_flags.h>

#include <CGAL/Algebraic_kernel_2/Root_box_1.h>

namespace CGAL {

namespace internal {

template < class IsolatingApproximation, class Refiner_ >
struct Well_seperating_interval {
  
  //! instance's first template parameter
  typedef IsolatingApproximation Isolating_approximation;
  
  //! instance's second template parameter
  typedef Refiner_ Refiner;
  
  //! bound type
  typedef typename Isolating_approximation::Bound Bound;

  typedef Root_box_1< Bound > Root_box;

  /*! 
   * t_test to decide whether other roots are far enough away
   * \param p polynomial (already scaled and shifted)
   * \param r radius of current interval
   * \param c scale factor
   * \return true if t test returns true
   */
  template < class Polynomial_1_ > 
  bool _t_test (const Polynomial_1_& p,
                const Bound& r, const Bound& c) const {
    Bisolve_telemetry_code(t_ttest.start();)
    Bound y(0);
    for (int i = 1; i <= CGAL::degree(p); i++) {
      Bound t =  CGAL::abs(p[i]) * CGAL::ipower(r, i);
      y += t;
    }
    bool res = (CGAL::abs(p[0]) > c * y);
    Bisolve_telemetry_code(t_ttest.stop();)
    return res;
  }
  
  Well_seperating_interval() {
  }
 
public:
  
  template< class InputIterator >
  void operator()(Isolating_approximation& iapp, 
                  const Bound& factor,
                  InputIterator begin, InputIterator end, InputIterator diff) {
    
    typedef typename InputIterator::value_type::first_type Polynomial_1; // TODO not first_type
    
    typedef typename CGAL::Polynomial_traits_d< Polynomial_1 >::Coefficient_type Coefficient;

    Root_box& rb = iapp.rb();

    Refiner refine;
    
    // execute lense filter only if there is a single factor with non-trivial multiplicity
    bool filter = ((std::distance(begin,end) == 1) && (diff->second > 1));
    bool result = true;
    bool t3over2 = false;

    do {

      std::pair< Bound, Bound > intv;
  
      result = true;
      
      // refine
      Bisolve_telemetry_code(t_refine_app.start();)

      if (!rb.is_exact()) {

        intv = refine(iapp);
        
        // mid
        rb.set_mid((intv.first + intv.second) / Bound(2));
	
        // rad
        rb.set_rad((intv.second - intv.first) / Bound(2)); 
        
        if (CGAL::sign(rb.rad()) == CGAL::ZERO) {
          rb.set_rad(Bound(1));
          rb.set_exact();
        }
	
      } else {
        rb.set_rad(rb.rad() / Bound(16)); // is 16 a good value
      }

      Bisolve_telemetry_code(t_refine_app.stop();)
      
      if (!t3over2) {
	
	Bound deg(CGAL::degree(diff->first));

	if (!rb.is_exact() && 
            diff->second == 1 && !CGAL_BISOLVE_DISABLE_COMBINATORIAL_CERTIFICATION) {
          //std::cout << "diffsecond == 1" << std::endl;
	  // do nothing

	} else {

	  Bisolve_telemetry_code(c_t3over2_total++;)
	  
	  // lense filter

#if !CGAL_BISOLVE_ENABLE_LENSE_FILTER_IN_T_TEST
          Bisolve_telemetry_code(t_filter.start();)

	  if (filter && !rb.is_exact()) {
	    
	    Bound a0 = intv.first, b0 = intv.second, w0 = b0 - a0;
	    
	    Bound a(a0);
	    Bound b(b0);
	    
	    while (true) {

	      intv = refine(iapp);
	      
	      a = intv.first;
	      b = intv.second;
	      Bound w = b-a;
	      Bound m = (b+a) / Bound(2);
	      
	      if (CGAL::sign(w) == CGAL::ZERO) {
	      	rb.set_mid(m);
                rb.set_rad(Bound(1));
    		rb.set_exact();
		
                break;
	      }

	      Bound delta = CGAL::min(CGAL::abs(a0 - a), CGAL::abs(b0 - b));
	      
	      if (CGAL::compare(w, delta/(6*deg+3)) == CGAL::SMALLER) {
		    // interval is safe
		
		    // mid
		    rb.set_mid(m);
                
		    // rad
		    rb.set_rad(w/Bound(2));
		
		    Bisolve_telemetry_code(c_lensefilter_pass++;)
                //std::cout << "FILTER PASS" << std::endl;
		    t3over2 = true;
		    break;
	      } 

	      // else 
	      
	      if (CGAL::compare(w, w0/(6*deg*(deg+3))) == CGAL::SMALLER) {
		    // do t-test
		    t3over2 = false;
		    filter = false; // no more filter
		    break;
	      }
	    }
	  }
	  Bisolve_telemetry_code(t_filter.stop();)
#endif // !CGAL_BISOLVE_ENABLE_LENSE_FILTER_IN_T_TEST
	  
	  if (!t3over2) {
	    // lense filter failed, now taylor shift
	    Bisolve_telemetry_code(t_tshift.start();)

#if !CGAL_BISOLVE_USE_RS_AK
            typedef CGAL::Fraction_traits< Bound > FT;
            typename FT::Numerator_type mn;
            typename FT::Denominator_type md;
            typename FT::Decompose decompose;
            decompose(rb.mid(), mn, md);
#else
            typedef CGAL::internal::Float_traits< Bound > FT;
            typename FT::Get_mantissa mantissa;
            typename FT::Get_exponent exponent;
            typedef CGAL::Coercion_traits< typename CGAL::Get_arithmetic_kernel< Bound >::Arithmetic_kernel::Integer, Coefficient > CTm;
            typename CTm::Cast castm;
            typedef CGAL::Coercion_traits< long, Coefficient > CTe;
            typename CTm::Type mn(castm(mantissa(rb.mid())));
            typename CTe::Type md(1);
            long e = exponent(rb.mid());
            if (e < 0) {
              md <<= -e;
            } else {
              mn <<= e;
            }
#endif
            
	    Polynomial_1 pm = 
	      CGAL::translate_homogeneous(CGAL::differentiate(diff->first), mn, md);
	    Bisolve_telemetry_code(t_tshift.stop();)
            bool test = _t_test(pm, factor * rb.rad(), Bound(3)/Bound(2)); 
            //std::cout << "T_3/2 predicate, multf=" << diff->second << ", res=" << test << std::endl;
	    t3over2 = test;
	    result &= t3over2; 
	    if (t3over2) {
	      rb.set_used_t_test();
	    }
	  }
	}
      }

#if !CGAL_BISOLVE_USE_RS_AK
      typedef CGAL::Fraction_traits< Bound > FT;
      typename FT::Numerator_type mn;
      typename FT::Denominator_type md;
      typename FT::Decompose decompose;
      decompose(rb.mid(), mn, md);
#else
      typedef CGAL::internal::Float_traits< Bound > FT;
      typename FT::Get_mantissa mantissa;
      typename FT::Get_exponent exponent;
      typedef CGAL::Coercion_traits< typename CGAL::Get_arithmetic_kernel< Bound >::Arithmetic_kernel::Integer, Coefficient > CTm;
      typename CTm::Cast castm;
      typedef CGAL::Coercion_traits< long, Coefficient > CTe;
      typename CTm::Type mn(castm(mantissa(rb.mid())));
      typename CTe::Type md(1);
      long e = exponent(rb.mid());
      if (e < 0) {
        md <<= -e;
      } else {
        mn <<= e;
      }
#endif

      // t_1 for all other factors
      for (InputIterator it = begin; result && (it != end); it++) {
        if (it != diff /* not current factor */) {
          Bisolve_telemetry_code(t_tshift.start();)
          Polynomial_1 pm = 
            CGAL::translate_homogeneous(it->first, mn, md);
          Bisolve_telemetry_code(t_tshift.stop();)
          bool test = _t_test(pm, factor * rb.rad(), Bound(1)); 
          //std::cout << "T_1 predicate, multf=" << it->second << ", res=" << test << std::endl;
          result &= test;
        }
      }
      
      //std::cout << "result: " << result << std::endl;
      
    } while (!result);
    
    if (rb.used_t_test()) {
#if CGAL_BISOLVE_RELAXED_RATIOS
      rb.set_rad((Bound(2)) * rb.rad());
#else
      rb.set_rad((Bound(2)) * rb.rad());
#endif
    }
    
  }

};

} // namespace internal

} // namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_2_WELL_SEPERATING_INTERVAL_H
// EOF
