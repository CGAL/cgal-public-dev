// Copyright (c) 2009, 2010, 2011, 2012 Max-Planck-Institut fuer Informatik (Germany).
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
// Author(s): Eric Berberich <eric.berberich@cgal.org>
//            Pavel Emeliyanenko <asm@mpi-inf.mpg.de>

#ifndef CGAL_ALGEBRAIC_KERNEL_2_BI_SOLVE_2_H
#define CGAL_ALGEBRAIC_KERNEL_2_BI_SOLVE_2_H

/*! \file
 * The header file for the Bi_solve_2 class.
 */

// TODO dependency graph of actions to take ... search for right design pattern

#include <CGAL/config.h>

#ifndef CGAL_BISOLVE_ENABLE_NTL_FACTORIZE
// TODO 2012 symbolic_exports.h???
#include <CGAL/symbolic_exports.h>
#endif

#include <CGAL/Algebraic_kernel_2/Bi_solve_2_flags.h>

#include <vector>
#include <list>
#include <set>

#include <boost/optional.hpp>

#include <CGAL/Handle_with_policy.h>

#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>

#include <CGAL/Algebraic_kernel_d/Bitstream_coefficient_kernel.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_rndl_tree_traits.h>

#include <CGAL/Algebraic_kernel_2/Bi_slice_certify.h>
#include <CGAL/Algebraic_kernel_2/Well_seperating_interval.h>

namespace CGAL {

namespace internal {

//! template rep class for new bi-solve
template < class CertifierTraits >
class Bi_solve_2_rep {

public:
    
  // this instance's template parameter
  typedef CertifierTraits Certifier_traits;
  
  //! type of AK_1
  typedef typename Certifier_traits::Algebraic_kernel_d_1 Algebraic_kernel_d_1;

  //! type of Coefficient
  typedef typename Algebraic_kernel_d_1::Coefficient Coefficient;

  //! type for univariate polynomials
  typedef typename Algebraic_kernel_d_1::Polynomial_1 Polynomial_1;

  //! type for algebraic reals
  typedef typename Algebraic_kernel_d_1::Algebraic_real_1 Algebraic_real_1;

  //! type for Multiplicity
  typedef typename Algebraic_kernel_d_1::Multiplicity_type Multiplicity_type;

  //! type for Bound
  typedef typename Algebraic_kernel_d_1::Bound Bound;

  //! type of univariate polynomial traits
  typedef CGAL::Polynomial_traits_d< Polynomial_1 > PT_1;

  //! type of bivariate polynomial
  typedef typename CGAL::Polynomial_type_generator< Coefficient, 2>::Type Polynomial_2;

  //! type of bivariate polynomial traits
  typedef CGAL::Polynomial_traits_d< Polynomial_2 > PT_2;

  //! type of bivariate algebraic real
  typedef typename Certifier_traits::Algebraic_real_2 Algebraic_real_2;

  //! type of certifier
  typedef internal::Bi_slice_certify< Certifier_traits > Certifier;

  //! isolating box of a root
  typedef typename Certifier::Root_box Root_box;

  //! bigfloat interval 
  typedef typename Certifier::BFI BFI;

  //! bigfloats
  typedef typename Certifier::BigFloat BigFloat;

  //! solution
  typedef typename Certifier::Solution_1 Solution_1;

  //! solutions
  typedef typename Certifier::Solution_1_container Solution_1_container;

  //! solutions iterator
  typedef typename Solution_1_container::const_iterator Solution_1_container_const_iterator;

  //! approximator
  typedef typename Certifier::Approximate_solution_1 Approximate_solution_1;

#if !CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS
  //! iterator to active intervals
  typedef typename Certifier::Active_intervals_set Active_intervals_set;
#endif

  //! type of a box
  typedef typename Algebraic_real_2::Box_2 Box_2;
  
  typedef typename Certifier::Algebraic_real_2_x_less_than Algebraic_real_2_x_less_than;
  typedef typename Certifier::Algebraic_real_2_y_less_than Algebraic_real_2_y_less_than;
  typedef typename Certifier::Algebraic_real_2_xy_less_than Algebraic_real_2_xy_less_than;
 
  struct Solution_1_less_than {
    bool operator()(const Solution_1& s1, const Solution_1& s2) const {
      
      return (CGAL::compare(s1.value(), s2.value()) == CGAL::SMALLER);
    }
  };

public:
 
  typedef typename Certifier::template Fiber< Algebraic_real_2_y_less_than > Fiber_x;
  
  typedef std::map< Solution_1, Fiber_x, Solution_1_less_than > Fiber_map_x;
  
  typedef typename Fiber_map_x::iterator Fiber_map_x_iterator;
  typedef typename Fiber_map_x::const_iterator Fiber_map_x_const_iterator;

  typedef typename Certifier::template Fiber< Algebraic_real_2_x_less_than > Fiber_y;
  
  typedef std::map< Solution_1, Fiber_y, Solution_1_less_than > Fiber_map_y;
  
  typedef typename Fiber_map_y::iterator Fiber_map_y_iterator;
  typedef typename Fiber_map_y::const_iterator Fiber_map_y_const_iterator;

public:

  Bi_solve_2_rep() :
    _m_achieved(0),
    _m_second_is_diff(false)
  {}

  Bi_solve_2_rep(const Polynomial_2& f) :
    _m_achieved(0),
    _m_f(f), _m_g(CGAL::differentiate(f)),
    _m_second_is_diff(true)
  {}

  Bi_solve_2_rep(const Polynomial_2& f, const Polynomial_2& g) :
    _m_achieved(0),
    _m_f(f), _m_g(g),
    _m_second_is_diff(false)
  {}

public:

  Algebraic_kernel_d_1 _m_ak1;

  mutable int _m_achieved;
  
  //! copy of f
  Polynomial_2 _m_f;

  //! copy of g
  Polynomial_2 _m_g;

  //! indicates whether g = f_y
  bool _m_second_is_diff;

  //! transposed f
  mutable Polynomial_2 _m_ft;

  //! transposed g
  mutable Polynomial_2 _m_gt;

  //! resultant(f,g,x)
  mutable boost::optional< Polynomial_1 > _m_res_in_y;
  
  //! and its factors
  mutable std::vector< std::pair< Polynomial_1, Multiplicity_type > > _m_res_in_y_factors;

  //! resultant(f,g,y)
  mutable boost::optional< Polynomial_1 > _m_res_in_x;

  //! and its factors
  mutable std::vector< std::pair< Polynomial_1, Multiplicity_type > > _m_res_in_x_factors;

  //! certifier_traits
  mutable boost::optional< Certifier_traits > _m_certifier_traits;

  //! upper bound 
  mutable boost::optional< long > _m_grid_upper_bound_y_log2_abs;

  //! map of fibers for x-coordinates
  mutable Fiber_map_x _m_fibers_x;

  //! map of fibers for y-coordinates
  mutable Fiber_map_y _m_fibers_y;

};

} // namespace internal
 

//! template class for new bi-solve
 template < class CertifierTraits, 
            class Rep_ = CGAL::internal::Bi_solve_2_rep< CertifierTraits > >
  class Bi_solve_2 : public CGAL::Handle_with_policy< Rep_ > {

public:
  
  enum Bi_solve_2_status {
    CGAL_BISOLVE_CONSTRUCTED = 0,

    CGAL_BISOLVE_RESULTANT_IN_X = 1,
    CGAL_BISOLVE_RESULTANT_IN_Y = 2,

    CGAL_BISOLVE_SQUARE_FREE_FACTORIZATION_IN_X = 4,
    CGAL_BISOLVE_SQUARE_FREE_FACTORIZATION_IN_Y = 8,

    CGAL_BISOLVE_ISOLATION_IN_X = 16,
    CGAL_BISOLVE_ISOLATION_IN_Y = 32,

    CGAL_BISOLVE_SEPERATION_IN_X =  64,
    CGAL_BISOLVE_SEPERATION_IN_Y = 128,
 
    CGAL_BISOLVE_ALL_BUT_LIFTING = 255,

    CGAL_BISOLVE_LIFTING = 256
  };


  //! this instance's first template parameter
  typedef CertifierTraits Certifier_traits;

  //! this instance's second template parameer;
  typedef Rep_ Rep;
 
  //! base class
  typedef CGAL::Handle_with_policy< Rep > Base;

  //! type of AK_1
  typedef typename Rep::Algebraic_kernel_d_1 Algebraic_kernel_d_1;

  //! coefficient type
  typedef typename Algebraic_kernel_d_1::Coefficient Coefficient;

  //! type for univariate polynomials
  typedef typename Algebraic_kernel_d_1::Polynomial_1 Polynomial_1;

  //! type for algebraic reals
  typedef typename Algebraic_kernel_d_1::Algebraic_real_1 Algebraic_real_1;

  //! type for Multiplicity
  typedef typename Algebraic_kernel_d_1::Multiplicity_type Multiplicity_type;

  //! type for Bound
  typedef typename Algebraic_kernel_d_1::Bound Bound;

  //! type of univariate polynomial traits
  typedef typename Rep::PT_1 PT_1;

  //! type of bivariate polynomial
  typedef typename Rep::Polynomial_2 Polynomial_2;

  //! type of bivariate polynomial traits
  typedef typename Rep::PT_2 PT_2;

  //! type of bivariate algebraic real
  typedef typename Rep::Algebraic_real_2 Algebraic_real_2;

  //! type of certifier
  typedef typename Rep::Certifier Certifier;

  //! isolating box of a root
  typedef typename Certifier::Root_box Root_box;

  //! bigfloat interval 
  typedef typename Certifier::BFI BFI;

  //! bigfloats
  typedef typename Certifier::BigFloat BigFloat;

  //! solution
  typedef typename Certifier::Solution_1 Solution_1;

  //! less than
  typedef typename Rep::Solution_1_less_than Solution_1_less_than;

  //! solutions
  typedef typename Certifier::Solution_1_container Solution_1_container;

  //! solutions iterator
  typedef typename Solution_1_container::const_iterator Solution_1_container_const_iterator;
 
public:

  //! approximator
  typedef typename Certifier::Approximate_solution_1 Approximate_solution_1;

#if !CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS
  //! iterator to active intervals
  typedef typename Certifier::Active_intervals_set Active_intervals_set;
#endif

  //! type of a box
  typedef typename Algebraic_real_2::Box_2 Box_2;
  
  //! type of fiber
  typedef typename Rep::Fiber_x Fiber_x;

  typedef typename Rep::Fiber_map_x Fiber_map_x;
  
  typedef typename Fiber_map_x::iterator Fiber_map_x_iterator;
  typedef typename Fiber_map_x::const_iterator Fiber_map_x_const_iterator;

  typedef typename Rep::Fiber_y Fiber_y;

  typedef typename Rep::Fiber_map_y Fiber_map_y;
  
  typedef typename Fiber_map_y::iterator Fiber_map_y_iterator;
  typedef typename Fiber_map_y::const_iterator Fiber_map_y_const_iterator;


  struct Fiber_map_inserter {
    
  public:

    Fiber_map_inserter(Fiber_map_x& fmx, Fiber_map_y& fmy) :
      _m_fmx(fmx),
      _m_fmy(fmy) {
    }

    template < class InputIterator >
    void operator()(InputIterator begin, InputIterator end) {
      
      for (InputIterator it = begin; it != end; it++) {
        // TODO let Algebraic_real_2 store two Solution_1
        Fiber_map_x_iterator fxit = _m_fmx.find(Solution_1(it->x(), it->mult_x()));
        CGAL_assertion(fxit != _m_fmx.end());
        fxit->second.roots.insert(*it);
        // TODO let Algebraic_real_2 store two Solution_1
        Fiber_map_y_iterator fyit = _m_fmy.find(Solution_1(it->y(), it->mult_y()));
        CGAL_assertion(fyit != _m_fmy.end());
        fyit->second.roots.insert(*it);
      }
    }
    
  private:
    Fiber_map_x& _m_fmx;
    Fiber_map_y& _m_fmy;
    
  };

  struct Refiner {
    
    std::pair< Bound, Bound > operator()(const Solution_1& s) {
      Approximate_solution_1 approx_sol;
      if (s.minprec() < 2) {
        s.increase_minprec(2);
      } else {
        s.double_minprec();
      }
      std::pair< Bound, Bound > ab = approx_sol(s,s.minprec());
      return ab;
    }
    
  };

public:

  //!\name Constuctors
  //!@{

  //! function-call operator to solve bivariate system of two polynomials
  Bi_solve_2() :
    Base()  {
  }


  //! function-call operator to solve bivariate system of two polynomials
   Bi_solve_2(const Polynomial_2& f, const Polynomial_2& g,
              int achieve = CGAL_BISOLVE_CONSTRUCTED) :
    Base(f,g) {
    _achieve(achieve);
  }

  //! function-call operator to solve bivariate system
  Bi_solve_2(const Polynomial_2& f,
            int achieve = CGAL_BISOLVE_CONSTRUCTED) :
    Base(f) {

    _achieve(achieve);
  }
  
  //!@}

  //!\name Access 
  //!@{

  //! the algebraic kernel used
  Algebraic_kernel_d_1* kernel() {
    return &(this->ptr()->_m_ak1);
  }
  
  //! the first polynomial
  const Polynomial_2 f() const {
    return this->ptr()->_m_f;
  }

  //! the second polynomial
  const Polynomial_2 g() const {
    return this->ptr()->_m_g;
  }

  //! the first polynomial
  const Polynomial_2 ft() const {
    return this->ptr()->_m_ft;
  }

  //! the second polynomial
  const Polynomial_2 gt() const {
    return this->ptr()->_m_gt;
  }

  //!@}

public:

  // for experts only
  void _achieve(int achieve) {
    
    if (achieve == CGAL_BISOLVE_CONSTRUCTED) {
      return;
    }

    // PHASE 1
    _resultants(achieve);
    _square_free_factorization(achieve);
    _isolation(achieve);

    // PHASE 2
    this->_seperation(achieve, 
                      this->ptr()->_m_fibers_x.begin(), this->ptr()->_m_fibers_x.end(),
                      this->ptr()->_m_fibers_y.begin(), this->ptr()->_m_fibers_y.end());
   
    // PHASE 3
    if (achieve & CGAL_BISOLVE_LIFTING) {
      _lifting(achieve, 
               this->ptr()->_m_fibers_x.begin(), this->ptr()->_m_fibers_x.end(),
               this->ptr()->_m_fibers_y.begin(), this->ptr()->_m_fibers_y.end());
    }
  }

protected:
    
  /*!
   * ensure that approximations are well seperated and minimal
   * required precisions for solutions are set
   */
  template < class InputIterator >
  void _increase_minprecs(InputIterator begin, InputIterator end) const {
     
    InputIterator it = begin;
    
    for (; it != end; it++) {
      
      bool sufficient = true;
      
      long prec = 2; // GMP need at least 2
      
      do {
       
       sufficient = true;
       
       // NOTE: Do not use Approximate_solution_1 here!
       typename Algebraic_kernel_d_1::Approximate_absolute_1 approx_abs;
       
       std::pair< Bound, Bound > app = approx_abs(it->value(), prec);
       
       if (it != begin) {
         std::pair< Bound, Bound > ppp = approx_abs(CGAL::cpp0x::prev(it)->value(), prec); 
         //std::cout << "ppps: "<< dbl(ppp.second) << std::endl;
         //std::cout << "apfi: "<< dbl(app.first) << std::endl;
         sufficient = sufficient && (CGAL::compare(app.first, ppp.second) == CGAL::LARGER);
       }
       if (sufficient && (CGAL::cpp0x::next(it) != end)) {
         std::pair< Bound, Bound > npp = approx_abs(CGAL::cpp0x::next(it)->value(), prec);
         //std::cout << "apse: "<< dbl(app.second) << std::endl;
         //std::cout << "nppf: "<< dbl(npp.second) << std::endl;
         sufficient = sufficient && (CGAL::compare(app.second, npp.first) == CGAL::SMALLER);
       } 
       
       if (sufficient) {
	    if (it != begin) {
              CGAL::cpp0x::prev(it)->increase_minprec(prec);
	    }
            it->increase_minprec(prec);
	    if (CGAL::cpp0x::next(it) != end) {
              CGAL::cpp0x::next(it)->increase_minprec(prec);
	    }
       } else {
         prec *= 2;
       }
       
      } while (!sufficient);
    }
    
    return;
  }


  //! isolatesunivariate real roots
  template < class FactorsInputIterator, class OutputIterator > 
  void _isolate(int v /* x or y */,
                FactorsInputIterator fac_begin,
                FactorsInputIterator fac_end,
                OutputIterator oi) const {
    
  
    // Isolate real roots of resultants (with multiplicities!)
    // -------------------------------------------------------
    
    // Refine until isolating intervals
    // --------------------------------

    Solution_1_container sols;
 
    typename Algebraic_kernel_d_1::Solve_1 solve;
    
    for (FactorsInputIterator fit = fac_begin; fit != fac_end; fit++) {
      
      std::list< Algebraic_real_1 > roots;
      Bisolve_telemetry_code(
        t_solve1[0].start();  
        t_solve1[v].start();
      )
      // we need to isolate all roots, to determine whether combinatorial filter can be used
      solve(fit->first, true, std::back_inserter(roots));
      Bisolve_telemetry_code(
        t_solve1[v].stop();  
        t_solve1[0].stop();
      )

      if (v == 1) {
        Bisolve_out("# x-roots (" << fit->second << "): "<< roots.size() << std::endl);
      } else {
        Bisolve_out("# y-roots (" << fit->second << "): "<< roots.size() << std::endl);
      }
      
      for (typename std::list< Algebraic_real_1 >::const_iterator 
             rit = roots.begin(); rit != roots.end(); rit++) {
        // handle-repped class!
        Solution_1 s(*rit, fit->second);
        sols.push_back(s);
      }
     
    }
    
    // ensure order
    Bisolve_telemetry_code(
      t_sort[0].start();
      t_sort[v].start();
    )
    std::sort(sols.begin(), sols.end(), Solution_1_less_than());
    Bisolve_telemetry_code(
      t_sort[v].stop();
      t_sort[0].stop();
    )

    _increase_minprecs(sols.begin(), sols.end());

    CGAL_assertion_code(
      // overlap check //
      {
        typename Solution_1_container::iterator sit = sols.begin();
        for (;sit != sols.end(); sit++) {
          if (sit != sols.begin()) {
            // TODO use Approximate_absolute_1
            CGAL_assertion_msg(CGAL::compare(CGAL::cpp0x::prev(sit)->value().high(), sit->value().low()) == CGAL::SMALLER, 
                               "Found left overlap after isolation");
          }
          if (CGAL::cpp0x::next(sit) != sols.end()) {
            // TODO use Approximate_absolute_1
            CGAL_assertion_msg(CGAL::compare(sit->value().high(), CGAL::cpp0x::next(sit)->value().low()) == CGAL::SMALLER, 
                               "Found right overlap after isolation");
          } 
        }
      }
    )

    std::copy(sols.begin(), sols.end(), oi);

  }


  //! isolate and seperate univariate real roots
  template < class FiberMapIterator, class FactorsInputIterator > 
  void _seperate(int v /* x or y */,
                 int degree_res,
                 FiberMapIterator sols_begin,
                 FiberMapIterator sols_end,
                 FactorsInputIterator fac_begin,
                 FactorsInputIterator fac_end) {
    
    // Refine until isolating intervals are well-seperating
    // ----------------------------------------------------
    
#if CGAL_BISOLVE_USE_RESULTANT_COFACTORS
#if CGAL_BISOLVE_RELAXED_RATIOS
    // Values: just 8
    Bound factor(8);
#else
    // Values: Choose n = 3 * degree of non-square-free resultant + 2
    Bound factor(Bound(3) * degree_res + Bound(2));
#endif
#else // !CGAL_BISOLVE_USE_RESULTANT_COFACTORS
      // Values: just 3
    Bound factor(3);
#endif // !CGAL_BISOLVE_USE_RESULTANT_COFACTORS

    typedef CGAL::internal::Well_seperating_interval< Solution_1, Refiner > Wellsep;
    
    Wellsep wellsep;
    
    for (FactorsInputIterator fit = fac_begin; fit != fac_end; fit++) {
      for (FiberMapIterator sit = sols_begin; sit != sols_end; sit++) {
        
        Solution_1 sol = sit->first; // handle rep?

        // filter out solutions not belonging to the current factors
#if CGAL_BISOLVE_ENABLE_NTL_FACTORIZE
        if (//sol.multiplicity() == fit->second ||
            kernel()->sign_at_1_object()(fit->first, sol.value()) != CGAL::ZERO) {
          continue;
        }
#else
        if (sol.multiplicity() != fit->second) {
          continue;
        }
#endif
        Bisolve_out("- Refine root (" << dbl(sol.value()) 
                    << ") to hold circle-tests, mult=" << sol.multiplicity() 
                    << ", minprec=" << sol.minprec() <<
                    std::endl)
        Bisolve_telemetry_code(
          t_wellsep[0].start();
          t_wellsep[v].start();
        )

        wellsep(sol, factor, fac_begin, fac_end, fit);

        Bisolve_telemetry_code(
          t_wellsep[v].stop();
          t_wellsep[0].stop();
        )
      }
    }
    
    CGAL_assertion_code(// overlap check //
     {
       for (FiberMapIterator sit = sols_begin; sit != sols_end; sit++) {
        
         if (sit != sols_begin) {
           //std::cout << "predr: " << CGAL::cpp0x::prev(sit)->first.rb().right() << std::endl;
           //std::cout << "thisl: " << sit->first.rb().left() << std::endl;
           CGAL_assertion_msg(CGAL::compare(CGAL::cpp0x::prev(sit)->first.rb().right(), sit->first.rb().left()) == CGAL::SMALLER, 
                              "Found left overlap after well-sep");
         }
         if (CGAL::cpp0x::next(sit) != sols_end) {
           //std::cout << "thisr: " << sit->first.rb().right() << std::endl;
           //std::cout << "succl: " << CGAL::cpp0x::next(sit)->first.rb().left() << std::endl;
           CGAL_assertion_msg(CGAL::compare(sit->first.rb().right(), CGAL::cpp0x::next(sit)->first.rb().left()) == CGAL::SMALLER, 
                              "Found right overlap after well-sep");
         } 
       }
     }
    )
  }


public:

  void set_resultant_in_x(const Polynomial_1& rx) {
    this->ptr()->_m_res_in_x = rx;
  }

  void set_resultant_in_y(const Polynomial_1& ry) {
    this->ptr()->_m_res_in_y = ry;
  }


public:

  Certifier_traits* certifier_traits() {
    
    //CGAL_precondition(this->ptr()->_m_achieved == CGAL_BISOLVE_ALL_BUT_LIFTING);

    if (!this->ptr()->_m_certifier_traits) {
      this->ptr()->_m_certifier_traits = 
        Certifier_traits(this->kernel(),
                         this->ptr()->_m_f, this->ptr()->_m_g, this->ptr()->_m_ft, this->ptr()->_m_gt, 
                         *this->ptr()->_m_res_in_x, *this->ptr()->_m_res_in_y);
      this->ptr()->_m_certifier_traits->_m_res_in_x_factors_ptr = &(this->ptr()->_m_res_in_x_factors);
      this->ptr()->_m_certifier_traits->_m_res_in_y_factors_ptr = &(this->ptr()->_m_res_in_y_factors);
    }
    
    return &(*this->ptr()->_m_certifier_traits);
  } 

  const Certifier_traits* certifier_traits() const {

    // TODO ensure that it's initialized
    CGAL_precondition(this->ptr()->_m_certifier_traits);

    return &(*this->ptr()->_m_certifier_traits);
  } 



protected:

  long upper_bound_log2_abs_front2back(const Solution_1& front, const Solution_1& back) const {

    // compute largest absolute boundary of the solutions' isolating intervals
    typedef CGAL::internal::Bitstream_coefficient_kernel< Bound > BCK_rat;
    typedef CGAL::internal::Bitstream_descartes_rndl_tree_traits< BCK_rat > Bs_traits_rat;

    Bs_traits_rat traits_rat;
    
    typename Bs_traits_rat::Upper_bound_log2_abs_approximator
    upper_bound_log2_abs_approximator = 
    traits_rat.upper_bound_log2_abs_approximator_object();
    
    Bound ll = front.rb().left();
    Bound hh = back.rb().right();
    
    if (CGAL::is_zero(front.rb().rad())) {
      ll -= Bound(1)/Bound(4);
    }
    if (CGAL::is_zero(back.rb().rad())) {
      hh += Bound(1)/Bound(4);
    }  
    Bound m = CGAL::max(CGAL::abs(ll), CGAL::abs(hh));
    
    Bisolve_out("Rational upper bound: " << dbl(m) << std::endl)

    bool is_zero;
    
    long ub;

    upper_bound_log2_abs_approximator.initial_upper_bound(m, ub, is_zero); // TODO do something with is_zero?
    
    return ub;
    
  }

  long grid_upper_bound_y_log2_abs() const {
    
    if (!this->ptr()->_m_grid_upper_bound_y_log2_abs) {
      
      long ub(0);

      if (!this->ptr()->_m_fibers_y.empty()) {

        Solution_1 front = (this->ptr()->_m_fibers_y.begin())->first;
        Solution_1 back = CGAL::cpp0x::prev(this->ptr()->_m_fibers_y.end())->first;
        
        ub = upper_bound_log2_abs_front2back(front, back);
      
      }
      
      this->ptr()->_m_grid_upper_bound_y_log2_abs = ub;
      
    }
    
    return *this->ptr()->_m_grid_upper_bound_y_log2_abs;
    
  }

protected:
  
  // PHASE 1: PROJECT
  // ================

  // Phase 1a: Compute resultants
  // ----------------------------

  void _resultants(int achieve) const {

    if (!(this->ptr()->_m_achieved & CGAL_BISOLVE_RESULTANT_IN_X) &&
        (achieve & CGAL_BISOLVE_RESULTANT_IN_X)) {
      
      Bisolve_out("Started computing resultant in x ... " << std::flush)
    
      if (!this->ptr()->_m_res_in_x) {
        Bisolve_telemetry_code(
          t_res[0].start();
          t_res[1].start();
        )
        this->ptr()->_m_res_in_x = 
          CGAL::resultant(this->ptr()->_m_f , this->ptr()->_m_g);

        if(!CGAL::is_zero(this->ptr()->_m_f) &&
            !CGAL::is_zero(this->ptr()->_m_g) &&
            CGAL::is_zero(*this->ptr()->_m_res_in_x)) {
            std::cerr << "Zero resultant in x..\n";
            throw internal::Zero_resultant_exception<Polynomial_2>(this->ptr()->_m_f);
        }
          
        Bisolve_telemetry_code(
          t_res[1].stop();
          t_res[0].stop();
        )
      } else {
        Bisolve_out(" already set ... " << std::flush)
      }

#if CGAL_BISOLVE_DEBUG
      std::stringstream ss;
      CGAL::set_pretty_mode(ss);
      ss << *this->ptr()->_m_res_in_x << "\n";
      std::ofstream oo("oo.oo"); // ofstream has buffer overflow so use
      oo << ss.str();            // strings first
      oo.close();
#endif

      // canonicalize should be called after saving to the file since otherwise
      // results will differ from that of computed by algorithm
      Bisolve_telemetry_code(
        t_res[0].start();
        t_res[1].start();
      )
      this->ptr()->_m_res_in_x = CGAL::canonicalize(*this->ptr()->_m_res_in_x);
      Bisolve_telemetry_code(
        t_res[1].stop();
        t_res[0].stop();
      )
   
      this->ptr()->_m_achieved |= CGAL_BISOLVE_RESULTANT_IN_X;  

      Bisolve_out(" finished")
      Bisolve_telemetry_code(Bisolve_out(" in time: " << t_res[1].time()))
      Bisolve_out("." << std::endl << std::endl)

    }
 
    if (!(this->ptr()->_m_achieved & CGAL_BISOLVE_RESULTANT_IN_Y) &&
        (achieve & CGAL_BISOLVE_RESULTANT_IN_Y)) {

      Bisolve_out("Started computing resultant in y ... " << std::flush)

      typename PT_2::Swap swap;
      this->ptr()->_m_ft = swap(this->ptr()->_m_f,0,1);
      this->ptr()->_m_gt = swap(this->ptr()->_m_g,0,1);
      
      if (!this->ptr()->_m_res_in_y) {
        Bisolve_telemetry_code(
          t_res[0].start();
          t_res[2].start();
        )
        this->ptr()->_m_res_in_y = 
            CGAL::resultant(this->ptr()->_m_ft , this->ptr()->_m_gt);

        if(!CGAL::is_zero(this->ptr()->_m_ft) &&
            !CGAL::is_zero(this->ptr()->_m_gt) && CGAL::is_zero(*this->ptr()->_m_res_in_y)) {
            std::cerr << "Zero resultant in y..\n";
            throw internal::Zero_resultant_exception<Polynomial_2>(this->ptr()->_m_f);
        }
            
        Bisolve_telemetry_code(
          t_res[2].stop();
          t_res[0].stop();
        )
      } else {
        Bisolve_out(" already set ... " << std::flush)
      }
      
#if CGAL_BISOLVE_DEBUG
      std::stringstream ss;
      CGAL::set_pretty_mode(ss);
      ss << *this->ptr()->_m_res_in_y;
      std::ofstream oo("oo.oo"); // ofstream has buffer overflow so use
      oo << ss.str();            // strings first
      oo.close();
#endif
      // canonicalize should be called after saving to the file since otherwise
      // results will differ from that of computed by algorithm
      Bisolve_telemetry_code(
        t_res[0].start();
        t_res[2].start();
      )
      this->ptr()->_m_res_in_y = CGAL::canonicalize(*this->ptr()->_m_res_in_y);
      Bisolve_telemetry_code(
        t_res[2].stop();
        t_res[0].stop();
      )

      this->ptr()->_m_achieved |= CGAL_BISOLVE_RESULTANT_IN_Y;  

      Bisolve_out(" finished")
      Bisolve_telemetry_code(Bisolve_out(" in time: " << t_res[2].time()))
      Bisolve_out("." << std::endl << std::endl)

    }
    

#if CGAL_BISOLVE_WRITE_ONLY_RESULTANTS // TODO implement in terms of achieve!
      // early exits if we're only interested in resultants
      //return oi;
#endif 
  }


  // Phase 1b: Make resultant square-free
  // ------------------------------------

  void _square_free_factorization(int achieve) const {
 
#if !CGAL_BISOLVE_ENABLE_NTL_FACTORIZE
    typename PT_1::Square_free_factorize_up_to_constant_factor square_free_factorize;
#endif 

    if (!(this->ptr()->_m_achieved & CGAL_BISOLVE_SQUARE_FREE_FACTORIZATION_IN_X) &&
        (achieve & CGAL_BISOLVE_SQUARE_FREE_FACTORIZATION_IN_X)) {
 
      Bisolve_out("Started square-free factorization for resultant in x ... " << std::flush)
 
      CGAL_assertion(this->ptr()->_m_achieved & CGAL_BISOLVE_RESULTANT_IN_X);
      
      Bisolve_telemetry_code(t_sqfree[0].start();)
      Bisolve_telemetry_code(t_sqfree[1].start();)
#if CGAL_BISOLVE_ENABLE_NTL_FACTORIZE
        // from symbolic_exports.h
      CGAL::factorize_NTL(*this->ptr()->_m_res_in_x, std::back_inserter(this->ptr()->_m_res_in_x_factors));
#else
      square_free_factorize(*this->ptr()->_m_res_in_x, std::back_inserter(this->ptr()->_m_res_in_x_factors));
#endif
      Bisolve_telemetry_code(t_sqfree[1].stop();)
      Bisolve_telemetry_code(t_sqfree[0].stop();)
        
      this->ptr()->_m_achieved |= CGAL_BISOLVE_SQUARE_FREE_FACTORIZATION_IN_X;  

      Bisolve_out(" finished")
      Bisolve_telemetry_code(Bisolve_out(" in time: " << t_sqfree[1].time()))
      Bisolve_out("." << std::endl << std::endl)

    }
    
    if (!(this->ptr()->_m_achieved & CGAL_BISOLVE_SQUARE_FREE_FACTORIZATION_IN_Y) &&
        (achieve & CGAL_BISOLVE_SQUARE_FREE_FACTORIZATION_IN_Y)) {
      
      Bisolve_out("Started square-free factorization for resultant in y ... " << std::flush)

      CGAL_assertion(this->ptr()->_m_achieved & CGAL_BISOLVE_RESULTANT_IN_Y);

      Bisolve_telemetry_code(t_sqfree[0].start();)
      Bisolve_telemetry_code(t_sqfree[2].start();)
#if CGAL_BISOLVE_ENABLE_NTL_FACTORIZE
        // from symbolic_exports.h
      CGAL::factorize_NTL(*this->ptr()->_m_res_in_y, std::back_inserter(this->ptr()->_m_res_in_y_factors));
#else
      square_free_factorize(*this->ptr()->_m_res_in_y, std::back_inserter(this->ptr()->_m_res_in_y_factors));
#endif
      Bisolve_telemetry_code(t_sqfree[2].stop();)
      Bisolve_telemetry_code(t_sqfree[0].stop();)

      this->ptr()->_m_achieved |= CGAL_BISOLVE_SQUARE_FREE_FACTORIZATION_IN_Y; 

      Bisolve_out(" finished")
      Bisolve_telemetry_code(Bisolve_out(" in time: " << t_sqfree[2].time()))
      Bisolve_out("." << std::endl << std::endl)
 
    }
    
  }


  // Phase 1c: Isolation
  // -------------------

  void _isolation(int achieve) {

    if (!(this->ptr()->_m_achieved & CGAL_BISOLVE_ISOLATION_IN_X) &&
        (achieve & CGAL_BISOLVE_ISOLATION_IN_X)) {
      
      Bisolve_out("Started isolation for resultant in x ... " << std::flush)

      CGAL_assertion(this->ptr()->_m_achieved & CGAL_BISOLVE_SQUARE_FREE_FACTORIZATION_IN_X);
      
      Solution_1_container sols_x;

      // TODO save Solve_1 for y if all x-roots are simple (not possible with active-intervals-sets)???
      _isolate(1 /* = x */, 
               this->ptr()->_m_res_in_x_factors.begin(), this->ptr()->_m_res_in_x_factors.end(),
               std::back_inserter(sols_x));
      
      for (Solution_1_container_const_iterator sit = sols_x.begin(); 
           sit != sols_x.end(); sit++) {
        
        // create fiber
        Fiber_x fiber;

        const Algebraic_real_1 alpha = sit->value();
        
        //std::pair< Fiber_map_x_iterator, bool > res = 
          this->ptr()->_m_fibers_x.insert(std::make_pair(*sit, fiber));
      }

      Bisolve_out("done. #CandidatesX: " << this->ptr()->_m_fibers_x.size() << std::endl << std::endl)
      
      this->ptr()->_m_achieved |= CGAL_BISOLVE_ISOLATION_IN_X;  
    }
    
    if (!(this->ptr()->_m_achieved & CGAL_BISOLVE_ISOLATION_IN_Y) &&
        (achieve & CGAL_BISOLVE_ISOLATION_IN_Y)) {
      
      Bisolve_out("Started isolation for resultant in y ... " << std::flush)

      CGAL_assertion(this->ptr()->_m_achieved & CGAL_BISOLVE_SQUARE_FREE_FACTORIZATION_IN_Y);

      Solution_1_container sols_y;

      // TODO save Solve_1 for y if all x-roots are simple (not possible with active-intervals-sets)???
      _isolate(2 /* = x */, 
               this->ptr()->_m_res_in_y_factors.begin(), this->ptr()->_m_res_in_y_factors.end(),
               std::back_inserter(sols_y));
      
      for (Solution_1_container_const_iterator sit = sols_y.begin(); 
           sit != sols_y.end(); sit++) {

        // create fiber
        Fiber_y fiber;

        const Algebraic_real_1 alpha = sit->value();

        //std::pair< Fiber_map_y_iterator, bool > res = 
          this->ptr()->_m_fibers_y.insert(std::make_pair(*sit, fiber));

      }

      Bisolve_out("done. #CandidatesY: " << this->ptr()->_m_fibers_y.size() << std::endl << std::endl)

      this->ptr()->_m_achieved |= CGAL_BISOLVE_ISOLATION_IN_Y;  
    }
    
  }


  // PHASE 2: Seperation
  // ===================

void _seperation(int achieve,
                 Fiber_map_x_iterator x_begin, 
                 Fiber_map_x_iterator x_beyond,
                 Fiber_map_y_iterator y_begin,
                 Fiber_map_y_iterator y_beyond) {
  
    if (!(this->ptr()->_m_achieved & CGAL_BISOLVE_SEPERATION_IN_X) &&
        (achieve & CGAL_BISOLVE_SEPERATION_IN_X)) {
      
      Bisolve_out("Started well-seperation for resultant in x ... " << std::flush)

      CGAL_assertion(this->ptr()->_m_achieved & CGAL_BISOLVE_ISOLATION_IN_X);

      _seperate(1 /* = x */, 
                CGAL::degree(*this->ptr()->_m_res_in_x), 
                x_begin, x_beyond,
                this->ptr()->_m_res_in_x_factors.begin(), this->ptr()->_m_res_in_x_factors.end());
      
      Bisolve_out("done." << std::endl)
      
      this->ptr()->_m_achieved |= CGAL_BISOLVE_SEPERATION_IN_X;  
    }

    if (!(this->ptr()->_m_achieved & CGAL_BISOLVE_SEPERATION_IN_Y) &&
        (achieve & CGAL_BISOLVE_SEPERATION_IN_Y)) {
      
      Bisolve_out("Started well-seperation for resultant in y ... " << std::flush)

      CGAL_assertion(this->ptr()->_m_achieved & CGAL_BISOLVE_ISOLATION_IN_Y);

      _seperate(2 /* = x */, 
                CGAL::degree(*this->ptr()->_m_res_in_y), 
                y_begin, y_beyond,
                this->ptr()->_m_res_in_y_factors.begin(), this->ptr()->_m_res_in_y_factors.end());
      
      Bisolve_out("done." << std::endl)
    
      this->ptr()->_m_achieved |= CGAL_BISOLVE_SEPERATION_IN_Y;  
    }
  }


  // PHASE 3: Lifting
  // ================

  void _lifting(int achieve, 
                Fiber_map_x_iterator x_begin, 
                Fiber_map_x_iterator x_beyond,
                Fiber_map_y_iterator y_begin,
                Fiber_map_y_iterator y_beyond) {

    // TODO add correct CGAL_preconditions
    if (!(achieve & CGAL_BISOLVE_LIFTING)) {
      return;
    }

    if (this->ptr()->_m_achieved & CGAL_BISOLVE_LIFTING) {
      return;
    }

    Bisolve_out("================== #GridCandidates: " <<
                (std::distance(x_begin, x_beyond) * std::distance(y_begin, y_beyond))  << " ==============\n")

    // for copying solutions to fibers define
    // inserter for fibers in both direction
    Fiber_map_inserter fiber_map_inserter(this->ptr()->_m_fibers_x, this->ptr()->_m_fibers_y);
    
    Bisolve_telemetry_code(t_validate.start();)
    
    bool full_x_range = (x_begin == this->ptr()->_m_fibers_x.begin() && x_beyond == this->ptr()->_m_fibers_x.end());
    bool full_y_range = (y_begin == this->ptr()->_m_fibers_y.begin() && y_beyond == this->ptr()->_m_fibers_y.end());
    Bisolve_out("full x-range: " << full_x_range << std::endl);
    Bisolve_out("full y-range: " << full_y_range << std::endl);

    // set traits
    Certifier certifier(this->certifier_traits());

    // copy ALL solutions, to not get problems with indices below
    certifier.fill_xy_solution_lists(this->ptr()->_m_fibers_x.begin(), this->ptr()->_m_fibers_x.end(),
                                     this->ptr()->_m_fibers_y.begin(), this->ptr()->_m_fibers_y.end());

    // =======================================================================
    // first run

    Bisolve_out("\n\n=============== first run ==============\n")
    Bisolve_out("processing " << std::distance(x_begin, x_beyond) << " x-fibers\n")
    
    // init

#if CGAL_BISOLVE_DISABLE_COMBINATORIAL_CERTIFICATION
    // disable combinatorial filter
    certifier.set_combinatorial(false);
#else
    // combinatorial only works if full y-range
    certifier.set_combinatorial(full_y_range);
#endif
    
#if !CGAL_BISOLVE_DISABLE_BIDIRECTIONAL_CERTIFICATION
    // Remark: bidirection only works if full x-range is given
    // for bidirection the first direction needs to be run dry
    certifier.set_dry_run(full_x_range); 
#endif         
    
    // array to index y-solutions
    std::vector< int > indices(this->ptr()->_m_fibers_y.size()); 
    
    // setup the list of y-intervals
    for (unsigned i = 0; i < this->ptr()->_m_fibers_y.size(); i++) {
      indices[i] = i;
    }
    
    // FIRST LOOP
    for (Fiber_map_x_iterator xit = x_begin; xit != x_beyond; xit++) {
      
      CGAL_assertion_code(Bound r = xit->first.rb().rad();)
      CGAL_assertion(CGAL::sign(r) != CGAL::ZERO);


      // choose right fiber
      certifier.set_fiber(&(xit->second));
#if CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS && CGAL_BISOLVE_CURVE_ANALYSIS
#error Curve analysis needs active intervals
#endif
      
      //! list of found solutions
      std::list< Algebraic_real_2 > fiber_sols;
      
      typename std::vector< int >::iterator ibegin = indices.begin();
      
      std::advance(ibegin, std::distance(this->ptr()->_m_fibers_y.begin(), y_begin));
      typename std::vector< int >::iterator iend = ibegin;
      std::advance(iend, std::distance(y_begin, y_beyond));
      
      // TODO use iterator instead of i
      int i = std::distance(this->ptr()->_m_fibers_x.begin(), xit);
      certifier(i, ibegin, iend,
                grid_upper_bound_y_log2_abs() + 1,
                std::back_inserter(fiber_sols));

      // copy solutions to fibers
      fiber_map_inserter(fiber_sols.begin(), fiber_sols.end());
    }
    
#if !CGAL_BISOLVE_DISABLE_BIDIRECTIONAL_CERTIFICATION
    if (full_x_range) { // enable only if full range

      // =======================================================================
      // second run

      Bisolve_out("\n\n=============== second run ==============\n")
      Bisolve_out("processing " << certifier.fiber_list().size() << " y-fibers\n")
      
      // init

      certifier.set_dry_run(false);
      certifier.set_reverse_dir(true);
      
#if CGAL_BISOLVE_DISABLE_COMBINATORIAL_CERTIFICATION
      certifier.set_combinatorial(false);
#else
      // combinatorial in other only works if full x-range
      certifier.set_combinatorial(full_x_range);
#endif
      
      // SECOND LOOP
      for (typename Certifier::Fiber_list::const_iterator yi = certifier.fiber_list().begin(); 
           yi != certifier.fiber_list().end(); yi++) {
        
        // compute index
        int y_index = yi->first;
        typename Certifier::Index_list ilist = yi->second;

        Fiber_map_y_iterator yit = this->ptr()->_m_fibers_y.begin();
        std::advance(yit, y_index);
        
        CGAL_assertion_code(Bound r = yit->first.rb().rad();)
        CGAL_assertion(CGAL::sign(r) != CGAL::ZERO);
        
        // determine upper bound
        Fiber_map_x_const_iterator xfrit = this->ptr()->_m_fibers_x.begin();
        std::advance(xfrit, ilist.front());
        Fiber_map_x_const_iterator xbait = this->ptr()->_m_fibers_x.begin();
        std::advance(xbait, ilist.back());
        
        Solution_1 front = xfrit->first;
        Solution_1 back = xbait->first;
        
        long ub = upper_bound_log2_abs_front2back(front, back);
        
#if !CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS
        // choose right fiber
        certifier.set_fiber(&(yit->second));
#endif
        std::list< Algebraic_real_2 > fiber_sols;
        
        // TODO use iterator instead of y_index
        certifier(y_index, ilist.begin(), ilist.end(),
                  ub + 1,
                  std::back_inserter(fiber_sols));
        
        // copy to fibers
        fiber_map_inserter(fiber_sols.begin(), fiber_sols.end());
      }
    }

#endif     
    
    Bisolve_out("_lifting end" << std::endl)

    // achieved if full ranges
    if (full_x_range && full_y_range) {
        this->ptr()->_m_achieved |= CGAL_BISOLVE_LIFTING;
    }

    
    Bisolve_telemetry_code(t_validate.stop();)

  }
  
public:

  //!\name Function call operators
  //!@{

  //! function-call operator to solve bivariate system
  template < class OutputIterator >
  OutputIterator operator()(OutputIterator oi) {
    
    // PHASE 1
    _resultants(CGAL_BISOLVE_RESULTANT_IN_X | CGAL_BISOLVE_RESULTANT_IN_Y);
    _square_free_factorization(CGAL_BISOLVE_SQUARE_FREE_FACTORIZATION_IN_X |
                               CGAL_BISOLVE_SQUARE_FREE_FACTORIZATION_IN_Y);
    _isolation(CGAL_BISOLVE_ISOLATION_IN_X | CGAL_BISOLVE_ISOLATION_IN_Y);

    Fiber_map_x_iterator x_begin = this->ptr()->_m_fibers_x.begin();
    Fiber_map_x_iterator x_beyond = this->ptr()->_m_fibers_x.end();
 
    Fiber_map_y_iterator y_begin = this->ptr()->_m_fibers_y.begin();
    Fiber_map_y_iterator y_beyond = this->ptr()->_m_fibers_y.end();
    
    // PHASE 2
    _seperation(CGAL_BISOLVE_SEPERATION_IN_X | CGAL_BISOLVE_SEPERATION_IN_Y, 
                x_begin, x_beyond, y_begin, y_beyond);
    
    // PHASE 3
    _lifting(CGAL_BISOLVE_LIFTING, x_begin, x_beyond, y_begin, y_beyond);
    
    // copy solutions stored in fibers to oi
    for (Fiber_map_x_const_iterator fit = this->ptr()->_m_fibers_x.begin();
         fit != this->ptr()->_m_fibers_x.end(); fit++) {
      std::copy(fit->second.roots.begin(), fit->second.roots.end(), oi);
    }

    return oi;
    
  }
  
  
#if 0
// TODO operator disabled as not compatible with achievement structure so far

  //! function-call operator to solve bivariate system in a box
  template < class OutputIterator >
    OutputIterator operator()(const Box_2& box, OutputIterator oi) {

    // PHASE 1
    _resultants(CGAL_BISOLVE_RESULTANT_IN_X | CGAL_BISOLVE_RESULTANT_IN_Y);
    _square_free_factorization(CGAL_BISOLVE_SQUARE_FREE_FACTORIZATION_IN_X |
                               CGAL_BISOLVE_SQUARE_FREE_FACTORIZATION_IN_Y);
    _isolation(CGAL_BISOLVE_ISOLATION_IN_X | CGAL_BISOLVE_ISOLATION_IN_Y);
    
    Fiber_map_x_iterator x_begin = this->ptr()->_m_fibers_x.lower_bound(Solution_1(Algebraic_real_1(box[0]), -1));
    Fiber_map_x_iterator x_beyond = this->ptr()->_m_fibers_x.upper_bound(Solution_1(Algebraic_real_1(box[1]), -1));
 
    Fiber_map_y_iterator y_begin = this->ptr()->_m_fibers_y.lower_bound(Solution_1(Algebraic_real_1(box[2]), -1));
    Fiber_map_y_iterator y_beyond = this->ptr()->_m_fibers_y.upper_bound(Solution_1(Algebraic_real_1(box[3]), -1));

    // PHASE 2
    _seperation(CGAL_BISOLVE_SEPERATION_IN_X | CGAL_BISOLVE_SEPERATION_IN_Y,
                x_begin, x_beyond, y_begin, y_beyond);
    
    // PHASE 3
    _lifting(CGAL_BISOLVE_LIFTING, x_begin, x_beyond, y_begin, y_beyond);
    
    // copy solutions stored in fibers to oi
    for (Fiber_map_x_const_iterator fit = this->ptr()->_m_fibers_x.begin();
         fit != this->ptr()->_m_fibers_x.end(); fit++) {
      std::copy(fit->second.roots.begin(), fit->second.roots.end(), oi);
    }

    // TODO maintain status, e.g., if some solutions are already computed then search for right output 
    //      only the right one, otherwise compute only the missing ones

  }
#endif

  //! function-call operator to solve bivariate system
  template < class OutputIterator >
OutputIterator operator()(const Algebraic_real_1& x, OutputIterator oi, bool local = false) {
    
    // PHASE 1
    _resultants(CGAL_BISOLVE_RESULTANT_IN_X | CGAL_BISOLVE_RESULTANT_IN_Y);
    _square_free_factorization(CGAL_BISOLVE_SQUARE_FREE_FACTORIZATION_IN_X |
                               CGAL_BISOLVE_SQUARE_FREE_FACTORIZATION_IN_Y);
    _isolation(CGAL_BISOLVE_ISOLATION_IN_X | CGAL_BISOLVE_ISOLATION_IN_Y);

    Fiber_map_x_iterator x_begin, x_beyond;

    if (CGAL_BISOLVE_DISABLE_LOCALIZED_OPERATOR || !local) {
      x_begin = this->ptr()->_m_fibers_x.begin();
      x_beyond = this->ptr()->_m_fibers_x.end();
    } else {
      x_begin = fiber_map_iterator(x);
      if (x_begin == this->ptr()->_m_fibers_x.end()) {
        return oi;
      }
      x_beyond = CGAL::cpp0x::next(x_begin);
    }
    Fiber_map_y_iterator y_begin = this->ptr()->_m_fibers_y.begin();
    Fiber_map_y_iterator y_beyond = this->ptr()->_m_fibers_y.end();

    // PHASE 2
    _seperation(CGAL_BISOLVE_SEPERATION_IN_X | CGAL_BISOLVE_SEPERATION_IN_Y,
                x_begin, x_beyond, y_begin, y_beyond);
    
    // PHASE 3
    _lifting(CGAL_BISOLVE_LIFTING, x_begin, x_beyond, y_begin, y_beyond);
  
    // copy solutions stored in fibers to oi
    Fiber_map_x_const_iterator fit = fiber_map_const_iterator(x);
    if (fit != this->ptr()->_m_fibers_x.end()) {
      std::copy(fit->second.roots.begin(), fit->second.roots.end(), oi);
    }
    
    // TODO maintain status, e.g., if some solutions are already computed then search for right output 
    //      only the right one, otherwise compute only the missing ones
    
    return oi;
  }

  //!@}

  //!\name Access members
  //!@{

  const Polynomial_1& res_y_f_fy() const {
    CGAL_precondition(this->ptr()->_m_achieved & CGAL_BISOLVE_RESULTANT_IN_X);
    return *this->ptr()->_m_res_in_x;
  }

  template < class OutputIterator >
  OutputIterator res_y_f_fy_factors(OutputIterator oi) const {
      CGAL_precondition(this->ptr()->_m_achieved & CGAL_BISOLVE_RESULTANT_IN_X);
      std::copy(this->ptr()->_m_res_in_x_factors.begin(), this->ptr()->_m_res_in_x_factors.end(), oi);
      return oi;
  }

  // REMARK iterators changed
  //! returns range of x-solutions
  std::pair< Fiber_map_x_const_iterator, Fiber_map_x_const_iterator > x_range() const {
    return std::make_pair(this->ptr()->_m_fibers_x.begin(), this->ptr()->_m_fibers_x.end());
  }

  //! returns range of y-solutions
  std::pair< Fiber_map_y_const_iterator, Fiber_map_y_const_iterator > y_range() const {
    return std::make_pair(this->ptr()->_m_fibers_y.begin(), this->ptr()->_m_fibers_y.end());
  }

  template < class OutputIterator > 
  OutputIterator copy_solutions_in_x(OutputIterator oi) const {
    
    for (Fiber_map_x_const_iterator it = this->ptr()->_m_fibers_x.begin(); 
         it != this->ptr()->_m_fibers_x.end(); it++) {
      *oi++ = std::make_pair(it->first.value(), it->first.multiplicity());
    }
    return oi;
  }

  Fiber_map_x_iterator fiber_map_iterator(const Algebraic_real_1& x) {
    return this->ptr()->_m_fibers_x.find(Solution_1(x, -1));
  }

  Fiber_map_x_const_iterator fiber_map_const_iterator(const Algebraic_real_1& x) const {
    return this->ptr()->_m_fibers_x.find(Solution_1(x, -1));
  }

  bool has_fiber(const Algebraic_real_1& x) const {
    return fiber_map_const_iterator(x) != this->ptr()->_m_fibers_x.end();
  }

  Fiber_x& fiber(const Algebraic_real_1& x) {
    CGAL_precondition(has_fiber(x));
    return fiber_map_iterator(x)->second;
  }

  const Fiber_x& fiber(const Algebraic_real_1& x) const {
    CGAL_precondition(has_fiber(x));
    return fiber_map_const_iterator(x)->second;
  }
  
  bool has_vert_line_f(const Algebraic_real_1& x) const {
    CGAL_assertion(has_fiber(x));
    // TODO add precond that certifier has been called for this fiber
    return fiber_map_const_iterator(x)->second.vert_line_f;
  }

  bool has_vert_line_g(const Algebraic_real_1& x) const {
    CGAL_assertion(has_fiber(x));
    // TODO add precond that certifier has been called for this fiber
    return fiber_map_const_iterator(x)->second.vert_line_g;
  }

  //@}

}; // Bi_solve_2

} // namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_2_BI_SOLVE_2_H
// EOF
