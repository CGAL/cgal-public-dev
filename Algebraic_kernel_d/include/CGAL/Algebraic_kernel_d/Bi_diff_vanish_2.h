// Copyright (c) 2010, 2011, 2012 Max-Planck-Institut fuer Informatik (Germany).
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
// Author(s): Eric Berberich  <eric.berberich@cgal.org>

#ifndef CGAL_BI_DIFF_VANISH_2_H
#define CGAL_BI_DIFF_VANISH_2_H

/*! \file
 * The header file for the Bi_diff_vanish_2 class.
 */

#include <CGAL/config.h>

#include <set>
#include <map>

#include <CGAL/Algebraic_kernel_2/Bi_solve_2_flags.h>
#include <CGAL/Algebraic_kernel_2/Bi_solve_2.h>

namespace CGAL {

namespace internal {
  
//! template rep class for new bi-solve
template < class CertifierTraits >
class Bi_diff_vanish_2_rep : public Bi_solve_2_rep< CertifierTraits > {

public:
    
  typedef CertifierTraits Certifier_traits;

  typedef Bi_solve_2_rep< Certifier_traits > Base;

  typedef typename Certifier_traits::Algebraic_kernel_d_1 Algebraic_kernel_d_1;

  typedef typename Algebraic_kernel_d_1::Coefficient Coefficient;

  typedef typename Base::Polynomial_2 Polynomial_2;

  typedef typename Base::Algebraic_real_1 Algebraic_real_1;
  
  typedef typename Base::Algebraic_real_2 Algebraic_real_2;
  
  typedef CGAL::Bi_solve_2< Certifier_traits > Pure_bi_solve_2;

  typedef typename Algebraic_kernel_d_1::Multiplicity_type Multiplicity_type;
  
  typedef std::set< Algebraic_real_2 > Solution_set;

  typedef std::map< std::list< Polynomial_2 >, Solution_set > Solution_set_map;

  typedef std::map< Algebraic_real_1, Solution_set_map > Solution_sets;

  typedef std::map< std::pair< Polynomial_2, Polynomial_2 >, Pure_bi_solve_2 > Bisolves_map;

  typedef std::map< Algebraic_real_2, int > Diff_map;

public:

  //! standard constructor
  Bi_diff_vanish_2_rep() :
    Base() {
  }

  //! standard constructor
  Bi_diff_vanish_2_rep(const Polynomial_2& f) :
    Base(f) {
  }

public:

  mutable std::vector< Polynomial_2 > _m_diffs;

  mutable Bisolves_map _m_bisolves;

  mutable Solution_sets _m_sets;
  
  mutable std::map< Algebraic_real_1, Diff_map > _m_diff_counters;

}; // Bi_diff_vanish_2_rep

} // namespace internal


//! template class for vanishing differenations
template < class CertifierTraits >
  class Bi_diff_vanish_2 : public Bi_solve_2< CertifierTraits, CGAL::internal::Bi_diff_vanish_2_rep< CertifierTraits > > {

public:
    
  //! this instance's first template parameter
  typedef CertifierTraits Certifier_traits;


  //! rep class
  typedef CGAL::internal::Bi_diff_vanish_2_rep< Certifier_traits > Rep;
  
  //! base class
  typedef Bi_solve_2< Certifier_traits, Rep > Base;

  //! type of AK_1
  typedef typename Certifier_traits::Algebraic_kernel_d_1 Algebraic_kernel_d_1;

  //! Pure_bi_solve_2
  typedef typename Rep::Pure_bi_solve_2 Pure_bi_solve_2;

  //! coefficient type
  typedef typename Algebraic_kernel_d_1::Coefficient Coefficient;

  //! type for algebraic reals
  typedef typename Rep::Algebraic_real_1 Algebraic_real_1;

  //! type of bivariate polynomial
  typedef typename Rep::Polynomial_2 Polynomial_2;

  typedef typename Rep::Polynomial_1 Polynomial_1;

  //! type for multiplicity
  typedef typename Rep::Multiplicity_type Multiplicity_type;

  //! type of bivariate algebraic real
  typedef typename Rep::Algebraic_real_2 Algebraic_real_2;

protected:

  typedef typename Rep::Solution_set Solution_set;

  typedef typename Rep::Solution_set_map Solution_set_map;
  
  typedef typename Rep::Solution_sets Solution_sets;

  typedef typename Rep::Bisolves_map Bisolves_map;

  typedef typename Rep::Diff_map Diff_map;

public:

  //!\name Constuctors
  //!@{

  //! default constructor
  Bi_diff_vanish_2() {
  }

  //! standard constructor
  Bi_diff_vanish_2(const Polynomial_2& f) :
  Base(f, 
       Bi_diff_vanish_2::CGAL_BISOLVE_RESULTANT_IN_X | 
       Bi_diff_vanish_2::CGAL_BISOLVE_SQUARE_FREE_FACTORIZATION_IN_X | 
       Bi_diff_vanish_2::CGAL_BISOLVE_ISOLATION_IN_X) {

    // some init
    this->ptr()->_m_diffs.push_back(Base::f());
    this->ptr()->_m_diffs.push_back(Base::g());
  }

  //!@}

public:

  //! function-call operator to compute which differentiations vanish
  template < class OutputIterator >
    OutputIterator operator()(const Algebraic_real_1& x, Multiplicity_type mult, OutputIterator oi, bool local = false) {

    // std::cout << "x: "<< CGAL::to_double(x) << std::endl;

    Bisolve_telemetry_code(t_bdv.start();)
    
    typename std::map< Algebraic_real_1, std::map< Algebraic_real_2, int > >::iterator it = 
      this->ptr()->_m_diff_counters.find(x);

    // no sets computed at x
    if (it == this->ptr()->_m_diff_counters.end()) {

      // map to store number of vanishing diff for each solution
      std::map< Algebraic_real_2, int > diffs;

      std::list< Polynomial_2 > chain;
      chain.push_back(Base::f());

      for (int i = 1; i <= CGAL::degree(Base::f()); i++) {
        
        chain.push_back(diff(i));

        Solution_set seti = solution_set(x, chain, local);
        if (seti.empty()) {
          break;
        }
   
        // increase diff counter by one for each element in nextset
        for (typename Solution_set::const_iterator sit = seti.begin();
             sit != seti.end(); sit++) {
          //std::cout << "BDV increase by one: " << *sit << std::endl;
          if (1 == i) {
            diffs[*sit] = 1;
          } else {
            diffs[*sit] += 1;
          }
        }
        
      }

      // store diffs
      it = this->ptr()->_m_diff_counters.insert(it, std::make_pair(x, diffs));
    }
    
    // copy all elements 
    for (typename Diff_map::const_iterator dit = it->second.begin();
         dit != it->second.end(); dit++) {
      // std::cout << "BDV rit: "<< dit->first << " mult: " << dit->second << std::endl;
      *oi++ = *dit;
    }

    Bisolve_telemetry_code(t_bdv.stop();)

    return oi;
    
  }

  //! number of diffs stored
  int max_diff() const {
    return this->ptr()->_m_diffs.size();
  }

  //! return a given diff
  const Polynomial_2& diff(size_t i) const {
    CGAL_precondition(i > 0);
    
    if (i >= this->ptr()->_m_diffs.size()) {
      this->ptr()->_m_diffs.resize(std::max(i+1, this->ptr()->_m_diffs.size()));
      this->ptr()->_m_diffs[i] =
            CGAL::canonicalize(CGAL::differentiate(diff(i-1)));
    }

    CGAL_precondition(i < this->ptr()->_m_diffs.size());
    return this->ptr()->_m_diffs[i];
  }


 protected:

  Solution_set solution_set(const Algebraic_real_1& x, std::list< Polynomial_2 > chain, bool local) {

    if (chain.size() < 2) {
      return Solution_set(); // empty;
    }
    
    typename Solution_sets::iterator it = this->ptr()->_m_sets.find(x);

    if (it == this->ptr()->_m_sets.end()) {
      
      CGAL_assertion(2 == chain.size());
      // chain[0] = f, chain[1] = f_y

      Solution_set fiber_roots;
      
      // initial set of solutions
      // std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BASE START" << std::endl;
      Bisolve_telemetry_code(t_bdv_bs0.start();)
      Base::operator()(x, std::insert_iterator< Solution_set >(fiber_roots, fiber_roots.begin()), local);
      Bisolve_telemetry_code(t_bdv_bs0.stop();)
      // std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BASE STOP" << std::endl;
     
      //std::cout << "#fb: " << fiber_roots.size() << std::endl;
     
      Solution_set_map smap;
      smap.insert(std::make_pair(chain, fiber_roots));

      it = this->ptr()->_m_sets.insert(it, std::make_pair(x,smap));
      
    }
    
    CGAL_assertion(it != this->ptr()->_m_sets.end());
    
    typename Solution_set_map::iterator sit = it->second.find(chain);

    if (sit == it->second.end()) {

      Solution_set next;
      
      std::list< Polynomial_2 > rchain;
      std::copy(chain.begin(), CGAL::cpp0x::prev(chain.end()), std::back_inserter(rchain));

      //std::cout << "rchainlength: " << rchain.size() << std::endl;
      
      typename std::list< Polynomial_2 >::const_iterator last, pred;
      last = CGAL::cpp0x::prev(chain.end());
      pred = CGAL::cpp0x::prev(last);
     
      Polynomial_1 gres_x(0), gres_y(0);
      bool failed_x = true;
      bool failed_y = true;

      Bisolve_telemetry_code(t_common.start();)

      Polynomial_2 h;
      h = CGAL::gcd_up_to_constant_factor(*pred, *last);

      Bisolve_telemetry_code(t_common.stop();)
      
      //std::cout << "BDV h: " << h << std::endl;
            
      if (h != Polynomial_2(Polynomial_1(1))) {

        std::list< Polynomial_2 > rrchain1;
        std::copy(chain.begin(), CGAL::cpp0x::prev(CGAL::cpp0x::prev(chain.end())), std::back_inserter(rrchain1));
        rrchain1.push_back(h);
        Solution_set next1 = solution_set(x, rrchain1, local);

        std::list< Polynomial_2 > rrchain2;
        std::copy(chain.begin(), CGAL::cpp0x::prev(CGAL::cpp0x::prev(chain.end())), std::back_inserter(rrchain1));
        rrchain2.push_back(*pred/h);
        rrchain2.push_back(*last/h);

        Solution_set next2 = solution_set(x, rrchain2, local);

        // union
        
        std::set_union(next1.begin(), next1.end(), next2.begin(), next2.end(),
                       std::insert_iterator< Solution_set >(next, next.begin()));

      } else {
        
        Pure_bi_solve_2 ibs = bisolve(*pred, *last);
        if (!failed_x) {
          ibs.set_resultant_in_x(gres_x);
        }
        if (!failed_y) {
          ibs.set_resultant_in_y(gres_y);
        }         

        Solution_set diff_roots;
        
        if (0 != CGAL::degree(ibs.g())) {
          // only if chance for further increase
          
          // std::cout << "BDV ibsf: " << ibs.f() << std::endl;
          // std::cout << "BDV ibsg: " << ibs.g() << std::endl;
          
          Bisolve_telemetry_code(t_bdv_bs.start();)
          ibs(x, std::insert_iterator< Solution_set >(diff_roots, diff_roots.begin()), local);
          Bisolve_telemetry_code(t_bdv_bs.stop();)
          //std::cout << "BDV diff_roots.size: " << diff_roots.size() << std::endl;
          
        }
        
        Solution_set prev = solution_set(x, rchain, local);
        
        //std::cout << "#prev: " << prev.size() << std::endl;
        for (typename Solution_set::const_iterator it = prev.begin();
             it != prev.end(); it++) {
          //std::cout << "prevEl: " << *it << std::endl;
        }
        for (typename Solution_set::const_iterator it = diff_roots.begin();
             it != diff_roots.end(); it++) {
          //std::cout << "thisEl: " << *it << std::endl;
        }
        
        std::set_intersection(prev.begin(), prev.end(),
                              diff_roots.begin(), diff_roots.end(),
                              std::insert_iterator< Solution_set >(next, next.begin()));
        
      }
      
      
      //std::cout << "#next: " << next.size() << std::endl;
      //std::cout << "#sets: " << sit->second.size() << std::endl;
      sit = it->second.insert(sit, std::make_pair(chain, next));

    }
    
    //std::cout << "set at x = " << CGAL::to_double(x) << " with chain = " << "TODO" << " has size = " << sit->second.size() << std::endl;
    
    // finally return set
    return sit->second;
    
  }

  Pure_bi_solve_2 bisolve(const Polynomial_2& f, const Polynomial_2& g) const {

    std::pair< Polynomial_2, Polynomial_2 > ppair(f,g);

    typename Bisolves_map::iterator it = this->ptr()->_m_bisolves.find(ppair);

    if (it == this->ptr()->_m_bisolves.end()) {
      Bisolve_telemetry_code(t_bdv_construct.start();)
      it = this->ptr()->_m_bisolves.insert(it, std::make_pair(ppair, Pure_bi_solve_2(f,g)));
      Bisolve_telemetry_code(t_bdv_construct.stop();)
    }
    
    return it->second;
    
  }

}; // Bi_diff_vanish_2
 
} // namespace CGAL

#endif // CGAL_BI_DIFF_VANISH_2_H
// EOF
