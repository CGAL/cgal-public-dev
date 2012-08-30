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
// Author(s): Eric Berberich <eric.berberich@cgal.org>
//            Pavel Emeliyanenko <asm@mpi-inf.mpg.de>

#ifndef CGAL_ALGEBRAIC_KERNEL_2_BI_SLICE_CERTIFY_H 
#define CGAL_ALGEBRAIC_KERNEL_2_BI_SLICE_CERTIFY_H 

/*! \file
 * The header file for the Bi_slice_certify class.
 */


#include <CGAL/config.h>
#include <CGAL/Algebraic_kernel_2/Bi_solve_2_flags.h>

#include <CGAL/Handle_with_policy.h>

#include <complex>

#include <iostream>

#include <boost/optional.hpp>
#include <CGAL/algorithm.h>

#include <CGAL/Polynomial_traits_d.h>

#if CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS
#include <CGAL/Algebraic_kernel_2/Range_analysis_2.h>
#endif

#include <CGAL/Algebraic_kernel_d/Interval_evaluate_1.h>
#include <CGAL/Algebraic_kernel_d/Interval_evaluate_2.h>
#include <CGAL/Algebraic_kernel_d/Float_traits.h>

#include <CGAL/Algebraic_kernel_2/Root_box_1.h>

#include <CGAL/Timer.h>

namespace CGAL {

namespace internal {

//! type of univariate solution
template < class AlgebraicReal_1, class Bound_, class MultiplicityType >
struct Solution_templ_1_rep {
    
  //! this instance's first template parameter
  typedef AlgebraicReal_1 Algebraic_real_1;
  
  //! this instance's second template parameter
  typedef Bound_ Bound;
  
  //! this instance's third template parameter
  typedef MultiplicityType Multiplicity_type;
  
  //! information on seperation
  typedef Root_box_1< Bound > Root_box;
  
  //! default constructor
  Solution_templ_1_rep() :
    multiplicity(0), 
    // for Teissier
    multiplicity_res_fx_fy(0),
    minprec(-1) {
    }
  
  //! construction from value and multiplicity
  Solution_templ_1_rep(const Algebraic_real_1& v, const Multiplicity_type& m) :
    value(v), 
    multiplicity(m), 
    // for Teissier  
    multiplicity_res_fx_fy(0),
    minprec(-1) {   
  } 

  //! increase precsion to new value
  void increase_minprec(long mp) const {
    if (mp > minprec) {
      minprec = mp;
    }
  }
  
  //! double precision
  void double_minprec() const {
    minprec *= 2;
  }

  //! the coordinate
  Algebraic_real_1 value;
  
  //! the multiplicity of the coordinate as root of res
  Multiplicity_type multiplicity;

  // for Teissier
  //! the multiplicity of the coordinate as root of res
  Multiplicity_type multiplicity_res_fx_fy;

  //! the minimal required precision for approximating
  mutable long minprec;
  
  //! information on seperation
  Root_box rb;
  
};


//! type of univariate solution
template < class AlgebraicReal_1, class Bound_, class MultiplicityType >
class Solution_templ_1 : 
 public CGAL::Handle_with_policy< Solution_templ_1_rep< AlgebraicReal_1, Bound_, MultiplicityType > > {
  
 public:  
  //! this instance's first template parameter
  typedef AlgebraicReal_1 Algebraic_real_1;
  
  //! this instance's second template parameter
  typedef Bound_ Bound;
  
  //! this instance's third template parameter
  typedef MultiplicityType Multiplicity_type;

  //! Rep class
  typedef Solution_templ_1_rep< Algebraic_real_1, Bound, Multiplicity_type > Rep;

  //! Base
  typedef CGAL::Handle_with_policy< Rep > Base;
  
  //! information on seperation
  typedef typename Rep::Root_box Root_box;
  
  //! default constructor
  Solution_templ_1() :
    Base() {
  }
  
  //! construction from value and multiplicity
  Solution_templ_1(const Algebraic_real_1& v, const Multiplicity_type& m) :
    Base(Rep(v,m)) {   
  } 

  Algebraic_real_1 value() const {
    return this->ptr()->value;
  }

  long minprec() const {
    return this->ptr()->minprec;
  }

  Multiplicity_type multiplicity() const {
    return this->ptr()->multiplicity;
  }

  // for Tessier
  Multiplicity_type multiplicity_res_fx_fy() const {
    return this->ptr()->multiplicity_res_fx_fy;
  }

  Root_box& rb() {
    return this->ptr()->rb;
  }

  Root_box rb() const {
    return this->ptr()->rb;
  }

  //! increase precsion to new value
  void increase_minprec(long mp) const {
    this->ptr()->increase_minprec(mp);
  }

  //! increase precsion to new value
  void double_minprec() const {
    this->ptr()->double_minprec();
  }

  //! output operator
  void write(std::ostream& os) const {
    os << "  -- value     : " << CGAL::to_double(this->ptr()->value) << std::endl;
    os << "  -- mult      : " << this->ptr()->multiplicity << std::endl;
    // for Teissier
    os << "  -- mult_fx_fy: " << this->ptr()->multiplicity_res_fx_fy << std::endl;
    os << "  -- minprec   : " << this->ptr()->minprec << std::endl;
    os << "  -- rb        : " << std::endl;
    os << this->ptr()->rb;
  } 
  
 
};


//!\brief Box to contain a root of a f(x,y)=g(x,y)=0 system
template < class CertifierTraits >
class Bi_slice_certify {
  
public:

  //! this instance's first template parameter
  typedef CertifierTraits Certifier_traits;
  
  //! type of solution  
  typedef typename Certifier_traits::Algebraic_real_2 Algebraic_real_2;
  
  //! type of AK
  typedef typename Certifier_traits::Algebraic_kernel_d_1 Algebraic_kernel_d_1;
  
  //! type of univariate polynomial
  typedef typename Certifier_traits::Polynomial_1 Polynomial_1;
  
  //! type of bivariate polynomial
  typedef typename Certifier_traits::Polynomial_2 Polynomial_2;
  
  //! type of algebraic real
  typedef typename Certifier_traits::Algebraic_real_1 Algebraic_real_1;
  
  //! type of Coefficient
  typedef typename Algebraic_kernel_d_1::Coefficient Coefficient;
  
  //! type of Bound
  typedef typename Certifier_traits::Bound Bound;
  
  //! isolating box of a root
  typedef typename Certifier_traits::Root_box Root_box;
  
  //! type of Multiplicity
  typedef typename Algebraic_kernel_d_1::Multiplicity_type Multiplicity_type;
  
  typedef Solution_templ_1< Algebraic_real_1, Bound, Multiplicity_type > Solution_1;
    
  //! functor to approximate a solution to certain precision, respecting a
  //! set minimal precision
  class Approximate_solution_1 {
    public:

      //! approximates a solution to given precision, respecting the minimal
      //! precision required
      std::pair< Bound, Bound > operator()(const Solution_1& sol, long prec) {
	    typename Algebraic_kernel_d_1::Approximate_absolute_1 approx_abs;
	
	    return approx_abs(sol.value(), std::max(prec, sol.minprec()));
      }
  };

  //! arithmetic kernel  
  typedef typename CGAL::Get_arithmetic_kernel< Bound >::Arithmetic_kernel AK;
  
  //! bigfloat interval 
  typedef typename AK::Bigfloat_interval BFI;
  
  //! our lovely bigfloats
  typedef typename CGAL::Bigfloat_interval_traits< BFI >::Bound BigFloat;  
  
  //! a vector of solutions
  typedef std::vector< Solution_1 > Solution_1_container;
  
  //! lists indices of solutions along an x- or y-fiber
  typedef std::vector< int > Index_list;
  
  //! maps index of y-fiber to set of x-solutions associated with it
  //! (or vise a versa)
  typedef std::map< int, Index_list > Fiber_list;
  
#if !CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS
  //! active intervals set
  typedef typename Certifier_traits::Active_intervals_set Active_intervals_set;
  
  //! iterator to set element
  typedef typename Certifier_traits::Active_intervals_set_iterator Active_intervals_set_iterator;
  
  //! iterator to set element (const version)
  typedef typename Certifier_traits::Active_intervals_set_const_iterator Active_intervals_set_const_iterator;
#endif  

private:
  
  //! additional properties to be attached to the candidate
  typedef typename Certifier_traits::Candidate_props Props_base;
  
private:

#if CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS
  typedef CGAL::Range_analysis_2< Polynomial_2, BFI, false > IA_real_bfi;
#endif

  //! private struct to store overlapping intervals
  struct Candidate_props : public Props_base {

    //! default constructors
    Candidate_props() :
        Props_base(), 
        vert_line_f(false),
        vert_line_g(false),
#if !CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS
        nf(0), ng(0),
#endif
        prec(2) {
    }

    int v_index;    // index of a solution along a fiber

    bool vert_line_f; // f has a vertical line at candidate

    bool vert_line_g; // g has a vertical line at candidate

    //! \c fb - the head of \c nf f-nodes associated with this box
#if !CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS
    mutable boost::optional< Active_intervals_set_iterator > fb;
    mutable int nf;

    //! \c gb - the head of \c ng g-nodes associated with this box
    mutable boost::optional< Active_intervals_set_iterator > gb;
    mutable int ng;

    //! has overlap?
    mutable boost::optional< bool > _m_has_overlap;

#endif

    //! store current prec for x- and y-approx
    mutable long prec;

    //! timer to later control number of bitstream subdivisions
    mutable CGAL::Timer t_norm;


    // TODO move to traits???
#if !CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS
    void refine_ais(Certifier_traits* traits, 
                    Active_intervals_set& f_ai_set, 
                    Active_intervals_set& g_ai_set) const {
            
      Bisolve_telemetry_code(t_ai_main_subdiv.start();)
      

      if (!vert_line_f) {
        Active_intervals_set_iterator ni, cb, ce;
        int n = 0;
        int i = 0;
        ni = *fb;
        for (; i < nf; i++, ni = ce) {
          if (!traits->is_real(f_ai_set, ni)) {
            continue;
          }
          std::pair< Active_intervals_set_iterator, Active_intervals_set_iterator > cbe = 
            traits->refine(f_ai_set, ni); 
          cb = cbe.first;
          ce = cbe.second;
          n += std::distance(cb, ce);
          if (i == 0) {
            fb = cb;
          }
        }
        nf = n;
      }
      
      if (!vert_line_g) {
        Active_intervals_set_iterator ni, cb, ce;
        int n = 0;
        int i = 0;
        ni = *gb;
        for (; i < ng; i++, ni = ce) {
          if (!traits->is_real(g_ai_set, ni)) {
            continue;
          }
          std::pair< Active_intervals_set_iterator, Active_intervals_set_iterator > cbe = 
            traits->refine(g_ai_set, ni); 
          cb = cbe.first;
          ce = cbe.second;
          n += std::distance(cb, ce);
          if (i == 0) {
            gb = cb;
          }
        }
        ng = n;
      }
      Bisolve_telemetry_code(t_ai_main_subdiv.stop();)
        
      // need to recompute overlaps
      _m_has_overlap = boost::none;
  
    }
    
    // TODO move to traits???
    bool has_overlap(Certifier_traits* traits, 
                     Active_intervals_set& f_ai_set, 
                     Active_intervals_set& g_ai_set) {

      if (!_m_has_overlap) {
        
        // store and only recompute if refine_ai has been called
        
        int n_overlaps = 0;
        
        if (this->vert_line_f) {
          
          n_overlaps = this->ng;
          
        } else if (this->vert_line_g) {
          
          n_overlaps = this->nf;
          
        } else {
          
          Active_intervals_set_iterator fi = *(this->fb), gi = *(this->gb);
          
          int nf = 0, ng = 0;
          
          while (nf < this->nf && ng != this->ng) {
            
            Bound fl = traits->lower(f_ai_set, fi), fh = traits->upper(f_ai_set, fi),
                  gl = traits->lower(g_ai_set, gi), gh = traits->upper(g_ai_set, gi);
            
            //std::cout << "fl: " << dbl(fl) << std::endl;
            //std::cout << "fh: " << dbl(fh) << std::endl;
            //std::cout << "gl: " << dbl(gl) << std::endl;
            //std::cout << "gh: " << dbl(gh) << std::endl;
            
            if (fl <= gh && gl <= fh) { // intervals overlap
              n_overlaps++;
              
              if (fh < gh) {
                fi++, nf++;
              } else {
                gi++, ng++;
              }
            } else {
              if (gl > fh) {
                fi++, nf++;
              } else {
                gi++, ng++;
              }
            }
          }
        }
      
        _m_has_overlap = (n_overlaps > 0);
  
      }
      
      return *_m_has_overlap;  
      
    }
#endif

  };

public:

  struct Algebraic_real_2_x_less_than {
    bool operator()(const Algebraic_real_2& s1, const Algebraic_real_2& s2) const {

      return (s1.compare_x(s2) == CGAL::SMALLER);
    }
  };

  struct Algebraic_real_2_y_less_than {
    bool operator()(const Algebraic_real_2& s1, const Algebraic_real_2& s2) const {

      return (s1.compare_y(s2) == CGAL::SMALLER);
    }
  };

  struct Algebraic_real_2_xy_less_than {
    bool operator()(const Algebraic_real_2& s1, const Algebraic_real_2& s2) const {

      return (s1.compare_xy(s2) == CGAL::SMALLER);
    }
  };

  
  struct Fiber_base {

#if !CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS
    //! active intervals
    Active_intervals_set f_ai_set;

    //! active intervals
    Active_intervals_set g_ai_set;
#endif

    //! vertical line for f
    bool vert_line_f;

    //! vertical line for g
    bool vert_line_g;

    //! f_end
    typename Polynomial_2::const_iterator f_end;
   
    //! g_end
    typename Polynomial_2::const_iterator g_end;

  };

  template< typename LessThan >
  struct Fiber : 
    public Fiber_base {
    
    typedef LessThan Less_than;

    Fiber() {
    };

    //! roots along fiber
    std::set< Algebraic_real_2, Less_than > roots;

  };


  enum {
     CGAL_BI_SLICE_UNDECIDED = 0,   // a candidate is yet undecided
     CGAL_BI_SLICE_EXIN_PREDICATES, // requires theoretical ex- and in-clusion predicate
     CGAL_BI_SLICE_REJECTED,        // a candidate being rejected
     CGAL_BI_SLICE_CERTIFIED        // a candidate being certified
  };

  
public:
  
  //!\name Constructors
  //!@{
  
  //! standard constructor for two polynomials
  Bi_slice_certify(Certifier_traits *traits) :
    _m_traits(traits),
    _m_combinatorial(true), 
    _m_dry_run(false), _m_reverse_dir(false) {
    CGAL_precondition(_m_traits != 0);
  }

   //!@}

public:
  //!\name setter / getters

  //! disabling combinatorial tests
  void set_combinatorial(bool combinatorial) {
    _m_combinatorial = combinatorial;
  }

  //! "dry run": certifier adds undecided candidates to \c fiber_list and
  //! returns without running the norm test (possibly with several tree 
  //! subdivisions)
  void set_dry_run(bool dry_run) {
      _m_dry_run = dry_run;
  }

  //! reverse_dir = true: certification proceeds along x-fiber
  //! otherwise along y-fiber
  void set_reverse_dir(bool reverse_dir) {
      
      if(_m_reverse_dir != reverse_dir) {  
          // swap polynomials & resultants
          _m_traits->swap_coordinates();
      }
      _m_reverse_dir = reverse_dir;
  }

  //! fills up a list of 1D solutions and respective root boxes along x- and
  //! y-directions separately. These lists are used to address candidates by
  //! x/y-indices
  template < typename SolutionXInputIterator, typename SolutionYInputIterator >
  void fill_xy_solution_lists(SolutionXInputIterator x_begin, SolutionXInputIterator x_beyond,
                              SolutionYInputIterator y_begin, SolutionYInputIterator y_beyond) {
    for (SolutionXInputIterator xit = x_begin; xit != x_beyond; xit++) {
      _m_sls_x.push_back(xit->first); // key,value-pair of fiber map
    }
    for (SolutionYInputIterator yit = y_begin; yit != y_beyond; yit++) {
      _m_sls_y.push_back(yit->first); // key,value-pair of fiber map
    }
  }

  //! sets up sets for active intervals
  void set_fiber(Fiber_base* fiber) {
    _m_fiber = fiber;
  }

  //! access the list of active y-fibers
  Fiber_list& fiber_list() {
      return _m_fiber_list;
  }


protected:

  //!\name Helpers
  //!@{
  
  bool _has_vertical_line_at_alpha(const Polynomial_2& p, 
                                   const Algebraic_real_1& alpha,
                                   typename Polynomial_2::const_iterator& p_end) {

    Bisolve_telemetry_code(t_vl_detect.start();)
    
    typename Polynomial_2::const_iterator prev;
    // ignore vanishing leading coefficients
    for (p_end = p.end(); p_end != p.begin(); p_end--) {
      prev = CGAL::cpp0x::prev(p_end);
      // Remark: sign_at is faster than iterative gcd of 
      //         prevs and Construct_polynomial()(alpha)
      //         esp. for alpha without vertical line
      //         for alpha WITH vertical line, penatly is
      //         still less than loss for non-vertical alpha
      CGAL_assertion(_m_traits != 0);
      CGAL_assertion(_m_traits->kernel() != 0);
      if (_m_traits->kernel()->sign_at_1_object()(*prev, alpha) != CGAL::ZERO) {
          break;
      }
    }
    Bisolve_telemetry_code(t_vl_detect.stop();)

    return (p_end == p.begin());

  }

  //!@}
    

  //!\name Certification
  //!@{

public:

  //! certifies (or rejects) candidates of an x- or y-fiber
  //! \c f_index - an index of a fiber along which candidates are to be
  //! processed (in x- or y- depending on whether \c _m_reverse_dir is set)
  //! \c [i_begin;i_end] - a set of candidate indices to process along this
  //! fiber
  //! \c upper_bound_log2_abs - defining the upper bound on the roots 
  //! \c oi - an output iterator of \c Algebraic_real_2
  template< class InputIterator, class OutputIterator >
  OutputIterator operator()(int f_index,
                            InputIterator i_begin,
                            InputIterator i_beyond,
                            long upper_bound_log2_abs,
                            OutputIterator oi) {

    const Solution_1_container& sls = (!_m_reverse_dir ? _m_sls_y : _m_sls_x);
    _m_sol = (!_m_reverse_dir ? _m_sls_x[f_index] : _m_sls_y[f_index]);
    
    _m_traits->init_fiber(_m_sol.rb(), _m_sol.multiplicity());
    
    Algebraic_real_1 alpha = _m_sol.value();
    
    int num_total = std::distance(i_begin, i_beyond);
    int num_decided = 0, num_certified = 0;

    Bisolve_out("\n\nCertify " << num_total << " candidates at "
               << dbl(alpha) << std::endl);

    const Polynomial_2& f = _m_traits->f(),  g = _m_traits->g();

   // Check for vertical lines
 
    Bisolve_out("- Prep: Detect vertical lines")

    bool vert_line_f = 
        _m_fiber->vert_line_f = _has_vertical_line_at_alpha(f, alpha, _m_fiber->f_end);
    
    bool vert_line_g = 
        _m_fiber->vert_line_g = _has_vertical_line_at_alpha(g, alpha, _m_fiber->g_end);
 
    if (vert_line_f) {
      Bisolve_out("\n-- found vertical line as component of f\n");
    }

    if (vert_line_g) {
      Bisolve_out("\n--found vertical line as component of g\n");
    }

    CGAL_assertion(!(vert_line_f && vert_line_g));

    if (i_begin == i_beyond) {
      return oi;
    }

    // Prepare candidates

    // a list of candidates that are yet undecided
    std::list< Candidate_props > candidates;

    InputIterator idx = i_begin;

#if !CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS
    Bisolve_telemetry_code(t_ai_init.start();)

    if (!vert_line_f) {
      _m_fiber->f_ai_set = 
        _m_traits->active_intervals_set_at(alpha, f.begin(), _m_fiber->f_end, upper_bound_log2_abs); 
    }
    if (!vert_line_g) {
      _m_fiber->g_ai_set = 
        _m_traits->active_intervals_set_at(alpha, g.begin(), _m_fiber->g_end, upper_bound_log2_abs); 
    }

    Bisolve_telemetry_code(t_ai_init.stop();)
#endif

#if !CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS

    // Evaluate trees, i.e., all nodes' intervals are contained in candidates;
    // a candidate stores nodes 
    Bisolve_out("- Prep: Refine active intervals wrt to candidates' intervals"
              << std::endl);
    

    Bisolve_telemetry_code(t_ai_prep.start();)
    { // scope for prep-phase
      
      // default # of subdivision rounds
      int j = 0, def_subdivs = 1;
      
      Active_intervals_set_iterator fb = _m_fiber->f_ai_set.begin();
      Active_intervals_set_iterator gb = _m_fiber->g_ai_set.begin();
      Active_intervals_set_iterator fi = fb;
      Active_intervals_set_iterator gi = gb;
      Active_intervals_set_iterator fe = _m_fiber->f_ai_set.end();
      Active_intervals_set_iterator ge = _m_fiber->g_ai_set.end();
      
      // we assume that boxes & tree nodes are sorted by y-coordinate
      while ((vert_line_f || fi != fe) && 
             (vert_line_g || gi != ge) && 
             (idx != i_beyond)) {
      
        if (!vert_line_f && !_m_traits->is_real(_m_fiber->f_ai_set, fi)) {
          fi++;
          continue;
        }
        if (!vert_line_g && !_m_traits->is_real(_m_fiber->g_ai_set, gi)) {
          gi++;
          continue;
        }
        
        Bound fl, fh, gl, gh; 
        
        if (!vert_line_f) {
          fl = _m_traits->lower(_m_fiber->f_ai_set, fi);
          fh = _m_traits->upper(_m_fiber->f_ai_set, fi);
        }
        if (!vert_line_g) {
          gl = _m_traits->lower(_m_fiber->g_ai_set, gi);
          gh = _m_traits->upper(_m_fiber->g_ai_set, gi);
        }
        
        //std::cout << "fl: " << dbl(fl) << std::endl;
        //std::cout << "fh: " << dbl(fh) << std::endl;
        //std::cout << "gl: " << dbl(gl) << std::endl;
        //std::cout << "gh: " << dbl(gh) << std::endl;
        
        CGAL_assertion(vert_line_f || fl <= fh);
        CGAL_assertion(vert_line_g || gl <= gh);
        
        const Root_box& rbox = sls[*idx].rb();
        Bound yl = rbox.left(), yh = rbox.right();
        
        //std::cout << "yl: " << dbl(yl) << std::endl;
        //std::cout << "yh: " << dbl(yh) << std::endl;
        
        // boxes are sorted w.r.t. y-coordinate: take the next interval
        if (!vert_line_f && fh < yl) {
          
          Bisolve_out("-- Skipping f-interval to the left: fh=" << dbl(fh) << ", yl=" << dbl(yl) << "\n");
          fi++;
          continue;
          
        } else if (!vert_line_g && gh < yl) {
          
          Bisolve_out("-- Skipping g-interval to the left: gh=" << dbl(gh) << ", yl=" << dbl(yl) << "\n");
          gi++;
          continue;
          
        } else if ((!vert_line_f && yh < fl) || (!vert_line_g && yh < gl)) {
          
          Bisolve_out("-- Box " << j << " (" << dbl(rbox.mid())
                     << ") rejected\n");
          num_decided++, j++, idx++;
          continue;
          
        }
        
        // [fl; fh] / [gl; gh] is included in [yl; yh]
        bool f_incl = vert_line_f || (yl <= fl && fh <= yh);
        bool g_incl = vert_line_g || (yl <= gl && gh <= yh);
        
        //std::cout << "f_incl: " << f_incl << std::endl;
        //std::cout << "g_incl: " << g_incl << std::endl;
        
        if (f_incl && g_incl && 
            (((vert_line_f || vert_line_g) || (fl <= gh && gl <= fh))) /* overlap if there is no vert line*/ ) {
          
          Candidate_props ps;
          
          if (!vert_line_f) {
            if (!ps.fb) {
              ps.fb = fi;
            }
            fi++;
            ps.nf += 1;
          } else {
            ps.nf = 0;
          }
          if (!vert_line_g) {
            if (!ps.gb) {
              ps.gb = gi;
            }
            gi++;
            ps.ng += 1;
          } else {
            ps.ng = 0;
          }
          
          ps.v_index = *idx;
          
          ps.vert_line_f = vert_line_f;
          ps.vert_line_g = vert_line_g;
          
          candidates.push_back(ps);
          
          Bisolve_out("- Box " << j << " (" << dbl(rbox.mid())
                     << ") is a candidate\n");
          
          idx++, j++;
          
        } else {
          
          // at this point we know that there is an overlap btw [yl; yh]
          // and [fl; fh], [gl; gh]
          // check if tree intervals are disjoint -> skip the one to the left
          if (!vert_line_f && !vert_line_g) {
            if (fh < gl) {
              fi++; 
              continue;
            } else if (gh < fl) {
              gi++;
              continue;
            }
          } 
          
          Bisolve_telemetry_code(t_ai_prep_subdiv.start();)
          Active_intervals_set_iterator cb, ce;
          if (!f_incl) {// list iterators don't get corrupt after insertion..
            int cc = def_subdivs;
            //std::cout << "f-subdiv" << std::endl;
            while (cc--) {
              std::pair< Active_intervals_set_iterator, Active_intervals_set_iterator > cbe = 
                _m_traits->refine(_m_fiber->f_ai_set, fi);
              cb = cbe.first;
              ce = cbe.second;
              fi = cb; // go over this interval once again
              if (cb == ce)
                break;
            }
          }
          
          if (!g_incl) {
            int cc = def_subdivs;
            //std::cout << "g-subdiv" << std::endl;
            while (cc--) {
              std::pair< Active_intervals_set_iterator, Active_intervals_set_iterator > cbe = 
                _m_traits->refine(_m_fiber->g_ai_set, gi);
              cb = cbe.first;
              ce = cbe.second;
              gi = cb; // go over this interval once again
              if (cb == ce)
                break;
            }
          }
          
          Bisolve_telemetry_code(t_ai_prep_subdiv.stop();)
        }

        //std::cout << "iend: " << (idx == i_beyond) << std::endl;
        //std::cout << "fend: " << (fi == fe) << std::endl;
        //std::cout << "gend: " << (gi == fe) << std::endl;

      }
      
      // all excess boxes are automatically rejected...
      int n_rejected = std::distance(idx, i_beyond);
      num_decided += n_rejected;
      Bisolve_out("-- The remaining " << n_rejected << " boxes are rejected\n");
      
    } // end scope prep-phase
    Bisolve_telemetry_code(t_ai_prep.stop();)
    
#else

   for (; idx != i_beyond; idx++) {
      // fill candidates
      Candidate_props ps;
      ps.v_index = *idx;
      candidates.push_back(ps);
   }

#endif // !CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS


    // MAINLOOP decide each candidate

   long round = 0;

#if CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS
   // only needed for IA-rejection // TODO add flag
   long ia_prec = 100;
#endif
   
    // THIS is the main loop!
    Bisolve_telemetry_code(t_main_loop.start();)

    while (num_decided < num_total) {
   
      Bisolve_out("\n- Mainloop - #decided:   " << num_decided << " out of " << num_total)
      Bisolve_out("\n- Mainloop - #certified: " << num_certified << "\n\n")

      typename std::list< Candidate_props >::iterator ips;
      
      int num_with_overlap = num_certified;
#if CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS


      // only needed for IA-rejection TODO: add flag
      Bisolve_out("######## setting ia_prec to: " << ia_prec << "\n");
      
      long oldp = CGAL::set_precision(BFI(), ia_prec);
      // precision change, no normalization
      _m_ia_real_f.set_polynomial(_m_traits->f(), true, true);
      _m_ia_real_g.set_polynomial(_m_traits->g(), true, true);
      CGAL::set_precision(BFI(), oldp); // precision only matters during 
      // conversion
      ia_prec = ia_prec * 2;

#else

      // find number of overlaps
      for (ips = candidates.begin(); ips != candidates.end(); ips++) {
        if (ips->has_overlap(_m_traits, _m_fiber->f_ai_set, _m_fiber->g_ai_set)) {
          num_with_overlap++;
        } 
      }
      
#endif

#if CGAL_BISOLVE_ENABLE_ARCAVOID 
      // numerical certification
      bool numerically_certified_single_overlap = false;
      Bisolve_telemetry_code(c_arca_cert_total++;)
      Bisolve_telemetry_code(t_arca_cert.start();)
      std::list< std::pair< Active_intervals_set_const_iterator, Active_intervals_set_const_iterator > >  overlaps;
      typename Active_intervals_set::Overlapping_regions overlapping_regions;
      overlapping_regions (_m_fiber->f_ai_set, _m_fiber->g_ai_set, std::back_inserter (overlaps));
      Bisolve_out("-- Numerical certifier found " << overlaps.size() << " overlaps." << std::endl)

      if (overlaps.size() == 1) {
        // single overlap ...
        if (overlaps.front().first->touch_real() && // the first should be sufficient
            overlaps.front().second->touch_real()) {
          // .. that is touching the real axis
          Bisolve_telemetry_code(c_arca_cert_success++;)
          numerically_certified_single_overlap = true;
        }
      }
      Bisolve_telemetry_code(t_arca_cert.stop();)
#endif

      Bisolve_out("-- #Candidates: " << candidates.size() << std::endl)


      // try to decide candidates in bottom-up (left-right, if reversed) fashio:
      for (ips = candidates.begin(); ips != candidates.end(); /* ips++ is at TWO positions in LOOP */) {
        
        Bisolve_out("-- Try to decide candidate #" << std::distance(candidates.begin(), ips) << std::endl)

        int flag = CGAL_BI_SLICE_UNDECIDED;

        //std::cout << "flag1: "<< flag << std::endl;

#if CGAL_BISOLVE_ENABLE_ARCAVOID 
        if (numerically_certified_single_overlap) {
          if (ips->has_overlap(_m_traits, _m_fiber->f_ai_set, _m_fiber->g_ai_set)) {
            Bisolve_out("-- Candidate certified: Numerical single complex overlap" << std::endl)
            flag = CGAL_BI_SLICE_CERTIFIED;
          } else {
            Bisolve_out("-- Candidate rejected: Numerical single complex overlap" << std::endl)
            flag = CGAL_BI_SLICE_REJECTED;
          }
        }
#endif

        //std::cout << "flag2: "<< flag << std::endl;

        // counting arguments
        if (flag == CGAL_BI_SLICE_UNDECIDED) {
          flag = _counting_decision(*ips, num_total, num_decided, num_certified);
        }
        
        //std::cout << "flag3: "<< flag << std::endl;

#if !CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS
        
        // check for overlaps of f- and g-intervals
        if (flag == CGAL_BI_SLICE_UNDECIDED && !ips->has_overlap(_m_traits, _m_fiber->f_ai_set, _m_fiber->g_ai_set)) {
          flag = CGAL_BI_SLICE_REJECTED;
          Bisolve_out("-- Candidate rejected: No f-g-overlap" << std::endl)
        }

        // combinatorial decision
        if (flag == CGAL_BI_SLICE_UNDECIDED && _m_combinatorial) {
          
          Bisolve_telemetry_code(t_combinatorial.start();)
          flag = _combinatorial_decision(*ips, num_total, num_decided,
                                         num_certified, num_with_overlap);
          Bisolve_telemetry_code(t_combinatorial.stop();)
            
        }
#endif

        //std::cout << "flag4: "<< flag << std::endl;

        if (flag != CGAL_BI_SLICE_REJECTED && flag != CGAL_BI_SLICE_CERTIFIED) {
          // run the final norm test only if there are not too many candidates are
          // on the slice
          if (_m_dry_run || num_with_overlap > _m_sol.multiplicity()) {
            if (!_m_dry_run) {
              Bisolve_out("-- Skip norm test, as too many active candidates left"
                          << std::endl)
                }
            flag = CGAL_BI_SLICE_UNDECIDED;
          } else {
            Bisolve_telemetry_code(t_exin.start();)
            // inclusion predicate, which is skiped during first round, i.e., make undecided
            flag = (!round ? CGAL_BI_SLICE_UNDECIDED : _exin_decision(*ips));
            Bisolve_telemetry_code(t_exin.stop();)
          }
        }
        
        //std::cout << "flag5: "<< flag << std::endl;

        if (flag == CGAL_BI_SLICE_UNDECIDED) { 
          if (_m_dry_run) {
          // prepare for reverse run
            _add_to_fiber(f_index, ips->v_index);
          }
          ips++;
          continue; // nope.. still undecided.. NEXT candidate
        }
        
        //std::cout << "flagD: "<< flag << std::endl;

        // reaching here, the candidate is decided!

        CGAL_precondition(flag == CGAL_BI_SLICE_REJECTED ||
                          flag == CGAL_BI_SLICE_CERTIFIED);
        
        // hurra!! decided!!
        num_decided++;
        
        if (flag == CGAL_BI_SLICE_CERTIFIED) {
          // hiphiphurra - even certified
          *oi++ = _construct(ips->v_index);
          num_certified++;
        }
        
        ips = candidates.erase(ips); // points to an element immediately after the one being deleted
      }
      
      // at this point we can either run several tree subdivisions or
      // save undecided candidates and return immediately
      if (_m_dry_run) { 
        Bisolve_telemetry_code(t_main_loop.stop();)
        return oi;
      }
          
#if !CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS
      
      // refine intervals contained in candidates, only if not first round
      if (round > 0) {
        for (ips = candidates.begin(); ips != candidates.end(); ips++) {
          
          int count = 0; 
          double ratio = 0.01; // TODO choose a good ratio
          CGAL::Timer t_ai_subdiv;
          
          while ((ips->t_norm.time() == 0 && t_ai_subdiv.time() == 0) || 
                 (t_ai_subdiv.time() < ratio * ips->t_norm.time())) {
            
            ++count;
            
            t_ai_subdiv.start(); // local timer, do not delete
            
            ips->refine_ais(_m_traits, _m_fiber->f_ai_set, _m_fiber->g_ai_set);
            
            t_ai_subdiv.stop(); // local timer, do not delete
            
          }
          
          ips->t_norm.reset();
          
          Bisolve_out("-- AI: refined interval " << count << " times" << std::endl);
          Bisolve_out("-- AI: increasing xy-prec to " << ips->prec << std::endl);
        }
      }
#endif // CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS

      ++round;
      
    } // while() // end of main loop
    
    Bisolve_telemetry_code(t_main_loop.stop();)
    
    if (num_decided != num_total) {
      Bisolve_out("- FATAL: not all " << num_total << " roots decided!!\n")
    }

    CGAL_postcondition(num_decided == num_total);

    Bisolve_out("- All " << num_total << " candidates decided, #certified: " << num_certified )
    
    return oi;
  }

  //!@}
private:
  //!\name Certification (private)
  //!@{
  
  //! constructs a 2D solution
  Algebraic_real_2 _construct(int v_index) const {
    
    const Solution_1& s = (!_m_reverse_dir ? _m_sls_y[v_index] :
                           _m_sls_x[v_index]);
    // TODO construct with solutions
    if(!_m_reverse_dir) { 
      return Algebraic_real_2(_m_sol.value(), _m_sol.multiplicity(),
                                 s.value(), s.multiplicity());
    } else {
      return Algebraic_real_2(s.value(), s.multiplicity(),
                                 _m_sol.value(), _m_sol.multiplicity());
    }
  }
  
  
  //! decides candidates based on counting arguments, i.e., the multiplicity
  //! returns \c CGAL_BI_SLICE_UNDECIDED , or \c CGAL_BI_SLICE_REJECTED ,
  int _counting_decision(const Candidate_props& ps,
                              int n_total, int n_decided, int n_certified) {
    
    // counting argument
    if (_m_sol.multiplicity() == n_certified) {
      Bisolve_out("- Candidate rejected: all mults used " << std::endl)
        return CGAL_BI_SLICE_REJECTED;
    }
    
#if !CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS

    // TODO move to props (has_overlap - i.e. overlap check only triggered if this holds)
    // not root contained test
    if ((!ps.vert_line_f && ps.nf == 0) || (!ps.vert_line_g && ps.ng == 0)) {
        Bisolve_out("- Candidate rejected: no f/g interval" << std::endl)
        return CGAL_BI_SLICE_REJECTED; 
    }

#endif    

    return CGAL_BI_SLICE_UNDECIDED;
    
  }

    //! decides candidates based on combinatorial argument, i.e., the number
  //! of overlaps
  //! returns \c CGAL_BI_SLICE_UNDECIDED or \c CGAL_BI_SLICE_CERTIFIED 
  int _combinatorial_decision(const Candidate_props& ps,
                              int n_total, int n_decided, int n_certified, 
                              int n_with_overlap) {
    
    // the following is not exectuted of numerical WITHOUT cominatorial is used
    if (n_with_overlap == 1 && n_certified == 0 && n_decided == n_total - 1) {
      // only one undecided candidate left (this one!)
      // and no other candidate has been certified
      if (_m_sol.multiplicity() % 2 == 1) {
        // odd multiplicity (there must be a real intersection)
        // there must be an intersection at this candidate
        Bisolve_out("-- Candidate certified: Single real intersection left"
                    << std::endl)
          return CGAL_BI_SLICE_CERTIFIED;
      }
    }
    
    return CGAL_BI_SLICE_UNDECIDED;
  }

  //! decides candidate based on the interval arithmetic rejection and norm test (expensive)
  //! returns \c CGAL_BI_SLICE_UNDECIDED , \c CGAL_BI_SLICE_REJECTED
  //! or \c CGAL_BI_SLICE_CERTIFIED
  // TODO seperate into ex-predicate & in-predicate?
  int _exin_decision(const Candidate_props& ps) {

    Bisolve_telemetry_code(t_approx.start();)
    Approximate_solution_1 approx_sol;
    
    ps.t_norm.reset();
    ps.t_norm.start();

    const Solution_1& s = (!_m_reverse_dir ? 
			   _m_sls_y[ps.v_index] :
			   _m_sls_x[ps.v_index]);
    
    Bisolve_telemetry_code(c_approx++;)
      Bisolve_telemetry_code(g_max_prec = std::max(g_max_prec, ps.prec);)
    
    std::pair< Bound, Bound >
            ax = approx_sol(_m_sol, ps.prec),
            ay = approx_sol(s, ps.prec);
    Bisolve_telemetry_code(t_approx.stop();)

    Bound 
      x0 = (ax.first + ax.second) / Bound(2),
      y0 = (ay.first + ay.second) / Bound(2);

    ps.t_norm.stop();
    
    // try to reject candidate with interval arithmic
    // TODO check whether arcovoid needs exact predicates to be sure (e.g., every 1000th iteration)
#if CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS
    Bisolve_telemetry_code(t_exia2.start();)

    bool rejected = false;
    
#if 0 // 1 for Bigfloat-IA, 0 for Bound-IA // TODO add flag
    
    // this is to get rid of compiler stupidness..
    // TODO check the precision for GMP number types
    long pp = std::max((long)16, ps.prec / 10);
    long oldp = CGAL::set_precision(BFI(), pp);

    // l, r, b, t
    typename IA_real_bfi::Quad q(
            CGAL::lower(CGAL::convert_to_bfi(ax.first)),
            CGAL::upper(CGAL::convert_to_bfi(ax.second)),
            CGAL::lower(CGAL::convert_to_bfi(ay.first)),
            CGAL::upper(CGAL::convert_to_bfi(ay.second)));

    typename IA_real_bfi::NT l, h;
    if (!ps.vert_line_f) {
      // _m_ia_real_f.eval_range_AF1_RD_2(q, l, h);
      _m_ia_real_f.eval_range_AF1_2(_m_ia_real_f.internal_poly_2(), q, l, h);
      rejected = (CGAL::sign(l) == CGAL::sign(h) && CGAL::sign(l) != CGAL::ZERO);
    }
    
    if (!rejected) {
      if (!ps.vert_line_g) {
        //  _m_ia_real_g.eval_range_AF1_RD_2(q, l, h);
        _m_ia_real_g.eval_range_AF1_2(_m_ia_real_g.internal_poly_2(), q, l, h);
        rejected = (CGAL::sign(l) == CGAL::sign(h) && CGAL::sign(l) != CGAL::ZERO);
      }
    }
    
    CGAL::set_precision(BFI(), oldp);

#else
    
    typename internal::Interval_evaluate_2< Polynomial_2, Bound > ieval_2;

    CGAL::cpp0x::array<Bound,4> box = 
      CGAL::make_array(ax.first, ax.second, ay.first, ay.second);
    
    if (!ps.vert_line_f) {
      std::pair< Bound, Bound > f_ai_set = ieval_2(_m_traits->f(), box);
      rejected = 
        (CGAL::sign(f_ai_set.first) == CGAL::sign(f_ai_set.second) && 
         CGAL::sign(f_ai_set.first) != CGAL::ZERO);
    }
    
    if (!rejected) {
      if (!ps.vert_line_g) {
        std::pair< Bound, Bound > g_ai_set = ieval_2(_m_traits->g(), box);
        rejected = 
          (CGAL::sign(g_ai_set.first) == CGAL::sign(g_ai_set.second) && 
           CGAL::sign(g_ai_set.first) != CGAL::ZERO);
      }
    }
    
#endif

    if (rejected) {
      // no zero contained in both intervals
      Bisolve_telemetry_code(
        t_exia2.stop();
        ++c_ieval2_reject;
      )
      Bisolve_out("-- Candidate rejected with IA_2" << std::endl);
      return CGAL_BI_SLICE_REJECTED;
    }

    Bisolve_telemetry_code(t_exia2.stop();)
#endif // CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS

    Bisolve_telemetry_code(t_innorm.start();)
    if (_m_traits->norm_test(ps, x0, y0, s.rb(), s.multiplicity(), ps.prec)) {
      Bisolve_telemetry_code(t_innorm.stop();)
      return CGAL_BI_SLICE_CERTIFIED;
    }
    Bisolve_telemetry_code(t_innorm.stop();)
 
    //! reaching here means that the candidate is not yet decided
    Bisolve_out("-- Candidate still undecided" << std::endl);
    return CGAL_BI_SLICE_UNDECIDED;
  }

    //! adds an x-solution with index \c f_index to y-fiber with index
    //! \c v_index
    void _add_to_fiber(int f_index, int v_index) {
        
        typename Fiber_list::iterator i = _m_fiber_list.find(v_index);
        
        if(i == _m_fiber_list.end()) {
            i = _m_fiber_list.insert(
                std::make_pair(v_index, Index_list())).first;
        }
        i->second.push_back(f_index);
    }
    
    //!@}
private:

    //! traits object
    Certifier_traits *_m_traits; // also provides "_m_ak"

    //! fiber
    Fiber_base *_m_fiber;

#if CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS
    IA_real_bfi _m_ia_real_f;

    IA_real_bfi _m_ia_real_g;
#endif

    //! a solution in x- or y- defining a fiber along which candidates are
    //! processed (defaults to solution in x)
    Solution_1 _m_sol;

    //! a set of solutions along x- and y-directions (coord + multiplicity)
    Solution_1_container _m_sls_x, _m_sls_y;
    
    //! support for combinatorial test
    bool _m_combinatorial;

    //! fiber list: keeps candidate indices to certify in different directions
    Fiber_list _m_fiber_list;

    //! support for bidirectional filtering
    bool _m_dry_run, _m_reverse_dir;

}; // Bi_slice_certifier


//! output operator
 template < class AlgebraicReal_1, class Bound_, class MultiplicityType >
std::ostream& operator<<(std::ostream& os, 
                         const Solution_templ_1< AlgebraicReal_1, Bound_, MultiplicityType > & sol) {
  sol.write(os);
  return os;
} 
 
} // namespace internal

} // namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_2_BI_SLICE_CERTIFY_H
// EOF
