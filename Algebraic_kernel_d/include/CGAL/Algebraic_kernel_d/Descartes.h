// Copyright (c) 2006-2009, 2012 Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     :  Michael Hemmer <hemmer@mpi-inf.mpg.de>
//                  Eric Berberich <eric.berberich@cgal.org>
//
// ============================================================================

/*! \file CGAL/Descartes.h
  \brief Defines class CGAL::Descartes. 
  
  Isolate real roots of polynomials.

  This file provides a class to isolate real roots of polynomials,
  using the algorithm based on the method of Descartes.

  The polynomial has to be a univariate polynomial over any number
  type which is contained in the real numbers.

*/

#ifndef CGAL_ALGEBRAIC_KERNEL_D_DESCARTES_H
#define CGAL_ALGEBRAIC_KERNEL_D_DESCARTES_H

#include <CGAL/config.h>

#include <CGAL/Algebraic_kernel_d/Generic_isolator.h>

#include <CGAL/Algebraic_kernel_d/univariate_polynomial_utils.h>
#include <CGAL/Algebraic_kernel_d/construct_binary.h>

#include <boost/optional.hpp>
#include <boost/none.hpp>

#define POLYNOMIAL_REBIND( coeff ) \
  typename CGAL::Polynomial_traits_d<Polynomial>::template \
  Rebind<coeff,1>::Other::Type

namespace CGAL {

namespace internal {


/*! \brief A model of concept RealRootIsolator. 
 */
template < typename Polynomial_, typename Bound_, typename HandlePolicy > // no default on policy, should be decided in higher level
class Descartes_isolator_rep :
  public Generic_isolator_rep< Polynomial_, Bound_, HandlePolicy, CGAL::Tag_false > {

public:
    
//!\name Tags
//!@{
    
public:
//! no "refine_interval" function implemented
typedef CGAL::Tag_false Refine_interval_provided_tag;

//!@} // Tags


//!\name Public types
//!@{

//! first template parameter
typedef Polynomial_ Polynomial;

//! second template parameter
typedef Bound_ Bound;

//! third template parameter
typedef HandlePolicy Handle_policy;

//! the base class
typedef Generic_isolator_rep< Polynomial, Bound, Handle_policy, Refine_interval_provided_tag > Base;

//! the class itself
typedef Descartes_isolator_rep< Polynomial, Bound, Handle_policy > Self;

//!@} // Public types

  private:

// Integer or Numerator/Denominator type of bound.
typedef typename CGAL::Fraction_traits<Bound>::Numerator_type Integer;

typedef CGAL::Fraction_traits<Polynomial> FT_poly;
typedef Fraction_traits<Bound> FT_rat;


//!\name Constructors
//!@{

#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_DEFAULT_CONSTRUCTIBLE
  protected:
    Descartes_isolator_rep();
#else
#endif

#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_ASSIGNABLE
  protected:
    // needs a cocy constructor
    Descartes_isolator_rep& operator=(const Self&); // = disable 
#else
  public:
    //! copy constructor
    Descartes_isolator_rep(const Self& drep) : 
        Base(drep),    
        _m_left(drep._m_left), 
        _m_scale(drep._m_scale),
        _m_denom(drep._m_denom),
        _m_is_strong(drep._m_is_strong), 
        _m_prec(drep._m_prec),
        _m_interval(drep._m_interval) {
        
        int d = CGAL::degree(this->polynomial());

        _m_numerator = new Integer[d]; 
        _m_denominator_exponent = new Integer[d];
        _m_is_exact = new bool[d]; 
        for(int i = 0; i < this->number_of_real_roots(); i++) {
          _m_numerator[i] = drep._m_numerator[i];
          _m_denominator_exponent[i] = drep._m_denominator_exponent[i];
          _m_is_exact[i] = drep._m_is_exact[i];
        }
   }
   
   //! assignement operator
   Descartes_isolator_rep& operator=(const Self& drep) {
        if (&drep == this) {
           return;
        }
        Base::operator=(drep);
        _m_left = drep._m_left; 
        _m_scale = drep._m_scale;
        _m_denom = drep._m_denom;
        _m_is_strong = drep._m_is_strong; 
        _m_prec = drep._m_prec;
        _m_interval = drep._m_interval;
        
        int d = CGAL::degree(this->polynomial());
        
        _m_numerator = new Integer[d]; 
        _m_denominator_exponent = new Integer[d];
        _m_is_exact = new bool[d]; 
        for(int i = 0; i < this->number_of_real_roots(); i++) {
          _m_numerator[i] = drep._m_numerator[i];
          _m_denominator_exponent[i] = drep._m_denominator_exponent[i];
          _m_is_exact[i] = drep._m_is_exact[i];
        }
   }
#endif

  public:
/*! \brief Constructor from univariate square free polynomial.
  
  The RealRootIsolator provides isolating intervals for the real 
  roots of the polynomial.
  \pre the polynomial is square free      
*/
  Descartes_isolator_rep(const Polynomial& p, 
                         bool is_strong,
                         int prec) : 
    Base(p), 
    _m_is_strong(is_strong), 
    _m_prec(prec),
    _m_interval(boost::none) {      
      // nothing to do beyond initialization
    }

  // constructor with given interval 
  // (experimental)
  Descartes_isolator_rep(const Polynomial& p,
                         const Bound& left,
                         const Bound& right,
                         bool is_strong = false,
                         int prec = 2) : 
    Base(p) , 
    _m_is_strong(is_strong), 
    _m_prec(prec),
    _m_interval(std::make_pair(left, right)) {
      // nothing to do beyond initialization
    }

//!@} // Constructors

//!\name Destructor
//!@{

// destructor
~Descartes_isolator_rep() {
  delete[] _m_numerator;
  delete[] _m_denominator_exponent;
  delete[] _m_is_exact;
}

//!@} // Desctructor

  public:

//!\name Access functions
//!@{

//! isolates the real roots
void isolate() {
  
  if (this->is_isolated()) {
    return;
  }

  int d = CGAL::degree(this->polynomial());
  
  _m_numerator = new Integer[d]; 
  _m_denominator_exponent = new Integer[d];
  _m_is_exact = new bool[d];
  _m_number_of_real_roots = 0;
  if (d == 0) { 
    if (this->polynomial().is_zero()) {
      _m_number_of_real_roots = -1;
    }
    Base::_m_is_isolated = true;
    return;
  }
  
  if (this->_m_interval) {
    typename FT_rat::Decompose decompose;
    typedef typename FT_rat::Numerator_type Numerator;
    typedef typename FT_rat::Denominator_type Denominator;
    
    Numerator numleft, numright;
    Denominator denleft, denright;
    
    decompose(this->_m_interval->first  /* left  boundary */, numleft,  denleft);
    decompose(this->_m_interval->second /* right boundary */, numright, denright);
    
    _m_left = numleft * denright;
    _m_scale = numright * denleft - _m_left;
    _m_denom = denleft * denright;
    typedef typename Polynomial::NT Coefficient;
    Base::_m_polynomial.scale_down(Coefficient(denleft*denright));
  }
  
  intern_decompose(Base::_m_polynomial,typename FT_poly::Is_fraction());
  
  Base::_m_is_isolated = true;

}

//! returns the number of real roots
int number_of_real_roots() const { 
  return _m_number_of_real_roots;
}

/*! \brief returns  \f${l_i}\f$ the left bound of the isolating interval 
  for root  \f$root_{i}\f$.
  
  In case is_exact_root(i) is true,  \f$l_i = root_{i}\f$,\n
  otherwise:  \f$l_i < root_{i}\f$. 
  
  If  \f$i-1>=0\f$, then  \f$l_i > root_{i-1}\f$. \n
  If  \f$i-1>=0\f$, then  \f$l_i >= r_{i-1}\f$, 
  the right bound of  \f$root_{i-1}\f$\n
  
  \pre 0 <= i < number_of_real_roots()
*/
Bound left_bound(size_t i) const { 
  Integer numerator_, denominator_;
  left_bound(i, numerator_, denominator_);
  return Bound(numerator_) / Bound(denominator_);
}

/*! \brief returns  \f${r_i}\f$ the right bound of the isolating interval 
  for root  \f$root_{i}\f$.
  
  In case is_exact_root(i) is true,  \f$r_i = root_{i}\f$,\n
  otherwise:  \f$r_i > root_{i}\f$. 
  
  If  \f$i+1< n \f$, then  \f$r_i < root_{i+1}\f$,
  where \f$n\f$ is number of real roots.\n
  If  \f$i+1< n \f$, then  \f$r_i <= l_{i+1}\f$, 
  the left bound of  \f$root_{i+1}\f$\n
  
  
  \pre 0 <= i < number_of_real_roots()
*/
Bound right_bound(size_t i) const { 
  Integer numerator_, denominator_;
        right_bound(i, numerator_, denominator_);
  return Bound(numerator_) / Bound(denominator_);
}  

/*! \brief returns true if the isolating interval is degenerated to a 
  single point.
  
  If is_exact_root(i) is true, 
  then left_bound(size_t i) equals  \f$root_i\f$. \n
  If is_exact_root(i) is true, 
  then right_bound(size_t i) equals  \f$root_i\f$. \n 
*/
bool is_exact_root(size_t i) const { 
  return _m_is_exact[i];
}

//! \brief Returns whether the \c i th root is definitely a simple root of the isolated polynomial
virtual bool is_certainly_simple_root(size_t /* i */) const {
  return true;
}
  
//!\brief Returns whether the \c i th root is definitely a multiple root of the isolated polynomial
virtual bool is_certainly_multiple_root(size_t /* i */) const {
  return false;
}
    
//! returns an upper bound of the \c i th root multiplicity                      
virtual int upper_bound_for_multiplicity(size_t /* i */) const {
  return 1;
}

//! returns the \c i th root multiplicity                      
virtual int multiplicity_of_root(size_t /* i */) const {
  return 1;
}

//! Refine the <tt>i</tt>th isolating interval
// virtual void refine_interval(size_t /* i */); // TODO 2012 not implemented

//!@} // Access functions

  private:

    void intern_decompose(const Polynomial& p, ::CGAL::Tag_true){
        typename FT_poly::Decompose decompose;
        typename FT_poly::Numerator_type nump;
        typename FT_poly::Denominator_type dummy;
 
        decompose(p, nump, dummy);
        init_with(nump);
    }

    void intern_decompose(const Polynomial& p, ::CGAL::Tag_false){ 
        init_with(p);
    }
  
    template<class Polynomial__> 
    void init_with(const Polynomial__& p){
      typedef typename Polynomial__::NT Coeff;
      if (!_m_interval) {
        _m_left = -weak_upper_root_bound<Coeff>(p);
        _m_scale = - _m_left * Integer(2);
        _m_denom = Integer(1);
      }
      Polynomial__ r = ::CGAL::translate(p, Coeff(_m_left));
      Polynomial__ q = ::CGAL::scale_up(r, Coeff(_m_scale));
      zero_one_descartes<Coeff>(q,0,0);
    }

    void left_bound(size_t i, Integer& numerator_, Integer& denominator_) const {
      CGAL_precondition(i >= 0);
      CGAL_precondition(static_cast<int>(i) < _m_number_of_real_roots);
      construct_binary(_m_denominator_exponent[i], denominator_);
      numerator_= _m_scale * _m_numerator[i] + _m_left * denominator_;
      denominator_ = denominator_ * _m_denom;
    }
      
    void right_bound(size_t i, Integer& numerator_, Integer& denominator_) const {
      CGAL_precondition(i >= 0);
      CGAL_precondition(static_cast<int>(i) < _m_number_of_real_roots);
      if (_m_is_exact[i]){
        return left_bound(i,numerator_,denominator_);
      } else {
        construct_binary(_m_denominator_exponent[i],denominator_);
        numerator_= _m_scale * (_m_numerator[i]+1) + _m_left * denominator_;
        denominator_ = denominator_ * _m_denom;
      }
    }
    
    //! returns the polynomial $(1 + x)^n P(1/(1 + x))$.
    template <class Coeff__>
    /*
    typename 
    CGAL::Polynomial_traits_d<Polynomial> 
    ::template Rebind<Coeff__,1>::Other::Type
    */
    POLYNOMIAL_REBIND(Coeff__) 
        variation_transformation(const POLYNOMIAL_REBIND(Coeff__)& p) { 
        POLYNOMIAL_REBIND(Coeff__) r = reversal(p);
        return translate_by_one(r); 
    }

    //! Returns an upper bound on the absolute value of all roots of $P$.
    /*! The upper bound is a power of two. Only works for univariate
     * polynomials. 
     */
    template <class Coeff__>
    Integer weak_upper_root_bound(const POLYNOMIAL_REBIND(Coeff__)& p) { 
  
        typename Real_embeddable_traits<Coeff__>::Abs abs;
        const int n = CGAL::degree(p);
        Integer r(1);  // return value
        Coeff__ x(1);  // needed to "evaluate" the polynomial
        Coeff__ val;
        for (;;) {
            val = -abs(p[n]);
            for (int i = n-1; i >= 0; i--) {
                val = val*x + abs(p[i]);
            }
            if (val < Coeff__(0)) return r;
            r *= Integer(2);
            x = Coeff__(r);
        }
    }

    //! tests if the polynomial has no root in the interval.
    template <class Coeff__>
    bool not_zero_in_interval(const POLYNOMIAL_REBIND(Coeff__)& p)
    {
        if(CGAL::degree(p) == 0) return true;
        if(internal::sign_variations(variation_transformation<Coeff__>(p)) != 0)
            return false;
        return (p[0] != Coeff__(0) && p.evaluate(Coeff__(1)) != Coeff__(0));
    }

    //! Descartes algoritm to determine isolating intervals for the roots 
    //! lying in the interval (0,1).
    // The parameters $(i,D)$ describe the interval $(i/2^D, (i+1)/2^D)$.
    // Here $0\leq i < 2^D$.
    template <class Coeff__>
    void zero_one_descartes(const POLYNOMIAL_REBIND(Coeff__)& p,
                            Integer i, Integer d) { 
        // Determine the number of sign variations of the transformed 
        // polynomial $(1+x)^nP(1/(1+x))$. This gives the number of 
        // roots of $p$ in $(0,1)$.

        POLYNOMIAL_REBIND(Coeff__) r = variation_transformation<Coeff__>(p);
        int sign_vars = sign_variations(r);
        
        // no root
        if ( sign_vars == 0 ) return;

        // exactly one root
        // Note the termination criterion $P(0)\neq 0$ and $P(1)\neq 0$.
        // This ensures that the given interval is an isolating interval.
        if ( sign_vars == 1 
                && p[0] != Coeff__(0) 
                && p.evaluate(Coeff__(1)) != Coeff__(0) ) { 
            if (_m_is_strong) {
                strong_zero_one_descartes<Coeff__>(p,i,d);
                return;
            }
            else {
                _m_numerator[_m_number_of_real_roots] = i;
                _m_denominator_exponent[_m_number_of_real_roots] = d;
                _m_is_exact[_m_number_of_real_roots] = false;
                _m_number_of_real_roots++;
                return; 
            } 
        }

        // more than one root
        // Refine the interval.
        i = 2*i; 
        d = d+1;

        // Transform the polynomial such that the first half of the interval
        // is mapped to the unit interval.
        POLYNOMIAL_REBIND(Coeff__) q = scale_down(p,Coeff__(2));
 
        // Consider the first half of the interval.
        zero_one_descartes<Coeff__>(q,i,d);
     
        // Test if the polynomial is zero at the midpoint of the interval 
        POLYNOMIAL_REBIND(Coeff__)  s = translate_by_one(q);
        if ( s[0] == Coeff__(0) ) { 
            _m_numerator[_m_number_of_real_roots] = i + 1;
            _m_denominator_exponent[_m_number_of_real_roots] = d;
            _m_is_exact[_m_number_of_real_roots] = true;
            _m_number_of_real_roots++;
        }
         
        // Consider the second half of the interval. 
        zero_one_descartes<Coeff__>(s,i+1,d); 
    }

    //! Strong Descartes algoritm to determine isolating intervals for the 
    //! roots lying in the interval (0,1), where the first
    //! derivative have no sign change. \pre $P$ has only one root in the 
    //! interval given by $(i,d)$.
    // The parameters $(i,d)$ describe the interval $(i/2^d, (i+1)/2^d)$.
    // Here $0\leq i < d$.
    template <class Coeff__>
    void strong_zero_one_descartes(const POLYNOMIAL_REBIND(Coeff__)& p,
                                   Integer i, Integer d) { 

        // Test if the polynomial P' has no roots in the
        // interval. For further use in Newton, the interval should be not
        // too large.

        // test if isolating interval is smaller than epsilon
        // [l,r]  ->  r-l < epsilon
        // l = (r-l) * i/2^d + l
        // r = (r-l) * (i+1)/2^d + l
        // r-l = (r-l) * 1/2^d
        // r-l < epsilon = 2^(-k)
        // <=> (r-l) * 1/2^d < 2^(-k)
        // <=> 2^d > (r-l) / 2^(-k)
        // <=> 2^d > (r-l) * 2^k
        // <=> 2^(d-k) > (r-l)

      POLYNOMIAL_REBIND(Coeff__) pp = CGAL::differentiate(p);
        if (not_zero_in_interval<Coeff__>(pp)) { // P' has no root over current interval
            Integer tmp;
            construct_binary(d-_m_prec, tmp);  // tmp = 2^{d-k}
            if(tmp * _m_denom > _m_scale ) {
                _m_numerator[_m_number_of_real_roots] = i;
                _m_denominator_exponent[_m_number_of_real_roots] = d;
                _m_is_exact[_m_number_of_real_roots] = false;
                _m_number_of_real_roots++;
                return; 
            }
        }

        // either $P'$ fails the test, 
        // or the interval is too large
        // Refine the interval.
        i = 2*i; 
        d = d+1;

        // Transform the polynomial such that the first half of the interval
        // is mapped to the unit interval.
        POLYNOMIAL_REBIND(Coeff__) q = scale_down(p,Coeff__(2));
 
        // Test if the polynomial is zero at the midpoint of the interval 
        POLYNOMIAL_REBIND(Coeff__)  s = translate_by_one(q);
        if (s[0] == Coeff__(0) ) { 
            _m_numerator[_m_number_of_real_roots] = i + 1;
            _m_denominator_exponent[_m_number_of_real_roots] = d;
            _m_is_exact[_m_number_of_real_roots] = true;
            _m_number_of_real_roots++;
            return;
        }

        // Consider the first half of the interval.
        if(sign_variations(variation_transformation<Coeff__>(q)) == 1) {
            strong_zero_one_descartes<Coeff__>(q,i,d);
            return;
        }
         
        // Consider the second half of the interval. 
        strong_zero_one_descartes<Coeff__>(s,i+1,d); 
        return;
    } 

  private:

    //! Needed for the referencing counting mechanism
    virtual CGAL::Reference_counted_hierarchy<>* clone() {
      return new Descartes_isolator_rep(*this);
    }

  protected:

    //!\name Private members
    //!@{

    //! stores the number of real roots
    int _m_number_of_real_roots;   

    //! stores the numerators of the intervals' left boundary
    Integer* _m_numerator;   

    //! stores the denominator of the interval's boundaries
    Integer* _m_denominator_exponent; 

    //! is a root exact?
    bool* _m_is_exact;

    //! parameter for interval creation
    Integer _m_left;

    //! parameter for interval creation
    Integer _m_scale;

    //! parameter for interval creation
    Integer _m_denom;

    //! strong refinement required
    bool _m_is_strong;

    //! require intervals to be smaller than $2^{_m_prec}$
    int _m_prec;

    //! if not boost::none defines x-range
    boost::optional< std::pair< Bound, Bound > > _m_interval;

    //!@}

}; // class Descartes_isolator_rep

// TODO 2012 keep Descartes in internal namespace?

/*\brief Class to isolate real roots of a square-free univariate polynomial
 * using the Descartes method
 */
template< typename Polynomial_, typename Bound_, typename HandlePolicy = CGAL::Handle_policy_no_union >
class Descartes : 
  public internal::Generic_isolator< Polynomial_, Bound_, HandlePolicy, internal::Descartes_isolator_rep< Polynomial_, Bound_, HandlePolicy > > {
  
 public:
  
  //!\name Public types
  //!@{
  
  //! first template parameter 
  typedef Polynomial_ Polynomial;
  
  //! second template parameter
  typedef Bound_ Bound;

  //! third template parameter
  typedef HandlePolicy Handle_policy;
  
  //! the Descartes rep class
  typedef internal::Descartes_isolator_rep< Polynomial, Bound, Handle_policy > Rep;

  //! the base class
  typedef internal::Generic_isolator< Polynomial, Bound, Handle_policy, Rep > Base;       
  
  //! the class itself
  typedef Descartes< Polynomial, Bound, Handle_policy > Self;

  //!@} // Public types

  //!\name Constructors
  //!@{
  
#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_DEFAULT_CONSTRUCTIBLE
 protected:
    Descartes(); // = disable
#else
  // no default constructor is needed, as the 
  // standard constructor has default values on all arguments
#endif
  
#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_ASSIGNABLE
  protected:
    Descartes(const Self&); // = disable
    Descartes& operator=(const Self&); // = disable
#endif

 public:

  //! standard constructor
  Descartes(const Polynomial& p = Polynomial(0),
            bool isolate = true,
            bool is_strong = false,
            int prec = 2) : 
    Base(new Rep(p, is_strong, prec)) {    
      
    if (isolate) {
      this->isolate();
    }
    
  }
  
  //! constructor for x-range
  Descartes(const Polynomial& p,
            const Bound& left,
            const Bound& right,
            bool isolate = true,
            bool is_strong = false,
            int prec = 2) : 
    Base(new Rep(p, left, right, is_strong, prec)) {    

    if (isolate) {
      this->isolate();
    }

  }
    
  //!@} // Constructors

}; // Descartes

} // namespace internal

} //namespace CGAL


#endif // CGAL_ALGEBRAIC_KERNEL_D_DESCARTES_H
// EOF
