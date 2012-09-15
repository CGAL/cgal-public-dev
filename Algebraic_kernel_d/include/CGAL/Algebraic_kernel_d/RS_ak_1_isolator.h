// Copyright (c) 2012 Max-Planck-Institute Saarbruecken (Germany).
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
// $URL: svn+ssh://hemmer@scm.gforge.inria.fr/svn/cgal/trunk/Polynomial/include/CGAL/Polynomial.h $
// $Id: Polynomial.h 47254 2008-12-06 21:18:27Z afabri $
// 
//
// Author(s)     :  Eric Berberich <eric.berberich@cgal.org>

/*! \file RS_ak_1_isolator.h
  \brief Defines class CGAL::internal::RS_ak_1_isolator.h
  
  Isolate real roots of polynomials with Fabrice Roullier's Rs.

  The polynomial has to be a univariate polynomial over any number
  type which is contained in the real numbers.

  The isolator uses internally representations of RS_AK_1's algebraic reals. 
  This way it can provide "refine_interval". The class RS/isolator.h does 
  not offer this function.
*/

#ifndef CGAL_ALGEBRAIC_KERNEL_D_RS_AK_1_ISOLATOR_H
#define CGAL_ALGEBRAIC_KERNEL_D_RS_AK_1_ISOLATOR_H 1

#include <CGAL/config.h>

#if !CGAL_USE_RS3
#error Please add RS3 as external library to include this file.
#endif

#include <CGAL/Algebraic_kernel_d/Generic_isolator.h>

#include <CGAL/Algebraic_kernel_rs_gmpz_d_1.h>

namespace CGAL {

namespace internal {

/*! \brief An RS model of concept RealRootIsolatorRep. 
 */
template < typename Polynomial_, typename Bound_, typename HandlePolicy = CGAL::Handle_policy_no_union > 
class RS_ak_1_isolator_rep : 
  public Generic_isolator_rep< Polynomial_, Bound_, HandlePolicy, CGAL::Tag_true > {


//!\name Tags
//!@{

public:
    //! "refine_interval" function is implemented
    typedef CGAL::Tag_true Refine_interval_provided_tag;

//!@} // Tags

public:

    //!\name Public types
    //!@{

    //! First template parameter
    typedef Polynomial_ Polynomial;

    //! Second template parameter
    typedef Bound_ Bound;

    //! Third template parameter
    typedef HandlePolicy Handle_policy;

    //! base class 
    typedef Generic_isolator_rep< Polynomial, Bound, Handle_policy, Refine_interval_provided_tag > Base;

    //! the class itself
    typedef RS_ak_1_isolator_rep< Polynomial, Bound, Handle_policy > Self;

    //!@} // Public types
  
protected:
    
    //! Coefficient type of polynomial
    typedef typename Polynomial::NT Coefficient;

    //! internal algebraic kernel
    typedef CGAL::Algebraic_kernel_rs_gmpz_d_1 AK_1;

    //! type of algebraic real
    typedef typename AK_1::Algebraic_real_1 Algebraic_real_1;
    
public:
    //!\name Constructors
    //!@{

#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_DEFAULT_CONSTRUCTIBLE
protected:
    RS_ak_1_isolator_rep(); // = disable
#else
#endif
    
#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_ASSIGNABLE
protected:
    RS_ak_1_isolator_rep(const Self&) {} // = default
    RS_ak_1_isolator_rep& operator=(const Self&); // = disable
#endif

public:

    /*! \brief Constructor from univariate square free polynomial.
    
      The RealRootIsolator provides isolating intervals for the real 
      roots of the polynomial.
      \pre the polynomial is square free      
    */
    RS_ak_1_isolator_rep(const Polynomial& p,
                         bool isolate = true) :
      Base(p)
    {
      if (isolate) {
        this->isolate();
      }
    }
    
    // TODO 2012 interval constructor?

    //!@} // Constructors

public:
    //!\name Access functions
    //!@{
    
    //! isolates the roots of the polynomial (if not yet done)
    virtual void isolate() { 
      
      if (this->is_isolated()) {
        return; // nothing to do
      }
      
      typename AK_1::Solve_1()(Base::_m_polynomial, true, std::back_inserter(_m_real_roots));
      if (_m_real_roots.size() > 1) {
        for (std::vector< Algebraic_real_1 >::size_type i = 0;
             i < _m_real_roots.size() - 1;
             ++i) {
          // Isolate them against each other
          typename AK_1::Compare_1()(_m_real_roots[i],_m_real_roots[i+1]);
        }
      }
      
      // TODO implement isolated function

      Base::_m_is_isolated = true;
    }

    //! returns the number of real roots
    int number_of_real_roots() const { 
      return _m_real_roots.size();
    }

    /*! \brief returns true if the isolating interval is degenerated to a 
      single point.
      
      If is_exact_root(i) is true, 
      then left_bound(int i) equals  \f$root_i\f$. \n
      If is_exact_root(i) is true, 
      then right_bound(int i) equals  \f$root_i\f$. \n 
    */
    bool is_exact_root(size_t i) const {
      CGAL_assertion(i >= 0);
      CGAL_assertion(static_cast<int>(i) < this->number_of_real_roots());
      CGAL_precondition(this->is_isolated());
      return _m_real_roots[i].is_point();
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
      CGAL_assertion(i >= 0);
      CGAL_assertion(static_cast<int>(i) < this->number_of_real_roots());
      CGAL_precondition(this->is_isolated());
      return Bound(_m_real_roots[i].left());
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
      CGAL_assertion(i >= 0);
      CGAL_assertion(static_cast<int>(i) < this->number_of_real_roots());
      CGAL_precondition(this->is_isolated());
      return Bound(_m_real_roots[i].right());
    }  
    
    //! \brief Returns whether the \c i th root is definitely a simple root of the isolated polynomial
    virtual bool is_certainly_simple_root(size_t /* i */) const {
      return true;
    };
  
    //!\brief Returns whether the \c i th root is definitely a multiple root of the isolated polynomial
    virtual bool is_certainly_multiple_root(size_t /* i */) const {
      return false;
    };
    
    //! returns an upper bound of the \c i th root multiplicity                      
    virtual int upper_bound_for_multiplicity(size_t /* i */) const { 
      return 1;
    };

    //! returns the \c i th root multiplicity                      
    virtual int multiplicity_of_root(size_t /* i */) const { 
      return 1;
    };

    //! refines i-th root
    virtual void refine_interval(size_t i) {
      //std::cout << "RS refine_interval" << std::endl;
      CGAL_assertion(i >= 0);
      CGAL_assertion(static_cast<int>(i) < this->number_of_real_roots());
      CGAL_precondition(this->is_isolated());
      // TODO 2012 correct refinement?
      typename AK_1::Approximate_absolute_1()(_m_real_roots[i], 2*mpfi_get_prec(_m_real_roots[i].mpfi()));
    }

    //!@} // Access functions

private:

    //! Needed for the referencing counting mechanism
    virtual CGAL::Reference_counted_hierarchy<>* clone() {
      return new RS_ak_1_isolator_rep(*this);
    }

protected:
    //!\name Protected members
    //!@{
    
    //! the solutions
    std::vector< Algebraic_real_1 > _m_real_roots;
   
    //!@}

}; // RS_ak_1_isolator_rep


// TODO 2012 document RS_ak_1_isolator
template< typename Polynomial_, typename Bound_, typename HandlePolicy = CGAL::Handle_policy_no_union  >
class RS_ak_1_isolator : 
  public internal::Generic_isolator< Polynomial_, Bound_, HandlePolicy, internal::RS_ak_1_isolator_rep< Polynomial_, Bound_, HandlePolicy > >  {
  
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
  typedef internal::RS_ak_1_isolator_rep< Polynomial, Bound, Handle_policy > Rep;

  //! the base class
  typedef internal::Generic_isolator< Polynomial, Bound, Handle_policy, Rep > Base;       
   
  //! the class itself
  typedef RS_ak_1_isolator< Polynomial, Bound, Handle_policy > Self;

  //!@} // Public types

  //!\name Constructors
  //!@{
  
#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_DEFAULT_CONSTRUCTIBLE
 protected:
  RS_ak_1_isolator(); // = disable
#else
  // no default needed as standard constructor has default values for all its arguments
#endif

#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_ASSIGNABLE
 protected:
  RS_ak_1_isolator(const Self&); // = disable
  RS_ak_1_isolator& operator=(const Self&); // = disable
#endif

  
 public:

  //! standard constructor
  RS_ak_1_isolator(const Polynomial& p = Polynomial(0)) :
    Base(new Rep(p)) {    
  }
  
  //! TODO 2012 implement constructor for x-range?
  RS_ak_1_isolator(const Polynomial& p,
                   const Bound& left,
                   const Bound& right);
  // : Base(new Rep(p, left, rightk)) {    
    
  // TODO 2012 copy constructor or is generated one sufficient?

  //!@} // Constructors

}; // RS_ak_1_isolator

} // namespace internal

} //namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_D_AK_1_ISOLATOR_H
// EOF
