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
// $URL$
// $Id$
// 
//
// Author(s)     :  Eric Berberich <eric.berberich@cgal.org>
//
// ============================================================================

/*! \file CGAL/Algebraic_kernel_d/Generic_solator.h
  \brief Defines class for generic isolator (as base).
  
  An Isolator isolates the (real or complext) roots of a polynomial.
  This file provides several Isolators.

*/

#ifndef CGAL_ALGEBRAIC_KERNEL_D_GENERIC_ISOLATOR_H
#define CGAL_ALGEBRAIC_KERNEL_D_GENERIC_ISOLATOR_H

#ifndef CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_DEFAULT_CONSTRUCTIBLE
#define CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_DEFAULT_CONSTRUCTIBLE 1
#endif

#ifndef CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_ASSIGNABLE
#define CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_ASSIGNABLE 1
#endif

#include <CGAL/config.h>

#include <CGAL/tags.h>

#include <CGAL/Handle_with_policy.h>

namespace CGAL {

namespace internal {


  template< typename Tag, bool refine = boost::is_same< Tag, CGAL::Tag_true >::value /* (TODO actually: has_refine_interval< Rep, mpl_::true_ >::value)*/ >
    struct Refineable_isolator_rep {

      //! Refine the <tt>i</tt>th isolating interval
      virtual void refine_interval(size_t /* i */) = 0;

  };

  template< typename Refine >
  struct Refineable_isolator_rep< Refine, false> {
    // empty on purpose
    // Comment: An isolator_rep for which Generic_isolator::refine_interval should be callable, must derived from 
    //          Refineable_isolator_rep< Tag, true >
  };  

  template < typename Bound_, typename HandlePolicy, typename RefineIntervalProvidedTag > // no default on policy, should be decided in higher level
  class Generic_isolator_rep_base : 
    // REMARK: Do we really want to have isolators policy-based ref-counted?
    public HandlePolicy::template Hierarchy_base<CGAL_ALLOCATOR(char) >::Type,
    public Refineable_isolator_rep< RefineIntervalProvidedTag > {

    // this is the root class of all reps

  public:
    //!\name Public types
    //!@{ 
  
    //! the bound type
    typedef Bound_ Bound;

    //! the policy type
    typedef HandlePolicy Handle_policy;

    //!@} // Public types

    //!\name Tags
    //!@{

    // TODO 2012 implement tag dispatches for other functions (some might not be needed)
    // TODO 2012 collect all tags in one class and provide a class under a common typename for each isolator-rep 
    //           (and following for each interface class)

    //! is a squarefree polynomial required
    typedef CGAL::Tag_true  Squarefree_polynomial_required_tag;
    
    //! are exact coefficients required (or bitstream coefficients needed)
    typedef CGAL::Tag_true  Exact_coefficients_required_tag;
    
    //! isolator provides constructor for range (instead of x-axis)
    typedef CGAL::Tag_false Range_constructor_provided_tag;
        
    //! isolator provides method to refine intervals
    typedef CGAL::Tag_false Refine_interval_provided_tag;

    //! is the isolated capable of reporting complex roots, too
    typedef CGAL::Tag_false Complex_solving_provided_tag;

    //! are roots isolated upon termination (or needs interaction with a callback function)
    typedef CGAL::Tag_false Isolation_upon_termination_guaranteed_tag;

    //!@}

  protected:
    //!\name Constructors
    //!@{

    //! default constructor
    Generic_isolator_rep_base() :
      _m_is_isolated(false) {
    }

    //!@} // Constructors

  public:
    //!\name Destructor
    //!@{

    virtual ~Generic_isolator_rep_base() {
    }

    //!@} // Destructor

  public:
    //!\name Virtual functions
    //!@{

    //! returns true if all roots are isolated
    virtual bool is_isolated() const {
      return _m_is_isolated;
    }

    //!@} Virtual functions

  public:
    //!\name Pure virtual functions
    //!@{

    //! isolates the roots of the polynomial (if not yet done)
    virtual void isolate() = 0;

    //! returns the number of real roots
    virtual int number_of_real_roots() const = 0;
    
    //! the lower bound of the \c i th root
    virtual Bound left_bound(size_t /* i */) const = 0;
  
    //! the upper bound of the \c i th root
    virtual Bound right_bound(size_t /* i */) const = 0;

    //! returns whether i-th root is exact
    virtual bool is_exact_root(size_t /* i */) const = 0;

    //! \brief Returns whether the \c i th root is definitely a simple root of the isolated polynomial
    virtual bool is_certainly_simple_root(size_t /* i */) const = 0;
  
    //!\brief Returns whether the \c i th root is definitely a multiple root of the isolated polynomial
    virtual bool is_certainly_multiple_root(size_t /* i */) const = 0;
    
    //! returns an upper bound of the \c i th root multiplicity                      
    virtual int upper_bound_for_multiplicity(size_t /* i */) const = 0;

    //! returns the \c i th root multiplicity                      
    virtual int multiplicity_of_root(size_t /* i */) const = 0;
  
    //!@} // Pure virtual functions

  protected:

    //!\name Protected members
    //!@{

    //! Has isolation already taken place
    // TODO 2012 be more distinguishing (ie for each interval)
    // TODO 2012 make _number_of_real_roots a boost::optional and get rid of this?
    mutable bool _m_is_isolated;

    //!@}
  
  }; // Generic_isolator_rep_base

  //! class of generic isolator representation
  template< typename Polynomial_, typename Bound_, typename HandlePolicy, typename RefineIntervalProvidedTag > // no default on policy, should be decided in higher level
  class Generic_isolator_rep : public Generic_isolator_rep_base< Bound_, HandlePolicy, RefineIntervalProvidedTag > {

  public:
  
    //!\name Public types
    //!@{
  
    //! first template parameter
    typedef Polynomial_ Polynomial;
    
    //! second template parameter
    typedef Bound_ Bound;

    //! third template parameter
    typedef HandlePolicy Handle_policy;

    //! the class itself
    typedef Generic_isolator_rep< Polynomial, Bound, Handle_policy, RefineIntervalProvidedTag > Self;

    //!@} // Public types
  
    //!\name Constructors
    //!@{

  protected:

    //! default constructor (is protected, as only derived classes are allowed to call it)
    Generic_isolator_rep() :
      _m_polynomial(0) {
    }

#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_ASSIGNABLE
  protected:
    Generic_isolator_rep(const Self&) {} // = default
    Generic_isolator_rep& operator=(const Self&); // = disable
#endif


  protected:

    //! standard constructor (is protected, as only derived classes are allowed to call it)
    Generic_isolator_rep(const Polynomial& p) :
      _m_polynomial(p) {
    }

    //!@} // Constructors

  public:
  
    //!\name Destructor
    //!@{

    //! Destructor (does nothing)
    virtual ~Generic_isolator_rep() {
    }
  
    //!@} // Destructor

    //!\name Access functions
    //!@{

    //! returns the polynomial
    virtual Polynomial polynomial() const {
      return _m_polynomial;
    }

    //!@} // Access functions

  protected:

    //!\name Protected members
    //!@{

    //! instance of polynomial
    Polynomial _m_polynomial;

    //!@} // Protected members

    // Comment: All derived reps are required to implement pure virtual functions of base rep!

  }; // Generic_isolator_rep

} // namespace internal

} // namespace CGAL

#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_isolator_reps.h>  

namespace CGAL {

namespace internal {
    
  // TODO 201x add more: Sturm, Continued Fraction, Backshear with Arcavoid

  // in the following CGAL::Null_tag is just used as a dummy class. struct Dummy {}; would also suffice here

  template < class Isolator, typename Size_t /* = size_t */ > 
  struct _Refine_interval {
    
    void operator()(Isolator& isolator, Size_t i) {
      CGAL_precondition(i >= 0);
      CGAL_precondition(static_cast<int>(i) < isolator.number_of_real_roots());
      CGAL_precondition(isolator.is_isolated());
      dynamic_cast< typename Isolator::RRep* >(isolator.ptr())->refine_interval(i);
    }
    
  };
  
  template < class Isolator > 
  struct _Refine_interval< Isolator, CGAL::Null_tag > {
    
    void operator()(Isolator& isolator, CGAL::Null_tag i) {} // to silent the compiler
  
  };

  
  //! generic isolator class (handle)
  template < typename Polynomial_, typename Bound_, typename HandlePolicy, typename IsolatorRep /* = Generic_isolator_rep_base< Bound_, HandlePolicy > */ > // no default on policy, should be decided in higher level
  class Generic_isolator : 
    public CGAL::Handle_with_policy< Generic_isolator_rep_base< Bound_, HandlePolicy, typename IsolatorRep::Refine_interval_provided_tag > > {

    //!\name Tags
    //!@{
    
  public:
    //! CGAL::Tag_true if "refine_interval" is provided, otherwise CGAL::Tag_false
    typedef typename IsolatorRep::Refine_interval_provided_tag Refine_interval_provided_tag;
    
    //!@} // Tags


  public:

    //!\name Public types
    //!@{

    //! first template parameter 
    typedef Polynomial_ Polynomial;

    //! second template parameter
    typedef Bound_ Bound;

    //! third template parameter
    typedef HandlePolicy Handle_policy;
    
    //! fourth template parameter
    typedef IsolatorRep Isolator_rep;

    // the generic representation class
    typedef Generic_isolator_rep_base< Bound, Handle_policy, Refine_interval_provided_tag > Rep;

    // the generic representation class
    typedef Generic_isolator_rep< Polynomial, Bound, Handle_policy, Refine_interval_provided_tag > PRep;
    
    //! refinable rep 
    typedef Refineable_isolator_rep< Refine_interval_provided_tag > RRep;

    // the base type
    typedef CGAL::Handle_with_policy< Rep > Base;
    
    //! the class itself
    typedef Generic_isolator< Polynomial, Bound, Handle_policy, Isolator_rep > Self;
    
    //!@} // Public types

    //!\name Constructors
    //!@{

#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_DEFAULT_CONSTRUCTIBLE
  protected:
    Generic_isolator(); // = disable
#else
   protected:
    //! default constructor
    Generic_isolator() : Base(typename Base::Use_with_initialize_with()) {} // = default
#endif

#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_ASSIGNABLE
  protected:
    //    Generic_isolator(const Self&); // = disable
    //Generic_isolator& operator=(const Self&); // = disable
#endif

  protected:
    
    //! special constructor to avoid activated default constructor in this class for derived classes that
    Generic_isolator(typename Base::Use_with_initialize_with uwiw) : Base(uwiw) {} 

  public:
    //! standard constructor from rep
    Generic_isolator(Rep *rep) :
      Base(rep) {
    }
    
    //!@} // Constructors

  public:
    //!\name Destructor
    //!@{

    //! empty virtual desctructor
    virtual ~Generic_isolator() {
    }

    //!@} // Destructor

    //!\name Access members
    //!@{

    //! returns the polynomial
    virtual Polynomial polynomial() const {
      return dynamic_cast<const PRep*>(this->ptr())->polynomial();
    }

    //! returns true if all roots have been isolated
    virtual bool is_isolated() const {
      return this->ptr()->is_isolated();
    }

    //! isolates the roots of the polynomial (if not yet done)
    virtual void isolate() {
      this->ptr()->isolate();
    }

    //! returns the number of real roots
    virtual int number_of_real_roots() const { 
      CGAL_precondition(this->is_isolated());
      return this->ptr()->number_of_real_roots();
    }
 
    //! the lower bound of the \c i th root
    virtual Bound left_bound(size_t i) const {
      CGAL_precondition(i >= 0);
      CGAL_precondition(static_cast<int>(i) < number_of_real_roots());
      CGAL_precondition(this->is_isolated());
      return this->ptr()->left_bound(i);
    } 
  
    //! the upper bound of the \c i th root
    virtual Bound right_bound(size_t i) const {
      CGAL_precondition(i >= 0);
      CGAL_precondition(static_cast<int>(i) < number_of_real_roots());
      CGAL_precondition(this->is_isolated());
      return this->ptr()->right_bound(i);
    }

    //! returns whether i-th root is exact
    virtual bool is_exact_root(size_t i) const { 
      CGAL_precondition(i >= 0);
      CGAL_precondition(static_cast<int>(i) < number_of_real_roots());
      CGAL_precondition(this->is_isolated());
      return this->ptr()->is_exact_root(i);
    }

    //! returns true if the <tt>i</tt>th root is known to be a simple root of the curve.
    virtual bool is_certainly_simple_root(size_t i) const {
      CGAL_precondition(i >= 0);
      CGAL_precondition(static_cast<int>(i) < number_of_real_roots());
      CGAL_precondition(this->is_isolated());
      return this->ptr()->is_certainly_simple_root(i);
    }
	
    //! returns true if the <tt>i</tt>th root is known to be a multiple root of the curve.
    virtual bool is_certainly_multiple_root(size_t i) const {
      CGAL_precondition(i >= 0);
      CGAL_precondition(static_cast<int>(i) < number_of_real_roots());
      CGAL_precondition(this->is_isolated());
      return this->ptr()->is_certainly_multiple_root(i);
    }
    
    //! returns the multiplicity of the root if know, otherwise -1
    virtual int multiplicity_of_root(size_t i) const {
      CGAL_precondition(i >= 0);
      CGAL_precondition(static_cast<int>(i) < number_of_real_roots());
      CGAL_precondition(this->is_isolated());
      return this->ptr()->multiplicity_of_root(i);
    }
    
    //! returns an upper bound for the multiplicity of the ith root
    virtual int upper_bound_for_multiplicity(size_t i) const {
      CGAL_precondition(i >= 0);
      CGAL_precondition(static_cast<int>(i) < number_of_real_roots());
      CGAL_precondition(this->is_isolated());
      return this->ptr()->upper_bound_for_multiplicity(i);
    }
    
  private:

    // some logic to disable "refine_interval(size_t)" by only providing "refine_interval(CGAL::Null_tag)", 
    // in case Refine_interval_provided_tag is not equal to CGAL::Tag_true:

    typedef size_t refine_interval_arg_type; // refine_intervals's argument type 
    typedef typename boost::mpl::if_< boost::is_same< CGAL::Tag_true, typename Isolator_rep::Refine_interval_provided_tag >, refine_interval_arg_type, CGAL::Null_tag >::type cond_true_f_arg_type; 
    
    // partially specialized functor must be non-class in C++03, but needs access to ptr() member, thus we grant friendship:
    friend class _Refine_interval< Self, cond_true_f_arg_type >;

  public:
    //! Refine the <tt>i</tt>th isolating interval
    virtual void refine_interval(cond_true_f_arg_type /* should be size_t, and calls to fail if it is CGAL::Null_type */ i) {
      // needed as a function specialization does not work in C++03
      _Refine_interval< Self, cond_true_f_arg_type >()(*this, i);
    }

    //!@} // Access members

  }; // Generic_isolator
  
} // namespace internal

} // namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_D_GENERIC_ISOLATOR_H
// EOF
