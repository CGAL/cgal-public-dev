// Copyright (c) 2006-2012 Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//                 Eric Berberich <eric.berberich@cgal.org>
//
// ============================================================================

#ifndef CGAL_ALGEBRAIC_KERNEL_D_BITSTREAM_DESCARTES_ISOLATOR_REPS_H
#define CGAL_ALGEBRAIC_KERNEL_D_BITSTREAM_DESCARTES_ISOLATOR_REPS_H 1

#include <CGAL/config.h>

#include <CGAL/Algebraic_kernel_d/Generic_isolator.h>

#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_rndl_tree_traits.h>

// NOTE: If this flag is set, you need EXACUS!
#if CGAL_ACK_BITSTREAM_USES_E08_TREE
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_E08_tree.h>
#else
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_rndl_tree.h>
#endif

#include <CGAL/Algebraic_kernel_d/exceptions.h>

namespace CGAL { 

namespace internal {


/*
 * \brief The base class for all variants of the Bitstream Descartes method.
 *
 */
template< class BitstreamCoefficientKernel, typename HandlePolicy > // no default on policy, should be decided in higher level
class Generic_bitstream_descartes_isolator_rep : 
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

  //!\name Public types
  //!@{

  //! first template parameter: Bitstream coefficient kernel
  typedef BitstreamCoefficientKernel Bitstream_coefficient_kernel;
  
  //! the policy class
  typedef HandlePolicy Handle_policy;

  //! The Coeeficient type of the input polynomial
  typedef typename Bitstream_coefficient_kernel::Coefficient Coefficient;

  //! The polynomial type
  typedef typename Bitstream_coefficient_kernel::Polynomial Polynomial;

  //! How the boundaries of the isolating intervals are represented
  typedef typename Bitstream_coefficient_kernel::Bound Bound;
  
  //! the base class
  typedef Generic_isolator_rep< Polynomial, Bound, Handle_policy, Refine_interval_provided_tag > Base;

  //! the class itself
  typedef Generic_bitstream_descartes_isolator_rep< Bitstream_coefficient_kernel, Handle_policy > Self;

  protected:
// TODO 2012 make traits protected? if nothing in AK_3 needs it!
  //! The traits class for approximations
  typedef CGAL::internal::Bitstream_descartes_rndl_tree_traits< Bitstream_coefficient_kernel >
  Bitstream_descartes_rndl_tree_traits;

  public:
//! The type of the used Bitstream Descartes tree
#if CGAL_ACK_BITSTREAM_USES_E08_TREE
  typedef CGAL::internal::Bitstream_descartes_E08_tree<Bitstream_descartes_rndl_tree_traits> Bitstream_tree;
#else
  typedef CGAL::internal::Bitstream_descartes_rndl_tree<Bitstream_descartes_rndl_tree_traits> Bitstream_tree;
#endif

  //! The used integer type
  typedef typename Bitstream_coefficient_kernel::Integer Integer;
  
  protected:
  //! The type for the iterator of the nodes of the bitstream tree
  typedef typename Bitstream_tree::Node_iterator Node_iterator;
  
  //! The same as constant iterator
  typedef typename Bitstream_tree::Node_const_iterator Node_const_iterator;

  //!@} // Public types

  //!\name Constructors
  //!@{

    //! Default constructor (does nothing)
#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_DEFAULT_CONSTRUCTIBLE
  protected:
#else
  public:
#endif
  Generic_bitstream_descartes_isolator_rep() {} // = default

  // needs no special assignable-implementation as no pointers 
  // and this is a Rep of a handle class (depending on Handle_policy!)
#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_ASSIGNABLE
  protected:
  // Generic_bitstream_descartes_isolator_rep(const Self&); // = default
  Generic_bitstream_descartes_isolator_rep& operator=(const Self&); // = disable
#endif

  protected:
  /*! 
   * Constructor computing an interval containing all real roots of \c f,
   * and initialising the Bitstream Descartes tree
   */
  Generic_bitstream_descartes_isolator_rep(const Polynomial& p,
                                           Bitstream_coefficient_kernel bck) : 
    Base(p), 
    _m_traits(Bitstream_descartes_rndl_tree_traits(bck)) {
    
    Integer lower,upper;
    long log_div;
    this->interval(p, lower, upper, log_div);
    //AcX_DSTREAM("f: " << f << std::endl);
    if (CGAL::degree(p) > 0) {
      _m_bitstream_tree 
#if CGAL_ACK_BITSTREAM_USES_E08_TREE
        = Bitstream_tree(-log_div,
                         p.begin(),
                         p.end(),
                         typename Bitstream_tree::Monomial_basis_tag(),
                         _m_traits);
#else
      = Bitstream_tree(lower,upper,log_div,
                       p.begin(),
                       p.end(),
                       typename Bitstream_tree::Monomial_basis_tag(),
                       _m_traits);
#endif
      
      if (_m_bitstream_tree.begin() == _m_bitstream_tree.end()) {
        _m_number_of_intervals = 0;
      } else {
        _m_number_of_intervals = 1;
      }
    } else {
      _m_number_of_intervals=0;
    }
  }
  
  /*! 
   * Constructor that copies the Bitstream tree given from outside
   * and initialising the Bitstream Descartes tree
   * The tree must "fit" to the polynomial
   */
  Generic_bitstream_descartes_isolator_rep(const Polynomial& p,
                                           Bitstream_tree tree,
                                           Bitstream_coefficient_kernel bck) : 
    Base(p),
    _m_bitstream_tree(tree),
    _m_traits(Bitstream_descartes_rndl_tree_traits(bck)) {
    
    tree.set_traits(_m_traits);
    
    _m_number_of_intervals = 0;
    
    for(Node_iterator curr = _m_bitstream_tree.begin();
        curr != _m_bitstream_tree.end();
        curr++) {
      _m_number_of_intervals++;
    }
  }
  
  //!@} // Constructors

  //!\name Destructor
  //!@{
  public:

  //! Destructor (does nothing)
  virtual ~Generic_bitstream_descartes_isolator_rep() {
  }
  
  //!@} // Desctructor

  public:

  //!\name Access functions
  //!@{
   
  // polynomial() is inherited from Generic_isolator_rep
 
  // is_isolated() const is inherited from Generic_isolator_rep

  /*! 
   * \brief isolates the root of \c f
   *
   * The mechanism is the following: The \c bitstream_tree member of the
   * object is transformed via subdivision until the 
   * \c termination_condition routine of the object returns true. When this
   * happens, the \c process_nodes routine of the object is called.
   */
  virtual void isolate() {
    
    if (this->is_isolated()) {
      return;
    }
    
    //AcX_DSTREAM("Starting isolation" << std::endl);
    
    Node_iterator curr = _m_bitstream_tree.begin(),dummy,new_curr;
    
    if(curr == _m_bitstream_tree.end()) {
      Base::_m_is_isolated = true;
      return;
    }
    
    int newly_created;
    
    while (!this->termination_condition()) {
      
      if (curr == _m_bitstream_tree.end()) {
        curr = _m_bitstream_tree.begin();
      }
      
      if (_m_bitstream_tree.max_var(curr) == 1) {
        ++curr;
      }
      else {
        //AcX_DSTREAM("Subdivision at " 
        //<< CGAL::to_double(_m_bitstream_tree.lower(curr)) << " " 
        //<< CGAL::to_double(_m_bitstream_tree.upper(curr)) << std::flush);
        newly_created = _m_bitstream_tree.subdivide(curr,dummy,new_curr);
        _m_number_of_intervals += newly_created-1;
        curr = new_curr;
        //AcX_DSTREAM("done" << std::endl);
      }
      
    }
    this->process_nodes();
    Base::_m_is_isolated = true;
  }
  
  //! returns the number of detected isolating intervals
  virtual int number_of_real_roots() const {
    return _m_number_of_intervals;
  }
  
  //! The lower bound of the \c i th root
  virtual Bound left_bound(size_t i) const  {
    CGAL_precondition(i >= 0);
    CGAL_precondition(i < _m_number_of_intervals);
    CGAL_precondition(this->is_isolated());
    Node_const_iterator curr = _m_bitstream_tree.begin();
    std::advance(curr,i);
    return _m_bitstream_tree.lower(curr);
  } 
  
  //! The upper bound of the \c i th root
  virtual Bound right_bound(size_t i) const {
    CGAL_precondition(i >= 0);
    CGAL_precondition(i < _m_number_of_intervals);
    CGAL_precondition(this->is_isolated());
    Node_const_iterator curr = _m_bitstream_tree.begin();
    std::advance(curr,i);
    return _m_bitstream_tree.upper(curr);
  } 

  //! returns whether root is exact 
  bool is_exact_root(size_t CGAL_precondition_code(i)) const { 
    CGAL_precondition(i >= 0);
    CGAL_precondition(i < _m_number_of_intervals);
    CGAL_precondition(this->is_isolated());
    return false; 
  }

  // is_certainly_simple_root(size_t /* i */) const must be implemented by derived classes
 
  // is_certainly_multiple_root(size_t /* i */) const must be implemented by derived classes

  //! returns an upper bound of the \c i th root multiplicity                      
  virtual int upper_bound_for_multiplicity(size_t i) const {
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < _m_number_of_intervals);
    CGAL_precondition(this->is_isolated());
    Node_const_iterator curr = _m_bitstream_tree.begin();
    std::advance(curr,i);
    return _m_bitstream_tree.min_var(curr);
  } 

  // int multiplicity_of_root( size_t /* i */) const must be implemented by derived classes

  //! Computes a better approximation of the \c i th root of the
  virtual void refine_interval(size_t i) {
    //std::cout << "Bitstream refine_interval" << std::endl;
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < _m_number_of_intervals);
    CGAL_precondition(this->is_isolated());
    Node_iterator curr = _m_bitstream_tree.begin(), begin, end, 
      new_begin, helper;
    std::advance(curr,i);
    int intervals = 1;
    end = curr;
    ++end;
    begin = curr;
    do {
      //std::cout << _m_bitstream_tree.lower(begin) << " " << _m_bitstream_tree.upper(begin) << std::endl;
      //std::cout << _m_bitstream_tree.min_var(begin) << " " << _m_bitstream_tree.max_var(begin) << std::endl;
      int new_intervals = _m_bitstream_tree.subdivide(begin,new_begin,helper);
      intervals += new_intervals-1;
      begin = new_begin;
      curr = helper;
      
      // Fixes the bug when a interval splits, and the leftmost subintervals
      // has no children with sign variation >=1
      if (intervals == 1) {
        break;
      }
      if (new_intervals == 0) {
        continue;
      }
      
      while(curr != end) {
        intervals += _m_bitstream_tree.subdivide(curr,new_begin,helper)-1;
        curr = helper;
      }
      
    }
    while (intervals != 1);
    //std::cout << "Refined " << left_bound(i) << " " << right_bound(i) << std::endl; 
    
  }

  // more bitstream specific code
    
  //! access to kernel
  Bitstream_coefficient_kernel bck() const {
    return _m_traits.bck();
  }

  protected:
  //! access to traits
  Bitstream_descartes_rndl_tree_traits traits() const {
    return _m_traits;
  }

  public:

  //! access to tree
  Bitstream_tree tree() const {
    return _m_bitstream_tree;
  }

  /*!
   * \brief Computes an interval containing all real roots of \c p,
   * using the Fujiwara root bound.
   *
   * So far, the \c log_div variable is always set to zero, this means
   * that <i>[lower,upper]</i> is the interval containing all real roots
   */
  virtual void interval(const Polynomial& p, Integer& lower, 
                        Integer& upper, long& log_div) {
    
    
    typename Bitstream_descartes_rndl_tree_traits::Lower_bound_log2_abs
      lower_bound_log2_abs = this->_m_traits.lower_bound_log2_abs_object();
    typename 
      Bitstream_descartes_rndl_tree_traits::Upper_bound_log2_abs_approximator
      upper_bound_log2_abs_approximator 
      = this->_m_traits.upper_bound_log2_abs_approximator_object();
    //AcX_DSTREAM("Fujiwara bound.." << p <<  std::endl);
#if CGAL_ACK_BITSTREAM_USES_E08_TREE
    log_div = -CGAL::internal::Fujiwara_root_bound_log
      (p.begin(),
       p.end(),
       lower_bound_log2_abs,
       upper_bound_log2_abs_approximator
      );
#else
    log_div = -CGAL::internal
      ::Fujiwara_root_bound_log
      (p.begin(),
       p.end(),
       lower_bound_log2_abs,
       upper_bound_log2_abs_approximator
      );
#endif
    
    //AcX_DSTREAM("Fujiwara returns " << log_div << std::endl);
    // To be sure
    log_div--;
    lower=Integer(-1);
    upper=Integer(1);
    return;
    
  }

  /*!
   * \brief Gives an opportunity to process the nodes after
   * the subdivision steps are finished
   *
   * This method must be specialised by derived classes, but can
   * remain empty in many cases.
   */ 
  virtual void process_nodes() = 0;

  /*! 
   * \brief When does the isolation algorithm terminate?
   *
   * This method must be specialised by derived classes
   */      
  virtual bool termination_condition() = 0;
   

private:

  //!\name Protected members
  //!@{

  //! type to distinguish used constructor
protected:
  
  //! The tree of the Bitstream Descartes method
  mutable Bitstream_tree _m_bitstream_tree;

  //! The traits class
  // TODO 2012 const ref?
  Bitstream_descartes_rndl_tree_traits _m_traits;
    
  //! The number of detected isolating intervals
  mutable size_t _m_number_of_intervals;
  
  //!@} // Protected members

}; // Generic_bitstream_descartes_isolator_rep

/*
 * \brief Representation for square free polynomials
 */
template< class BitstreamCoefficientKernel, typename HandlePolicy = CGAL::Handle_policy_no_union>
class Square_free_bitstream_descartes_isolator_rep :
  public Generic_bitstream_descartes_isolator_rep< BitstreamCoefficientKernel, HandlePolicy > {
  
  
public:
 
  //!\name Public types
  //!@{
  
  //! first template parameter: Bitstream coefficient kernel
  typedef BitstreamCoefficientKernel Bitstream_coefficient_kernel;

  //! the policy class
  typedef HandlePolicy Handle_policy;

  //! The traits class for approximations
  typedef CGAL::internal::Bitstream_descartes_rndl_tree_traits< Bitstream_coefficient_kernel >
  Bitstream_descartes_rndl_tree_traits;

  //! The generic representation
  typedef Generic_bitstream_descartes_isolator_rep< Bitstream_coefficient_kernel, Handle_policy > Base;
  
  //! the class itself
  typedef Square_free_bitstream_descartes_isolator_rep< Bitstream_coefficient_kernel, Handle_policy > Self;

  //! Polynomial type
  typedef typename Base::Polynomial Polynomial;
  
  //! Iterator for the leaves in the bitstream tree
  typedef typename Base::Node_iterator Node_iterator;
  
  //! The type of the tree that controls the Bitstream instance
  typedef typename Base::Bitstream_tree Bitstream_tree;
  
  //!@}

  //!\name Tags
  //!@{
    
public:
  //! "refine_interval" function is implemented
  typedef typename Base::Refine_interval_provided_tag Refine_interval_provided_tag;

  //!@} // Tags

  //!\name Constructors
  //!@{

#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_DEFAULT_CONSTRUCTIBLE
  protected:
  Square_free_bitstream_descartes_isolator_rep();
#else
  public:
    //! Default constructor (does nothing)
  Square_free_bitstream_descartes_isolator_rep() :
    Base()
  {
  };
#endif

  // needs no special assignable-implementation as no pointers 
  // and this is a Rep of a handle class (depending on Handle_policy!)
#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_ASSIGNABLE
  protected:
    Square_free_bitstream_descartes_isolator_rep(const Self&) {} // = default
    Square_free_bitstream_descartes_isolator_rep& operator=(const Self&); /* = disable */
#endif

  public:
  /*! 
   * \brief Constructor with the square free polynomial <tt>f<tt>.
   */
  Square_free_bitstream_descartes_isolator_rep(const Polynomial& p,
                                               Bitstream_coefficient_kernel bck) :
    Base(p, bck) {
      // nothing beyond initialization
  }

  /*! 
   * \brief Constructor with the square free polynomial <tt>f<tt>.
   */
  Square_free_bitstream_descartes_isolator_rep(const Polynomial& p,
                                               Bitstream_tree tree,
                                               Bitstream_coefficient_kernel bck) :
     Base(p, tree, bck) {
      // nothing beyond initialization
  }

  //!@} // Constructors

  //!\name Access functions
  //!@{

  //! Always true
  virtual bool is_certainly_simple_root(size_t /* i */) const {
    return true;
  }
	
  //! Always false
  virtual bool is_certainly_multiple_root(size_t /* i */) const {
    return false;
  }

  //! Always 1
  virtual int multiplicity_of_root(size_t /* i */) const {
    return 1;
  }
      
  //! nothing to do here	
  virtual void process_nodes() {
    return;
  }

  /*!
   * \brief Terminates when all detected roots are simple
   */
  virtual bool termination_condition() {
    for(Node_iterator curr=Base::_m_bitstream_tree.begin();
        curr != Base::_m_bitstream_tree.end();curr++) {
      if(Base::_m_bitstream_tree.max_var(curr)!=1) {
        return false;
      }
    }
    return true;
  }

  //!@}

  private:

  //! Needed for reference counting
  virtual CGAL::Reference_counted_hierarchy<>* clone() {
    return new Square_free_bitstream_descartes_isolator_rep(*this);
  }
	
};

} // namespace internal

} // namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_D_BITSTREAM_DESCARTES_REP_H 1
