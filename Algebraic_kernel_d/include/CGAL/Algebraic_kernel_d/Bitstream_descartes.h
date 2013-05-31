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
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//                 Eric Berberich <eric.berberich@cgal.org>
//
// ============================================================================

#ifndef CGAL_BITSTREAM_DESCARTES_H
#define CGAL_BITSTREAM_DESCARTES_H 1

#if !CGAL_AK_USE_OLD_BITSTREAM_DESCARTES

#include <CGAL/config.h>

#include <CGAL/Algebraic_kernel_d/Generic_isolator.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_isolator_reps.h>

namespace CGAL {

namespace internal {

// TODO 2012 keep Bitstream_descartes in internal namespace? (include BS_isolator_rep.h here?)

// TODO 2012 document BitstreamCoefficientKernel concept
/*!
 * \brief Class for the Bitstream Descartes method
 *
 * Class for the real root isolation of polynomials, using the Bitstream
 * Descartes method. The polynomials coefficient type is arbitrary, the 
 * approximations of the coefficient type are obtained with the 
 * \c BitstreamCoefficientKernel parameter. For the requirements
 * of this traits class, see the documentation of the 
 * BitstreamCoefficientKernel concept, which basically needs 
 * a few types, Is_zero and Convert_to_bfi).
 *
 * Internally, an instance of CGAL::Bitstream_descartes_rndl_tree is explored
 * in a specific way. That exploration strategy depends on the constructor
 * that is used to create the object. A tag is passed that defines the
 * variant of the Bitstream Descartes method: The Square_free_descartes_tag
 * starts the usual Bitstream method for square free integer polynomials.
 * With the M_k_descartes tag, it is able to handle one multiple root in 
 * favourable situations, the Backshear_descartes_tag allows to isolate
 * even more complicated polynomials, if the multiple roots with even
 * multiplicity can be refined from outside. See the corresponding
 * constructors for more information.
 * 
 */
template < class BitstreamCoefficientKernel, typename HandlePolicy = CGAL::Handle_policy_no_union >
class Bitstream_descartes : 
  public internal::Generic_isolator< typename BitstreamCoefficientKernel::Polynomial,
                                     typename BitstreamCoefficientKernel::Bound, 
                                     HandlePolicy,
                                     internal::Square_free_bitstream_descartes_isolator_rep< BitstreamCoefficientKernel, HandlePolicy > > {

public:

  //!\name Tags
  //!@{

  //! "refine_interval" function is implemented
  typedef CGAL::Tag_true Refine_interval_provided_tag;

  //!@} // Tags

 public:

  //!\name Public types
  //!@{

  //! first template parameter: Bitstream coefficient kernel
  typedef BitstreamCoefficientKernel Bitstream_coefficient_kernel;

  //! second template parameter
  typedef HandlePolicy Handle_policy;

  //! The polynomial type
  typedef typename Bitstream_coefficient_kernel::Polynomial Polynomial;

  //! How the boundaries of the isolating intervals are represented
  typedef typename Bitstream_coefficient_kernel::Bound Bound;

  //! the representation class
  typedef internal::Square_free_bitstream_descartes_isolator_rep< Bitstream_coefficient_kernel, Handle_policy > Rep;

  //! the base type
  typedef internal::Generic_isolator< Polynomial, Bound, Handle_policy, Rep > Base;

  //! the class itself
  typedef Bitstream_descartes< Bitstream_coefficient_kernel, Handle_policy > Self;

  //! type of bitstream tree
  typedef typename Rep::Bitstream_tree Bitstream_tree;

  //!@} // Public types

  protected:

  //! type of bitstream tree
  typedef typename Rep::Bitstream_descartes_rndl_tree_traits Bitstream_descartes_rndl_tree_traits;

  //!\name Constructors
  //!@{

#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_DEFAULT_CONSTRUCTIBLE
 protected:
  Bitstream_descartes(); // = disable
#else // CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_DEFAULT_CONSTRUCTIBLE
 public:
  //! Default constructor
 Bitstream_descartes() : Base(new Rep()) {}
#endif // CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_DEFAULT_CONSTRUCTIBLE

#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_ASSIGNABLE
 protected:
  Bitstream_descartes(const Self&); // = disable
  Bitstream_descartes& operator=(const Self&); // = disable
#endif // CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_ASSIGNABLE

 public:
  /*! 
   * \brief Constructor for a polynomial \c f
   */
  Bitstream_descartes(const Polynomial& p,
                      // TODO 2012 default bck? works only if no additional information is needed
                      Bitstream_coefficient_kernel bck = Bitstream_coefficient_kernel(),
                      bool isolate = true)
    : Base(new Rep(p, bck))
  {
    if (isolate) {
      this->isolate();
    }
  }

  /*! 
   * \brief Constructor for the square free Descartes method,
   * using a precomputed tree
   *
   * The polynomial \c f must not have multiple real roots. The 
   * Bitstream Descartes tree is traversed in a bfs manner until
   * all leaves have sign variation zero or one.
   * The tree must be adequate for the polynomial. 
   * Use that constructor only if you know what you're doing!
   */
  Bitstream_descartes(const Polynomial& p,
                      Bitstream_tree tree,
                      // TODO 2012 default bck?
                      Bitstream_coefficient_kernel bck = Bitstream_coefficient_kernel(),
                      bool isolate = true)
    : Base(new Rep(p, tree, bck))
  {
    if (isolate) {
      this->isolate();
    }
  }

  // TODO 2012 constructor for [a,b] range?

  //!@}

  //!\name Access members
  //!@{

  // the ones beyond those inherited from Base()

  //! Returns the kernel class
  Bitstream_coefficient_kernel bck() const {
    return this->ptr()->bck();
  } 

  //! returns the underlying tree
  Bitstream_tree tree() const {
    return this->ptr()->tree();
  }

  protected:

  //! Returns the traits class
  Bitstream_descartes_rndl_tree_traits traits() const {
    return this->ptr()->traits();
  }

  //!@} // Access members 

}; // Bitstream_descartes

} // namespace internal

} // namespace CGAL

#else // CGAL_AK_USE_OLD_BITSTREAM_DESCARTES

#include <CGAL/basic.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_type_generator.h>

#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_rndl_tree_traits.h>

#if CGAL_ACK_BITSTREAM_USES_E08_TREE
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_E08_tree.h>
#else // CGAL_ACK_BITSTREAM_USES_E08_TREE
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_rndl_tree.h>
#endif // CGAL_ACK_BITSTREAM_USES_E08_TREE

#include <CGAL/Algebraic_kernel_d/exceptions.h>

namespace CGAL { 

namespace internal {

//! enum to distinguish between different descartes instances
enum Bitstream_descartes_type {
  GENERIC_DESCARTES = 0, 
  SQUARE_FREE_DESCARTES = 1, //!< uses Square_free_descartes_tag constructor
  M_K_DESCARTES = 2, //!< uses M_k_descartes_tag constructor
  BACKSHEAR_DESCARTES = 3, //!< uses Backshear_descartes_tag constructor
  VERT_LINE_ADAPTER_DESCARTES = 4 // ! < uses Vert_line_adapter_descartes
};


//! forward declaration
template<typename BitstreamCoefficientKernel>
class Bitstream_descartes;

//! Tag for the square free Bitstream Descartes method
struct Square_free_descartes_tag {};

//! Tag for  the Bitstream m-k-Descartes method
struct M_k_descartes_tag {};

//! Tag for the Backshear Descartes method
struct Backshear_descartes_tag {};

//! Tag for the Exchange-Descartes method
struct Vert_line_adapter_descartes_tag {};


/*
 * \brief Thrown whenever a non-specialised virtual member function is called
 */
class Virtual_method_exception {};


/*
 * \brief The base class for all variants of the Bitstream Descartes method.
 *
 */
template<typename BitstreamCoefficientKernel,
      typename Policy=CGAL::Handle_policy_no_union>
class Generic_descartes_rep 
  : public Policy::template Hierarchy_base<CGAL_ALLOCATOR(char) >::Type {
  
public:

  //! first template parameter: Bitstream coefficient kernel
  typedef BitstreamCoefficientKernel Bitstream_coefficient_kernel;

  //! The traits class for approximations
  typedef CGAL::internal::Bitstream_descartes_rndl_tree_traits< Bitstream_coefficient_kernel >
  Bitstream_descartes_rndl_tree_traits;
  
  //! The Handle class
  typedef Bitstream_descartes<Bitstream_coefficient_kernel> Handle;
  
  //! The Coeeficient type of the input polynomial
  typedef typename Bitstream_coefficient_kernel::Coefficient Coefficient;
  
  //! The polynomial type
  typedef typename Bitstream_descartes_rndl_tree_traits::POLY Polynomial;
  
  typedef Generic_descartes_rep<Bitstream_coefficient_kernel> Self;
  
  //! The type of the used Bitstream Descartes tree
#if CGAL_ACK_BITSTREAM_USES_E08_TREE
  typedef CGAL::internal::Bitstream_descartes_E08_tree
  <Bitstream_descartes_rndl_tree_traits> 
  Bitstream_tree;
#else // CGAL_ACK_BITSTREAM_USES_E08_TREE
  typedef CGAL::internal::Bitstream_descartes_rndl_tree
  <Bitstream_descartes_rndl_tree_traits> 
  Bitstream_tree;
#endif // CGAL_ACK_BITSTREAM_USES_E08_TREE
  
  
  //! The used integer type
  typedef typename Bitstream_coefficient_kernel::Integer Integer;
  
  //! The type for the iterator of the nodes of the bitstream tree
  typedef typename Bitstream_tree::Node_iterator Node_iterator;
  
  //! The same as constant iterator
  typedef typename Bitstream_tree::Node_const_iterator Node_const_iterator;
  
  //! How the boundaries of the isolating intervals are represented
  typedef typename Bitstream_coefficient_kernel::Bound Bound;
  
  //! Default constructor (does nothing)
  Generic_descartes_rep(Bitstream_descartes_type type = GENERIC_DESCARTES) :
    type_(type) {
  };
  
  /*! 
   * Constructor computing an interval containing all real roots of \c f,
   * and initialising the Bitstream Descartes tree
   */
  Generic_descartes_rep(Bitstream_descartes_type type, 
                        Polynomial f,
                        Bitstream_coefficient_kernel bck) : 
    type_(type), 
    _m_f(f), 
    _m_traits(Bitstream_descartes_rndl_tree_traits(bck)), 
    _m_is_isolated(false) {
    
    Integer lower,upper;
    long log_div;
    this->interval(f, lower, upper, log_div);
    //AcX_DSTREAM("f: " << f << std::endl);
    if (CGAL::degree(f) > 0) {
      _m_bitstream_tree 
#if CGAL_ACK_BITSTREAM_USES_E08_TREE
        = Bitstream_tree(-log_div,
                         f.begin(),
                         f.end(),
                         typename Bitstream_tree::Monomial_basis_tag(),
                         _m_traits);
#else // CGAL_ACK_BITSTREAM_USES_E08_TREE
      = Bitstream_tree(lower,upper,log_div,
                       f.begin(),
                       f.end(),
                       typename Bitstream_tree::Monomial_basis_tag(),
                       _m_traits);
#endif // CGAL_ACK_BITSTREAM_USES_E08_TREE
      
      if (_m_bitstream_tree.begin() == _m_bitstream_tree.end()) {
        number_of_intervals = 0;
      } else {
        number_of_intervals = 1;
      }
    } else {
      number_of_intervals=0;
    }
  }
  
  /*! 
   * Constructor that copies the Bitstream tree given from outside
   * and initialising the Bitstream Descartes tree
   * The tree must "fit" to the polynomial
   */
  Generic_descartes_rep(Bitstream_descartes_type type, 
                        Polynomial f,
                        Bitstream_tree tree,
                        Bitstream_coefficient_kernel bck) : 
    type_(type), 
    _m_f(f), 
    _m_bitstream_tree(tree),
    _m_traits(Bitstream_descartes_rndl_tree_traits(bck)), 
    _m_is_isolated(false) {
    
    tree.set_traits(_m_traits);
    
    number_of_intervals = 0;
    
    for(Node_iterator curr = _m_bitstream_tree.begin();
        curr != _m_bitstream_tree.end();
        curr++) {
      number_of_intervals++;
    }
  }
  
  //! Destructor (does nothing)
  virtual ~Generic_descartes_rep() {
  }
  
  //! Needed for the referencing counting mechanism
  virtual CGAL::Reference_counted_hierarchy<>* clone() {
    return new Generic_descartes_rep(*this);
  }
  
  /*! 
   * \brief Computes a better approximation of the \c i th root of the
   * polynomial
   */
  virtual void refine_interval(int i) const {
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < number_of_intervals);
    Node_iterator curr = _m_bitstream_tree.begin(), begin, end, 
      new_begin, helper;
    std::advance(curr,i);
    int intervals = 1;
    end = curr;
    ++end;
    begin=curr;
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
  
  /*! 
   * \brief isolates the root of \c f
   *
   * The mechanism is the following: The \c bitstream_tree member of the
   * object is transformed via subdivision until the 
   * \c termination_condition routine of the object returns true. When this
   * happens, the \c process_nodes routine of the object is called.
   */
  virtual void isolate() {
    
    //AcX_DSTREAM("Starting isolation" << std::endl);
    
    Node_iterator curr = _m_bitstream_tree.begin(),dummy,new_curr;
    
    if(curr == _m_bitstream_tree.end()) {
      _m_is_isolated = true;
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
        number_of_intervals += newly_created-1;
        curr = new_curr;
        //AcX_DSTREAM("done" << std::endl);
      }
      
    }
    this->process_nodes();
    _m_is_isolated = true;
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
#else // CGAL_ACK_BITSTREAM_USES_E08_TREE
    log_div = -CGAL::internal
      ::Fujiwara_root_bound_log
      (p.begin(),
       p.end(),
       lower_bound_log2_abs,
       upper_bound_log2_abs_approximator
      );
#endif // CGAL_ACK_BITSTREAM_USES_E08_TREE
    
    //AcX_DSTREAM("Fujiwara returns " << log_div << std::endl);
    // To be sure
    log_div--;
    lower=Integer(-1);
    upper=Integer(1);
    return;
    
  }
  
  //! returns the number of detected isolating intervals
  virtual int number_of_real_roots() const {
    return number_of_intervals;
  }
  
  //! The lower bound of the \c i th root
  virtual Bound left_bound(int i) const  {
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < number_of_intervals);
    Node_const_iterator curr = _m_bitstream_tree.begin();
    std::advance(curr,i);
    return _m_bitstream_tree.lower(curr);
  } 
  
  //! The upper bound of the \c i th root
  virtual Bound right_bound(int i) const {
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < number_of_intervals);
    Node_const_iterator curr = _m_bitstream_tree.begin();
    std::advance(curr,i);
    return _m_bitstream_tree.upper(curr);
  } 
  
  //! Returns the polynomial which is isolated
  Polynomial polynomial() const {
    return _m_f;
  }
  
  /*! 
   * \brief When does the isolation algorithm terminate?
   *
   * This method must be specialised by derived classes
   */      
  virtual bool termination_condition() {
    throw Virtual_method_exception();
    return false;
  }
  
  /*!
   * \brief Gives an opportunity to process the nodes after
   * the subdivision steps are finished
   *
   * This method must be specialised by derived classes, but can
   * remain empty in many cases.
   */ 
  virtual void process_nodes() {
    throw Virtual_method_exception();
    return;
  }
  
  /*! \brief Returns whether the \c i th root is definitely a simple root
   * of the isolated polynomial
   *
   * Must be specialised by derived class
   */
  virtual bool is_certainly_simple_root(int) const {
    throw Virtual_method_exception();
    return false;
  }
  
  /*! \brief Returns whether the \c i th root is definitely a multiple root
   * of the isolated polynomial
   *
   * Must be specialised by derived class
   */
  virtual bool is_certainly_multiple_root(int) const {
    throw Virtual_method_exception();
    return false;
  }
  
  
  virtual int multiplicity_of_root(int CGAL_assertion_code(i)) const {
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < number_of_intervals);
    return -1;
  }
  
  virtual int upper_bound_for_multiplicity(int i) const {
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < number_of_intervals);
    Node_const_iterator curr = _m_bitstream_tree.begin();
    std::advance(curr,i);
    return _m_bitstream_tree.min_var(curr);
  } 
  
  //! Must be specialized by the derived class
  virtual int degree_of_gcd() const {
    throw Virtual_method_exception();
    return -1;
  }
  
  //! Must be specialized by the derived class
  virtual Polynomial square_free_part() const {
    throw Virtual_method_exception();
    return Polynomial();
  }
  
  //! Must be specialized by the derived class
  virtual Handle inverse_transform_isolator() const {
    throw Virtual_method_exception();
    return Handle();
  }

  //! returns whether roots are already isolated
  bool is_isolated() const {
    return _m_is_isolated;
  }
  
  //! access to tree
  Bitstream_tree tree() const {
    return _m_bitstream_tree;
  }

  //! access to kernel
  Bitstream_coefficient_kernel bck() const {
    return _m_traits.bck();
  }

  //! access to traits
  Bitstream_descartes_rndl_tree_traits traits() const {
    return _m_traits;
  }

  //! type to distinguish used constructor
  Bitstream_descartes_type type_;
   
protected:
  
 //! Polynomial which is isolated
  Polynomial _m_f;
  
  //! The tree of the Bitstream Descartes method
  mutable Bitstream_tree _m_bitstream_tree;

  //! The traits class
  Bitstream_descartes_rndl_tree_traits _m_traits;
    
  //! The number of detected isolating intervals
  int number_of_intervals;
  
  //! Has isolation already taken place
  mutable bool _m_is_isolated;
  
};

/*
 * \brief Representation for square free polynomials
 */
template<typename BitstreamCoefficientKernel,
      typename Policy=CGAL::Handle_policy_no_union>
class Square_free_descartes_rep 
  : public Generic_descartes_rep<BitstreamCoefficientKernel> {
  
  
public:
  
  //! first template parameter: Bitstream coefficient kernel
  typedef BitstreamCoefficientKernel Bitstream_coefficient_kernel;

  //! The traits class for approximations
  typedef CGAL::internal::Bitstream_descartes_rndl_tree_traits< Bitstream_coefficient_kernel >
  Bitstream_descartes_rndl_tree_traits;

  //! The generic representation
  typedef Generic_descartes_rep< Bitstream_coefficient_kernel, Policy > Base;
  
  //! Polynomial type
  typedef typename Base::Polynomial Polynomial;
  
  //! Iterator for the leaves in the bitstream tree
  typedef typename Base::Node_iterator Node_iterator;
  
  //! The type of the tree that controls the Bitstream instance
  typedef typename Base::Bitstream_tree Bitstream_tree;
  
  /*! 
   * \brief Constructor with the square free polynomial <tt>f<tt>.
   */
  Square_free_descartes_rep(
      Polynomial f,
      Bitstream_coefficient_kernel bck) :
    Base(SQUARE_FREE_DESCARTES, f, bck) {
  }

  /*! 
   * \brief Constructor with the square free polynomial <tt>f<tt>.
   */
  Square_free_descartes_rep(
      Polynomial f,
      Bitstream_tree tree,
      Bitstream_coefficient_kernel bck) :
    Base(SQUARE_FREE_DESCARTES, f, tree, bck) {
  }

  //! Needed for reference counting
  virtual CGAL::Reference_counted_hierarchy<>* clone() {
    return new Square_free_descartes_rep(*this);
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
	
  //! nothing to do here	
  virtual void process_nodes() {
    return;
  }
	
  //! Polynomial is square free, so gcd is 1
  virtual int degree_of_gcd() const {
    return 0;
  }
        
  //! Polynomial is square free
  virtual Polynomial square_free_part() const {
    return this->_m_f;
  }

  //! Always true
  virtual bool is_certainly_simple_root(int ) const {
    return true;
  }
	
  //! Always false
  virtual bool is_certainly_multiple_root(int ) const {
    return false;
  }
      
};

/*
 * \brief Representation for polynomials with at most one multiple root
 */
template<typename BitstreamCoefficientKernel,
      typename Policy=CGAL::Handle_policy_no_union>
class M_k_descartes_rep 
  : public Generic_descartes_rep<BitstreamCoefficientKernel> {
	
	
public:
	
  //! first template parameter: Bitstream coefficient kernel
  typedef BitstreamCoefficientKernel Bitstream_coefficient_kernel;

  //! The traits class for approximations
  typedef CGAL::internal::Bitstream_descartes_rndl_tree_traits< Bitstream_coefficient_kernel >
  Bitstream_descartes_rndl_tree_traits;

  //! The generic representation
  typedef Generic_descartes_rep< Bitstream_coefficient_kernel, Policy > Base;

  //! Polynomial type
  typedef typename Base::Polynomial Polynomial;

  //! Iterator for the leaves of the Bitstream Descartes tree
  typedef typename Base::Node_iterator Node_iterator;

  //! Constant iterator for the leaves
  typedef typename Base::Node_const_iterator Node_const_iterator;

  //! The interval boundaries are represented in this type
  typedef typename  Bitstream_coefficient_kernel::Bound Bound;

  //! The type of the tree that controls the Bitstream instance
  typedef typename Base::Bitstream_tree Bitstream_tree;

  /*!
   * \brief Constructor for a polynomial <tt>f<tt>, not necessarily square
   * free
   *
   * The values <tt>m</tt>
   * and <tt>k</tt> need to be the exact number of real roots of <tt>f</tt>
   * counted without multiplicity and the degree of the greatest common
   * divisor of <tt>f</tt> with its partial derivative, respectively.
   */ 
  M_k_descartes_rep(Polynomial f,int m, int k,
                    Bitstream_coefficient_kernel bck) :
    Base(M_K_DESCARTES, f, bck),
    number_of_roots(m),
    gcd_degree(k),
    index_of_multiple(-1) {
  }

  M_k_descartes_rep(Polynomial f,int m, int k,
                    Bitstream_tree tree,
                    Bitstream_coefficient_kernel bck) :
    Base(M_K_DESCARTES, f, tree, bck),
    number_of_roots(m),
    gcd_degree(k),
    index_of_multiple(-1) {
  }

  //! Default constructor
  M_k_descartes_rep() { 
  }

  //! Needed for reference counting
  virtual CGAL::Reference_counted_hierarchy<>* clone() {
    return new M_k_descartes_rep(*this);
  }

  /*!
   * \brief Termination condition
   *
   * If <tt>m-1</tt> simple and one more leaf is detected, the Bitstream
   * Descartes method is stopped. If the minimal sign
   * variation drops under <tt>k</tt> in each leaf, a
   * \c Non_generic_position_exception is thrown.
   */
  virtual bool termination_condition() {
    int counted_simple_roots=0;
    int max_max_var = 0; 
    for(Node_iterator curr=Base::_m_bitstream_tree.begin();
        curr != Base::_m_bitstream_tree.end(); curr++) {
      int max_var = Base::_m_bitstream_tree.max_var(curr);
      if(max_var > max_max_var) {
        max_max_var = max_var;
      }
      if(max_var == 1) { // && Base::_m_bitstream_tree.max_var(curr)==1) {
        ++counted_simple_roots;
      }
    }      
    //AcX_DSTREAM("Situation: " << this->number_of_intervals << " intervals " << this->number_of_roots << " are expected" << std::endl);
    if (this->number_of_intervals == this->number_of_roots 
        && counted_simple_roots >= number_of_roots-1) {
      return true;
    }
    if (max_max_var <= gcd_degree) {
      throw CGAL::internal::Non_generic_position_exception();
    }
        
    return false;

  }

  //! The index of the (possibly) multiple root is computed here.
  virtual void process_nodes() {
    int i = 0;
    for (Node_iterator curr=Base::_m_bitstream_tree.begin();
        curr != Base::_m_bitstream_tree.end(); curr++) {
      if(Base::_m_bitstream_tree.max_var(curr) > 1 ) {
        index_of_multiple = i;
        return;
      } else {
        ++i;
      }
    }
    return;
  }

  //! Returns k
  virtual int degree_of_gcd() const {
    return gcd_degree;
  }
        
  //! True for all roots except for the candidate
  virtual bool is_certainly_simple_root(int i) const {
    return (i!=index_of_multiple);
  }
	
  //! Always false
  virtual bool is_certainly_multiple_root(int) const {
    return false;
  }
      

protected:
    
  //! The "m"
  int number_of_roots;

  //! The "k"
  int gcd_degree;

  //! The candidate's index
  int index_of_multiple;

};

    
template<typename BitstreamCoefficientKernel,
             typename EventRefinement,
             typename Policy=CGAL::Handle_policy_no_union>
class Backshear_descartes_rep 
  : public Generic_descartes_rep<BitstreamCoefficientKernel> {
                                   

public:

  //! first template parameter: Bitstream coefficient kernel
  typedef BitstreamCoefficientKernel Bitstream_coefficient_kernel;

  //! second template parameter: Event_refinement
  typedef EventRefinement Event_refinement;

  //! The traits class for approximations
  typedef CGAL::internal::Bitstream_descartes_rndl_tree_traits< Bitstream_coefficient_kernel >
  Bitstream_descartes_rndl_tree_traits;

  //! The generic representation
  typedef Generic_descartes_rep< Bitstream_coefficient_kernel, Policy > Base;

  typedef typename Base::Polynomial Polynomial;

  typedef typename Base::Node_iterator Node_iterator;

  typedef std::list<int>::iterator Marking_iterator;

  typedef std::list<int>::const_iterator Marking_const_iterator;

  typedef typename Base::Node_const_iterator Node_const_iterator;

  typedef typename Base::Bound Bound;

  Backshear_descartes_rep(
      Polynomial f,
      int number_of_non_event_points,
      int number_of_events,
      Event_refinement event_refinement,
      Bitstream_coefficient_kernel bck) :
    Base(BACKSHEAR_DESCARTES, f, bck), 
    number_of_non_event_points(number_of_non_event_points),
    number_of_events(number_of_events),
    event_refinement(event_refinement) {
  }
          
  Backshear_descartes_rep() {
  }

  virtual CGAL::Reference_counted_hierarchy<>* clone() {
    return new Backshear_descartes_rep(*this);
  }

  virtual void isolate() {
      
    Node_iterator curr = Base::_m_bitstream_tree.begin(),sub_begin,new_curr;

    if(curr == Base::_m_bitstream_tree.end()) {
      this->_m_is_isolated = true;
      return;
    }
    markings.clear();
    markings.push_back(this->check_marking(curr));

    Marking_iterator curr_mark = markings.begin(),mark_helper;

    int newly_created;

    while(! this->termination_condition()) {
      //AcX_DSTREAM("Subdivision..." << number_of_intervals << std::endl);
      if (curr == Base::_m_bitstream_tree.end()) {
        curr = Base::_m_bitstream_tree.begin();
        CGAL_assertion(curr_mark == markings.end());
        curr_mark = markings.begin();
      }
      if(Base::_m_bitstream_tree.max_var(curr) == 1) {
        ++curr;
        ++curr_mark;
        //AcX_DSTREAM("nothing happend" << std::endl);
      }
      else {
        newly_created = 
          Base::_m_bitstream_tree.subdivide(curr,sub_begin,new_curr);
        mark_helper = markings.erase(curr_mark);
        curr_mark = mark_helper;
        for(Node_iterator tmp_curr = sub_begin;
            tmp_curr != new_curr;
            tmp_curr++) {
          markings.insert(curr_mark,check_marking(tmp_curr));
        }
        Base::number_of_intervals += newly_created-1;
        curr = new_curr;
        //AcX_DSTREAM(newly_created << " new intervals, marking size: " << markings.size() << std::endl);
        
      }
    }
    this->process_nodes();
    this->_m_is_isolated = true;
  }



  virtual bool termination_condition() {
    int marked_intervals = 0;
    int unmarked_odd_intervals = 0;
    Node_iterator curr = Base::_m_bitstream_tree.begin();
    Marking_iterator curr_mark = markings.begin();
    for(;curr != Base::_m_bitstream_tree.end(); curr++) {
      if((*curr_mark) >= 0) {
        ++marked_intervals;
      }
      else {
        if (Base::_m_bitstream_tree.min_var(curr) % 2 == 1) { // odd
          ++unmarked_odd_intervals;
        }
      }
      ++curr_mark;
    }
    CGAL_assertion(curr_mark == markings.end());
    return ((marked_intervals == number_of_events) 
            && (unmarked_odd_intervals == number_of_non_event_points));
  }

  virtual void process_nodes() {
    Node_iterator curr=Base::_m_bitstream_tree.begin(),curr_helper;
    Marking_iterator curr_mark = markings.begin();
    while(curr!=Base::_m_bitstream_tree.end()) {
      if(((*curr_mark) == -1) && 
         (Base::_m_bitstream_tree.min_var(curr) % 2 == 0)) {
        ++curr;
        curr_helper = curr;
        curr_helper--;
        Base::_m_bitstream_tree.erase(curr_helper);
        curr_mark = markings.erase(curr_mark);
        Base::number_of_intervals--;
      } else {
        ++curr_mark;
        ++curr;
      }
    }
    CGAL_assertion(curr_mark == markings.end());

    //AcX_DSTREAM(markings.size() << " " << number_of_non_event_points << " " << number_of_events << std::endl);
    CGAL_assertion(static_cast<int>(markings.size())
                   ==number_of_non_event_points + number_of_events);
    return;
  }

  virtual bool is_certainly_simple_root(int i) const {
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < Base::number_of_intervals);
    Node_const_iterator curr=Base::_m_bitstream_tree.begin();
    std::advance(curr,i);
    return (Base::_m_bitstream_tree.max_var(curr) == 1);
  }
	
  virtual bool is_certainly_multiple_root(int i) const {
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < Base::number_of_intervals);
    Marking_const_iterator curr = markings.begin();
    std::advance(curr,i);
    return (*curr>=0);
  }
          

protected:

  int number_of_non_event_points;

  int number_of_events;

  Event_refinement event_refinement;

  std::list<int> markings;

protected:

  int check_marking(Node_iterator node) {
    Bound lower = Base::_m_bitstream_tree.lower(node),
      upper = Base::_m_bitstream_tree.upper(node);
    for(int i = 0; i < number_of_events; i++) {
      while(true) {
        if(CGAL::compare(event_refinement.lower_bound(i),lower)
           !=CGAL::NEGATIVE
           && 
           CGAL::compare(event_refinement.upper_bound(i),upper)
           !=CGAL::POSITIVE) {
          //Event inside the interval
          return i;
        }
        if(CGAL::compare(event_refinement.lower_bound(i),upper)
           ==CGAL::POSITIVE
           ||
           CGAL::compare(event_refinement.upper_bound(i),lower)
           ==CGAL::NEGATIVE) {
          //This event is outside
          break;
        }
        event_refinement.refine(i);
              
      }
    }
    return -1;
  }
          
};

/*
 * \brief Adaptor for roots of a vert line 
 * (needed as dummy in surface analysis)
 *
 */
template<typename BitstreamCoefficientKernel,
      typename VertLine,
      typename Policy=CGAL::Handle_policy_no_union>
class Vert_line_adapter_descartes_rep 
  : public Generic_descartes_rep<BitstreamCoefficientKernel> {
        
public:
        
  //! first template parameter: Bitstream coefficient kernel
  typedef BitstreamCoefficientKernel Bitstream_coefficient_kernel;

  //! second template parameter: Vert_line
  typedef VertLine Vert_line;

  //! The traits class for approximations
  typedef CGAL::internal::Bitstream_descartes_rndl_tree_traits< Bitstream_coefficient_kernel >
  Bitstream_descartes_rndl_tree_traits;

  //! The generic representation
  typedef Generic_descartes_rep< Bitstream_coefficient_kernel, Policy > Base;
        
  //! type of Curve_analysis_2
  typedef typename Vert_line::Curve_analysis_2 Curve_analysis_2;

  //! type of Curve_kernel_2;
  typedef typename Curve_analysis_2::Algebraic_kernel_with_analysis_2 
  Curve_kernel_2;

  //! The Coeeficient type of the input polynomial
  typedef typename Bitstream_coefficient_kernel::Coefficient Coefficient;
        
  //! The polynomial type
  typedef typename CGAL::Polynomial_type_generator<Coefficient,1>::Type 
  Polynomial;
        
  typedef Vert_line_adapter_descartes_rep
  <Bitstream_coefficient_kernel, Vert_line, Policy> Self;
        
  //! The used integer type
  typedef typename Bitstream_coefficient_kernel::Integer Integer;
        
  //! How the boundaries of the isolating intervals are represented
  typedef typename Bitstream_coefficient_kernel::Bound Bound;
        
  //! The type for the inverse isolator
  typedef typename Base::Handle Handle;
        
  /*! 
   * \brief Constructor
   */
  template<typename InputIterator>
  Vert_line_adapter_descartes_rep(InputIterator begin,
                                  InputIterator end,
                                  Bitstream_coefficient_kernel bck)
    : Base(VERT_LINE_ADAPTER_DESCARTES) 
  {
    for (InputIterator it = begin; it != end; it++) {
      root_vec.push_back(std::make_pair(*it, 4));
    }
	    
    this->_m_is_isolated = true;
    this->_m_traits = Bitstream_descartes_rndl_tree_traits(bck);
    this->_m_f  = Polynomial(0);
    this->number_of_intervals 
      = static_cast<int>(root_vec.size());
    // Isolate all real roots until intervals are disjoint:
    for (int i = 1; i < this->number_of_real_roots(); i++ ){
      while(left_bound(i) < right_bound(i-1) ) {
        if (right_bound(i)-left_bound(i) < 
            right_bound(i-1) - left_bound(i-1) ) {
          refine_interval(i-1);
        } else {
          refine_interval(i);
        }
      }

    }
  }
        
        
  //! Destructor (does nothing)
  virtual ~Vert_line_adapter_descartes_rep() {
  }

  //! Needed for the referencing counting mechanism
  virtual CGAL::Reference_counted_hierarchy<>* clone() {
    return new Vert_line_adapter_descartes_rep(*this);
  }
        
  virtual void refine_interval(int i) const {
    root_vec[i] = std::make_pair(root_vec[i].first, root_vec[i].second * 2);
  }
        
  virtual void isolate() const {
  }

 
  //! The lower bound of the \c i th root
  virtual Bound left_bound(int i) const  {
    typename Curve_kernel_2::Approximate_absolute_y_2
      approx_y = this->bck().point().xy().kernel()
      ->approximate_absolute_y_2_object();
    return approx_y(root_vec[i].first.first.algebraic_real_2
                    (root_vec[i].first.second),
                    root_vec[i].second).first;
  }
    
  //! The upper bound of the \c i th root
  virtual Bound right_bound(int i) const {
    typename Curve_kernel_2::Approximate_absolute_y_2
      approx_y = this->bck().point().xy().kernel()
      ->approximate_absolute_y_2_object();
    return approx_y(root_vec[i].first.first.algebraic_real_2
                    (root_vec[i].first.second),
                    root_vec[i].second).second;
  }


  /*! \brief Returns whether the \c i th root is definitely a simple root
   * of the isolated polynomial
   *
   */
  virtual bool is_certainly_simple_root(int /* i */) const {
    return false;
  }
        
  /*! \brief Returns whether the \c i th root is definitely 
   * a multiple root
   * of the isolated polynomial
   *
   */
  virtual bool is_certainly_multiple_root(int /* i */) const {
    return false;
  }

protected:

  //! Roots stored as pair of a AcX::Vert_line and an integer denoting the
  //! index. Also, current precision of each root is stored
  mutable std::vector<std::pair<std::pair<Vert_line, int>,int> > root_vec;


};

/*!
 * \brief Class for the Bitstream Descartes method
 *
 * Class for the real root isolation of polynomials, using the Bitstream
 * Descartes method. The polynomials coefficient type is arbitrary, the 
 * approximations of the coefficient type are obtained with the 
 * \c BitstreamCoefficientKernel parameter. For the requirements
 * of this traits class, see the documentation of the traits.

 *
 * Internally, an instance of CGAL::Bitstream_descartes_rndl_tree is explored
 * in a specific way. That exploration strategy depends on the constructor
 * that is used to create the object. A tag is passed that defines the
 * variant of the Bitstream Descartes method: The Square_free_descartes_tag
 * starts the usual Bitstream method for square free integer polynomials.
 * With the M_k_descartes tag, it is able to handle one multiple root in 
 * favourable situations, the Backshear_descartes_tag allows to isolate
 * even more complicated polynomials, if the multiple roots with even
 * multiplicity can be refined from outside. See the corresponding
 * constructors for more information.
 * 
 */
template<typename BitstreamCoefficientKernel>
class Bitstream_descartes 
  : ::CGAL::Handle_with_policy<
    CGAL::internal::Generic_descartes_rep<BitstreamCoefficientKernel> > {
    
public:
    
  //! Bitstream coefficient kernel
  typedef BitstreamCoefficientKernel Bitstream_coefficient_kernel;

  //! the traits class
  typedef CGAL::internal::Bitstream_descartes_rndl_tree_traits< Bitstream_coefficient_kernel >
  Bitstream_descartes_rndl_tree_traits;
      
  // The generic representation class
  typedef 
  CGAL::internal::Generic_descartes_rep< Bitstream_coefficient_kernel > Rep;

  // The Handle type
  typedef ::CGAL::Handle_with_policy<Rep> Base;
    
  //! The coefficients of the polynomial
  typedef typename Bitstream_coefficient_kernel::Coefficient Coefficient;
    
  //! The polynomial's type
  typedef typename CGAL::Polynomial_type_generator<Coefficient,1>::Type 
  Polynomial;
    
  typedef Bitstream_descartes<Bitstream_coefficient_kernel> Self;
    
  // Type for the Bitstream Descartes tree
  typedef typename Rep::Bitstream_tree Bitstream_tree;

  //! Type for Integers
  typedef typename Bitstream_coefficient_kernel::Integer Integer;

  //! Type for the interval boundaries of the isolating intervals
  typedef typename Bitstream_coefficient_kernel::Bound Bound;

  //! Iterator type for the leaves of the Descartes tree
  typedef typename Bitstream_tree::Node_iterator Node_iterator;

  //! Const iterator for the leaves
  typedef typename Bitstream_tree::Node_const_iterator Node_const_iterator;

  //! Default constructor
  Bitstream_descartes() : Base(new Rep()) {}

  //! Copy constructor
  Bitstream_descartes(const Self& other) : Base(static_cast<const Base&>(other))
  {}
  
  /*! 
   * \brief Constructor for a polynomial \c f
   *
   * See the documentation of the constrctor 
   * with \c Square_free_descartes_tag
   */
  Bitstream_descartes(Polynomial f,
                      Bitstream_coefficient_kernel bck = Bitstream_coefficient_kernel(),
                      bool isolate = true)
    : Base(new CGAL::internal::Square_free_descartes_rep
           <Bitstream_coefficient_kernel>(f, bck))
  {
    if (isolate) {
      this->isolate();
    }
  }

  /*! 
   * \brief Constructor for the square free Descartes method
   *
   * The polynomial \c f must not have multiple real roots. The 
   * Bitstream Descartes tree is traversed in a bfs manner until
   * all leaves have sign variation zero or one.
   */
  Bitstream_descartes(Square_free_descartes_tag ,
                      Polynomial f,
                      Bitstream_coefficient_kernel bck = Bitstream_coefficient_kernel(),
                      bool isolate = true)
    : Base(new CGAL::internal::Square_free_descartes_rep
           <Bitstream_coefficient_kernel>(f, bck))
  {
    if (isolate) {
      this->isolate();
    }
  }

  /*! 
   * \brief Constructor for the square free Descartes method,
   * using a precomputed tree
   *
   * The polynomial \c f must not have multiple real roots. The 
   * Bitstream Descartes tree is traversed in a bfs manner until
   * all leaves have sign variation zero or one.
   * The tree must be adequate for the polynomial. 
   * Use that constructor only if you know what you're doing!
   */
  Bitstream_descartes(Square_free_descartes_tag ,
                      Polynomial f,
                      Bitstream_tree tree,
                      Bitstream_coefficient_kernel bck = Bitstream_coefficient_kernel(),
                      bool isolate = true)
    : Base(new CGAL::internal::Square_free_descartes_rep
           <Bitstream_coefficient_kernel>(f, tree, bck))
  {
    if (isolate) {
      this->isolate();
    }
  }

  /*! 
   * \brief Constructor for the m-k-Descartes method
   *
   * The polynomial \c f must have exactly \c m real roots, counted without
   * multiplicity, and the degree of <tt>gcd(f,f')</tt> must be \c k. In this
   * case, the constructor either isolates the real roots of \c f sucessfully
   * or a Non_generic_position_exception is thrown. Such an exception
   * certainly occurs if \c f has more than one multiple real root. If \c f
   * has at most one multiple root over the complex numbers, the roots are
   * certainly isolated with success.
   */
  Bitstream_descartes(M_k_descartes_tag /* t */,
                      Polynomial f,int m,int k,
                      Bitstream_coefficient_kernel bck = Bitstream_coefficient_kernel(),
                      bool isolate = true)
    : Base(new CGAL::internal::M_k_descartes_rep
           <Bitstream_coefficient_kernel>(f, m, k, bck))
  {
    if (isolate) {
      this->isolate();
    }
  }


  Bitstream_descartes(M_k_descartes_tag /* t */,
                      Polynomial f,int m,int k,
                      Bitstream_tree tree,
                      Bitstream_coefficient_kernel bck = Bitstream_coefficient_kernel(),
                      bool isolate = true)
    : Base(new CGAL::internal::M_k_descartes_rep
           <Bitstream_coefficient_kernel>(f, m, k, tree, bck))
  {
    if (isolate) {
      this->isolate();
    }
  }
  

  /*! 
   * \brief Constructor for the Backshear-Decartes method
   * 
   * The polynomial \c f must have exactly \c number_of_real_roots
   * many real roots, counted without multiplicity. Additionally, a set of
   * \c number_of_events root can be refined to arbitrary precision with the
   * \c event_refinement object. This must support three operations
   * for each <tt>0<=i<number_of_events</tt>:
   * <ul><li>lower_bound(i), upper_bound(i) gives an interval (not
   * necessarily isolating) of some root of \c f</li>
   * <li>refine(i) refines the corresponding interval</li></ul>
   * Note that the roots in \c event_refinement need not be sorted. All roots
   * which are not covered by \c event_refinement must have odd multiplicity.
   */
  template<typename EventRefinement>
  Bitstream_descartes(Backshear_descartes_tag ,
                      Polynomial f,
                      int number_of_real_roots,
                      int number_of_events,
                      EventRefinement event_refinement,
                      Bitstream_coefficient_kernel bck = Bitstream_coefficient_kernel(),
                      bool isolate = true)
    : Base(new 
           CGAL::internal::Backshear_descartes_rep
           <Bitstream_coefficient_kernel, EventRefinement>
           (f, number_of_real_roots-number_of_events,
            number_of_events, event_refinement, bck))
  {
    if (isolate) {
      this->isolate();
    }
  }


  /*! 
   * \brief Constructor for the Vert-line-adapter-Descartes method
   *
   */
  template<typename InputIterator>
  Bitstream_descartes(Vert_line_adapter_descartes_tag /* t */,
                      InputIterator begin,
                      InputIterator end,
                      Bitstream_coefficient_kernel bck)
    : Base(new CGAL::internal::Vert_line_adapter_descartes_rep
           <Bitstream_coefficient_kernel,
           typename InputIterator::value_type::first_type>
           (begin, end, bck) )
  {
    // No isolation necessary
  }

  //! return the type of the used descartes method
  Bitstream_descartes_type type() const {
    return this->ptr()->type_;
  }

  //! Return the polynomial
  Polynomial polynomial() const {
    CGAL_assertion(is_isolated());
    return this->ptr()->polynomial();
  }

  //! Returns the kernel class
  Bitstream_coefficient_kernel bck() const {
    return this->ptr()->bck();
  } 

  //! Returns the traits class
  Bitstream_descartes_rndl_tree_traits traits() const {
    return this->ptr()->traits();
  } 

  //! Number of real roots of the polynomial
  int number_of_real_roots() const {
    CGAL_assertion(is_isolated());
    return this->ptr()->number_of_real_roots();
  }

  //! Refine the <tt>i</tt>th isolating interval
  void refine_interval(int i) const {
    CGAL_assertion(is_isolated());
    this->ptr()->refine_interval(i);
  }

  //! The left bound of the <tt>i</tt>th isolating interval
  Bound left_bound(int i) const  {
    CGAL_assertion(is_isolated());
    return this->ptr()->left_bound(i);
  }

  //! The left bound of the <tt>i</tt>th isolating interval
  void left_bound(int i, 
                  Integer& numerator, 
                  Integer& denominator) const {
    typedef CGAL::Fraction_traits<Bound> Fraction_traits; 
    typename Fraction_traits::Decompose decompose;
    decompose(left_bound(i),numerator,denominator);
  }

  //! The right bound of the <tt>i</tt>th isolating interval
  Bound right_bound(int i) const  {
    CGAL_assertion(is_isolated());
    return this->ptr()->right_bound(i);
  }

  //! The right bound of the <tt>i</tt>th isolating interval
  void right_bound(int i, 
                   Integer& numerator, 
                   Integer& denominator) const {
    typedef CGAL::Fraction_traits<Bound> Fraction_traits; 
    typename Fraction_traits::Decompose decompose;
    decompose(right_bound(i),numerator,denominator);
  }

  //! The length of the <tt>i</tt>th isolating interval
  Bound length(int i) const {
    CGAL_assertion(is_isolated());
    return (this->ptr()->right_bound(i) - 
            this->ptr()->left_bound(i));
  }

  bool is_exact_root(int) const { return false; }

  /*! 
   * \brief Returns true if the <tt>i</tt>th root is known to be a simple 
   * root of the curve.
   */
  bool is_certainly_simple_root(int i) const {
    CGAL_assertion(is_isolated());
    return this->ptr()->is_certainly_simple_root(i);
  }
	
  /*! 
   * \brief Returns true if the <tt>i</tt>th root is known to be a multiple 
   * root of the curve.
   */
  bool is_certainly_multiple_root(int i) const {
    CGAL_assertion(is_isolated());
    return this->ptr()->is_certainly_multiple_root(i);
  }

        
  /*! 
   * \brief Returns the multiplicity of the root if know, otherwise -1
   */
  int multiplicity_of_root(int i) const {
    CGAL_assertion(is_isolated());
    return this->ptr()->multiplicity_of_root(i);
  }

  /*!
   * Returns an upper bound for the multiplicity of the ith root
   */
  int upper_bound_for_multiplicity(int i) const {
    CGAL_assertion(is_isolated());
    return this->ptr()->upper_bound_for_multiplicity(i);
  }

  /*!
   * \brief Returns the isolator of the polynomial f(1/x + q), if known
   */
  Self inverse_transform_isolator() const {
    return this->ptr()->inverse_transform_isolator();
  }


public:
  
  //! Starts the isolation of the real roots.
  void isolate() {
    CGAL_assertion(!is_isolated());
    this->ptr()->isolate();
  }
  
  //! returns whether all roots are isolated
  bool is_isolated() const {
    return this->ptr()->is_isolated();
  }

  //! returns the underlying tree
  Bitstream_tree tree() const {
    return this->ptr()->tree();
  }
  
  //! returns the degree of the gcd of f and its derivative, if known
  int degree_of_gcd() const {
    return this->ptr()->degree_of_gcd();
  }
  
  //! returns the square free part of f, if known
  Polynomial square_free_part() const {
    return this->ptr()->square_free_part();
  }
  
};


} // namespace internal

} // namespace CGAL

#endif // CGAL_AK_USE_OLD_BITSTREAM_DESCARTES

#endif // CGAL_BITSTREAM_DESCARTES_H
// EOF
