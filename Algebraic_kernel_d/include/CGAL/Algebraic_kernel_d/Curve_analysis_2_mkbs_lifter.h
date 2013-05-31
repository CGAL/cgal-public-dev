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
// Author(s)     :  Michael Kerber <mkerber@mpi-inf.mpg.de>
//                  Eric Berberich <eric.berberich@cgal.org>
//
// ============================================================================

/*! \file CGAL/Curve_analysis_2_mkbs_lifter.h
  \brief Defines isolator for lifting x-coordinates in curve analyses

  mkbs can be explained in the following 4 ways: 
  "mk" 
    - m number of real roots, k degree of local gcd
    - Michael Kerber (and mainly his implementation)
  "bs"
    - backshear (as it requires backshearing)
    - bitstream
*/

#ifndef CGAL_ALGEBRAIC_KERNEL_D_CURVE_ANALYSIS_2_MKBS_LIFTER_H
#define CGAL_ALGEBRAIC_KERNEL_D_CURVE_ANALYSIS_2_MKBS_LIFTER_H

#include <CGAL/config.h>

#include <CGAL/tags.h>

#include <CGAL/Algebraic_kernel_d/Generic_isolator.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_isolator_reps.h>

namespace CGAL {

namespace internal {

/*
 * \brief Representation for polynomials with at most one multiple root
 */
template< class BitstreamCoefficientKernel, typename HandlePolicy > // no default on policy, should be decided in higher level
class Mk_bitstream_descartes_isolator_rep : 
  public Generic_bitstream_descartes_isolator_rep< BitstreamCoefficientKernel, HandlePolicy > {
	
    // tags are inherited

public:
	
  //! first template parameter: Bitstream coefficient kernel
  typedef BitstreamCoefficientKernel Bitstream_coefficient_kernel;

  //! second template parameter: Handle_policy
  typedef HandlePolicy Handle_policy;

  //! The generic representation
  typedef Generic_bitstream_descartes_isolator_rep< Bitstream_coefficient_kernel, Handle_policy > Base;

  //! the class itself
  typedef Mk_bitstream_descartes_isolator_rep< Bitstream_coefficient_kernel, Handle_policy > Self;

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

  //!\name Constructors
  //!@{

#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_DEFAULT_CONSTRUCTIBLE
  private:
    Mk_bitstream_descartes_isolator_rep();
#else
  public:
    //! Default constructor
    Mk_bitstream_descartes_isolator_rep() : Base() { 
    }
#endif

#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_ASSIGNABLE
  protected:
    // Mk_bitstream_descartes_isolator_rep(const Self&) // = default
    Mk_bitstream_descartes_isolator_rep& operator=(const Self&); // = disable
#endif

  public:
  /*!
   * \brief Constructor for a polynomial <tt>f<tt>, not necessarily square
   * free
   *
   * The values <tt>m</tt>
   * and <tt>k</tt> need to be the exact number of real roots of <tt>f</tt>
   * counted without multiplicity and the degree of the greatest common
   * divisor of <tt>f</tt> with its partial derivative, respectively.
   */ 
  Mk_bitstream_descartes_isolator_rep(const Polynomial& p, int m, int k,
                                      Bitstream_coefficient_kernel bck) :
    Base(p, bck),
    _m_number_of_roots(m),
    _m_gcd_degree(k),
    _m_index_of_multiple(-1) {
  }

#if 1 // TODO 2012 CA/CPA does not need this constructor (only _test_mkbs_lifter seems to call it)
  /*!
   * \brief Constructor for a polynomial <tt>f<tt>, not necessarily square
   * free, and a given bitstream tree
   *
   * The values <tt>m</tt>
   * and <tt>k</tt> need to be the exact number of real roots of <tt>f</tt>
   * counted without multiplicity and the degree of the greatest common
   * divisor of <tt>f</tt> with its partial derivative, respectively.
   */ 
  Mk_bitstream_descartes_isolator_rep(const Polynomial& p,int m, int k,
                                      Bitstream_tree tree,
                                      Bitstream_coefficient_kernel bck) :
    Base(p, tree, bck),
    _m_number_of_roots(m),
    _m_gcd_degree(k),
    _m_index_of_multiple(-1) {
  }
#endif

  //!@} // Constructors

  //!\name Access functions
  //!@{

  /*!
   * \brief Termination condition
   *
   * If <tt>m-1</tt> simple and one more leaf is detected, the Bitstream
   * Descartes method is stopped. If the minimal sign
   * variation drops under <tt>k</tt> in each leaf, a
   * \c Non_generic_position_exception is thrown.
   */
  virtual bool termination_condition() {
    size_t counted_simple_roots = 0;
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
    //AcX_DSTREAM("Situation: " << this->number_of_intervals << " intervals " << this->_m_number_of_roots << " are expected" << std::endl);
    if (this->_m_number_of_intervals == this->_m_number_of_roots 
        && counted_simple_roots >= _m_number_of_roots-1) {
      return true;
    }
    if (max_max_var <= _m_gcd_degree) {
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
        _m_index_of_multiple = i;
        return;
      } else {
        ++i;
      }
    }
    return;
  }

  //! Returns k
  virtual int degree_of_gcd() const {
    return _m_gcd_degree;
  } 
        
  //! True for all roots except for the candidate
  virtual bool is_certainly_simple_root(size_t i) const {
    return (static_cast<int>(i) != _m_index_of_multiple);
  }
	
  //! Always false
  virtual bool is_certainly_multiple_root(size_t /* i */) const {
    return false;
  }

  //! 
  virtual int multiplicity_of_root(size_t /* i */) const {
    bool _WARNING_DO_NOT_USE_multiplicity_of_root_which_is_not_yet_implemented;
    CGAL_error(); // TODO 2012 how to forbid access to this
    return -1;
  }


  //!@} // Access functions
      
  private:

  //! Needed for reference counting
  virtual CGAL::Reference_counted_hierarchy<>* clone() {
    return new Mk_bitstream_descartes_isolator_rep(*this);
  }

protected:
    
  //!\name Private members
  //!@{

  //! The "m"
  size_t _m_number_of_roots;

  //! The "k"
  int _m_gcd_degree;

  //! The candidate's index
  mutable int _m_index_of_multiple;
   
  //!@} // Private members      
}; // Mk_bitstream_descartes_isolator_rep

  // TASK replace BitstreamDescartes in Backshear_descartes with an "active interval" to allow arvavoid as alternative, 
  //      ie get rid of Generic_bitstream_descartes_rep-"base"
template< class BitstreamCoefficientKernel, typename EventRefinement, typename HandlePolicy > // no default on policy, should be decided in higher level
class Backshear_isolator_rep : 
  public Generic_bitstream_descartes_isolator_rep< BitstreamCoefficientKernel, HandlePolicy > {

    // tags are inherited

public:

  //!\name Public types
  //!@{

  //! first template parameter: Bitstream coefficient kernel
  typedef BitstreamCoefficientKernel Bitstream_coefficient_kernel;

  //! second template parameter: Event_refinement
  typedef EventRefinement Event_refinement;

  //! third template parameter: Handle_policy
  typedef HandlePolicy Handle_policy;

  //! The generic representation
  typedef Generic_bitstream_descartes_isolator_rep< Bitstream_coefficient_kernel, Handle_policy > Base;

  //! the class itself
typedef Backshear_isolator_rep< Bitstream_coefficient_kernel, Event_refinement, Handle_policy > Self;

  typedef typename Base::Polynomial Polynomial;

  typedef typename Base::Node_iterator Node_iterator;

  typedef std::list<int>::iterator Marking_iterator;

  typedef std::list<int>::const_iterator Marking_const_iterator;

  typedef typename Base::Node_const_iterator Node_const_iterator;

  typedef typename Base::Bound Bound;

  //!@} // Public types

  //!\name Constructors
  //!@{

#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_DEFAULT_CONSTRUCTIBLE
  private:
    Backshear_isolator_rep();
#else
  public:
    //! Default constructor
    Backshear_isolator_rep() : Base() { 
    }
#endif

#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_ASSIGNABLE
  protected:
    //Backshear_isolator_rep(const Self&) // = default
    Backshear_isolator_rep& operator=(const Self&); // = disable
#endif

  public:
  Backshear_isolator_rep(const Polynomial& p,
                         int number_of_non_event_points,
                         int number_of_events,
                         Event_refinement event_refinement,
                         Bitstream_coefficient_kernel bck) :
    Base(p, bck), 
    _m_number_of_non_event_points(number_of_non_event_points),
    _m_number_of_events(number_of_events),
    _m_event_refinement(event_refinement) {
  }
          
  //!@} // Constructors

  //!\name Access functions
  //!@{

  virtual void isolate() {
      
    Node_iterator curr = Base::_m_bitstream_tree.begin(),sub_begin,new_curr;

    if(curr == Base::_m_bitstream_tree.end()) {
      this->_m_is_isolated = true;
      return;
    }
    _m_markings.clear();
    _m_markings.push_back(this->check_marking(curr));

    Marking_iterator curr_mark = _m_markings.begin(),mark_helper;

    int newly_created;

    while(! this->termination_condition()) {
      //AcX_DSTREAM("Subdivision..." << Base::_m_number_of_intervals << std::endl);
      if (curr == Base::_m_bitstream_tree.end()) {
        curr = Base::_m_bitstream_tree.begin();
        CGAL_assertion(curr_mark == _m_markings.end());
        curr_mark = _m_markings.begin();
      }
      if(Base::_m_bitstream_tree.max_var(curr) == 1) {
        ++curr;
        ++curr_mark;
        //AcX_DSTREAM("nothing happend" << std::endl);
      }
      else {
        newly_created = 
          Base::_m_bitstream_tree.subdivide(curr,sub_begin,new_curr);
        mark_helper = _m_markings.erase(curr_mark);
        curr_mark = mark_helper;
        for(Node_iterator tmp_curr = sub_begin;
            tmp_curr != new_curr;
            tmp_curr++) {
          _m_markings.insert(curr_mark,check_marking(tmp_curr));
        }
        Base::_m_number_of_intervals += newly_created-1;
        curr = new_curr;
        //AcX_DSTREAM(newly_created << " new intervals, marking size: " << _m_markings.size() << std::endl);
        
      }
    }
    this->process_nodes();
    this->_m_is_isolated = true;
  }

  virtual bool is_certainly_simple_root(size_t i) const {
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < Base::_m_number_of_intervals);
    Node_const_iterator curr=Base::_m_bitstream_tree.begin();
    std::advance(curr,i);
    return (Base::_m_bitstream_tree.max_var(curr) == 1);
  }
	
  virtual bool is_certainly_multiple_root(size_t i) const {
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < Base::_m_number_of_intervals);
    Marking_const_iterator curr = _m_markings.begin();
    std::advance(curr,i);
    return (*curr>=0);
  }

  virtual int multiplicity_of_root(size_t i) const {
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < Base::_m_number_of_intervals);
    bool _WARNING_DO_NOT_USE_multiplicity_of_root_which_is_not_yet_implemented;
    CGAL_error(); // TODO 2012 avoid call to multiplicity of root
    return -1;
  }

  // bitstream functions

  virtual bool termination_condition() {
    int marked_intervals = 0;
    int unmarked_odd_intervals = 0;
    Node_iterator curr = Base::_m_bitstream_tree.begin();
    Marking_const_iterator curr_mark = _m_markings.begin();
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
    CGAL_assertion(curr_mark == _m_markings.end());
    return ((marked_intervals == _m_number_of_events) 
            && (unmarked_odd_intervals == _m_number_of_non_event_points));
  }

  virtual void process_nodes() {
    Node_iterator curr=Base::_m_bitstream_tree.begin(),curr_helper;
    Marking_iterator curr_mark = _m_markings.begin();
    while(curr!=Base::_m_bitstream_tree.end()) {
      if(((*curr_mark) == -1) && 
         (Base::_m_bitstream_tree.min_var(curr) % 2 == 0)) {
        ++curr;
        curr_helper = curr;
        curr_helper--;
        Base::_m_bitstream_tree.erase(curr_helper);
        curr_mark = _m_markings.erase(curr_mark);
        Base::_m_number_of_intervals--;
      } else {
        ++curr_mark;
        ++curr;
      }
    }
    CGAL_assertion(curr_mark == _m_markings.end());

    //AcX_DSTREAM(_m_markings.size() << " " << _m_number_of_non_event_points << " " << _m_number_of_events << std::endl);
    CGAL_assertion(static_cast<int>(_m_markings.size())
                   ==_m_number_of_non_event_points + _m_number_of_events);
    return;
  }
          
protected:

  int check_marking(Node_iterator node) {
    Bound lower = Base::_m_bitstream_tree.lower(node),
      upper = Base::_m_bitstream_tree.upper(node);
    for(int i = 0; i < _m_number_of_events; i++) {
      while(true) {
        if(CGAL::compare(_m_event_refinement.lower_bound(i),lower)
           !=CGAL::NEGATIVE
           && 
           CGAL::compare(_m_event_refinement.upper_bound(i),upper)
           !=CGAL::POSITIVE) {
          //Event inside the interval
          return i;
        }
        if(CGAL::compare(_m_event_refinement.lower_bound(i),upper)
           ==CGAL::POSITIVE
           ||
           CGAL::compare(_m_event_refinement.upper_bound(i),lower)
           ==CGAL::NEGATIVE) {
          //This event is outside
          break;
        }
        _m_event_refinement.refine(i);
              
      }
    }
    return -1;
  }

//!@} // Access functions

private:

  // needed for Handle_with_policy
  virtual CGAL::Reference_counted_hierarchy<>* clone() {
    return new Backshear_isolator_rep(*this);
  }

private:

  //!\name Private members
  //!@{

  int _m_number_of_non_event_points;

  int _m_number_of_events;

  Event_refinement _m_event_refinement;

  std::list<int> _m_markings;
  
  //!@} // Private members
        
}; // Backshear_isolator_rep

/*
 * \brief Thrown whenever a non-specialised virtual member function is called
 */
class Virtual_method_exception {};

//! will be used for lifting over events and intervals in Curve_analysis
template< class BitstreamCoefficientKernel, typename HandlePolicy = CGAL::Handle_policy_no_union >
class Mkbs_lifter : 
  public internal::Generic_isolator< typename BitstreamCoefficientKernel::Polynomial,
                                     typename BitstreamCoefficientKernel::Bound,
                                     HandlePolicy,
                                     Generic_bitstream_descartes_isolator_rep< BitstreamCoefficientKernel, HandlePolicy > > {
  // combines Descartes/Bitstream-Decartes/RS, m-k-Descartes, Backshear

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

  //! second template parameter: Handle_policy
  typedef HandlePolicy Handle_policy;

  //! The polynomial type
  typedef typename Bitstream_coefficient_kernel::Polynomial Polynomial;

  //! How the boundaries of the isolating intervals are represented
  typedef typename Bitstream_coefficient_kernel::Bound Bound;

  //! the generic representation
  typedef Generic_bitstream_descartes_isolator_rep< Bitstream_coefficient_kernel, Handle_policy > Generic_bs_rep;

  //! the representation class for square-free case1
  // TODO 2012 pick right square-free rep or make square-free even a parameter
  typedef Square_free_bitstream_descartes_isolator_rep< Bitstream_coefficient_kernel, Handle_policy > Square_free_rep;

  //! the representation class for mk case
  typedef Mk_bitstream_descartes_isolator_rep <Bitstream_coefficient_kernel, Handle_policy > Mk_rep;

  //! the type of the tree
  typedef typename Mk_rep::Bitstream_tree Bitstream_tree;

  // ! the representation class for backshear case
  // Remark: the second template parameter is only given in the constructor, thus the following type 
  //         is not given here
  // typedef internal::Backshear_isolator_rep<Bitstream_coefficient_kernel, Event_refinement, Handle_policy > Backshear_rep;

  //! the base type
 typedef internal::Generic_isolator< Polynomial, Bound, Handle_policy, Generic_bitstream_descartes_isolator_rep< Bitstream_coefficient_kernel, Handle_policy > /* in fact, provide all of the used reps! */ > Base;

  //! the class itself
  typedef Mkbs_lifter< Bitstream_coefficient_kernel, Handle_policy > Self;

  //!@} // Public types

public:

  //!\name Constructors 
  //!@{

#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_DEFAULT_CONSTRUCTIBLE
  private:
    Mkbs_lifter();
#else
  public:
    //! Default constructor
    Mkbs_lifter() : Base(new Square_free_rep()) { 
    }
#endif

#if 0
#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_ASSIGNABLE
  protected:
    Mkbs_lifter(const Self&); // = disable
    Mkbs_lifter& operator=(const Self&); // = disable
#endif
#endif

  public:

  // TODO 2012 replace all constructors by one ctor for event and one ctor for internal
  // and let this class decide which is the best isolator (e.g. with help of modular arithmetic)

  /*! 
   * \brief Constructor for a polynomial \c f
   */
  Mkbs_lifter(const Polynomial& p,
              Bitstream_coefficient_kernel bck = Bitstream_coefficient_kernel(),
              bool isolate = true) : 
    Base(new Square_free_rep(p, bck))
  {
    if (isolate) {
      this->isolate();
    }
  }

#if 1 // TODO 2012 CA/CPA does not need this constructor (only _test_mkbs_lifter seems to call it)
  /*! 
   * \brief Constructor for a polynomial \c f and \c tree
   */
  Mkbs_lifter(const Polynomial& p,
              Bitstream_tree tree,
              Bitstream_coefficient_kernel bck = Bitstream_coefficient_kernel(),
              bool isolate = true) : 
  Base(new Square_free_rep(p, tree, bck))
  {
    if (isolate) {
      this->isolate();
    }
  }
#endif

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
  Mkbs_lifter(const Polynomial& p,
              int m, int k,
              Bitstream_coefficient_kernel bck = Bitstream_coefficient_kernel(),
              bool isolate = true) : 
    Base(new Mk_rep(p, m, k, bck))
  {
    if (isolate) {
      this->isolate();
    }
  }

  // TODO 2012 CA/CPA does not need this constructor (only _test_mkbs_lifter seems to call it)
  /*! 
   * \brief Constructor for the m-k-Descartes method and a given bitstream tree
   *
   * The polynomial \c f must have exactly \c m real roots, counted without
   * multiplicity, and the degree of <tt>gcd(f,f')</tt> must be \c k. In this
   * case, the constructor either isolates the real roots of \c f sucessfully
   * or a Non_generic_position_exception is thrown. Such an exception
   * certainly occurs if \c f has more than one multiple real root. If \c f
   * has at most one multiple root over the complex numbers, the roots are
   * certainly isolated with success.
   */
  Mkbs_lifter(const Polynomial& p,
              int m, int k,
              Bitstream_tree tree,
              Bitstream_coefficient_kernel bck = Bitstream_coefficient_kernel(),
              bool isolate = true) : 
    Base(new Mk_rep(p, m, k, tree, bck))
  {
    if (isolate) {
      this->isolate();
    }
  }

  /*! 
   * \brief Constructor for the Backshear method
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
  template< class EventRefinement >
  Mkbs_lifter(const Polynomial& p,
              int number_of_real_roots,
              int number_of_events,
              EventRefinement event_refinement,
              Bitstream_coefficient_kernel bck = Bitstream_coefficient_kernel(),
              bool isolate = true) : 
 Base(new Backshear_isolator_rep< Bitstream_coefficient_kernel, EventRefinement, Handle_policy >(p, number_of_real_roots-number_of_events, number_of_events, event_refinement, bck))
  {
    if (isolate) {
      this->isolate();
    }
  }
  
  //!@} // Constructors

  //! Access members
  //!@{

  // beyond the ones inherited from Base()

  //! returns the degree of the gcd of f and its derivative, if known
  int degree_of_gcd() const {
    if (dynamic_cast< const Square_free_rep* >(this->ptr())) {
      return 0;
    }
    if (dynamic_cast< const Mk_rep* >(this->ptr())) {
      return  dynamic_cast< const Mk_rep* >(this->ptr())->degree_of_gcd();
    }
    // else
    throw Virtual_method_exception();
  }
  
  //! returns the square free part of f, if known
  Polynomial square_free_part() const {
    if (dynamic_cast< const Square_free_rep* >(this->ptr())) {
      return dynamic_cast< const Square_free_rep* >(this->ptr())->polynomial();
    }
    // else
    throw Virtual_method_exception();
  }

  // TODO 2012 CA/CPA does not need tree() function (only _test_mkbs_lifter seems to call it)
  //! returns the tree
  Bitstream_tree tree() const {
    return dynamic_cast<const Generic_bs_rep*>(this->ptr())->tree();
  }

  //!@} Access members

};



} // namespace internal

} // namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_D_CURVE_ANALYSIS_2_MKBS_LIFTER_H
// EOF
