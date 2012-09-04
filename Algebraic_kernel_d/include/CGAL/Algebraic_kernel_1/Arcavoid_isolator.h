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
// $URL:  $
// $Id:   $
//
//
// Author(s)     :  Alexander Kobel <akobel@mpi-inf.mpg.de>
//
// ============================================================================

/*! \file CGAL/Arcavoid/Arcavoid_isolator.h
  \brief GenericIsolator interface for Arcavoid
*/

#ifndef CGAL_ARCAVOID_ISOLATOR_H
#define CGAL_ARCAVOID_ISOLATOR_H

#include <CGAL/config.h>

#include <CGAL/Handle_with_policy.h>
#include <CGAL/Algebraic_kernel_d/Generic_isolator.h>

#include <CGAL/Polynomial_traits_d.h>

#include <CGAL/Algebraic_kernel_1/Arcavoid_root_isolator.h>

namespace CGAL {

namespace internal {

// Currently only forwarding to Arcavoid< BCK, Arcavoid_real_root_isolator_tag >
// TODO: directly derive the latter from Generic_isolator_rep_base

template< class BitstreamCoefficientKernel, typename HandlePolicy >
class Arcavoid_isolator_rep :
  public Generic_isolator_rep_base< typename BitstreamCoefficientKernel::Bound, HandlePolicy, CGAL::Tag_true >,
  public Arcavoid< BitstreamCoefficientKernel, Arcavoid_real_root_isolator_tag >
{

  //!\name Tags
  //!@{

public:
  //! "refine_interval" function is implemented
  typedef CGAL::Tag_true Refine_interval_provided_tag;

  //!@} // Tags

public:

  //!\name Public typedefs
  //!@{

  //! this instance's template parameter
  typedef BitstreamCoefficientKernel Bitstream_coefficient_kernel;

  //! the policy class
  typedef HandlePolicy Handle_policy;

  //! The polynomial type
  typedef typename Bitstream_coefficient_kernel::Polynomial Polynomial;

  //! type of Bound
  typedef typename Bitstream_coefficient_kernel::Bound Bound;

  //! the first base class
  typedef Generic_isolator_rep_base< Bound, Handle_policy, Refine_interval_provided_tag > Base;

  //! the class itself
  typedef Arcavoid_isolator_rep< Bitstream_coefficient_kernel, Handle_policy > Self;

  //! the second base class
  typedef Arcavoid< BitstreamCoefficientKernel, Arcavoid_real_root_isolator_tag > Arcavoid_base;

  typedef typename Arcavoid_base::Algebraic_complex_1_iterator Algebraic_complex_1_iterator;
  typedef typename Arcavoid_base::Algebraic_complex_1_const_iterator Algebraic_complex_1_const_iterator;

  //!@} // Public typedefs

public:

  //!\name Constructors
  //!@{

#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_DEFAULT_CONSTRUCTIBLE
protected:
  Arcavoid_isolator_rep(); // = disable
#else
public:
  //! Default constructor
  Arcavoid_isolator_rep() : Base(new Rep()) {}
#endif

  // needs no special assignable-implementation as no pointers
  // and this is a Rep of a handle class
#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_ASSIGNABLE
protected:
  Arcavoid_isolator_rep(const Self &rhs) : Base (rhs), Arcavoid_base (rhs) {} // = default
  Self& operator=(const Self&); // = disable
#endif


public:
  //! standard constructor
  template< class CoefficientInputIterator >
  Arcavoid_isolator_rep (CoefficientInputIterator first,
                         CoefficientInputIterator beyond,
                         const Bitstream_coefficient_kernel &bck) :
    Base(), Arcavoid_base (bck, first, beyond) {
    this->_m_is_isolated = true;
  }

  //! polynomial constructor
  Arcavoid_isolator_rep(const Polynomial& f,
                        const Bitstream_coefficient_kernel &bck) :
    Base(), Arcavoid_base(bck, f) {
    this->_m_is_isolated = true;
  }

  //!@} // Constructors

public:

  //!\name Access functions
  //!@{

  // Currently only forwarding to Arcavoid< BCK, Arcavoid_real_root_isolator_tag >
  // TODO: directly derive the latter from Generic_isolator_rep_base


  Polynomial polynomial () const { return Arcavoid_base::polynomial(); }
  void isolate() {}
  int number_of_real_roots () const { return static_cast< int > (Arcavoid_base::number_of_real_roots()); }
  Bound left_bound (size_t i) const { return Arcavoid_base::left_bound (i); }
  Bound right_bound (size_t i) const { return Arcavoid_base::right_bound (i); }
  bool is_exact_root (size_t i) const { return Arcavoid_base::is_exact_root (i); }
  bool is_certainly_simple_root (size_t i) const { return Arcavoid_base::is_certainly_simple_root (i); }
  bool is_certainly_multiple_root (size_t i) const { return Arcavoid_base::is_certainly_multiple_root (i); }
  bool get_upper_bound_for_multiplicity (int i) const { return Arcavoid_base::get_upper_bound_for_multiplicity (i); }
  int upper_bound_for_multiplicity (size_t i) const { return Arcavoid_base::get_upper_bound_for_multiplicity (i); }
  int multiplicity_of_root (size_t i) const { return Arcavoid_base::get_upper_bound_for_multiplicity (i); }
  void refine_interval (size_t i) { Arcavoid_base::refine_interval (i); }

  const Algebraic_complex_1_iterator real_roots_begin () { return Arcavoid_base::real_roots_begin(); }
  const Algebraic_complex_1_iterator real_roots_end () { return Arcavoid_base::real_roots_end(); }
  const Algebraic_complex_1_const_iterator real_roots_begin () const { return Arcavoid_base::real_roots_begin(); }
  const Algebraic_complex_1_const_iterator real_roots_end () const { return Arcavoid_base::real_roots_end(); }

  //!@} // Access functions

protected:

  // real root iterator begin
  Algebraic_complex_1_const_iterator begin() const { return Arcavoid_base::real_roots_begin(); }

  // real root iterator end
  Algebraic_complex_1_const_iterator end() const { return Arcavoid_base::real_roots_end(); }

  //! Needed for reference counting
  virtual CGAL::Reference_counted_hierarchy<>* clone() {
    return new Arcavoid_isolator_rep(*this);
  }
};

template < class BitstreamCoefficientKernel, typename HandlePolicy = CGAL::Handle_policy_no_union >
class Arcavoid_isolator :
  public internal::Generic_isolator< typename BitstreamCoefficientKernel::Polynomial,
                                     typename BitstreamCoefficientKernel::Bound,
                                     HandlePolicy,
                                     internal::Arcavoid_isolator_rep< BitstreamCoefficientKernel, HandlePolicy > > {
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
  typedef internal::Arcavoid_isolator_rep< Bitstream_coefficient_kernel, Handle_policy > Rep;

  //! the base type
  typedef internal::Generic_isolator< Polynomial, Bound, Handle_policy, Rep > Base;

  //! the class itself
  typedef Arcavoid_isolator< Bitstream_coefficient_kernel, Handle_policy > Self;

  //! type underlying Arcavoid real root isolator
  typedef typename Rep::Arcavoid_base Arcavoid_base;

  typedef typename Rep::Algebraic_complex_1_iterator Algebraic_complex_1_iterator;
  typedef typename Rep::Algebraic_complex_1_const_iterator Algebraic_complex_1_const_iterator;

  //!@} // Public types

  //!\name Constructors
  //!@{

#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_DEFAULT_CONSTRUCTIBLE
 protected:
  Arcavoid_isolator(); // = disable
#else
 public:
  //! Default constructor
 Arcavoid_isolator() : Base(new Rep()) {}
#endif

#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_ASSIGNABLE
 protected:
  Arcavoid_isolator(const Self&); // = disable
  Arcavoid_isolator& operator=(const Self&); // = disable
#endif

 public:
  /*!
   * \brief Constructor for a polynomial \c f
   */
  Arcavoid_isolator(const Polynomial& p,
                    // TODO 2012 default bck?
                    Bitstream_coefficient_kernel bck = Bitstream_coefficient_kernel(),
                    bool isolate = true)
    : Base(new Rep(p, bck))
  {
    if (isolate) {
      this->isolate();
    }
  }

  //! standard constructor
  template< class CoefficientInputIterator >
  Arcavoid_isolator(CoefficientInputIterator first, CoefficientInputIterator beyond,
                    Bitstream_coefficient_kernel bck = Bitstream_coefficient_kernel(),
                    bool isolate = true)
    : Base(new Rep (first, beyond, bck))
  {
    if (isolate) {
      this->isolate();
    }
  }

  //!@}

  //!\name Access members
  //!@{

  // the ones beyond those inherited from Base()

  //! returns the polynomial
  Polynomial polynomial() const {
    return dynamic_cast< const Rep* >(this->ptr())->polynomial();
  }

  //! Returns the kernel class
  Bitstream_coefficient_kernel bck() const {
    return this->ptr()->bck();
  }

  //!@} // Access members

}; // Arcavoid_isolator

} // namespace CGAL

} // namespace internal

#endif // CGAL_ARCAVOID_ISOLATOR_H
// EOF
