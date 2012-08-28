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
// Author(s)     :  Michael Kerber <mkerber@mpi-inf.mpg.de>
//                  Eric Berberich <eric.berberich@cgal.org>
//
// ============================================================================

/*! \file CGAL/Algebraic_kernel_3/Algebraic_surface_3_lifter.h
  \brief Defines various isolators for polynomials.
*/

#ifndef CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_KERNEL_3_LIFTER_H
#define CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_KERNEL_3_LIFTER_H 1

#include <CGAL/config.h>

#include <CGAL/tags.h>

#include <CGAL/Handle_with_policy.h>

#include <CGAL/Algebraic_kernel_d/Generic_isolator.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_isolator_reps.h> // for Generic_bitstream_descartes_isolator_rep
#include <CGAL/Algebraic_kernel_d/Curve_analysis_2_mkbs_lifter.h> // contains Mk_bitstream_descartes_isolator_rep

namespace CGAL {

namespace internal {

/*
 * \brief Adaptor for roots of a vert line 
 * (needed as dummy in surface analysis)
 *
 */
template< typename BitstreamCoefficientKernel, typename VerticalLine, typename HandlePolicy > // no default on policy, should be decided in higher level
class Vertical_line_adapter_rep : 
  public Generic_bitstream_descartes_isolator_rep< BitstreamCoefficientKernel, HandlePolicy > { // TODO 2012 derived from new bck-only-rep
       
  // tags are inherited
    
public:
        
  //!@\name Public typedefs
  //!@{

  //! first template parameter: Bitstream coefficient kernel
  typedef BitstreamCoefficientKernel Bitstream_coefficient_kernel;

  //! second template parameter: Vertical_line
  typedef VerticalLine Vertical_line;

  //! third template parameter: Handle_policy
  typedef HandlePolicy Handle_policy;

  //! The generic representation
  typedef Generic_bitstream_descartes_isolator_rep< Bitstream_coefficient_kernel, Handle_policy > Base;

  //! the class itself
  typedef Vertical_line_adapter_rep< Bitstream_coefficient_kernel, Vertical_line, Handle_policy > Self;
      
  //! The polynomial type
  typedef typename Bitstream_coefficient_kernel::Polynomial Polynomial;

  //! How the boundaries of the isolating intervals are represented
  typedef typename Bitstream_coefficient_kernel::Bound Bound;

  //! type of Curve_analysis_2
  typedef typename Vertical_line::Curve_analysis_2 Curve_analysis_2;

  //! type of Curve_kernel_2;
  typedef typename Curve_analysis_2::Algebraic_kernel_with_analysis_2 Curve_kernel_2;

  //!@} // Public typedefs
       
  //!@\name Constructors
  //!@{

  //! Default constructor (does nothing)
#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_DEFAULT_CONSTRUCTIBLE
  protected:
#else
  public:
#endif
  Vertical_line_adapter_rep() {} // = default

  // needs no special assignable-implementation as no pointers 
  // and this is a Rep of a handle class (depending on Handle_policy!)
#if CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_ASSIGNABLE
  protected:
  // Vertical_line_adapter_rep(const Self&); // = default
  Vertical_line_adapter_rep& operator=(const Self&); // = disable
#endif

  public:

  /*! 
   * \brief Constructor
   */
  template<typename InputIterator>
  Vertical_line_adapter_rep(InputIterator begin, InputIterator end,
                            Bitstream_coefficient_kernel bck)
    : Base(Polynomial(0), bck) 
  {
    for (InputIterator it = begin; it != end; it++) {
      _m_root_vec.push_back(std::make_pair(*it, 4));
    }
    
    // TODO 2012 move to isolate()

    this->_m_number_of_intervals = static_cast<int>(_m_root_vec.size());
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

    this->_m_is_isolated = true;

  }
  
  //!@} // Constructors

  //!\name Destructor
  //!@{
        
  //! Destructor (does nothing)
  virtual ~Vertical_line_adapter_rep() {
  }

  //!@}

  
public:

  //!\name Acces members
  //!@{

  virtual void isolate() {
    // TODO implement

  }
  
  // number of real_roots is derived from base class

  //! The lower bound of the \c i th root
  virtual Bound left_bound(size_t i) const  {
    typename Curve_kernel_2::Approximate_absolute_y_2
      approx_y = this->bck().point().xy().kernel()
      ->approximate_absolute_y_2_object();
    return approx_y(_m_root_vec[i].first.first.algebraic_real_2
                    (_m_root_vec[i].first.second),
                    _m_root_vec[i].second).first;
  }
    
  //! The upper bound of the \c i th root
  virtual Bound right_bound(size_t i) const {
    typename Curve_kernel_2::Approximate_absolute_y_2
      approx_y = this->bck().point().xy().kernel()
      ->approximate_absolute_y_2_object();
    return approx_y(_m_root_vec[i].first.first.algebraic_real_2
                    (_m_root_vec[i].first.second),
                    _m_root_vec[i].second).second;
  }

  // is_exact_root is derived from base class

  /*! \brief Returns whether the \c i th root is definitely a simple root
   * of the isolated polynomial
   *
   */
  virtual bool is_certainly_simple_root(size_t /* i */) const {
    return false;
  }
        
  /*! \brief Returns whether the \c i th root is definitely 
   * a multiple root
   * of the isolated polynomial
   *
   */
  virtual bool is_certainly_multiple_root(size_t /* i */) const {
    return false;
  }

  // upper_bound_for_multiplicity is derived from base class

  // TODO 2012 dispatch for multiplicity of root
  int multiplicity_of_root( size_t /* i */) const {
    bool _WARNING_DO_NOT_USE_multiplicity_of_root_which_is_not_yet_implemented;
    CGAL_error(); // TODO 2012 how to forbid access to this
    return -1;
  }

  //! Computes a better approximation of the \c i th root of the
  virtual void refine_interval(size_t i) const {
    _m_root_vec[i].second *= 2;
  }

  private:

  void process_nodes() { // TODO 2012 get rid of empty implementation by introducing non-bitstream-bck-rep
  }

  bool termination_condition() { // TODO 2012 get rid of empty implementation by introducing non-bitstream-bck-rep
    return false;
  }

  private:

  //! Needed for the referencing counting mechanism
  virtual CGAL::Reference_counted_hierarchy<>* clone() {
    return new Vertical_line_adapter_rep(*this);
  }

  //!@} // Access members

protected:

  //! Roots stored as pair of a Vertical_line and an integer denoting the
  //! index. Also, current precision of each root is stored
  mutable std::vector<std::pair<std::pair<Vertical_line, int>, int> > _m_root_vec;

}; // Vertical_line_adapter_rep


//! will be used for lifting over faces, edges, vertices in Surface_analysis
template < class BitstreamCoefficientKernel, typename HandlePolicy = CGAL::Handle_policy_no_union >
class Algebraic_surface_3_lifter :
  public internal::Generic_isolator< typename BitstreamCoefficientKernel::Polynomial,
                                     typename BitstreamCoefficientKernel::Bound, 
                                     HandlePolicy, 
 Square_free_bitstream_descartes_isolator_rep< BitstreamCoefficientKernel, HandlePolicy >
    /*Generic_bitstream_descartes_isolator_rep< BitstreamCoefficientKernel, HandlePolicy > */
    > {


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

  //! second template parameter: Handle_policy
  typedef HandlePolicy Handle_policy;

  //! The polynomial type
  typedef typename Bitstream_coefficient_kernel::Polynomial Polynomial;

  //! How the boundaries of the isolating intervals are represented
  typedef typename Bitstream_coefficient_kernel::Bound Bound;

  //! generic rep that serves as base class for all reps
  typedef Generic_bitstream_descartes_isolator_rep< Bitstream_coefficient_kernel, Handle_policy > BRep;

  //! the representation class for square-free case1
  // TODO 2012 pick right square-free rep or make square-free even a parameter
  typedef Square_free_bitstream_descartes_isolator_rep< Bitstream_coefficient_kernel, Handle_policy > Square_free_rep;

  //! the representation class for mk case
  typedef Mk_bitstream_descartes_isolator_rep <Bitstream_coefficient_kernel, Handle_policy > Mk_rep;

  //! the following needs a parameter given only for the corresponding constructor, thus it's commented out here:
  // typedef Vertical_line_adapter_rep<Bitstream_coefficient_kernel, typename InputIterator::value_type::first_type, Handle_policy > Vertical_line_rep;

  //! the base type
  typedef internal::Generic_isolator< Polynomial, Bound, Handle_policy, Square_free_rep /* in fact, provide all of the used reps! */ > Base;

  //! the class itself
  typedef Algebraic_surface_3_lifter< Bitstream_coefficient_kernel, Handle_policy > Self;

  //!@} // Public types

 public:

  //!\name Constructors 
  //!@{

  //! default constructor // TODO 2012 remove?
  Algebraic_surface_3_lifter() : Base(typename Base::Use_with_initialize_with()) {
  }

  // Comment: Copy and assign constructor are default!

  // TODO 2012 replace all constructors by one ctor for event and one ctor for internal
  // and let this class decide which is the best isolator (e.g. with help of modular arithmetic)

  //! bitstream constructor (gets a bivariate bck with x = alpha stored)
  Algebraic_surface_3_lifter(const Polynomial& p,
                             const Bitstream_coefficient_kernel& bck,
                             bool isolate = true) :
    // use square-free bitstream descartes
    Base(new Square_free_rep(p, bck))
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
  Algebraic_surface_3_lifter(const Polynomial& p,
                             int m, int k,
                             Bitstream_coefficient_kernel bck = Bitstream_coefficient_kernel(),
                             bool isolate = true) : 
    Base(new Mk_rep(p, m, k, bck))
  {
    if (isolate) {
      this->isolate();
    }
  }

  template<typename InputIterator>
  Algebraic_surface_3_lifter(InputIterator begin,
                             InputIterator end,
                             Bitstream_coefficient_kernel bck,
                             bool isolate = true) : 
    Base(new Vertical_line_adapter_rep<Bitstream_coefficient_kernel, typename InputIterator::value_type::first_type, Handle_policy >(begin, end, bck))
  {
    if (isolate) {
      this->isolate();
    }
  }

  //!@} // Constructors

  //! Access members
  //!@{
  //! access to kernel

  //! access to kernel
  Bitstream_coefficient_kernel bck() const {
    return dynamic_cast< const BRep* >(this->ptr())->bck();
  }

  //! The length of the <tt>i</tt>th isolating interval
  // added for the purpose
  Bound length(size_t i) const {
    CGAL_precondition(i >= 0);
    CGAL_precondition(static_cast<int>(i) < this->number_of_real_roots());
    CGAL_precondition(this->is_isolated());
    return (this->right_bound(i) - this->left_bound(i));
  }

  enum Rep_type {
    SF = 0,
    MK = 1,
    VL = 2
  };

  Rep_type type() const {
    if (dynamic_cast< const Square_free_rep* >(this->ptr())) {
      return SF;
    }
    if (dynamic_cast< const Mk_rep* >(this->ptr())) {
      return MK;
    }
    // else
    return VL;
  }

  //!@} Access members

}; // Algebraic_kernel_3_lifter

} // namespace internal

} // namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_KERNEL_3_LIFTER_H
// EOF
