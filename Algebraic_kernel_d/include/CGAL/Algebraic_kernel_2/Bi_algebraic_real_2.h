// Copyright (c) 2009, 2010, 2011, 2012 Max-Planck-Institut fuer Informatik (Germany).
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


#ifndef CGAL_ALGEBRAIC_KERNEL_2_BI_ALGEBRAIC_REAL_2_H
#define CGAL_ALGEBRAIC_KERNEL_2_BI_ALGEBRAIC_REAL_2_H

/*! \file
 * The header file for the Bisolve_2 class.
 */

#include <CGAL/config.h>
#include <CGAL/Handle_with_policy.h>

namespace CGAL {

template < class AlgebraicKernel_1 >
class Bi_algebraic_real_2;

namespace internal {


//! Template class to represent a pair of algerbaic reals
template < class AlgebraicKernel_1 >
class Bi_algebraic_real_2_rep {

public:

  //! the first template parameter
  typedef AlgebraicKernel_1 Algebraic_kernel_1;

  //! the real type
  typedef typename Algebraic_kernel_1::Algebraic_real_1 Algebraic_real_1;

  //! Multiplicity type
  typedef typename Algebraic_kernel_1::Multiplicity_type Multiplicity_type;

  //! the bound type
  typedef typename Algebraic_kernel_1::Bound Bound;

  //!\name Constructors
  //!@{

  //! default constructor
  Bi_algebraic_real_2_rep() {
  }
  
  //! constructor from two coordinates
  Bi_algebraic_real_2_rep(
      const Algebraic_real_1& x, const Multiplicity_type& mx,
      const Algebraic_real_1& y, const Multiplicity_type& my) :
    _m_x(x),
    _m_mx(mx),
    _m_y(y),
    _m_my(my) {
  }

  //!@}

private:

  //!\name Private members
  //!@{

  //! x-coordinate
  Algebraic_real_1 _m_x;

  //! multiplicity of x-coordinate
  Multiplicity_type _m_mx;

  //! y-coordinate
  Algebraic_real_1 _m_y;

  //! multiplicity of y-coordinate
  Multiplicity_type _m_my;

  friend class Bi_algebraic_real_2< AlgebraicKernel_1 >;

  //!@}

}; // Bi_algebraic_real_2_rep

} // namespace internal

//! Template class to represent a pair of algerbaic reals
template < class AlgebraicKernel_1 >
class Bi_algebraic_real_2 : 
    public CGAL::Handle_with_policy< internal::Bi_algebraic_real_2_rep< 
      AlgebraicKernel_1 > > {

public:

  //! the first template parameter
  typedef AlgebraicKernel_1 Algebraic_kernel_1;

  //! rep type
  typedef internal::Bi_algebraic_real_2_rep< Algebraic_kernel_1 > Rep;
  
  //! base type
  typedef CGAL::Handle_with_policy< Rep > Base;
  
  //! self
  typedef Bi_algebraic_real_2< Algebraic_kernel_1 > Self;

  //! the real type
  typedef typename Rep::Algebraic_real_1 Algebraic_real_1;

  //! Multiplicity type
  typedef typename Rep::Multiplicity_type Multiplicity_type;

  //! the bound type
  typedef typename Rep::Bound Bound;

  //! type of a box
  typedef CGAL::cpp0x::array<Bound,4> Box_2;

  //!\name Constructors
  //!@{

  //! default constructor
  Bi_algebraic_real_2() {
  }
  
  //! constructor from two coordinates
  Bi_algebraic_real_2(const Algebraic_real_1& x, const Multiplicity_type& mx,
                      const Algebraic_real_1& y, const Multiplicity_type& my) :
    Base(Rep(x,mx,y,my)) {
  }

  //!@}
  //!\name Getters
  //!@{

  // TODO replace x(), mult_x() by x_with_mult return pair< AR_1, Mult_type >?
  // TODO replace x(), mult_x() by x_with_mult return pair< AR_1, Mult_type >?

  //! the x-coordinate
  Algebraic_real_1 x() const {
    return this->ptr()->_m_x;
  }

  //! the y-coordinate
  Algebraic_real_1 y() const {
    return this->ptr()->_m_y;
  }

  //! the x-coordinate's multiplicity
  Multiplicity_type mult_x() const {
    return this->ptr()->_m_mx;
  }

  //! the y-coordinate's multiplicity
  Multiplicity_type mult_y() const {
    return this->ptr()->_m_my;
  }

  //!@}
public:
  //!\name Comparisons
  //!@{

  //!\brief Compares x-coordinates of two instances
  CGAL::Comparison_result compare_x(const Bi_algebraic_real_2& bar) const {
    return CGAL::compare(this->x(), bar.x());
  }

  //!\brief Compares y-coordinates of two instances
  CGAL::Comparison_result compare_y(const Bi_algebraic_real_2& bar) const {
    return CGAL::compare(this->y(), bar.y());
  }
  
  //!\brief Compares coordinates lexicographically (first x, then y)
  CGAL::Comparison_result compare_xy(const Bi_algebraic_real_2& bar) const {
    CGAL::Comparison_result cmp_x = this->compare_x(bar);
    if (cmp_x != CGAL::EQUAL) {
      return cmp_x;
    }
    // else
    return this->compare_y(bar);
  }
  
  //!\brief Compares coordinates lexicographically (first y, then x)
  CGAL::Comparison_result compare_yx(const Bi_algebraic_real_2& bar) const {
    CGAL::Comparison_result cmp_y = this->compare_y(bar);
    if (cmp_y != CGAL::EQUAL) {
      return cmp_y;
    }
    // else
    return this->compare_x(bar);
  }

  //!\brief Are two objects equal
  bool equal(const Bi_algebraic_real_2& bar) {
    return this->compare_xy(bar) == CGAL::EQUAL;
  }
  
  //! lexicographic <
  bool operator<(const Bi_algebraic_real_2& bar) const {
    return compare_xy(bar) == CGAL::SMALLER;
  }

  //!@}

  //!\name Containment
  //!@{

  //! tests whether bi-algebraic real is contained in a box
  bool is_contained_in(const Box_2& box) const {
    if (CGAL::compare(box[0], this->x()) == CGAL::LARGER) {
      return false;
    }
    if (CGAL::compare(box[1], this->x()) == CGAL::SMALLER) {
      return false;
    }
    if (CGAL::compare(box[2], this->y()) == CGAL::LARGER) {
      return false;
    }
    if (CGAL::compare(box[3], this->y()) == CGAL::SMALLER) {
      return false;
    }

    return true;
  }

  //!@}
  
  //!\name IO
  //!@{

  void write(std::ostream& os) const {

    typename Algebraic_kernel_1::Approximate_relative_1 approx_rel;
    std::pair< Bound, Bound > interval_x = approx_rel(x(),0);
    std::pair< Bound, Bound > interval_y = approx_rel(y(),0);

    os << "x: " << CGAL::to_double(x())
    //<< CGAL::lower(CGAL::convert_to_bfi(x()))
//        << " [" << interval_x.first << "; " << interval_x.second << "]"
       << " (" << mult_x() << ")";
    os << ", y: " << CGAL::to_double(y())
    //<< CGAL::lower(CGAL::convert_to_bfi(y()))
       //<< " [" << interval_y.first << "; " << interval_y.second << "]"
       << " (" << mult_y() << ")";
  }

  //!@}

}; // Bi_algebraic_real_2

//! output
template < class AlgebraicKernel_1 >
std::ostream& operator<<(std::ostream& os,
        const Bi_algebraic_real_2< AlgebraicKernel_1 >& v) {

    v.write(os);
    return os;
}

} // namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_2_BI_ALGEBRAIC_REAL_2
// EOF
