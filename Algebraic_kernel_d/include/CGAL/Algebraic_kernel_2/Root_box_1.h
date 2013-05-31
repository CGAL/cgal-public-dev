// Copyright (c) 2010, 2012 Max-Planck-Institut fuer Informatik (Germany).
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
// $URL: svn+ssh://asm@scm.gforge.inria.fr/svn/cgal/branches/experimental-packages/Algebraic_kernel_2/include/CGAL/Bi_slice_certify.h $
// $Id: Bi_slice_certify.h 57108 2010-06-25 13:54:53Z eric $
//
//
// Author(s): Eric Berberich <eric.berberich@cgal.org>
//            Pavel Emeliyanenko <asm@mpi-inf.mpg.de>


#ifndef CGAL_ALGEBRAIC_KERNEL_2_ROOT_BOX_1_H
#define CGAL_ALGEBRAIC_KERNEL_2_ROOT_BOX_1_H

/*! \file Root_box_1.h
 * Isolating box of a real root of a univariate polynomial with additional info
 */

#include <CGAL/config.h>

#include <iostream>

#include <boost/optional.hpp>

namespace CGAL {

namespace internal {

template < class Bound_ >
class Root_box_1 {

 public:
  //! instance parameter
  typedef Bound_ Bound;
  
  Root_box_1() :
    _m_is_exact(false), 
    _m_used_t_test(false) {
  }

  void set_mid(const Bound& m) {
    _m_left == boost::none;
    _m_right == boost::none;
    _m_mid = m;
  }
  
  Bound mid() const {
    return _m_mid;
  }

  void set_rad(const Bound& r) {
    _m_rad = r;
  }
  
  Bound rad() const {
    _m_left == boost::none;
    _m_right == boost::none;
    return _m_rad;
  }

  void set_exact() {
    _m_is_exact = true;
  }
  
  bool is_exact() const {
    return _m_is_exact;
  }
  
  void set_used_t_test() {
    _m_used_t_test = true;
  }
  
  bool used_t_test() const {
    return _m_used_t_test;
  }

  Bound left() const {
    if (!_m_left) {
      _m_left = _m_mid - _m_rad;
    }
    return *_m_left;
  }

  Bound right() const {
    if (!_m_right) {
      _m_right = _m_mid + _m_rad;
    }
    return *_m_right;
  }
  
  void write(std::ostream& os) {
    os << "    * m=" << CGAL::to_double(this->mid()) << " (" << this->mid()
            << ")" << std::endl;
    os << "    * r=" << CGAL::to_double(this->rad()) << " (" << this->rad()
            << ")" << std::endl;
    os << "    * exact=" << this->is_exact() << std::endl;
    os << "    * 3/2-test=" << this->used_t_test() << std::endl;
  }

 private:
  
  Bound _m_mid;
  Bound _m_rad;
  bool _m_is_exact;
  bool _m_used_t_test;
  
  mutable boost::optional< Bound > _m_left;
  mutable boost::optional< Bound > _m_right;
  
};

template < class Bound_ >
std::ostream& operator<<(std::ostream& os, Root_box_1< Bound_ > rb) {
  rb.write(os);
  return os;
}
 
} // namespace internal

} // namespace CGAL

#endif // CGAL_ROOT_BOX_1_H
// EOF

