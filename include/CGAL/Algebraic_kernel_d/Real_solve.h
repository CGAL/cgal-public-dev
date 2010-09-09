// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     :  Eric Berberich <eric@mpi-inf.mpg.de>
//                  Michael Hemmer <hemmer@mpi-inf.mpg.de>
//
// ============================================================================

/*! \file Real_solve.h
  \brief Defines class CGAL::Real_solve
  
  Isolate real roots of polynomials with Fabrice Roullier's Rs.

  The polynomial has to be a univariat polynomial over any number
  type which is contained in the real numbers.

*/

#ifndef CGAL_ALGEBRAIC_KERNEL_D_REAL_SOLVE_H
#define CGAL_ALGEBRAIC_KERNEL_D_REAL_SOLVE_H

#include <CGAL/config.h>

#include <CGAL/Algebraic_kernel_rs_gmpz_d_1.h>

namespace CGAL {

namespace internal {

/*! \brief A model of concept RealRootIsolator. 
 */
template <class Polynomial_, class Rational_> 
class Real_solve {

public:
    //! First template parameter
    typedef Polynomial_ Polynomial;

    //! Second template parameter
    typedef Rational_ Rational;

    //! Bound type of the isolating intervals
    typedef Rational_ Bound;
private:
    
    //! Coefficient type of polynomial
    typedef typename Polynomial::NT Coefficient;

    //! internal algebraic kernel
    typedef CGAL::Algebraic_kernel_rs_gmpz_d_1 AK_1;

    //! type of algebraic real
    typedef typename AK_1::Algebraic_real_1 Algebraic_real_1;
    
public:
    /*! \brief Constructor from univariate square free polynomial.
    
    The RealRootIsolator provides isolating intervals for the real 
    roots of the polynomial.
    \pre the polynomial is square free      
    */
   Real_solve(const Polynomial& p = Polynomial(Coefficient(0))) :
      _m_polynomial(p)
      //, _m_interval_given(false)
    {
      typename AK_1::Solve_1()(_m_polynomial, true, std::back_inserter(_m_real_roots));
      if (_m_real_roots.size() > 1) {
        for (std::vector< Algebraic_real_1 >::size_type i = 0;
             i < _m_real_roots.size() - 1;
             ++i) {
          // Isolate them against each other
          typename AK_1::Compare_1()(_m_real_roots[i],_m_real_roots[i+1]);
        }
      }
    }
    
    // TODO interval constructor?

public: // functions
    
    /*! \brief returns the defining polynomial*/ 
    Polynomial polynomial() const { 
      return _m_polynomial;
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
    bool is_exact_root(int i) const {
      return _m_real_roots[i].is_point();
    }
      
public:   
  
    /*! \brief returns  \f${l_i}\f$ the left bound of the isolating interval 
      for root  \f$root_{i}\f$.
      
      In case is_exact_root(i) is true,  \f$l_i = root_{i}\f$,\n
      otherwise:  \f$l_i < root_{i}\f$. 
         
      If  \f$i-1>=0\f$, then  \f$l_i > root_{i-1}\f$. \n
      If  \f$i-1>=0\f$, then  \f$l_i >= r_{i-1}\f$, 
      the right bound of  \f$root_{i-1}\f$\n

      \pre 0 <= i < number_of_real_roots()
    */
    Bound left_bound(int i) const { 
      CGAL_assertion(i >= 0);
      CGAL_assertion(i < this->number_of_real_roots());
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
    Bound right_bound(int i) const { 
      CGAL_assertion(i >= 0);
      CGAL_assertion(i < this->number_of_real_roots());
      return Bound(_m_real_roots[i].right());
    }  
    
private:
    
    //! the squarefree polynomial
    Polynomial _m_polynomial;
    
    //! the solutions
    std::vector<  Algebraic_real_1 > _m_real_roots;
   
    //! restricted interval?
    // TODO bool _m_interval_given;
    
};

} // namespace internal

} //namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_D_REAL_SOLVE_H
