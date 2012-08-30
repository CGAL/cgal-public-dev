// Copyright (c) 2010 Max-Planck-Institut fuer Informatik (Germany).
// All rights reserved.
//
// LGPL?
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
// Author(s): Eric Berberich    <eric@mpi-inf.mpg.de>
//            Pavel Emeliyanenko <asm@mpi-inf.mpg.de>
// ============================================================================


#ifndef CGAL_ROUNDING_AK_D_1_H
#define CGAL_ROUNDING_AK_D_1_H

/*! \file Rounding_ak_d_1.h
 * Wrapper for univariate algebraic kernel that rounds approximations
 */

#include <CGAL/convert_to_bfi.h>

namespace CGAL {

namespace internal {

template < class AlgebraicKernel_d_1 >
class Rounding_ak_d_1 : public AlgebraicKernel_d_1 {

  public:

    //! instance parameter
    typedef AlgebraicKernel_d_1 Algebraic_kernel_d_1;

    typedef Rounding_ak_d_1< Algebraic_kernel_d_1 > Rounding_algebraic_kernel_d_1;

    typedef typename Algebraic_kernel_d_1::Bound Bound;

    typedef typename Algebraic_kernel_d_1::Algebraic_real_1 Algebraic_real_1;

  private:

    static std::pair< Bound, Bound > _round(const std::pair< Bound, Bound >& intv, int prec) {
      
      Bound l, u; 
      
      typedef typename CGAL::Get_arithmetic_kernel< Bound >::Arithmetic_kernel::Bigfloat_interval BFI;
      long old_prec = CGAL::set_precision(BFI(), prec);

      prec = std::max(prec, 2);
      CGAL::set_precision(BFI(), prec);
      l = CGAL::lower(CGAL::convert_to_bfi(intv.first));
      u = CGAL::upper(CGAL::convert_to_bfi(intv.second));
      CGAL::set_precision(BFI(), old_prec);
      
      return std::make_pair(l,u);
    }

  public:
    
    struct Approximate_absolute_1:
      public std::binary_function<Algebraic_real_1,int,std::pair<Bound,Bound> >{
   
      std::pair<Bound,Bound> 
	operator()(const Algebraic_real_1& x, int prec) const {
	
	CGAL_precondition(prec >= 0);
	
	typename Algebraic_kernel_d_1::Approximate_absolute_1 approx;
	
	// max(2,prec) for GMP types
	prec = std::max(2,prec);
	return Rounding_algebraic_kernel_d_1::_round(approx(x, prec), prec);
      }
    };

    struct Approximate_relative_1:
      public std::binary_function<Algebraic_real_1,int,std::pair<Bound,Bound> > {
     
      std::pair<Bound,Bound> 
	operator()(const Algebraic_real_1& x, int prec) const {
	
	CGAL_precondition(prec >= 0);

	typename Algebraic_kernel_d_1::Approximate_relative_1 approx;
	
	// max(2,prec) for GMP types
	prec = std::max(2,prec);
	return Rounding_algebraic_kernel_d_1::_round(approx(x, prec), prec);
      }
    };


#define CGAL_ALGEBRAIC_KERNEL_1_PRED(Y,Z) Y Z() const { return Y(); }

  CGAL_ALGEBRAIC_KERNEL_1_PRED(Approximate_absolute_1,
      approximate_absolute_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Approximate_relative_1,
      approximate_relative_1_object);

#undef CGAL_ALGEBRAIC_KERNEL_1_PRED

};

} // namespace internal

} // namespace CGAL

#endif // CGAL_ROUNDING_AK_D_1
