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


#ifndef CGAL_ALGEBRAIC_KERNEL_2_ROUNDING_AK_D_1_H
#define CGAL_ALGEBRAIC_KERNEL_2_ROUNDING_AK_D_1_H

/*! \file Rounding_ak_d_1.h
 * Wrapper for univariate algebraic kernel that rounds approximations
 */

#include <CGAL/config.h>
#include <CGAL/convert_to_bfi.h>

#include <CGAL/Fraction_traits.h>
#include <CGAL/Algebraic_kernel_d/Real_embeddable_extension.h>

namespace CGAL {

namespace internal {

// TODO move specializations to other file

#if CGAL_USE_CORE

// Specialization for CORE::BigRat
template<>
class Real_embeddable_extension< CORE::BigRat > {
public:
  typedef CORE::BigRat Type;
#if 0

  struct Floor_log2_abs
    : public std::unary_function< CORE::BigRat, long > {
    long operator()( const CORE::BigRat& x ) const {
      return CORE::floorLg(x);
    }            
  };
     
#endif
   
  struct Ceil_log2_abs
    : public std::unary_function< CORE::BigRat, long > {
    long operator()( const CORE::BigRat& x ) const {
      return CORE::ceilLg(x);
    }
  };

#if 0
  struct Floor
    : public std::unary_function< CORE::BigRat, CORE::BigRat > {
    CORE::BigRat operator() (const CORE::BigRat& x) const { 
      return x;
    }
  };
  struct Ceil
    : public std::unary_function< CORE::BigRat, CORE::BigRat > {
    CORE::BigRat operator() (const CORE::BigRat& x) const { 
      return x;
    }
  };
#endif

};

#endif // CGAL_USE_CORE

#if CGAL_USE_GMP

// Specialization for Gmpq
template<>
class Real_embeddable_extension< Gmpq > {
public:
  typedef Gmpq Type;

#if 0
  struct Floor_log2_abs
    : public std::unary_function< Gmpq, long > {
    long operator()( const Gmpq& x ) const {
      CGAL_precondition(!CGAL::is_zero(x));
      return mpz_sizeinbase(x.mpz(),2)-1;
    }            
  };
#endif
        
  struct Ceil_log2_abs
    : public std::unary_function< Gmpq, long > {
    long operator()( const Gmpq& x ) const {

      typedef Fraction_traits< Gmpq > FT;
       FT::Numerator_type num;
       FT::Denominator_type denom;
       FT::Decompose decomp;
      
      decomp(x, num, denom);
       Real_embeddable_extension<  FT::Numerator_type >::Ceil_log2_abs log2_abs_num;
       Real_embeddable_extension<  FT::Denominator_type >::Floor_log2_abs log2_abs_den;
      
      return log2_abs_num(CGAL::abs(num)) - log2_abs_den(CGAL::abs(denom));
    }
  };

#if 0
  struct Floor
    : public std::unary_function< Gmpq, Gmpq > {
    Gmpq operator() (const Gmpq& x) const { 
      return x;
    }
  };
  struct Ceil
    : public std::unary_function< Gmpq, Gmpq > {
    Gmpq operator() (const Gmpq& x) const { 
      return x;
    }
  };
#endif
};


#endif // CGAL_USE_GMP



template < class AlgebraicKernel_d_1 >
class Rounding_ak_d_1 : public AlgebraicKernel_d_1 {

  public:

    //! instance parameter
    typedef AlgebraicKernel_d_1 Algebraic_kernel_d_1;

    typedef Rounding_ak_d_1< Algebraic_kernel_d_1 > Rounding_algebraic_kernel_d_1;

    typedef typename Algebraic_kernel_d_1::Bound Bound;
    typedef typename Algebraic_kernel_d_1::Coefficient Coefficient;

    typedef typename Algebraic_kernel_d_1::Algebraic_real_1 Algebraic_real_1;

  private:

    static std::pair< Bound, Bound > _round(const std::pair< Bound, Bound >& intv, int prec) {
      
      Bound l, u; 
      
      typedef typename CGAL::Get_arithmetic_kernel< Bound >::Arithmetic_kernel::Bigfloat_interval BFI;
      prec = std::max(prec, 2);
      long old_prec = CGAL::set_precision(BFI(), prec);
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
        
        typename Algebraic_kernel_d_1::Approximate_absolute_1 approx;
        
        std::pair<Bound, Bound> r = approx(x, prec);
        if (prec > 0) {
        
          long a_prec = std::max(2,prec); // for GMP types
          typename Real_embeddable_extension< Bound >::Ceil_log2_abs  log2_abs;
          a_prec = std::max(a_prec, 
                            a_prec + log2_abs(std::max(CGAL::abs(r.first), CGAL::abs(r.second))));
          return Rounding_algebraic_kernel_d_1::_round(approx(x, prec), a_prec);

        }
        
        // else
        return r;
      }
    };


    struct Approximate_relative_1:
      public std::binary_function<Algebraic_real_1,int,std::pair<Bound,Bound> > {

      std::pair<Bound,Bound>
	operator()(const Algebraic_real_1& x, int prec) const {
        
	typename Algebraic_kernel_d_1::Approximate_relative_1 approx;
        
        if (prec > 0) {
        
          prec = std::max(2,prec); // for GMP types
          return Rounding_algebraic_kernel_d_1::_round(approx(x, prec), prec);

        }
        
        // else
        return approx(x, prec);
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

#endif // CGAL_ALGEBRAIC_KERNEL_2_ROUNDING_AK_D_1
// EOF
