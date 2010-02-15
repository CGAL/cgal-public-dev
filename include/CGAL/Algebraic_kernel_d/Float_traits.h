// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : 
//
// ============================================================================

// TODO: Some comments are original EXACUS comments and aren't adapted. So
//         they may be wrong now.

// TODO: should exponent type be long or Integer ? 

#ifndef CGAL_ALGEBRAIC_KERNEL_D_FLOAT_TRAITS_H
#define CGAL_ALGEBRAIC_KERNEL_D_FLOAT_TRAITS_H

#include <CGAL/basic.h>

#if CGAL_USE_LEDA
#include <CGAL/leda_bigfloat.h>
#endif 

#if CGAL_USE_CORE
#include <CGAL/CORE_BigFloat.h>
#endif 

#if CGAL_USE_GMP
#include <CGAL/Gmpfr.h>
#endif

#include <CGAL/ipower.h>


CGAL_BEGIN_NAMESPACE

namespace internal {
    
// Don't define default, results in more convinient compiler messages 
template< class Type > class Float_traits;
// {
// public:        
//   typedef Null_functor    Get_mantissa;
//   typedef Null_functor    Get_exponent;  
//   typedef Null_functor    Mul_by_pow_of_2;
// };
    
#ifdef CGAL_USE_LEDA

// Specialization for leda_bigfloat
template<>
class Float_traits< leda_bigfloat > {
public:
  struct Get_mantissa
    : public std::unary_function< leda_bigfloat, leda_integer > {
    leda_integer operator()( const leda_bigfloat& x ) const {
      //std::cout << x.get_significant() << std::endl;
      return x.get_significant();                
    }
  };
        
  struct Get_exponent
    : public std::unary_function< leda_bigfloat, long > {
    long operator()( const leda_bigfloat& x ) const {
      return x.get_exponent().to_long();                
    }
  };

  struct Mul_by_pow_of_2
    : public std::binary_function< leda_bigfloat, long, leda_bigfloat> {
    leda_bigfloat operator()( const leda_bigfloat& a, long e ) const {
      return leda_bigfloat(a.get_significant(), a.get_exponent()+e);
    }
  };
};

#endif    
    
#ifdef CGAL_USE_CORE

// Specialization for CORE::BigFloat
template<>
class Float_traits< CORE::BigFloat > {
public:
      
  struct Get_mantissa
    : public std::unary_function< CORE::BigFloat, CORE::BigInt > {
    CORE::BigInt operator()( const CORE::BigFloat& x ) const { 
      return x.m();
    }
  };
        
  struct Get_exponent
    : public std::unary_function< CORE::BigFloat, long > {
    long operator()( const CORE::BigFloat& x ) const {
      return 14*x.exp(); // The basis is 8092                 
    }
  };

  struct Mul_by_pow_of_2
    : public std::binary_function
    < CORE::BigFloat, long , CORE::BigFloat> {
    CORE::BigFloat operator()( const CORE::BigFloat& a, long e ) const {
      return a*CORE::BigFloat::exp2(e);
    }
  };

};
#endif    


#if CGAL_USE_GMP
template<> class Float_traits< Gmpfr > {
  
  struct Get_mantissa_exponent
    : public std::unary_function< Gmpfr, std::pair<Gmpz,long> > {
    
    std::pair<Gmpz,long> operator()( const Gmpfr& x ) const {
      
      if(CGAL::is_zero(x)) 
        return std::make_pair(Gmpz(0),long(0));
      
      Gmpz z;
      long e=mpfr_get_z_exp(z.mpz(),x.fr());
      
      long zeros = mpz_scan1(z.mpz(),0);
      z >>= zeros;
      e +=  zeros;

      CGAL_postcondition(z % 2 != 0);
      CGAL_postcondition_code(if (e >= 0))
        CGAL_postcondition( x == (Gmpfr(z)) * CGAL::ipower(Gmpfr(2), e));
      CGAL_postcondition_code(else)
        CGAL_postcondition( x == (Gmpfr(z)) / CGAL::ipower(Gmpfr(2),-e));
      
      return std::make_pair(z,e);
    }
  };
public:  
  struct Get_mantissa
    : public std::unary_function< Gmpfr, Gmpz > {
    Gmpz operator()( const Gmpfr& x ) const {
      return Get_mantissa_exponent()(x).first;      
    }
  };
  
  struct Get_exponent
    : public std::unary_function< Gmpfr, long > {
    long operator()( const Gmpfr& x ) const { 
      return Get_mantissa_exponent()(x).second;      
    }
  };
    
struct Mul_by_pow_of_2
  : public std::binary_function< Gmpfr, Gmpz, Gmpfr> {
  Gmpfr operator()( const Gmpfr& a, long e ) const {

    //std::cout << "Mul_by_pow_of_2" <<std::endl;
 
    // std::pair<Gmpz,long> p = a.to_integer_exp();
    // p.second += e; 
    // Gmpfr result(p); 
    

    Gmpfr result; 
    if (e >= 0 ){
      mpfr_mul_2si (result.fr(), a.fr(), e, mpfr_get_default_rounding_mode());
      CGAL_postcondition(a * CGAL::ipower(Gmpfr(2),e) == result);
    }
    else{
      mpfr_div_2si (result.fr(), a.fr(), -e, mpfr_get_default_rounding_mode());
      CGAL_postcondition(a / CGAL::ipower(Gmpfr(2),-e) != result);
    }
    return result;
  }
};
};
} //namespace internal

#endif 

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_KERNEL_D_FLOAT_TRAITS_H
