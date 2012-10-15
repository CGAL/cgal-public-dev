// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
//
//
// Author(s)     :  Sarah Schaeffer
//
// ============================================================================


#ifndef CGAL_BITSTREAM_ISOLATORS_BITSTREAM_COEFFICIENT_H
#define CGAL_BITSTREAM_ISOLATORS_BITSTREAM_COEFFICIENT_H 1

#include <CGAL/basic.h>

#include <CGAL/Algebraic_kernel_d/Real_embeddable_extension.h>

namespace CGAL {

namespace internal {

// TODO add functor_object functions

template < class NT_ >
class Bitstream_coefficient {

public:

  //! this instance's template parameter
  typedef NT_ NT;

  //! the Arithmetic kernel
  typedef typename Get_arithmetic_kernel< NT >::Arithmetic_kernel 
  Arithmetic_kernel;

  //! the Rational number type
  typedef typename Arithmetic_kernel::Rational Rational;
  
  //! the type of a bigfloat interval
  typedef typename Arithmetic_kernel::Bigfloat_interval Bigfloat_interval;
  

  //! Functor to approximate bitstream as rational interval
  typedef CGAL::Null_functor Rational_approx;
/* e.g.:
  class Rational_approx {
    
  public:
    
    //! operator
    std::pair<Rational,Rational> 
    operator()(const NT& v, int precision_exponent) {
      
      CGAL_error(false);
      return std::make_pair(0,0);
    }

  };
*/
  
  //! Functor to approximate bitstream as bigfloat interval
  typedef CGAL::Null_functor Bigfloat_approx;
/* e.g.:
  class Bigfloat_approx {
    
  public:

    Bigfloat_interval operator()(const NT& v, int precision) {
      CGAL_error(false);
      return Bigfloat_interval(0);
    }
    
  };
*/

  //! Functor to decide whether a coefficient is zero
  typedef Null_functor Is_zero;
/* e.g.:
  class Is_zero {

  public:

    bool operator()(const NT& v) {
      return false;
    }
  };
*/

};


// Partial specializations

} //namespace internal

} //namespace CGAL

#ifdef CGAL_USE_LEDA

#include <CGAL/leda_integer.h>

namespace CGAL {

namespace internal {

template <>
class Bitstream_coefficient< leda_integer > {
  
public:

  //! the instance's number type
  typedef leda_integer NT;

  //! the Arithmetic kernel
  typedef Get_arithmetic_kernel< NT >::Arithmetic_kernel 
  Arithmetic_kernel;

  //! the Rational number type
  typedef Arithmetic_kernel::Rational Rational;
  
  //! the type of a bigfloat interval
  typedef Arithmetic_kernel::Bigfloat_interval Bigfloat_interval;

  
  //! Functor to approximate bitstream as rational interval  
  class Rational_approx {
    
  public:
    
    std::pair< Rational,Rational> 
    operator()(const NT& v, int precision_exponent) {
      
      NT shifted_v = v;
      NT last_shifted = v;
      Rational precision =
            power<Rational, int>(Rational(1) / Rational(2), precision_exponent);
      bool shift_sufficient = false;
      Rational relative_error;
    
      while (!shift_sufficient) {
        if (v == 0) {
          shift_sufficient = true;
          break;
        }
        shifted_v >>= 1;
        relative_error = Rational(CGAL::abs(v - shifted_v)) /
                         Rational(CGAL::abs(v));
        if (CGAL::compare(relative_error, precision) == CGAL::LARGER) {
          shifted_v = last_shifted;
          shift_sufficient = true;
        } 
        last_shifted = shifted_v;
      }
      std::pair<Rational,Rational> result = std::make_pair(v,v);
      return result;

    }

private:

  template <class NTBase, class NTExp>
  /*!
   *  computes base^e
   * \param base base
   * \param e exponent
   * \return base^e
   */
  NTBase power (const NTBase base, NTExp e) {
    NTBase result = base;
    if (e == NTExp(0)) {
      return NTBase(1);
    }
    for (NTExp i = NTExp(0); i < e - NTExp(1); i++) {
      result *= base;
    }
    return result;
  }

  };
  
  //! Is_zero functor
  typedef CGAL::Real_embeddable_traits< NT >::Is_zero Is_zero;
  
};

} //namespace internal

} //namespace CGAL



#include <CGAL/leda_rational.h>

namespace CGAL {

namespace internal {

template <>
class Bitstream_coefficient< leda_rational > {
  
public:

  //! the instance's number type
  typedef leda_rational NT;

  //! the Arithmetic kernel
  typedef Get_arithmetic_kernel< NT >::Arithmetic_kernel 
  Arithmetic_kernel;

  //! the Rational number type
  typedef Arithmetic_kernel::Rational Rational;
  
  //! the type of a bigfloat interval
  typedef Arithmetic_kernel::Bigfloat_interval Bigfloat_interval;

  
  //! Functor to approximate bitstream as rational interval  
  class Rational_approx {
    
  public:
    
    std::pair< Rational,Rational> 
    operator()(const NT& v, int precision_exponent) {
      
      // TODO respect precision_exponent
      std::pair<Rational,Rational> result = std::make_pair(v,v);
      return result;

    }

  };
  
  //! Is_zero functor
  typedef CGAL::Real_embeddable_traits< NT >::Is_zero Is_zero;
  
};

} //namespace internal

} //namespace CGAL

#endif // CGAL_USE_LEDA

#ifdef CGAL_USE_GMP

#include <CGAL/Gmpz.h>

namespace CGAL {

namespace internal {

template <>
class Bitstream_coefficient< Gmpz > {
  
public:

  //! the instance's number type
  typedef Gmpz NT;

  //! the Arithmetic kernel
  typedef Get_arithmetic_kernel< NT >::Arithmetic_kernel 
  Arithmetic_kernel;

  //! the Rational number type
  typedef Arithmetic_kernel::Rational Rational;
  
  //! the type of a bigfloat interval
  typedef Arithmetic_kernel::Bigfloat_interval Bigfloat_interval;

  
  //! Functor to approximate bitstream as rational interval  
  class Rational_approx {
    
  public:
    
    std::pair< Rational,Rational> 
    operator()(const NT& v, int precision_exponent) {
      
//      std::cout << "v: " << v << std::endl;

      /*
      long len = ceil_log2_abs(v);
      std::cout << "l: " << len << std::endl;
      // TODO respect precision_exponent
      NT shift = v;
      std::cout << "diff: " << len-precision_exponent << std::endl;

      long diff = len - precision_exponent;
      if ((diff - 1) > 0) {
        shift >>= (diff - 1);
      }

      std::cout << "shift: " << shift << std::endl;
      std::cout << std::endl;
      */

      NT shifted_v = v;
      NT last_shifted = v;
      Rational precision =
            power<Rational, int>(Rational(1) / Rational(2), precision_exponent);
      bool shift_sufficient = false;
      Rational relative_error;
    
      while (!shift_sufficient) {
        if (v == 0) {
          shift_sufficient = true;
          break;
        }
        shifted_v >>= 1;
        relative_error = Rational(CGAL::abs(v - shifted_v)) /
                         Rational(CGAL::abs(v));
        if (CGAL::compare(relative_error, precision) == CGAL::LARGER) {
          shifted_v = last_shifted;
          shift_sufficient = true;
        } 
        last_shifted = shifted_v;
      }     
//      std::cout << "precision: " << precision << std::endl;
//      std::cout << "result: " << shifted_v << std::endl;

      std::pair<Rational,Rational> result =
                                           std::make_pair(shifted_v, shifted_v);
      return result;

    }

private:

  template <class NTBase, class NTExp>
  /*!
   *  computes base^e
   * \param base base
   * \param e exponent
   * \return base^e
   */
  NTBase power (const NTBase base, NTExp e) {
    NTBase result = base;
    if (e == NTExp(0)) {
      return NTBase(1);
    }
    for (NTExp i = NTExp(0); i < e - NTExp(1); i++) {
      result *= base;
    }
    return result;
  }

  };
  
  //! Is_zero functor
  typedef CGAL::Real_embeddable_traits< NT >::Is_zero Is_zero;
  
};

} //namespace internal

} //namespace CGAL



#include <CGAL/Gmpq.h>

namespace CGAL {

namespace internal {

template <>
class Bitstream_coefficient< Gmpq > {
  
public:

  //! the instance's number type
  typedef Gmpq NT;

  //! the Arithmetic kernel
  typedef Get_arithmetic_kernel< NT >::Arithmetic_kernel 
  Arithmetic_kernel;

  //! the Rational number type
  typedef Arithmetic_kernel::Rational Rational;
  
  //! the type of a bigfloat interval
  typedef Arithmetic_kernel::Bigfloat_interval Bigfloat_interval;

  
  //! Functor to approximate bitstream as rational interval  
  class Rational_approx {
    
  public:
    
    std::pair< Rational,Rational> 
    operator()(const NT& v, int precision_exponent) {
      
      // TODO respect precision_exponent
      std::pair<Rational,Rational> result = std::make_pair(v,v);
      return result;

    }

  };
  
  //! Is_zero functor
  typedef CGAL::Real_embeddable_traits< NT >::Is_zero Is_zero;
  
};

} //namespace internal

} //namespace CGAL

#endif // CGAL_USE_GMP


#include <CGAL/CORE_BigInt.h>

namespace CGAL {

namespace internal {

template <>
class Bitstream_coefficient< CORE::BigInt > {

public:

  //! the instance's number type
  typedef CORE::BigInt NT;

  //! the Arithmetic kernel
  typedef Get_arithmetic_kernel< NT >::Arithmetic_kernel Arithmetic_kernel;

  //! the Rational number type
  typedef Arithmetic_kernel::Rational Rational;
  
  //! the type of a bigfloat interval
  typedef Arithmetic_kernel::Bigfloat_interval Bigfloat_interval;

  //! Functor to approximate bitstream as rational interval
  class Rational_approx {
    
  public:
    
    std::pair< Rational,Rational> 
    operator()(const NT& v, int precision_exponent) {

      NT shifted_v = v;
      NT last_shifted = v;
      Rational precision =
            power<Rational, int>(Rational(1) / Rational(2), precision_exponent);
      bool shift_sufficient = false;
      Rational relative_error;
    
      while (!shift_sufficient) {
        if (v == 0) {
          shift_sufficient = true;
          break;
        }
        shifted_v >>= 1;
        relative_error = Rational(CGAL::abs(v - shifted_v)) /
                         Rational(CGAL::abs(v));
        if (CGAL::compare(relative_error, precision) == CGAL::LARGER) {
          shifted_v = last_shifted;
          shift_sufficient = true;
        } 
        last_shifted = shifted_v;
      }
      std::pair<Rational,Rational> result = std::make_pair(v,v);
      return result;
    }

private:

  template <class NTBase, class NTExp>
  /*!
   *  computes base^e
   * \param base base
   * \param e exponent
   * \return base^e
   */
  NTBase power (const NTBase base, NTExp e) {
    NTBase result = base;
    if (e == NTExp(0)) {
      return NTBase(1);
    }
    for (NTExp i = NTExp(0); i < e - NTExp(1); i++) {
      result *= base;
    }
    return result;
  }

  };
  
  //! Is_zero functor
  typedef CGAL::Real_embeddable_traits< NT >::Is_zero Is_zero;
  
};

} //namespace internal

} //namespace CGAL


#include <CGAL/CORE_BigRat.h>

namespace CGAL {

namespace internal {

template <>
class Bitstream_coefficient< CORE::BigRat > {

public:

  //! the instance's number type
  typedef CORE::BigRat NT;

  //! the Arithmetic kernel
  typedef Get_arithmetic_kernel< NT >::Arithmetic_kernel Arithmetic_kernel;

  //! the Rational number type
  typedef Arithmetic_kernel::Rational Rational;
  
  //! the type of a bigfloat interval
  typedef Arithmetic_kernel::Bigfloat_interval Bigfloat_interval;

  //! Functor to approximate bitstream as rational interval
  class Rational_approx {
    
  public:
    
    std::pair< Rational,Rational> 
    operator()(const NT& v, int precision_exponent) {
      std::pair<Rational,Rational> result = std::make_pair(v,v);
      return result;
    }
  };
  
  //! Is_zero functor
  typedef CGAL::Real_embeddable_traits< NT >::Is_zero Is_zero;
  
};

} //namespace internal

} //namespace CGAL


#include <CGAL/Algebraic_kernel_d/Algebraic_real_d_1.h>

namespace CGAL {

namespace internal {

template< class Coefficient_, class Bound_, 
          class HandlePolicy, class RepClass >
class Bitstream_coefficient< 
 internal::Algebraic_real_d_1< Coefficient_, Bound_, HandlePolicy, RepClass > 
> {

public:
  
  typedef Coefficient_  Coefficient;
  typedef Bound_        Bound;
  typedef HandlePolicy  Handle_policy;
  typedef RepClass      Rep_class;
  
  //! the instance's number type
  typedef internal::Algebraic_real_d_1< Coefficient, Bound, 
                                         Handle_policy, Rep_class > NT;
  
  //! the Arithmetic kernel
  typedef typename Get_arithmetic_kernel< Bound >::Arithmetic_kernel 
  Arithmetic_kernel;
  
  //! the Rational number type
  typedef typename Arithmetic_kernel::Rational Rational;
  
  //! the type of a bigfloat interval
  typedef typename Arithmetic_kernel::Bigfloat_interval Bigfloat_interval;
  
private:

  // TODO introduce 'Get_algebraic_kernel'?
  //! the algebraic kernel
  typedef typename CGAL::Algebraic_kernel_d_1_generator<Coefficient, Rational>::
  Algebraic_kernel_with_qir_and_descartes_1  Algebraic_kernel_1;
  
  //! the approximations
  typedef typename Algebraic_kernel_1::Approximate_relative_1 
  Approximate_relative_1;

public:
  
  //! Functor to approximate bitstream as rational interval
  class Rational_approx {
    
  public:

    std::pair< Rational, Rational > 
    operator()(const NT& v, int precision_exponent) {
      
      Approximate_relative_1 approximation;
      // TODO fix it, as: only works because Rational is equal type as Bound 
      return approximation(v,precision_exponent);
    }

  };
  
  //! Is_zero functor
  typedef typename CGAL::Real_embeddable_traits< NT >::Is_zero Is_zero;

};



} //namespace internal

} //namespace CGAL

#endif // CGAL_BITSTREAM_ISOLATORS_BITSTREAM_COEFFICIENT_H
