/*
** Complex_embeddable_traits.h
** 
** Made by Alexander Kobel
** Login   <perpeduumimmobile@lidinoid>
** 
** Started on  Sun Jan  3 12:19:22 2010 Alexander Kobel
** Last update Sun Jan  3 12:19:22 2010 Alexander Kobel
*/

/****************************************************************
 * TODO:
 * - what should the equivalent of Quadrant be for Interval/uncertain
 *   types of T?
 * - should we offer Is_east(), Is_south(), ... functors for each
 *   direction, like Is_positive in RET?
 * - best representation of values in Quadrant (three options shown)?
 * - names:
 *   - directions in Quadrant (I'm fine with the compass names.)
 *   - Modulus vs. Abs vs. Radius vs. Norm
 *   - Argument/Arg vs. Phase vs. Angle vs. Azimuth
 *     (I deprecate Argument_type to avoid confusion with
 *     argument_type in std::unary_function< argument_type,
 *     return_type >. I like Azimuth)
 *   - Norm vs. Squared_abs (std::complex uses norm())
 * - What is expected from Unit_part?
 *   -> Unit_part of Gaussian integers is unclear
 *   -> z / |z| implicitly uses abs()...
 * - Should we offer Sign(z) = z / |z|?
 *   (cf. http://mathworld.wolfram.com/Sign.html )  
 * - Check/implement default implementations of Gcd, Div, ...
 * - Should we design CET s.t. Real_embeddable_traits can inherit
 *   from CET?
 ****************************************************************/

#ifndef   	CGAL_COMPLEX_EMBEDDABLE_TRAITS_H_
# define   	CGAL_COMPLEX_EMBEDDABLE_TRAITS_H_

#include <CGAL/basic.h>
#include <CGAL/Interval_nt.h>
#include <boost/detail/select_type.hpp>

namespace CGAL {

/*
enum Quadrant {
  // Enumerate directions s.t. -d is on the opposite side of d
  ORIGIN     = 0,
  EAST       =  1,
  WEST       = -1,
  NORTH_EAST =  2,
  SOUTH_WEST = -2,
  NORTH      =  3,
  SOUTH      = -3,
  NORTH_WEST =  4,
  SOUTH_EAST = -4
};
*/
// ORIGIN conflicts with CGAL/Origin.h
/*
enum Quadrant {
  // Enumerate directions in ccw order
  ORIGIN     = 0,
  EAST       = 1,
  NORTH_EAST = 2,
  NORTH      = 3,
  NORTH_WEST = 4,
  WEST       = 5,
  SOUTH_WEST = 6,
  SOUTH      = 7,
  SOUTH_EAST = 8
};
*/
/*
enum Quadrant {
  // Two least significant bits [1..0] correspond to the real part,
  // bits [3..2] correspond to the imaginary part.
  // Assumes representation of Sign and Quadrant as signed integers,
  // which is the common way for enums with negative values,
  // but not guaranteed/specified for arbitrary compilers by the C++
  // standard.
  ORIGIN     =  0, // 00:00
  EAST       =  1, // 00:01
  NORTH_EAST =  5, // 01:01
  NORTH      =  4, // 01:00
  NORTH_WEST =  7, // 01:11
  WEST       =  3, // 00:11
  SOUTH_WEST = -1, // 11:11
  SOUTH      = -4, // 11:00
  SOUTH_EAST = -3  // 11:01
};
*/

struct Represents_complex_tag {};
struct Represents_cartesian_complex_tag : public Represents_complex_tag {};
struct Represents_polar_complex_tag : public Represents_complex_tag {};

namespace INTERN_CET {

template< class T, class AST_is_zero >
struct Is_zero_selector { typedef AST_is_zero Type; };

template< class T >
struct Is_zero_selector< T, Null_functor > {
  struct Type : public std::unary_function< T, bool > {
    const bool operator() (const T& x) const {
      return x == T(0);
    }
  };
};

template< class Type_, class Is_complex_embeddable_, class Representation_ >
class Complex_embeddable_traits_base {
public:
  typedef Type_                         Type;
  typedef Is_complex_embeddable_        Is_complex_embeddable;
  typedef Representation_               Representation;

  typedef Null_tag                      Cartesian_type;
  typedef Null_tag                      Modulus_type;
  typedef Null_tag                      Azimuth_type;

  typedef Null_tag                      Boolean;
  typedef Null_tag                      Sign;

  typedef Null_functor Real;
  typedef Null_functor Imag;
  typedef Null_functor Arg;

  typedef Null_functor Abs;
  typedef Null_functor Norm;

  typedef Null_functor Conjugate;

  typedef Null_functor Real_sign;
  typedef Null_functor Imag_sign;
  typedef Null_functor Quadrant;
                        
  typedef Null_functor Real_is_finite;
  typedef Null_functor Imag_is_finite;
  typedef Null_functor Is_finite;
                       
  typedef Null_functor Real_is_positive;
  typedef Null_functor Real_is_negative;
  typedef Null_functor Real_is_zero;
  typedef Null_functor Imag_is_positive;
  typedef Null_functor Imag_is_negative;
  typedef Null_functor Imag_is_zero;
  typedef Null_functor Is_zero;
                       
  typedef Null_functor To_double;
  typedef Null_functor To_interval;
};

template< class Type_ >
class Complex_embeddable_traits_base< Type_, CGAL::Tag_true,
                                      Represents_cartesian_complex_tag > {
public:
  typedef Type_                         Type;
  typedef Tag_true                      Is_complex_embeddable;
  typedef
  Represents_cartesian_complex_tag      Representation;
  
  typedef
  typename Type::Cartesian_type         Cartesian_type;

private:
  typedef Algebraic_structure_traits< Cartesian_type >  Cartesian_AST;
  typedef Real_embeddable_traits< Cartesian_type >      Cartesian_RET;

public:
  typedef
  typename boost::detail::if_true <
    boost::is_base_of< Field_with_sqrt_tag,
                       typename Cartesian_AST::Algebraic_category >::value >
  ::template then< Cartesian_type,
                   Null_tag >::type     Modulus_type;
  // typedef Null_tag                      Modulus_type;

  typedef Null_tag                      Azimuth_type;

  typedef typename Cartesian_AST::Boolean       Boolean;
  typedef typename Cartesian_RET::Sign          Sign;

  struct Real : public std::unary_function< Type, Cartesian_type > {
    Cartesian_type & operator() (Type &z) const             { return z.real(); }
    const Cartesian_type & operator() (const Type &z) const { return z.real(); }
  };
  struct Imag : public std::unary_function< Type, Cartesian_type > {
    Cartesian_type & operator() (Type &z) const             { return z.imag(); }
    const Cartesian_type & operator() (const Type &z) const { return z.imag(); }
  };
  typedef Null_functor                  Arg;

  struct Abs : public std::unary_function< Type, Cartesian_type > {
    const Modulus_type operator() (const Type &z) const {
      return z.abs();
    }
  };
  struct Norm : public std::unary_function< Type, Cartesian_type > {
    const Cartesian_type operator() (const Type &z) const {
      return z.norm();
    }
  };

  struct Conjugate : public std::unary_function< Type, Type > {
    const Type operator() (const Type &z) const { return z.conj(); }
  };

  struct Real_sign : public std::unary_function< Type, Sign > {
    const Sign operator() (const Type &z) const {
      typename Cartesian_RET::Sign sign;
      return sign (Real() (z));
    }
  };
  struct Imag_sign : public std::unary_function< Type, Sign > {
    const Sign operator() (const Type &z) const {
      typename Cartesian_RET::Sign sign;
      return sign (Imag() (z));
    }
  };
  /*
  struct Quadrant : public std::unary_function< Type, ::CGAL::Quadrant > {
    const ::CGAL::Quadrant operator() (const Type &z) const {
      BOOST_STATIC_ASSERT ((boost::is_same< Sign, ::CGAL::Sign >::value));

      const Sign rs = Real_sign() (z);
      const Sign is = Imag_sign() (z);

      switch (rs) {
      case NEGATIVE:
        switch (is) {
        case NEGATIVE: return SOUTH_WEST;
        case ZERO:     return WEST;
        default:       CGAL_expensive_assertion (is == POSITIVE);
                       return NORTH_WEST; }
      case ZERO:
        switch (is) {
        case NEGATIVE: return SOUTH;
        case ZERO:     return ORIGIN;
        default:       CGAL_expensive_assertion (is == POSITIVE);
                       return NORTH; }
      default:
        CGAL_expensive_assertion (rs == POSITIVE);
        switch (is) {
        case NEGATIVE: return SOUTH_EAST;
        case ZERO:     return EAST;
        default:       CGAL_expensive_assertion (is == POSITIVE);
                       return NORTH_EAST; }
      }
    }
  };
  */

  struct Real_is_finite : public std::unary_function< Type, Boolean > {
    const Boolean operator() (const Type &z) const
    { return true; }
  };
  struct Complex_is_finite : public std::unary_function< Type, Boolean > {
    const Boolean operator() (const Type &z) const
    { return true; }
  };
  struct Is_finite : public std::unary_function< Type, Boolean > {
    const Boolean operator() (const Type &z) const
    { return true; }
  };

  struct Real_is_positive : public std::unary_function< Type, Boolean > {
    const Boolean operator() (const Type &z) const
    { return is_positive (Real() (z)); }
  };
  struct Real_is_negative : public std::unary_function< Type, Boolean > {
    const Boolean operator() (const Type &z) const
    { return is_negative (Real() (z)); }
  };
  struct Real_is_zero : public std::unary_function< Type, Boolean > {
    const Boolean operator() (const Type &z) const
    { return is_zero (Real() (z)); }
  };
  struct Imag_is_positive : public std::unary_function< Type, Boolean > {
    const Boolean operator() (const Type &z) const
    { return is_positive (Imag() (z)); }
  };
  struct Imag_is_negative : public std::unary_function< Type, Boolean > {
    const Boolean operator() (const Type &z) const
    { return is_negative (Imag() (z)); }
  };
  struct Imag_is_zero : public std::unary_function< Type, Boolean > {
    const Boolean operator() (const Type &z) const
    { return is_zero (Imag() (z)); }
  };
  struct Is_zero : public std::unary_function< Type, Boolean > {
    const Boolean operator() (const Type &z) const {
      const Boolean real_is_zero = Real_is_zero() (z);
      if (CGAL::certainly_not (real_is_zero))
        return false;
      return real_is_zero & Imag_is_zero() (z);
    }
  };

  struct To_double
    : public std::unary_function< Type, std::complex< double > > {
    const std::complex< double > operator() (const Type &z) const {
      return std::complex< double > (to_double (z.real()),
                                     to_double (z.imag()));
    }
  };
  struct To_interval
    : public std::unary_function< Type, std::complex< std::pair< double, double > > > {
    const std::complex< std::pair< double, double > > operator() (const Type &z) const {
      return std::complex< std::pair< double, double > > (to_interval (z.real()),
                                                          to_interval (z.imag()));
    }
  };
};

} // namespace INTERN_CET

template< class Type >
class Complex_embeddable_traits
  : public INTERN_CET::Complex_embeddable_traits_base< Type,
                                                       Tag_false,
                                                       Null_tag > {};

} // namespace CGAL

#endif 	    /* !CGAL_COMPLEX_EMBEDDABLE_TRAITS_H_ */
