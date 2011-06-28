// Copyright (c) 2006-2007 Max-Planck-Institute Saarbruecken (Germany).
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
// $URL$
// $Id$
//
//
// Author(s)     : Michael Hemmer    <hemmer@mpi-inf.mpg.de>
//
// =============================================================================

/*! \file NiX/Coercion_traits.h
 *  \brief Defines class NiX::Coercion_traits. 
 * 
 *  Provides the general definition of the \c Coercion_traits<A,B> class, with
 *  specializations for the builtin number types.
 */

#ifndef CGAL_COERCION_TRAITS_H
#define CGAL_COERCION_TRAITS_H 1

#include <CGAL/number_type_basic.h>


#include <iterator>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/type_traits/is_same.hpp>

#include <CGAL/compiler_config.h>

#if !defined(CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE) && !defined(CGAL_CFG_NO_CPP0X_DECLTYPE)
#define CGAL_NEW_COERCION_TRAITS 1
#include <type_traits>
#include <boost/mpl/if.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/or.hpp>
#endif

#if defined(CGAL_NEW_COERCION_TRAITS)
namespace CGAL {

// Conversion_traits work automatically for all types that are
// compatible with std::common_type
template<typename A, typename B>
struct Coercion_traits
{
  typedef typename std::common_type<A, B>::type Type;

  // If common_type works, the types are not necessarily implicitly
  // convertible to each other, as there could be a partial
  // specialization of common_type thus we have to check if the
  // implicit conversion is true. Are_explicit_interoperable is true in any
  // case, we check to be consistent though.

  // note: we don't care if A and B are convertible from and to each
  // other, but if they are convertible to the common_type.
  typedef typename boost::mpl::if_< 
    typename boost::mpl::and_< typename std::is_convertible<A, Type>::type,
                               typename std::is_convertible<B, Type>::type >::type,
    CGAL::Tag_true,
    CGAL::Tag_false >::type Are_implicit_interoperable;

  typedef typename boost::mpl::if_< 
    typename boost::mpl::and_< typename std::is_explicitly_convertible<A, Type>::type,
                               typename std::is_explicitly_convertible<B, Type>::type >::type,
    CGAL::Tag_true,
    CGAL::Tag_false >::type Are_explicit_interoperable;

  // if is_same<A,B>::value == true, we must only provide one operator()
  template<typename T, typename U = CGAL::Tag_false>
  struct Cast_i {
    typedef Type result_type;
    Type operator()(T&& x) { return Type(std::forward<T>(x)); }
    Type operator()(U&& x) { return Type(std::forward<U>(x)); }
    Type operator()(const T& x) { return Type(x); }
    Type operator()(const U& x) { return Type(x); }
  };

  typedef typename boost::mpl::if_<
    typename std::is_same<A, B>::type,
    Cast_i<A>,
    Cast_i<A, B> >::type Cast;
};
} // namespace CGAL
#endif // CGAL_NEW_COERCION_TRAITS

// Makro to define an additional operator for binary functors which takes
// two number types as parameters that are interoperable with the
// number type
#define CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT( NT, Result_type  ) \
  template < class CT_Type_1, class CT_Type_2 >                         \
  Result_type operator()( const CT_Type_1& x, const CT_Type_2& y ) const { \
    BOOST_STATIC_ASSERT((::boost::is_same<                              \
            typename Coercion_traits< CT_Type_1, CT_Type_2 >::Type, NT  \
            >::value));                                                 \
                                                                        \
    typename Coercion_traits< CT_Type_1, CT_Type_2 >::Cast cast;        \
    return operator()( cast(x), cast(y) );                              \
  }

#define CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( NT ) \
CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT( NT, NT )

#if !defined(CGAL_NEW_COERCION_TRAITS)
#define CGAL_DEFINE_COERCION_TRAITS_FROM_TO(FROM,TO)    \
  template <>                                           \
  struct Coercion_traits< FROM , TO >{                  \
    typedef Tag_true  Are_explicit_interoperable;       \
    typedef Tag_true  Are_implicit_interoperable;       \
    typedef TO Type;                                    \
    struct Cast{                                        \
      typedef Type result_type;                         \
      Type operator()(const TO& x)   const { return x;} \
      Type operator()(const FROM& x) const {            \
        return Type(x);}                                \
    };                                                  \
  };                                                    \
  template <>                                           \
  struct Coercion_traits< TO , FROM >{                  \
    typedef Tag_true  Are_explicit_interoperable;       \
    typedef Tag_true  Are_implicit_interoperable;       \
    typedef TO Type;                                    \
    struct Cast{                                        \
      typedef Type result_type;                         \
      Type operator()(const TO& x)   const { return x;} \
      Type operator()(const FROM& x) const {            \
        return Type(x);}                                \
    };                                                  \
  };      
#else
#define CGAL_DEFINE_COERCION_TRAITS_FROM_TO(FROM,TO)
#endif

#if !defined(CGAL_NEW_COERCION_TRAITS)
#define CGAL_DEFINE_COERCION_TRAITS_FROM_TO_TEM(FROM,TO,TEM)            \
    template <TEM>                                                      \
    struct Coercion_traits< FROM , TO >{                                \
        typedef Tag_true  Are_explicit_interoperable;                   \
        typedef Tag_true  Are_implicit_interoperable;                   \
        typedef TO Type;                                       \
        struct Cast{                                                    \
            typedef Type result_type;                          \
            Type operator()(const TO& x)   const { return x;}  \
            Type operator()(const FROM& x) const {             \
                return Type(x);}                               \
        };                                                              \
    };                                                                  \
    template <TEM>                                                      \
    struct Coercion_traits< TO , FROM >{                                \
        typedef Tag_true  Are_explicit_interoperable;                   \
        typedef Tag_true  Are_implicit_interoperable;                   \
        typedef TO Type;                                       \
        struct Cast{                                                    \
            typedef Type result_type;                          \
            Type operator()(const TO& x)   const { return x;}  \
            Type operator()(const FROM& x) const {             \
                return Type(x);}                               \
        };                                                              \
    };   
#else
#define CGAL_DEFINE_COERCION_TRAITS_FROM_TO_TEM(FROM,TO,TEM)
#endif
                                                 
#if !defined(CGAL_NEW_COERCION_TRAITS)
#define CGAL_DEFINE_COERCION_TRAITS_FOR_SELF(A)                         \
    template <>                                                         \
    struct Coercion_traits< A , A >{                                    \
        typedef Tag_true  Are_explicit_interoperable;                   \
        typedef Tag_true  Are_implicit_interoperable;                   \
        typedef A Type;                                        \
        struct Cast{                                                    \
            typedef Type result_type;                          \
            Type operator()(const A& x) const { return x;}     \
        };                                                              \
    };    
#else
#define CGAL_DEFINE_COERCION_TRAITS_FOR_SELF(A)
#endif 

#if !defined(CGAL_NEW_COERCION_TRAITS)
#define CGAL_DEFINE_COERCION_TRAITS_FOR_SELF_TEM(A,TEM)                 \
    template <TEM>                                                      \
    struct Coercion_traits< A , A >{                                    \
        typedef Tag_true  Are_explicit_interoperable;                   \
        typedef Tag_true  Are_implicit_interoperable;                   \
        typedef A Type;                                                 \
        struct Cast{                                                    \
            typedef Type result_type;                          \
            Type operator()(const A& x) const {return x;}      \
        };                                                              \
    };
#else    
#define CGAL_DEFINE_COERCION_TRAITS_FOR_SELF_TEM(A,TEM)
#endif


namespace CGAL {


namespace INTERN_CT{ 
template< class FROM, class TO >struct Cast_from_to{
    typedef TO result_type;
    TO operator()(const TO& x) const {return x;}
    TO operator()(const FROM& x) const {return TO(x);}
};
template< class TO>
struct Cast_from_to<TO,TO>{
    typedef TO result_type;
    TO operator()(const TO& x) const {return x;}
};
}


template<class A , class B> struct Coercion_traits;
template<class A , class B, int > struct Coercion_traits_for_level;
    


CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short,int)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short,long)
#ifdef CGAL_USE_LONG_LONG
  CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short,long long)
#endif
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short,float)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short,double)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short,long double)
        
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int,long)
#ifdef CGAL_USE_LONG_LONG
  CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int,long long)
#endif
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int,float)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int,double)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int,long double)

#ifdef CGAL_USE_LONG_LONG
  CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long,long long)
#endif
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long,float)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long,double)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long,long double)

#ifdef CGAL_USE_LONG_LONG
  CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long long,float)
  CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long long,double)
  CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long long,long double)
#endif

CGAL_DEFINE_COERCION_TRAITS_FROM_TO(float,double)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(float,long double)
      
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(double,long double)

//! Specialization for equal types.
#if !defined(CGAL_NEW_COERCION_TRAITS)
template <class A>    
struct Coercion_traits<A,A>{ 
    typedef Tag_true Are_explicit_interoperable;
    typedef Tag_true Are_implicit_interoperable;
    typedef A Type; 
    struct Cast{                                        
        typedef Type result_type;                             
        Type inline operator()(const A& x) const { 
            return x;
        }       
    };
};
#endif
    
CGAL_DEFINE_COERCION_TRAITS_FOR_SELF(short)
CGAL_DEFINE_COERCION_TRAITS_FOR_SELF(int)  
CGAL_DEFINE_COERCION_TRAITS_FOR_SELF(long)
#ifdef CGAL_USE_LONG_LONG
  CGAL_DEFINE_COERCION_TRAITS_FOR_SELF(long long)
#endif
CGAL_DEFINE_COERCION_TRAITS_FOR_SELF(float)
CGAL_DEFINE_COERCION_TRAITS_FOR_SELF(double)
CGAL_DEFINE_COERCION_TRAITS_FOR_SELF(long double)

enum COERCION_TRAITS_LEVEL {
    CTL_TOP          = 4,
    CTL_POLYNOMIAL   = 4, 
    CTL_COMPLEX      = 3,
    CTL_INTERVAL     = 2,
    CTL_SQRT_EXT     = 1 
};

template <class A, class B, int i > 
struct Coercion_traits_for_level: public Coercion_traits_for_level<A,B,i-1>{};

template <class A, class B> 
struct Coercion_traits_for_level<A,B,0> {
    typedef Tag_false Are_explicit_interoperable;
    typedef Tag_false Are_implicit_interoperable;
//    typedef Null_type               Type;
    typedef Null_functor Cast;
};

#if !defined(CGAL_NEW_COERCION_TRAITS)
template<class A , class B> 
struct Coercion_traits :public Coercion_traits_for_level<A,B,CTL_TOP>{};
#endif
 
} //namespace CGAL

#endif //NiX_COERCION_TRAITS_H
