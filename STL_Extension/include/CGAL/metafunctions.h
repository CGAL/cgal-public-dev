#ifndef PROPERGY_GENERATOR_METAFUNCTIONS_H
#define PROPERGY_GENERATOR_METAFUNCTIONS_H

#include <type_traits>
#include <boost/mpl/if.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/type_traits.hpp>
#include <CGAL/result_of.h>

namespace CGAL
{
namespace internal
{

/*----------------------------------------------------------------------------*/

// Utility meta-function to help selectively eliminate cases.
template<typename T>
struct to_void { typedef void type; };

/*----------------------------------------------------------------------------*/

// If the specialisations fail, we assume false. Therefore false negatives
// are possible - which is fine for our useage.
template<typename F, typename T, typename enable = void>
struct is_callable_with : boost::mpl::false_{};

/*----------------------------------------------------------------------------*/

// If C++11 is found, we can use decltype and declval.
#if 1

  // This specialisation is selected if the function return type can be deduced
  // when given a paramter of type T.
  template<typename F, typename T>
  struct is_callable_with<
    F,
    T,
    typename to_void<
      decltype( std::declval< F >() ( std::declval<T>() ) ) 
    >::type
  > : boost::mpl::true_{};

/*----------------------------------------------------------------------------*/

// Without C++11, we allow the user to specify the type 'argument_type', since
// it cannot be deduced automatically.
#else

/*----------------------------------------------------------------------------*/

  // If argument type exists, this is enabled and does the check. Otherwise we
  // assume false.
  template<typename F, typename T>
  struct is_callable_with<
    F,
    T,
    typename to_void<typename F::argument_type>::type
  > : boost::is_convertible< T, typename F::argument_type>::type {};

/*----------------------------------------------------------------------------*/

#endif

/*----------------------------------------------------------------------------*/

// We dereference iff the function is only callable with the iterator when 
// applied directly without dereferencing.
template <typename F, typename T>
struct do_dereference 
  :
    // If not( callable directly and not (callable after de-referencing) )
    boost::mpl::not_<
      boost::mpl::and_<        
        is_callable_with<F, T>,
        boost::mpl::not_<
          // Iterator traits doesn't work unless we remove cv qualifiers
          // and reference.
          is_callable_with<F, typename std::iterator_traits<
              typename boost::remove_cv<
                typename boost::remove_reference< T >::type 
              >::type
            >::value_type>
          >
        >
      >
{};

/*----------------------------------------------------------------------------*/

template<typename F, typename I>
struct get_result_type
{
  typedef 
    typename boost::mpl::if_<
      do_dereference<F,I>,
      typename cpp11::result_of<
        F(typename std::iterator_traits<I>::value_type) 
      >::type,
      typename cpp11::result_of<
        F(I)
      >::type
    >::type 
  type;
};

/*----------------------------------------------------------------------------*/

} // internal
} // CGAL

#endif
