/*
** Bounds.h
** 
** Made by Alexander Kobel
** Login   <perpeduumimmobile@lidinoid>
** 
** Started on  Wed Apr  7 23:19:49 2010 Alexander Kobel
** Last update Wed Apr  7 23:19:49 2010 Alexander Kobel
*/

#ifndef   	BOUNDS_H_
# define   	BOUNDS_H_

#include <climits>
#include <iterator>

#include <CGAL/basic.h>
#include <CGAL/number_utils.h>
#include <CGAL/Algebraic_kernel_d/Real_embeddable_extension.h>

namespace CGAL { 

namespace internal {

#ifdef CGAL_USE_LEDA

template<>
class Real_embeddable_extension< leda_rational > {
public:
  typedef leda_rational Type;
	typedef leda_integer Int;
	typedef Real_embeddable_extension< Int > REE_Int;

  struct Ceil_log2_abs
    : public std::unary_function< leda_rational, long > {
    long operator()( const leda_rational& x ) const {
			REE_Int::Ceil_log2_abs int_ceil_log2_abs;
			REE_Int::Floor_log2_abs int_floor_log2_abs;
			return int_ceil_log2_abs (leda::ceil (abs (x)));
			// return int_ceil_log2_abs (x.numerator())
			// 	- int_floor_log2_abs (x.denominator());
    }
  };

  struct Floor_log2_abs
    : public std::unary_function< leda_rational, long > {
    long operator()( const leda_rational& x ) const {
			REE_Int::Ceil_log2_abs int_ceil_log2_abs;
			REE_Int::Floor_log2_abs int_floor_log2_abs;
			return int_floor_log2_abs (leda::floor (abs (x)));
			// return int_floor_log2_abs (x.numerator())
			// 	- int_ceil_log2_abs (x.denominator());
    }
  };

  struct Floor
    : public std::unary_function< leda_rational, leda_rational > {
    leda_rational operator() (const leda_rational& x) const { return x;}
  };
  struct Ceil
    : public std::unary_function< leda_rational, leda_rational > {
    leda_rational operator() (const leda_rational& x) const { return x;}
  };
};

#endif // CGAL_USE_LEDA

struct Fujiwara_upper_bound_log2 {
	typedef const long value_type;

	template< class RandomAccessIterator >
	value_type operator() (RandomAccessIterator first,
												 const RandomAccessIterator beyond) const;

	template< class Polynomial >
	value_type operator() (const Polynomial &f) const {
		return operator() (f.begin(), f.end());
	}
};

struct Fujiwara_lower_bound_log2 {
	typedef const long value_type;

	template< class RandomAccessIterator >
	value_type operator() (RandomAccessIterator first,
												 const RandomAccessIterator beyond) const;

	template< class Polynomial >
	value_type operator() (const Polynomial &f) const {
		return operator() (f.begin(), f.end());
	}
};

// IMPLEMENTATION

template< class RandomAccessIterator >
Fujiwara_upper_bound_log2::value_type
Fujiwara_upper_bound_log2::operator() (RandomAccessIterator first,
																			 RandomAccessIterator beyond) const {
	typedef typename std::iterator_traits< RandomAccessIterator >::value_type NT;

	// using namespace CGALi; // FIXME: remove this (deprecated, should be internal)
	typename Real_embeddable_extension< NT >::Ceil_log2_abs ceil_log2_abs;
	typename Real_embeddable_extension< NT >::Floor_log2_abs floor_log2_abs;
	
	// TODO:
	// Is it clever to cut off trailing zeroes in the polynomial?
// 	while (first != beyond && CGAL::is_zero (*first))
// 		++first;
	while (first != beyond && CGAL::is_zero (*(beyond-1)))
		--beyond;

	if (std::distance (first, beyond) < 2)
		return 0;

	const size_t n = beyond - first - 1; // degree
	const long an_log2 = floor_log2_abs (*(first+n));
	
	long fb_log2 = std::numeric_limits< long >::min();
	
	for (size_t i = 1; i < n; ++i) {
		if (CGAL::is_zero (*(first+n-i)))
			continue;
		
		fb_log2 = CGAL::max (fb_log2,
												 (ceil_log2_abs (*(first+n-i)) - an_log2) / static_cast< long > (i));
	}
	
	if (! CGAL::is_zero (*first))
		fb_log2 = CGAL::max (fb_log2,
												 (ceil_log2_abs (*first) - an_log2) / static_cast< long > (n) + 1);
	
	return fb_log2 + 2;
}

template< class RandomAccessIterator >
Fujiwara_lower_bound_log2::value_type
Fujiwara_lower_bound_log2::operator() (const RandomAccessIterator first,
																			 const RandomAccessIterator beyond) const {
	typedef std::reverse_iterator< RandomAccessIterator > reverse_iterator;
	return Fujiwara_upper_bound_log2() (reverse_iterator (beyond),
																			reverse_iterator (first));
}

} // namespace internal

} // namespace CGAL

#endif 	    /* !BOUNDS_H_ */
