//  (C) Copyright Herve Bronnimann 2004.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <random>
#include <utility>
#include <functional>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <vector>
#include <list>
#include <set>
#include <iostream>
#include <iomanip>

#include <boost/timer.hpp>
#include <boost/algorithm/minmax_element.hpp>

#include <unistd.h>
#include <cstdlib>
#include <time.h>


namespace cgal {
  template < class ForwardIterator >
  std::pair< ForwardIterator, ForwardIterator >
  min_max_element(ForwardIterator first, ForwardIterator last)
  {
    typedef std::pair< ForwardIterator, ForwardIterator > FP;
    FP result(first, first);
    if (first != last)
      while (++first != last) {
	if (*first < *result.first)
	  result.first = first;
	if (*result.second < *first)
	  result.second = first;
      }
    return result;
  }

  template < class ForwardIterator, class CompareMin, class CompareMax >
  std::pair< ForwardIterator, ForwardIterator >
  min_max_element(ForwardIterator first,
		  ForwardIterator last,
		  CompareMin comp_min,
		  CompareMax comp_max)
  {
    typedef std::pair< ForwardIterator, ForwardIterator > FP;
    FP result(first, first);
    if (first != last)
      while (++first != last) {
	if (comp_min(*first, *result.first))
	  result.first = first;
	if (comp_max(*result.second, *first))
	  result.second = first;
      }
    return result;
  }
}

template <typename Value>
struct long_cmp
{
  bool operator()(Value const& a, Value const& b) const {
    usleep(30);
    return std::less<Value>()(a,b);
  }
};

template <class Value>
struct less_count : std::less<Value> {
  less_count(less_count<Value> const& lc) : _M_counter(lc._M_counter) {}
  less_count(int& counter) : _M_counter(counter) {}
  bool operator()(Value const& a, Value const& b) const {
    ++_M_counter;
    return std::less<Value>::operator()(a,b);
  }
  void reset() {
    _M_counter = 0;
  }
private:
  int& _M_counter;
};

inline int opt_minmax_count(int n) {
  if (n < 2) return 0;
  if (n == 2) return 1;
  return (n%2 == 0) ? 3*(n/2)-1 : 3*(n/2)+1;
}

int repeats = 10;

#define TIMER( n, cmd , cmdname ) \
  t.restart(); \
  for (int i=0; i<repeats; ++i) { cmd ; } \
  std::cout << "    " << std::setprecision(4) \
            << (double)n*repeats/t.elapsed()/1.0E6 \
            << "M items/sec  " << cmdname << "\n"

#define CTIMER( n, cmd , cmdname, count, opt )	\
  t.restart(); lc.reset(); \
  for (int i=0; i<repeats; ++i) { cmd ; } \
  std::cout << "    " << std::setprecision(4) \
            << (double)n*repeats/t.elapsed()/1.0E6 \
            << "M items/sec  " << cmdname \
  << " ("<< (count)/repeats << " vs " << opt << ")\n"

template <class CIterator>
void test_minmax_element(CIterator first, CIterator last, int n, const char* name)
{
  typedef typename std::iterator_traits<CIterator>::value_type vtype;

  boost::timer t;

  std::cout << "  ON " << name << " WITH OPERATOR<()\n";
  TIMER( n, boost::minmax_element(first, last),
  	 "boost::minmax_element" << name << "    ");
  TIMER( n, std::minmax_element(first, last),
  	 "std::minmax_element" << name << "    ");
  TIMER( n, cgal::min_max_element(first, last),
  	 "cgal::min_max_element" << name << "    ");

  // long_cmp<vtype> lcmp;

  // std::cout << "  ON " << name << " WITH long_cmp\n";
  // TIMER( n, boost::minmax_element(first, last, lcmp),
  // 	 "boost::minmax_element" << name << "    ");
  // TIMER( n, std::minmax_element(first, last, lcmp),
  // 	 "std::minmax_element" << name << "    ");
  // TIMER( n, cgal::min_max_element(first, last, lcmp, lcmp),
  // 	 "cgal::min_max_element" << name << "    ");


  std::cout << "  ON " << name << " WITH COUNTING OPERATOR<()\n";
  int i = 0;
  less_count<vtype> lc(i);

  CTIMER( n, boost::minmax_element(first, last, lc),
  	  "boost::minmax_element" << name << "    ",
  	  i, opt_minmax_count(n));
  CTIMER( n, std::minmax_element(first, last, lc),
  	  "std::minmax_element" << name << "    ",
    	  i, opt_minmax_count(n));
  CTIMER( n, cgal::min_max_element(first, last, lc, lc),
  	  "cgal::min_max_element" << name << "    ",
  	  i, opt_minmax_count(n));
}

template <class Container, class Iterator, class Value>
void test_container(Iterator first, Iterator last, int n, const char* name)
{
  Container c(first, last);
  typename Container::iterator fit(c.begin()), lit(c.end());
  test_minmax_element(fit, lit, n, name);
}

template <class Iterator>
void test_range(Iterator first, Iterator last, int n)
{
  typedef typename std::iterator_traits<Iterator>::value_type Value;
  // Test various containers with these values
  test_container< std::vector<Value>, Iterator, Value >(first, last, n, "<vector>");
  // test_container< std::list<Value>,   Iterator, Value >(first, last, n, "<list>  ");
  // test_container< std::multiset<Value>,    Iterator, Value >(first, last, n, "<set>   ");
}

template <class Value>
void test(int n)
{
  // Populate test vector with identical values
  std::cout << "IDENTICAL VALUES...   \n";
  std::vector<Value> test_vector(n, Value(1));
  typename std::vector<Value>::iterator first( test_vector.begin() );
  typename std::vector<Value>::iterator last( test_vector.end() );
  test_range(first, last, n);

  // Populate test vector with two values
  std::cout << "TWO DISTINCT VALUES...\n";
  typename std::vector<Value>::iterator middle( first + n/2 );
  std::fill(middle, last, Value(2));
  test_range(first, last, n);

  // Populate test vector with increasing values
  std::cout << "INCREASING VALUES...  \n";
  test_vector.clear();
  test_vector.reserve(n);
  for(int i = 0; i < n; ++i)
    test_vector.push_back(i);
  first = test_vector.begin();
  last = test_vector.end();
  test_range(first, last, n);
  
  // Populate test vector with decreasing values
  std::cout << "DECREASING VALUES...  \n";
  std::reverse(first, last);
  test_range(first, last, n);

  // Populate test vector with random values
  std::cout << "RANDOM VALUES...      \n";

  // std::uniform_int_distribution<Value> distribution(Value(0), Value(99));
  // std::mt19937 engine;

  test_vector.clear();
  test_vector.reserve(n);
  // std::generate_n(std::back_inserter(test_vector), n, std::bind(distribution, engine));
  
  srand(time(NULL));

  for(auto i = 0; i < n; ++i) {
    test_vector.push_back(rand());
  }

  first = test_vector.begin();
  last = test_vector.end();
  test_range(first, last, n);
}

int
main(int argc, char** argv)
{
  int n = 100;
  if (argc > 1) n = atoi(argv[1]);
  if (argc > 2) repeats = atoi(argv[2]);

  test<int>(n);

  return 0;
}
