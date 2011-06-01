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
#include <cassert>

#include <boost/timer.hpp>
#include <boost/algorithm/minmax_element.hpp>

#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#endif

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
    usleep(1);
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

int repeats = 10;

#define TIMER( n, cmd , cmdname ) \
  t.restart(); \
  for (int i=0; i<repeats; ++i) { cmd ; } \
  std::cout << "    " << std::setprecision(4) \
            << (double)n*repeats/t.elapsed()/1.0E6 \
            << "M items/sec  " << cmdname << "\n"

#define CTIMER( n, cmd , cmdname, count, opt )	\
  t.restart(); lc.reset(); \
  for (int i=0; i<repeats; ++i) { cmd; } \
  std::cout << "    " << std::setprecision(4) \
            << (double)n*repeats/t.elapsed()/1.0E6 \
            << "M items/sec"  << " ("<< (count)/repeats << ") "	\
            << cmdname << "\n"


std::string tests;

template <class CIterator>
void test_minmax_element(CIterator first, CIterator last, int n, const char* name)
{
  typedef typename std::iterator_traits<CIterator>::value_type vtype;
 
  boost::timer t;
  std::pair<CIterator,CIterator> res1;
  std::pair<CIterator,CIterator> res2;
  std::pair<CIterator,CIterator> res3;
 
  if(tests == "operator<") {
    std::cout << "  ON " << name << " WITH OPERATOR<()\n";
  
    TIMER( n,res1=boost::minmax_element(first, last),
	   "boost::minmax_element" << name << "    ");

    #ifndef CGAL_CFG_NO_CPP0X_MINMAX_ELEMENT
    TIMER( n, res2=std::minmax_element(first, last),
	   "std::minmax_element" << name << "    ");
    #endif

    TIMER( n, res3=cgal::min_max_element(first, last),
	   "cgal::min_max_element" << name << "    ");

    if(!(*(res1.first) == *(res2.first) && *(res1.second) == *(res2.second)
    	                                && 
    	 *(res1.first) == *(res3.first) && *(res1.second) == *(res3.second)
                                        &&
    	 *(res2.first) == *(res3.first) && *(res2.second) == *(res3.second)))
      std::cerr << "Different results";
  } else if(tests == "long_cmp") {
    long_cmp<vtype> lcmp;
 
    std::cout << "  ON " << name << " WITH long_cmp\n";

    TIMER( n, res1=boost::minmax_element(first, last, lcmp),
	   "boost::minmax_element" << name << "    ");

    #ifndef CGAL_CFG_NO_CPP0X_MINMAX_ELEMENT
    TIMER( n, res2=std::minmax_element(first, last, lcmp),
	   "std::minmax_element" << name << "    ");
    #endif

    TIMER( n, res3=cgal::min_max_element(first, last, lcmp, lcmp),
	   "cgal::min_max_element" << name << "    ");

    if(!(*(res1.first) == *(res2.first) && *(res1.second) == *(res2.second)
    	                                && 
    	 *(res1.first) == *(res3.first) && *(res1.second) == *(res3.second)
    	                                &&
    	 *(res2.first) == *(res3.first) && *(res2.second) == *(res3.second)))
      std::cerr << "Different results";
  } else if(tests == "counting_cmp") {
    std::cout << "  ON " << name << " WITH COUNTING OPERATOR<()\n";
    int i = 0;
    less_count<vtype> lc(i);
 
    CTIMER( n, res1=boost::minmax_element(first, last, lc),
	    "boost::minmax_element" << name << "    ",
	    i, opt_minmax_count(n));

    #ifndef CGAL_CFG_NO_CPP0X_MINMAX_ELEMENT
    CTIMER( n, res2=std::minmax_element(first, last, lc),
	    "std::minmax_element" << name << "    ",
	    i, opt_minmax_count(n));
    #endif

    CTIMER( n, res3=cgal::min_max_element(first, last, lc, lc),
	    "cgal::min_max_element" << name << "    ",
	    i, opt_minmax_count(n));
          
    if(!(*(res1.first) == *(res2.first) && *(res1.second) == *(res2.second)
    	 && 
    	 *(res1.first) == *(res3.first) && *(res1.second) == *(res3.second)
    	 &&
    	 *(res2.first) == *(res3.first) && *(res2.second) == *(res3.second)))
      std::cerr << "Different results";
  }
}

template <class Container, class Iterator, class Value>
void test_container(Iterator first, Iterator last, int n, const char* name)
{
  Container c(first, last);
  typename Container::iterator fit(c.begin()), lit(c.end());
  test_minmax_element(fit, lit, n, name);
}

std::string testContainer = "vector";

template <class Iterator>

void test_range(Iterator first, Iterator last, int n)
{
  typedef typename std::iterator_traits<Iterator>::value_type Value;
  // Test various containers with these values
  if(testContainer == "vector") {
    test_container< std::vector<Value>, Iterator, Value >(first, last, n, "<vector>");
  } else if(testContainer == "list") {
    test_container< std::list<Value>,   Iterator, Value >(first, last, n, "<list>  ");
  } else if(testContainer == "set") {
    test_container< std::multiset<Value>,    Iterator, Value >(first, last, n, "<set>   ");
  }
}

template <class Value>
void test(int n)
{
  // Populate test vector with identical values
  std::cout << "IDENTICAL VALUES   \n";
  std::vector<Value> test_vector(n, Value(1));
  typename std::vector<Value>::iterator first( test_vector.begin() );
  typename std::vector<Value>::iterator last( test_vector.end() );
  test_range(first, last, n);

  // Populate test vector with two values
  std::cout << "TWO DISTINCT VALUES\n";
  typename std::vector<Value>::iterator middle( first + n/2 );
  std::fill(middle, last, Value(2));
  test_range(first, last, n);

  // Populate test vector with increasing values
  std::cout << "INCREASING VALUES  \n";
  test_vector.clear();
  test_vector.reserve(n);
  
  std::generate_n(std::back_inserter(test_vector), n, 
  		  []() { static int i = 0; return ++i; });
  first = test_vector.begin();
  last = test_vector.end();
  test_range(first, last, n);
  
  // Populate test vector with decreasing values
  std::cout << "DECREASING VALUES  \n";
  std::reverse(first, last);
  test_range(first, last, n);

  // Populate test vector with random values
  std::cout << "RANDOM VALUES      \n";
  std::uniform_int_distribution<Value> distribution(Value(0), Value(n));
  std::mt19937 engine;

  test_vector.clear();
  test_vector.reserve(n);
  std::generate_n(std::back_inserter(test_vector), n, std::bind(distribution, engine));
  
  first = test_vector.begin();
  last = test_vector.end();
  test_range(first, last, n);
}

int
main(int argc, char** argv)
{
  int n;
  
  po::options_description desc("Allowed options for reload");
  desc.add_options()
    ("help", "produce help message")
    ("n", po::value<int>(&n)->default_value(100), "container size, default 100")
    ("r", po::value<int>(&repeats)->default_value(10), "repetitions, default 10")
    ("test", po::value<std::string>(&tests)->default_value("operator<"), 
     "which test: operator<, long_cmp, counting_cmp; default operator<")
    ("container", po::value<std::string>(&testContainer)->default_value("vector"),
     "which container: vector, list, set; default: vector")
    ;
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);    

  if(vm.count("help")) {
    std::cout << desc << std::endl;
    return EXIT_SUCCESS;
  }

  if(tests != "operator<"
     && tests != "long_cmp"
     && tests != "counting_cmp") {
    std::cout << "unknown tests option" << std::endl;
    return EXIT_FAILURE;
  }

  if(testContainer != "vector"
     && testContainer != "set"
     && testContainer != "list") {
    std::cout << "unknown testContainer option" << std::endl;
    return EXIT_FAILURE;
  }

  test<int>(n);

  return 0;
}
