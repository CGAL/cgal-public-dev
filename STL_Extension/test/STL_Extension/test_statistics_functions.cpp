#include <cmath>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <boost/iterator/counting_iterator.hpp>
#include <CGAL/algorithm.h>

bool approx_eq(double a, double b) {
  static double epsilon = 1e-3;
  std::cout << "Got: " << a << ", expected: " << b << std::endl;
  return std::fabs(a - b) <= epsilon;
}

template <typename T>
bool exact_eq(T a, T b) {
  std::cout << "Got: " << a << ", expected: " << b << std::endl;
  return a == b;
}

struct F {
  typedef double result_type;
  double operator()(double x) { return x; }
};

struct sqr {
  // typedef double result_type;
  double operator()(double x) { return x * x; }
};

int main() {

  std::cout << "[Testing statistics functions]" << std::endl;

  std::vector<double> vals;
  std::vector<int> vals_ints;

  std::copy(boost::counting_iterator<std::size_t>(0),
            boost::counting_iterator<std::size_t>(100),
            std::back_inserter(vals));

  std::copy(boost::counting_iterator<std::size_t>(0),
            boost::counting_iterator<std::size_t>(100),
            std::back_inserter(vals_ints));

  std::vector<double>::iterator b = vals.begin(), e = vals.end();
  std::vector<int>::iterator bi = vals_ints.begin(), ei = vals_ints.end();

  std::cout << "Testing 'mean_result'" << std::endl;
  assert(approx_eq(CGAL::mean_result(b, e, F()), 49.5));

  std::cout << "Testing 'max_result'" << std::endl;
  assert(approx_eq(CGAL::max_result(b, e, F()), 99));
  assert(exact_eq<int>(CGAL::max_result(bi, ei, F()), 99));

  std::cout << "Testing 'min_result'" << std::endl;
  assert(approx_eq(CGAL::min_result(b, e, F()), 0));
  assert(exact_eq<int>(CGAL::min_result(bi, ei, F()), 0));

  std::cout << "Testing 'count_result_in_interval'" << std::endl;
  assert(approx_eq(CGAL::count_result_in_interval(b, e, F(), 1, 5), 4));
  assert(approx_eq(CGAL::count_result_in_interval(bi, ei, F(), 0, 10), 10));
  assert(exact_eq<int>(CGAL::count_result_in_interval(bi, ei, F(), 1, 10), 9));
  assert(exact_eq<int>(CGAL::count_result_in_interval(bi, ei, F(), 0, 9), 9));

  std::cout << "Testing 'pearson'" << std::endl;
  assert(approx_eq(CGAL::pearson(vals.begin(), vals.begin() + 10, F(), sqr()),
                   0.9627));

  std::cout << "[Success]" << std::endl;
}