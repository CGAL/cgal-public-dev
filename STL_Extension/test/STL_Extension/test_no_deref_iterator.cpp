#include <cassert>
#include <list>
#include <vector>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <CGAL/iterator.h>

// Function taking iterator.
struct f {
  int operator()(std::vector<double>::iterator) { return 1729; }
};

int main() {
  std::cout << "[Testing 'No_deref_iterator']" << std::endl;

  // Create a vector of dou
  std::vector<double> vals_vec;
  std::list<double> vals_list;

  std::vector<std::vector<double>::iterator> iters_vec;
  std::vector<std::list<double>::iterator> iters_list;

  for (int i = 0; i < 100; ++i) {
    vals_vec.push_back(i);
    vals_list.push_back(i);
  }

  std::cout << "Testing copy vector" << std::endl;

  typedef std::vector<std::vector<double>::iterator>::iterator Iterator_vec;

  std::copy(CGAL::make_no_deref_iterator(vals_vec.begin()),
            CGAL::make_no_deref_iterator(vals_vec.end()),
            std::back_inserter(iters_vec));

  int n = 0;
  for (Iterator_vec i = iters_vec.begin(); i != iters_vec.end(); ++i)
    assert(n++ == *(*i));

  std::cout << "Testing copy list" << std::endl;

  typedef std::vector<std::list<double>::iterator>::iterator Iterator_list;

  std::copy(CGAL::make_no_deref_iterator(vals_list.begin()),
            CGAL::make_no_deref_iterator(vals_list.end()),
            std::back_inserter(iters_list));

  n = 0;
  for (Iterator_list i = iters_list.begin(); i != iters_list.end(); ++i)
    assert(n++ == *(*i));

  std::cout << "Testing transform" << std::endl;

  for (CGAL::No_deref_iterator<std::vector<double>::iterator> i =
           CGAL::No_deref_iterator<std::vector<double>::iterator>(
               vals_vec.begin());
       i != CGAL::No_deref_iterator<std::vector<double>::iterator>(
                vals_vec.end());
       ++i) {
    f foo;
    assert(foo(*i) == 1729);
  }

  CGAL::No_deref_iterator<std::list<double>::iterator> i1;
  CGAL::No_deref_iterator<std::vector<double>::iterator> i2;

  std::cout << "Testing assign" << std::endl;

  i1 = CGAL::No_deref_iterator<std::list<double>::iterator>(vals_list.begin());
  i2 = CGAL::No_deref_iterator<std::vector<double>::iterator>(vals_vec.begin());

  std::cout << "[Success]" << std::endl;
}