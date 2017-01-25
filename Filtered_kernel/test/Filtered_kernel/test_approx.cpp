#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>

int main() {
  std::cerr << typeid(CGAL::Epick::Do_intersect_3::Approximate_predicate).name() << std::endl;
  return 0;
}
