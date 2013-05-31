#include <CGAL/GMP_arithmetic_kernel.h>
#include <iostream>
#include <boost/lexical_cast.hpp>

int main (int argc, char **argv) {
  int max_iter = 1000;
  if (argc > 1)
    max_iter = boost::lexical_cast< int > (argv[1]);

  typedef CGAL::GMP_arithmetic_kernel AK;
  typedef AK::Bigfloat                BF;  // Gmpfr
  typedef AK::Bigfloat_interval       BFI; // Gmpfi

  CGAL::set_precision (BFI(), 256);

  BFI A (BF (-1.), BF (1.));
  BFI B (BF (-1.), BF (1.));
  BFI C;

  for (int i = 0; i < max_iter; ++i) {
    A *= B;
    C = A/B;
    C = A+B;
  }

  std::cout << "C = " << C << std::endl;

  mpfr_free_cache();

  return 0;
}
