#include <CGAL/basic.h>
#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Gmpz.h>

#ifndef CGAL_USE_CORE
#include <iostream>
int main ()
{
  std::cout << "Sorry, this example needs CORE ..." << std::endl; 
  return 0;
}
#else

#include "dummy_Bezier_traits.h"


typedef CGAL::Algebraic_kernel_d_1<CGAL::Gmpz>		AK;
typedef CGAL::Arr_Bezier_curve_traits_2<AK>		Bezier_traits_2;
typedef Bezier_traits_2::Bezier_cache			Bez_cache;
typedef Bez_cache::Polynomial_1				Polynomial_1;
typedef Bez_cache::Polynomial_traits_1			Poly_traits_1;
typedef Bez_cache::Algebraic_kernel_d_1			Alg_kernel_1;
typedef Bez_cache::Integer				Integer;

int main (int argc, char *argv[])
{

	Bez_cache	bez_cache;
	Polynomial_1	poly (2,2);
	Polynomial_1	poly1 (3,6);
	Alg_kernel_1	alg_kernel_1;
	bez_cache.get_vertical_tangencies (0, poly, 1);

	bool b = true;
	bez_cache.get_intersections (0, poly, 1, poly1, 1, 0, poly, 1, poly1, 1, b);

	bez_cache.mark_as_overlapping (0,0);
/*
  // Get the name of the input file from the command line, or use the default
  // Bezier.dat file if no command-line parameters are given.
  const char   *filename = (argc > 1) ? argv[1] : "Bezier.dat";

  // Open the input file.
  std::ifstream   in_file (filename);

  if (! in_file.is_open()) {
    std::cerr << "Failed to open " << filename << std::endl;
    return 1;
  }
*/
return 0;
}

#endif
