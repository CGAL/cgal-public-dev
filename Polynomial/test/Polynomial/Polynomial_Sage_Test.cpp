#include <iostream>
#include <cassert>

#include <CGAL/basic.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>

#include <CGAL/Test/_test_polynomial_traits_d.h>
#include <CGAL/Polynomial/Polynomial_Sage_type.h>
int main()
{
  typedef CGAL::Arithmetic_kernel AK; 
  typedef AK::Integer Integer;

  typedef CGAL::Polynomial_Sage<Integer> Poly;
  typedef CGAL::Polynomial_traits_d<Poly> PT; 

  Poly p;
  //std::cout << Poly.degree_dummy() << std::endl;

  return 0;
}

