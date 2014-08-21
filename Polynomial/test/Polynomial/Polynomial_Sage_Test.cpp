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
  CGAL::set_pretty_mode(std::cout);

  typedef CGAL::Arithmetic_kernel AK; 
  typedef AK::Integer Integer;
  typedef CGAL::Polynomial_Sage<int> Poly;
  typedef CGAL::Polynomial_traits_d<Poly> PT; 

  //CGAL::Test_Pol::test_multiple_dimensions(PT());  
  Poly k(-5,-2,3,-6,-7,3,-2,4); 
  Poly g(0,-1,5,-7,5);

  Poly f(-5);

  std::cout << "Return dummy degree: " << f.degree_dummy() << std::endl;

  return 0;
}

