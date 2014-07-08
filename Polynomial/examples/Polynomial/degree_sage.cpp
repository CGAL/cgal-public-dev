#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial/Polynomial_Sage_type.h>

int main(){  
  CGAL::set_pretty_mode(std::cout);
  //typedef Sage_Polynomial<int,2>::Type Poly_2;
  typedef CGAL::Polynomial_Sage<int> Rep;

  Rep two(1);
  //std::cout << two.degree(); 
}


