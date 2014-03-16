#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/Polynomial_Sage_type_generator.h>
#include <CGAL/Polynomial/Sage/Polynomial_Sage.h>

int main(){
  CGAL::set_pretty_mode(std::cout);
  typedef CGAL::Polynomial_Sage_type_generator<int,2>::Type Poly_2;
  typedef CGAL::Polynomial_traits_d<Poly_2>            PT_2;
  typedef PT_2::Coefficient_type                       Poly_1;
  typedef PT_2::Innermost_coefficient_type             Integer; 
   
  PT_2::Construct_polynomial construct_polynomial;
  
  // constructing a constant polynomial from int
  Poly_2 two(2); // = 2 
  std::cout << "A constant polynomial: " << two << std::endl;

  
  // construction from an iterator range of univariate polynomials
  
  std::list<Poly_1> univariate_coeffs;
  univariate_coeffs.push_back(Poly_1(3));
  univariate_coeffs.push_back(Poly_1(0));
  univariate_coeffs.push_back(Poly_1(5));
  Poly_2 F = // 5*y^2 + 3
    construct_polynomial(univariate_coeffs.begin(),univariate_coeffs.end());
  std::cout << "The bivariate polynomial F: " << F << std::endl;
  
    std::cout << "The degree of F with respect to y: "<< degree(F)       // = 4 
            << std::endl;

    //  std::cout << "The degree of F with respect to y from Sage: "<< degree_sage(F)       // = 4 
    //      << std::endl;

    
}


