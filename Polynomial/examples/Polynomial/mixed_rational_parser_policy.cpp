
#include <CGAL/config.h>
#include <CGAL/Polynomial/Polynomial_parser_d.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_type_generator.h>

#include <CGAL/Gmpq.h>

int main() {
   typedef CGAL::Polynomial_type_generator<CGAL::Gmpq,2>::Type 
      Polynomial_2;
    
      // parse bivariate polynomials with large integer and rational coefficients
   CGAL::Polynomial_parser_d<Polynomial_2,
            CGAL::Mixed_rational_parser_policy< Polynomial_2 > > parser;
  
   std::string str="x^5*8345345/3234234234*y^17-5/7*(-111345345y-x^2)+10";
    
   Polynomial_2 p2;
   parser(str, p2);
   
   std::cout << "p2: " << p2 << "\n\n";
   return 1;
}

