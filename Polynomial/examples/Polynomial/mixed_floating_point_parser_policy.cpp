
#include <CGAL/config.h>
#include <CGAL/Polynomial/Polynomial_parser_d.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_type_generator.h>

#include <CGAL/CORE_BigRat.h>

int main() {
   typedef CGAL::Polynomial_type_generator<CORE::BigRat,2>::Type 
      Polynomial_2;
    
      // parse bivariate polynomials with large integer, rational and fp coefficients
   CGAL::Polynomial_parser_d<Polynomial_2,
            CGAL::Mixed_floating_point_parser_policy< Polynomial_2 > > parser;
  
   std::string str="-13.003274*(x+y)^2-123.123e-10123x+112/233*(y+x^2-2323)^3";
    
   Polynomial_2 p2;
   parser(str, p2);
   
   std::cout << "p2: " << p2 << "\n\n";
   return 1;
}

