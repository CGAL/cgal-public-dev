
#include <CGAL/basic.h>
#include <CGAL/Polynomial/Polynomial_parser_d.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_type_generator.h>

int main() {

   typedefs CGAL::Polynomial_type_generator< double,2 >::Type
      Poly_double_2;
   
      // parse bivariate polynomial with double-precision coefficients
   Poly_double_2 p2;
   CGAL::Polynomial_parser_d< Poly_double_2,
        CGAL::Default_parser_policy< Poly_double_2 > > goofy_2;
   
   std::string str="87623.28374872*xy^10-x*y*8e+8\n\r\t\n+1e-11x+(y-x)^7";
   goofy_2(str, p2);

   std::cout << "p2: " << p2 << "\n\n";
   return 1;
}

