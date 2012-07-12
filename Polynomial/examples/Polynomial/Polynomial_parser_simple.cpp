#include <CGAL/config.h>

#include <CGAL/Polynomial/Polynomial_parser_d.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_type_generator.h>


int main(int argc, char** argv) {

  CGAL::set_pretty_mode(std::cout);

  { // simple parser: parses polynomials over doubles
    typedef typename CGAL::Polynomial_type_generator<double,1>::Type
      Poly_double_1;
      
      Poly_double_1 p1;
      CGAL::Polynomial_parser_d<Poly_double_1> goofy_1;

      // NOTE: escape symbols are deleted automatically
    std::string str="87623.28374872*x^10-x*x*8e+8\n\r\t\n+1e-11x+(1-x)^7";
    goofy_1(str, p1);

    std::cout << "p1: " << p1 << "\n\n";
        
  }
  {
      // trivariate parser with default policy
    typedef typename CGAL::Polynomial_type_generator<CGAL::Gmpz,3>::Type
      Polynomial_3;

    Polynomial_3 p3;
    CGAL::Polynomial_parser_d<Polynomial_3> goofy_3;

    // NOTE: variables can be named both in lower and upper cases
    std::string str="(234524523(X+y+z)^10-yYzZxX+13*x)^2-111";
    goofy_3(str, p3);

    std::cout << "p3: " << p3 << "\n\n";
      
  }

  return 0;

}

