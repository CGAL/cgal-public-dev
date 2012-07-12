// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Pavel Emeliyanenko    <asm@mpi-inf.mpg.de>
//
// ============================================================================

#include <CGAL/basic.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial/Polynomial_parser_d.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/Fraction_traits.h>

#include <CGAL/GMP_arithmetic_kernel.h>

// rename vars policy to parse 5-variate polynomials (OMG!)
template < class Poly_d_ >
struct Rename_vars_5 : public CGAL::Default_parser_policy < Poly_d_ > {

    //! template argument type
    typedef Poly_d_ Poly_d;
    //! base class
    typedef CGAL::Default_parser_policy < Poly_d_ > Base;
    //! type of polynomial coefficient
    typedef typename Base::Coeff Coeff;

    static const int n_var_names = 5;
    static const char *var_names_lower;
    static const char *var_names_upper;
};

template < class Poly_d_ >
const char * Rename_vars_5< Poly_d_ >::
        var_names_lower = "abcde";

template < class Poly_d_ >
const char * Rename_vars_5< Poly_d_ >::
        var_names_upper = "ABCDE";


template<typename ArithmeticKernel>
void test_routine_simple() {

  typedef ArithmeticKernel Arithmetic_kernel;
  typedef typename Arithmetic_kernel::Integer Integer;
  typedef typename Arithmetic_kernel::Rational Rational;

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
    typedef typename CGAL::Polynomial_type_generator<Integer,3>::Type
      Polynomial_3;

    Polynomial_3 p3;
    CGAL::Polynomial_parser_d<Polynomial_3> goofy_3;

    // NOTE: variables can be named both in lower and upper cases
    std::string str="(234524523(X+y+z)^10-yYzZxX+13*x)^2-111";
    goofy_3(str, p3);

    std::cout << "p3: " << p3 << "\n\n";
      
  }
  {
    typedef typename CGAL::Polynomial_type_generator<Integer,5>::Type
      Polynomial_5;

      Polynomial_5 p5;
     CGAL::Polynomial_parser_d< Polynomial_5, Rename_vars_5< Polynomial_5 > >
        goofy_5;

    // NOTE: abc is equivalent to a*b*c, i.e. '*' can be omitted
    std::string str="abcededeBABCBDEabcbdbdebdea+11";
    goofy_5(str, p5);

    std::cout << "p5: " << p5 << "\n\n";

    // NOTE: '=' is treated as LHS - RHS
    str="33(a+b-c+d-E)^4=ab2382734ee-897";
    goofy_5(str, p5);

    std::cout << "p5: " << p5 << "\n\n";
    

  }
}

int main(int argc, char** argv) {

  CGAL::set_pretty_mode(std::cout);

#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
    std::cout << "test for GMP" << std::endl;
    test_routine_simple<CGAL::GMP_arithmetic_kernel>();
#else
    std::cerr << "GMP tests skipped" << std::endl;
#endif
    return 0;
}

