// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Michael Kerber    <mkerber@mpi-inf.mpg.de> 
//
// ============================================================================

#include <CGAL/basic.h>
#include <CGAL/Algebraic_kernel_d/Polynomial_parser_d.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/Fraction_traits.h>

// namespace boost { // dirty hack...
// namespace detail {
// tss::~tss() {
//  
// }
// } }

template<typename ArithmeticKernel>
void test_routine() {

  typedef ArithmeticKernel Arithmetic_kernel;
  typedef typename Arithmetic_kernel::Integer Integer;
  typedef typename Arithmetic_kernel::Rational Rational;
  typename CGAL::Fraction_traits<Rational>::Compose compose;


  {
    typedef typename CGAL::Polynomial_type_generator<Integer,1>::Type 
      Polynomial_1;
    
    typedef CGAL::Polynomial_parser_d<Polynomial_1> Parser;
    
    Polynomial_1 t = CGAL::shift(Polynomial_1(1),1,0);
  
    std::string pol_str="x^17-(x*(x^3-x^2)+(x-1)*(x-2))";
    Polynomial_1 p1 = 
      CGAL::ipower(t,17)-(t*(CGAL::ipower(t,3)-CGAL::ipower(t,2))+(t-1)*(t-2));
    
    Polynomial_1 p2;
    Parser() (pol_str,p2);

    CGAL_assertion(p1==p2);
 

  }
  {
    typedef typename CGAL::Polynomial_type_generator<Integer,2>::Type 
      Polynomial_2;
    
    typedef CGAL::Polynomial_parser_d<Polynomial_2 > Parser;
    
    Polynomial_2 x = CGAL::shift(Polynomial_2(1),1,0);
    Polynomial_2 y = CGAL::shift(Polynomial_2(1),1,1);
  
    std::string pol_str="x^5*y^17-5*(y-x^2)+1";
    Polynomial_2 p1 = 
      CGAL::ipower(x,5)*CGAL::ipower(y,17)-5*(y-CGAL::ipower(x,2))+1;
    
    Polynomial_2 p2;
    Parser() (pol_str,p2);
    
    CGAL_assertion(p1==p2);
    
  }
  {
    typedef typename CGAL::Polynomial_type_generator<Integer,3>::Type 
      Polynomial_3;
    
    typedef CGAL::Polynomial_parser_d<Polynomial_3> Parser;
    
    Polynomial_3 x = CGAL::shift(Polynomial_3(1),1,0);
    Polynomial_3 y = CGAL::shift(Polynomial_3(1),1,1);
    Polynomial_3 z = CGAL::shift(Polynomial_3(1),1,2);
  
    std::string pol_str="x^5*y^10*(x-y-z)+z^0+1";
    Polynomial_3 p1 = 
      CGAL::ipower(x,5)*CGAL::ipower(y,10)*(x-y-z)+CGAL::ipower(z,0)+1;
    
    Polynomial_3 p2;
    Parser() (pol_str,p2);
    
    CGAL_assertion(p1==p2);
    
  }
  {
    typedef typename CGAL::Polynomial_type_generator<Rational,1>::Type 
      Polynomial_1;
    
    typedef CGAL::Polynomial_parser_d<Polynomial_1> Parser;
    
    Polynomial_1 t = CGAL::shift(Polynomial_1(1),1,0);
  
    std::string pol_str="1/7*x^17-(x*(x^3-x^2)+(x-1/2)*(x-2/5))";
    Polynomial_1 p1 = 
      compose(Integer(1),Integer(7))*				   \
      CGAL::ipower(t,17)-(t*(CGAL::ipower(t,3)-CGAL::ipower(t,2))+ \
      (t-compose(Integer(1),Integer(2)))* 	
      (t-compose(Integer(2),Integer(5))));
    
    Polynomial_1 p2;
    Parser() (pol_str,p2);

    CGAL_assertion(p1==p2);
 

  }
  
  {
    typedef typename CGAL::Polynomial_type_generator<Rational,2>::Type 
      Polynomial_2;
    
    typedef CGAL::Polynomial_parser_d<Polynomial_2,
            CGAL::Mixed_rational_parser_policy< Polynomial_2 > > Parser;
    
    Polynomial_2 x = CGAL::shift(Polynomial_2(1),1,0);
    Polynomial_2 y = CGAL::shift(Polynomial_2(1),1,1);
  
    std::string pol_str="x^5*8/3*y^17-5/7*(-111y-x^2)+10";
    Polynomial_2 p1 = 
    CGAL::ipower(x,5)*compose(Integer(8),Integer(3))*CGAL::ipower(y,17)-\
    compose(Integer(5),Integer(7))*(compose(Integer(-111),Integer(1))*y-CGAL::ipower(x,2))+\
    compose(Integer(10),Integer(1));
    
    Polynomial_2 p2;
    Parser() (pol_str,p2);
    
    CGAL_assertion(p1==p2);
    
  }
  {
    typedef typename CGAL::Polynomial_type_generator<Rational,3>::Type 
      Polynomial_3;
    
    typedef CGAL::Polynomial_parser_d<Polynomial_3> Parser;
    
    Polynomial_3 x = CGAL::shift(Polynomial_3(1),1,0);
    Polynomial_3 y = CGAL::shift(Polynomial_3(1),1,1);
    Polynomial_3 z = CGAL::shift(Polynomial_3(1),1,2);
  
    std::string pol_str="-5/4*x^5*y^10*(x-y-z)+z^0+1/1";
    Polynomial_3 p1 = 
      compose(Integer(-5),Integer(4))*CGAL::ipower(x,5)*CGAL::ipower(y,10)* \
	      (x-y-z)+CGAL::ipower(z,0)+compose(Integer(1),Integer(1));
    
    Polynomial_3 p2;
    Parser() (pol_str,p2);
    
    CGAL_assertion(p1==p2);
    
  }
  
}



int main(int argc, char** argv) {

  CGAL::set_pretty_mode(std::cout);

#ifdef CGAL_HAS_LEDA_ARITHMETIC_KERNEL
  std::cout << "test for LEDA" << std::endl;
  test_routine<CGAL::LEDA_arithmetic_kernel>();
#else
    std::cerr << "LEDA tests skipped" << std::endl;
#endif
#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL
    std::cout << "test for CORE" << std::endl;  
    test_routine<CGAL::CORE_arithmetic_kernel>();
#else
    std::cerr << "CORE tests skipped" << std::endl;
#endif
#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
    std::cout << "test for GMP" << std::endl;
    test_routine<CGAL::GMP_arithmetic_kernel>();
#else
    std::cerr << "GMP tests skipped" << std::endl;
#endif
    return 0;
}

