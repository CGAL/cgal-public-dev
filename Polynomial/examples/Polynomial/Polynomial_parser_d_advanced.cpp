#include <CGAL/config.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial/Polynomial_parser_d.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/Fraction_traits.h>
#include <CGAL/convert_to_bfi.h>

#include <CGAL/LEDA_arithmetic_kernel.h>
#include <CGAL/CORE_arithmetic_kernel.h>
#include <CGAL/GMP_arithmetic_kernel.h>

// maximal polynomial degree
#define PARSER_MAX_POLY_DEGREE 30
// bits to approximate bigfloat coefficients
#define PARSER_FLOAT_APPROX_BITS 53 

#define __STILL_ALIVE__ printf("line: %d\n", __LINE__);

//! example of custom parser policy to parse polynomials with
//! integer, rational and bigfloat coefficients (with supported rounding)
//! and degree check
template < class Poly_d_ >
struct Custom_parser_policy :
        public CGAL::Mixed_floating_point_parser_policy< Poly_d_ > {

    //! template argument type
    typedef Poly_d_ Poly_d;
    //! base class
    typedef CGAL::Mixed_floating_point_parser_policy< Poly_d > Base;
    //! type of polynomial coefficient
    typedef typename Base::Innermost_coefficient_type Innermost_coefficient_type;

    typedef typename CGAL::Get_arithmetic_kernel< Innermost_coefficient_type >::
             Arithmetic_kernel AK;
    //! integer number type
    typedef typename Base::Integer Integer;
    //! rational number type
    typedef typename Base::Rational Rational;
    //! BFI type
    typedef typename AK::Bigfloat_interval BFI;
    //! BigFloat type
    typedef typename CGAL::Bigfloat_interval_traits< BFI >::Bound BigFloat;
    //! input coefficient types
    typedef typename Base::CoefficientTypeID CoefficientTypeID;

    virtual Innermost_coefficient_type read_coeff_proxy(std::istream& is,
          CoefficientTypeID type) const {

        if(type == Base::COEFF_RATIONAL) {
            Integer num, denom;

            is >> CGAL::iformat(num); // read numerator
            is.get(); // skip '/'
            is >> CGAL::iformat(denom);
//             std::cout << "rational: " << num << "/" << denom << "\n";
            if(CGAL::is_zero(denom))
                throw CGAL::internal::Parser_exception("Zero divisor!");

            typedef CGAL::Fraction_traits< Rational > FT;
            return typename FT::Compose()(num, denom);

        } else if(type == Base::COEFF_FLOAT) {
            double ld;
            is >> CGAL::iformat(ld);
            BigFloat bf(ld);

            long prec = CGAL::set_precision(BFI(), PARSER_FLOAT_APPROX_BITS);
            BFI bfi = CGAL::convert_to_bfi(bf);
            CGAL::set_precision(BFI(), prec);

            return CGAL::lower(bfi);

        } else {
            
            return Base::read_coeff_proxy(is, type);
            
        }
    }

    //! checking for degree overflow: can be used in real-time applications
    virtual bool exponent_check(unsigned e) const {
        if(e > PARSER_MAX_POLY_DEGREE)
            return false;
        return true;
    }

protected:

};


// rename vars policy to parse 5-variate polynomials (OMG!)
template < class Poly_d_ >
struct Rename_vars_5 : public CGAL::Default_parser_policy < Poly_d_ > {

    //! template argument type
    typedef Poly_d_ Poly_d;
    //! base class
    typedef CGAL::Default_parser_policy < Poly_d_ > Base;
    //! type of polynomial coefficient
    typedef typename Base::Innermost_coefficient_type Innermost_coefficient_type;
    
    Rename_vars_5() :
        Base("abcde") {
    }
};

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

    typedef CGAL::Polynomial_parser_d<Polynomial_2,
        Custom_parser_policy< Polynomial_2 > > Custom_parser;
    // polynomial with invalid coefficients
    pol_str="123.3453x^5*8/000*y^17-5/0*(-111y-x^2)+10/2134234";

    if(!Custom_parser() (pol_str,p2)) {
        printf("don't worry: this was intended..\n");
    }
    // example polynomial with mixed integer/rational and fp-coefficients
    pol_str="234523.4(x-y)^5-y^17-5/32*(-111y-345.332342345x^2)+10/34234";

    if(!Custom_parser() (pol_str,p2)) {
        printf("wrong poly coefficients..\n");
    }
    std::cout << "constructed polynomial: " << p2 << "\n";
    
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
  {
    typedef typename CGAL::Polynomial_type_generator<CGAL::Gmpz,5>::Type
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

