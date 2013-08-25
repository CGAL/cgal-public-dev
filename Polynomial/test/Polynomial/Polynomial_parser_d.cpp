// Copyright (c) 2010, 2012 Max-Planck-Institut fuer Informatik (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Michael Kerber    <mkerber@mpi-inf.mpg.de> 

#include <CGAL/config.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/Fraction_traits.h>


#include <CGAL/Polynomial/Polynomial_parser_d.h>

template < class Poly_d_ >
struct My_policy :
        public CGAL::Mixed_rational_parser_policy< Poly_d_ > {

    //! template argument type
    typedef Poly_d_ Poly_d;
    //! base class
    typedef CGAL::Mixed_rational_parser_policy< Poly_d > Base;
    //! type of polynomial coefficient
    typedef typename Base::Coeff Coeff;

    typedef typename CGAL::Get_arithmetic_kernel< Coeff >::
             Arithmetic_kernel AK;
    //! integer number type
    typedef typename AK::Integer Integer;
    //! rational number type
    typedef typename AK::Rational Rational;
    //! input coefficient types
    typedef typename Base::CoeffType CoeffType;

    virtual Coeff read_coeff_proxy(std::istream& is,
          CoeffType type) const {

        if(type == Base::COEFF_RATIONAL) {
            Integer num, denom;

            is >> CGAL::iformat(num); // read numerator
            is.get(); // skip '/'
            is >> CGAL::iformat(denom);

            std::cout << "rational: " << num << "/" << denom << "\n";
            throw CGAL::internal::Parser_exception("error!");
            return Coeff(0);
        } else
            return Base::read_coeff_proxy(is, type);
    }

    //! checking for degree overflow: can be used in real-time applications
    virtual bool exponent_check(unsigned e) const {
        std::cout << "exponent_check: " << e << "\n";
        return true;
    }

protected:


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
        My_policy< Polynomial_2 > > Custom_parser;

    pol_str="x^5*8/000*y^17-5/0*(-111y-x^2)+10/2134234";

    if(!Custom_parser() (pol_str,p2)) {
        printf("smth bogus happend..\n");
    }
    std::cout << "\ndone\n";
    
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

