// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     :  
//
// ============================================================================


/*! \file Extended_descartes.C
 This is the test file for the class Extended_descartes.h
*/

#include <CGAL/basic.h>

#ifndef CGAL_ISOLATOR_USES_INPUT_RANGE
#define CGAL_ISOLATOR_USES_INPUT_RANGE 1
#endif

#ifndef CGAL_DESCARTES_VERBOSE
#define CGAL_DESCARTES_VERBOSE 1
#endif

#ifndef CGAL_DESCARTES_EXTENDED_VERBOSE
#define CGAL_DESCARTES_EXTENDED_VERBOSE 1
#endif

// include these traits here by 'hand', since not in release 3.3
#include <CGAL/Algebraic_extension_traits.h>
#include <CGAL/Scalar_factor_traits.h>

#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/Polynomial_traits_d.h>
#include <_test_real_root_isolator.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d/Real_bitstream_descartes.h>


template < class ArithmeticKernel >
void test_real_bitstream_descartes(){
  typedef ArithmeticKernel Arithmetic_kernel;
  typedef typename Arithmetic_kernel::Integer Integer;
  typedef typename Arithmetic_kernel::Rational Rational;
  {
    // general test of concept RealRootIsolator for Integer coefficients
    typedef ::CGAL::internal::Real_bitstream_descartes< Integer, Rational > 
      Isolator;
    CGAL::internal::test_real_root_isolator<Isolator>();
  }
  {
    // general test of concept RealRootIsolator for Rational coefficients
    typedef ::CGAL::internal::Real_bitstream_descartes< Rational, Rational > 
      Isolator;
    CGAL::internal::test_real_root_isolator<Isolator>();
  }    
}

template < class ArithmeticKernel >
void simple_test_real_bitstream_descartes(){
  
  typedef ArithmeticKernel Arithmetic_kernel;
  typedef typename Arithmetic_kernel::Rational Rational;
  typedef typename Arithmetic_kernel::Integer Integer;

  std::cout << "Rational polynomial1	:" << std::endl;

  std::list<Rational> coeffs;
  coeffs.push_back(2);
  coeffs.push_back(-3);
  coeffs.push_back(1);
  
  typedef ::CGAL::internal::
          Real_bitstream_descartes< Rational, Rational > Isolator;
  
  Isolator isolator(coeffs.begin(), coeffs.end());
  
  std::cout << "#Roots: " <<  isolator.number_of_real_roots() << std::endl;
  std::cout << "1st root: [" << isolator.left_bound(0) << "," 
            << isolator.right_bound(0) << "]" << std::endl;
  std::cout << "2nd root: [" << isolator.left_bound(1) << "," 
            << isolator.right_bound(1) << "]" << std::endl;


std::cout << "Mignotte-like polynomial:" << std::endl;

  std::list<Rational> coeff_mig;
  coeff_mig.push_back(-2);
  coeff_mig.push_back(800);
  coeff_mig.push_back(-80000);
  for (int i = 0; i < 37; i++) {
    coeff_mig.push_back(0);
  }
  coeff_mig.push_back(1048576);
  
  typedef ::CGAL::internal::Real_bitstream_descartes<Rational,Rational> 
                                                                      Isolator;
  
  Isolator isolator_mignotte(coeff_mig.begin(), coeff_mig.end());
  std::cout << "#Roots: " <<  isolator_mignotte.number_of_real_roots() <<
                std::endl;
  std::cout << "1st root: [" << isolator_mignotte.left_bound(0) << "," 
            << isolator_mignotte.right_bound(0) << "]" << std::endl;
  std::cout << "2nd root: [" << isolator_mignotte.left_bound(1) << "," 
            << isolator_mignotte.right_bound(1) << "]" << std::endl;
  std::cout << "3nd root: [" << isolator_mignotte.left_bound(2) << "," 
            << isolator_mignotte.right_bound(2) << "]" << std::endl;
  std::cout << "4nd root: [" << isolator_mignotte.left_bound(3) << "," 
            << isolator_mignotte.right_bound(3) << "]" << std::endl;
  


  std::cout << "Algebraic polynomial:" << std::endl;

  typedef typename CGAL::Algebraic_kernel_d_1_generator< Integer, Rational >::
    Algebraic_kernel_with_qir_and_descartes_1  Algebraic_kernel_1;
  typedef typename Algebraic_kernel_1::Algebraic_real_1 Algebraic_real_1;
  typedef std::list< std::pair< Algebraic_real_1, unsigned int > > Roots;

  std::list<Algebraic_real_1> coeff_real;
  Roots roots;
  typename Algebraic_kernel_1::Solve_1 solve;

  typedef typename Algebraic_kernel_1::Polynomial_1 Polynomial_1;

  typedef CGAL::Polynomial_traits_d<Polynomial_1> PT_1;
  Polynomial_1 x = typename PT_1::Shift()(Polynomial_1(1),1);
  Polynomial_1 polynomial = x*x - 2;
  solve(polynomial, std::back_inserter(roots));

  Algebraic_real_1 a = roots.back().first;    
  coeff_real.push_back(a);
  coeff_real.push_back(Algebraic_real_1(1));
  coeff_real.push_back(Algebraic_real_1(-2));

  typedef ::CGAL::internal::
        Real_bitstream_descartes< Algebraic_real_1, Rational > Isolator_real;
  Isolator_real isolator_real(coeff_real.begin(), coeff_real.end());
  std::cout << "#Roots: " <<  isolator_real.number_of_real_roots() << std::endl;
  std::cout << "1st root: [" << isolator_real.left_bound(0) << "," 
            << isolator_real.right_bound(0) << "]" << std::endl;
  std::cout << "2nd root: [" << isolator_real.left_bound(1) << "," 
            << isolator_real.right_bound(1) << "]" << std::endl;
  
}


int main(){
  
  CGAL::set_pretty_mode(std::cout);
  CGAL::set_pretty_mode(std::cerr);
  CGAL::set_pretty_mode(std::clog);

#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
  simple_test_real_bitstream_descartes< CGAL::GMP_arithmetic_kernel >();
//  test_real_bitstream_descartes<CGAL::GMP_arithmetic_kernel>();
#endif
#ifdef CGAL_HAS_LEDA_ARITHMETIC_KERNEL
//  simple_test_real_bitstream_descartes< CGAL::LEDA_arithmetic_kernel >();
//  test_real_bitstream_descartes<CGAL::LEDA_arithmetic_kernel>();
#endif
#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL
//  simple_test_real_bitstream_descartes< CGAL::CORE_arithmetic_kernel >();
//  test_real_bitstream_descartes<CGAL::CORE_arithmetic_kernel>();
#endif

  return EXIT_SUCCESS;
}
// EOF
