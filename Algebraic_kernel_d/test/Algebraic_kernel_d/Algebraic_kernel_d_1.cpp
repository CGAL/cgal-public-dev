// Copyright (c) 2006-2012 Max-Planck-Institute Saarbruecken (Germany).
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
// $URL: svn+ssh://eric@scm.gforge.inria.fr/svn/cgal/branches/unsorted-branches/eric/Numerical_algebraic_kernel_d/include/CGAL/Algebraic_kernel_d/Algebraic_curve_kernel_2.h $
// $Id: Algebraic_curve_kernel_2.h 70264 2012-07-04 13:01:47Z eric $
// 
//
// Author(s)     : Sebastian Limbach <slimbach@mpi-inf.mpg.de>
//                 Michael Hemmer    <hemmer@mpi-inf.mpg.de>
//
// ============================================================================

// Test of Algebraic_kernel

#define CGAL_TEST_ALL_AK_VARIANTS 1

#include <CGAL/config.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_rep_bfi.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_rep.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_quadratic_refinement_rep_bfi.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_coefficient_kernel.h>
#include <CGAL/Algebraic_kernel_d/Descartes.h>

#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Algebraic_kernel_2/Rounding_ak_d_1.h>

#include <CGAL/_test_algebraic_kernel_1.h>



template< class Coefficient_, class Bound_, class RepClass >
void test_algebraic_kernel_coeff_bound_rep() {
  typedef Coefficient_  Coefficient;
  typedef Bound_        Bound;
  typedef RepClass      Rep_class;

  typedef typename CGAL::Polynomial_type_generator<Coefficient,1>::Type
    Polynomial_1;
  typedef CGAL::internal::Algebraic_real_d_1
    < Coefficient, Bound, CGAL::Handle_policy_no_union, Rep_class >   Algebraic_real_1;

  typedef CGAL::internal::Descartes< Polynomial_1, Bound >               Descartes;
  typedef CGAL::internal::Bitstream_descartes<
    CGAL::internal::Bitstream_coefficient_kernel<Coefficient> >   BDescartes;

  typedef CGAL::Algebraic_kernel_d_1< Coefficient, Bound, Rep_class , Descartes>
    Kernel_Descartes;
  typedef CGAL::Algebraic_kernel_d_1< Coefficient, Bound, Rep_class , BDescartes>
    Kernel_BDescartes;

  typedef CGAL::internal::Rounding_ak_d_1< CGAL::Algebraic_kernel_d_1< Coefficient, Bound, Rep_class , Descartes> >
    Rounding_Kernel_Descartes;
  typedef CGAL::internal::Rounding_ak_d_1< CGAL::Algebraic_kernel_d_1< Coefficient, Bound, Rep_class , BDescartes> >
    Rounding_Kernel_BDescartes;

  std::cout << "Testing Kernel_Descartes" << std::endl;
  CGAL::test_algebraic_kernel_1(Kernel_Descartes());
  std::cout << "Testing Rounding_Kernel_Descartes" << std::endl;
  CGAL::test_algebraic_kernel_1(Rounding_Kernel_Descartes());
#if CGAL_TEST_ALL_AK_VARIANTS
  std::cout << "Testing Kernel_BDescartes" << std::endl;
  CGAL::test_algebraic_kernel_1(Kernel_BDescartes());
  std::cout << "Testing Rounding_Kernel_BDescartes" << std::endl;
  CGAL::test_algebraic_kernel_1(Rounding_Kernel_BDescartes());
#endif
}

template< class Coeff, class Bound >
void test_algebraic_kernel_coeff_bound() {
  test_algebraic_kernel_coeff_bound_rep<Coeff,Bound,
    CGAL::internal::Algebraic_real_rep< Coeff, Bound > > ();
#if CGAL_TEST_ALL_AK_VARIANTS
  test_algebraic_kernel_coeff_bound_rep<Coeff,Bound,
    CGAL::internal::Algebraic_real_rep_bfi< Coeff, Bound > > ();
  test_algebraic_kernel_coeff_bound_rep<Coeff,Bound,
    CGAL::internal::Algebraic_real_quadratic_refinement_rep_bfi< Coeff, Bound > > ();
#endif
}


template< class ArithmeticKernel >
void test_algebraic_kernel() {
  typedef ArithmeticKernel AK;
  typedef typename AK::Integer Integer;
  typedef typename AK::Rational Rational;

  test_algebraic_kernel_coeff_bound<Integer, Rational>();
#if CGAL_TEST_ALL_AK_VARIANTS
  test_algebraic_kernel_coeff_bound<Rational, Rational>();
  test_algebraic_kernel_coeff_bound
    <CGAL::Sqrt_extension< Integer, Integer>, Rational>();
  test_algebraic_kernel_coeff_bound
    <CGAL::Sqrt_extension< Rational, Integer>, Rational>();
  test_algebraic_kernel_coeff_bound
    <CGAL::Sqrt_extension< Rational, Rational>, Rational>();
#endif
}

int main() {
#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL
  std::cout << " TEST AK1 USING CORE " << std::endl;
  test_algebraic_kernel< CGAL::CORE_arithmetic_kernel >();
#else
  std::cout << " SKIPPED CORE " << std::endl;
#endif
#ifdef CGAL_HAS_LEDA_ARITHMETIC_KERNEL
  std::cout << " TEST AK1 USING LEDA " << std::endl;
  test_algebraic_kernel< CGAL::LEDA_arithmetic_kernel >();
#else
  std::cout << " SKIPPED LEDA " << std::endl;
#endif
#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
  std::cout << " TEST AK1 USING GMP " << std::endl;
  test_algebraic_kernel< CGAL::GMP_arithmetic_kernel >();
#else
  std::cout << " SKIPPED GMP " << std::endl;
#endif


  return 0;
}
