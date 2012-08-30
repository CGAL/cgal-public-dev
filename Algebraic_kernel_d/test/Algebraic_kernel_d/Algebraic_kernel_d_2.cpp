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
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

// code coverage test for Algebraic_curve_kernel_2

// #define CGAL_ACK_DEBUG_FLAG 1

//#define CGAL_AK_ENABLE_DEPRECATED_INTERFACE 1

#ifndef CGAL_ACK_CURVE_ANALYSES_USE_BISOLVE 
#define CGAL_ACK_CURVE_ANALYSES_USE_BISOLVE 1
#endif
#ifndef CGAL_ACK_CURVE_ANALYSES_BISOLVE_USE_TEISSIER
#define CGAL_ACK_CURVE_ANALYSES_BISOLVE_USE_TEISSIER 1
#endif

#include <CGAL/Algebraic_kernel_d/flags.h>
#include <CGAL/basic.h>

#include <CGAL/Arithmetic_kernel.h>

#include <CGAL/Algebraic_kernel_d_2.h>

#include <CGAL/Sqrt_extension.h>

#include <CGAL/_test_algebraic_kernel_2.h>


template< typename Coefficient >
void test_algebraic_kernel() {
    typedef CGAL::Algebraic_kernel_d_2<Coefficient> Algebraic_kernel_d_2;
    Algebraic_kernel_d_2 ak_2;
    CGAL::test_algebraic_kernel_2<Algebraic_kernel_d_2>(ak_2);
}

int main() {


#ifdef CGAL_HAS_LEDA_ARITHMETIC_KERNEL

#if CGAL_ACK_DEBUG_FLAG
    CGAL_ACK_DEBUG_PRINT << "TESTING LEDA" << std::endl;
#endif    
    {
      
      typedef CGAL::LEDA_arithmetic_kernel AK;
      test_algebraic_kernel<AK::Integer>();
      /*
      test_algebraic_kernel<AK::Rational>();
      test_algebraic_kernel<CGAL::Sqrt_extension<AK::Integer,AK::Integer> >();
      
      test_algebraic_kernel
	<CGAL::Sqrt_extension<AK::Rational,AK::Rational> >();
      */
    }
#else
    std::cerr << "LEDA tests skipped" << std::endl;
#endif


#if defined(CGAL_HAS_CORE_ARITHMETIC_KERNEL) && !CGAL_ACK_CURVE_ANALYSES_USE_BISOLVE
#if CGAL_ACK_DEBUG_FLAG
    CGAL_ACK_DEBUG_PRINT << "TESTING CORE" << std::endl;
#endif          
    {
      typedef CGAL::CORE_arithmetic_kernel AK;
      //      test_algebraic_kernel<AK::Integer>();
      /*
      test_algebraic_kernel<AK::Rational>();
      test_algebraic_kernel<CGAL::Sqrt_extension<AK::Integer,AK::Integer> >();
      
      test_algebraic_kernel
	<CGAL::Sqrt_extension<AK::Rational,AK::Rational> >();
      */
    }
#else
    std::cerr << "CORE tests skipped" << std::endl;
#endif


#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
#if CGAL_ACK_DEBUG_FLAG
    CGAL_ACK_DEBUG_PRINT << "TESTING GMP" << std::endl;
#endif       
    {
      typedef CGAL::GMP_arithmetic_kernel AK;
      test_algebraic_kernel<AK::Integer>();
      /*
      test_algebraic_kernel<AK::Rational>();
      test_algebraic_kernel<CGAL::Sqrt_extension<AK::Integer,AK::Integer> >();
      
      test_algebraic_kernel
	<CGAL::Sqrt_extension<AK::Rational,AK::Rational> >();
      */
    }
#else
    std::cerr << "GMP tests skipped" << std::endl;
#endif


    return 0;
}
