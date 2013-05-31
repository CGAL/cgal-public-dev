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
// Author(s)     :  
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

/*! \file Descartes.cpp
 This is the test file for the class CGAL::Descartes.
*/

#define CGAL_ALGEBRAIC_KERNEL_D_ISOLATORS_DISABLE_ASSIGNABLE 0
#define CGAL_TEST_REAL_ROOT_INTERVAL_ISOLATOR 1

#include <CGAL/basic.h>

// include these traits here by 'hand', since not in release 3.3
#include <CGAL/Algebraic_extension_traits.h>
#include <CGAL/Scalar_factor_traits.h>

#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/_test_real_root_isolator.h>

#include <CGAL/Algebraic_kernel_d/Descartes.h>
#include <CGAL/Arithmetic_kernel.h>

template <class AT>
void test_descartes(){
    typedef typename AT::Integer Integer;
    typedef typename AT::Rational Rational;
    {
        typedef typename CGAL::Polynomial_type_generator<Integer,1>::Type 
            Polynomial;
        typedef ::CGAL::internal::Descartes<Polynomial,Rational> Isolator;
        
        // general test of concept RealRootIsolator
        CGAL::internal::test_real_root_isolator<Isolator>();
    }{
        typedef typename CGAL::Polynomial_type_generator<Rational,1>::Type 
            Polynomial;
        typedef ::CGAL::internal::Descartes<Polynomial,Rational> Isolator;
        // general test of concept RealRootIsolator
        CGAL::internal::test_real_root_isolator<Isolator>();
    }    
}
    
int main(){
#ifdef CGAL_HAS_LEDA_ARITHMETIC_KERNEL
  std::cout << " TEST AK1 USING LEDA " << std::endl;
  test_descartes< CGAL::LEDA_arithmetic_kernel >();
#endif
#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL
  std::cout << " TEST AK1 USING CORE " << std::endl;
  test_descartes< CGAL::CORE_arithmetic_kernel >();
#endif
#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
  std::cout << " TEST AK1 USING GMP " << std::endl;
  test_descartes< CGAL::GMP_arithmetic_kernel >();
#endif

    return EXIT_SUCCESS;
}
// EOF
