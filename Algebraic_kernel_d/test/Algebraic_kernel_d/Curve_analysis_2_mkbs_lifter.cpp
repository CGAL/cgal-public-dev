// Copyright (c) 2012 Max-Planck-Institute Saarbruecken (Germany).
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
// $URL: svn+ssh://eric@scm.gforge.inria.fr/svn/cgal/branches/unsorted-branches/eric/Numerical_algebraic_kernel_d/include/CGAL/Algebraic_kernel_d/Bitstream_descartes.h $
// $Id: Bitstream_descartes.h 70669 2012-07-23 08:45:50Z eric $
//
// Author(s)     : Eric Berberich <eric.berberich@cgal.org>
//                 
//
// ============================================================================

/*! \file Curve_analysis_2_mkbs_lifter
 This is the test file for the class CGAL::internal::Mkbs_lifter

*/

#include <CGAL/config.h>

#include <CGAL/Arithmetic_kernel.h>

#include <CGAL/_test_mkbs_lifter.h>
   
int main() {
#ifdef CGAL_HAS_LEDA_ARITHMETIC_KERNEL  
  CGAL::internal::test_mkbs_lifter<CGAL::LEDA_arithmetic_kernel>();
#endif

#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL
  CGAL::internal::test_mkbs_lifter<CGAL::CORE_arithmetic_kernel>();
#endif

#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
  CGAL::internal::test_mkbs_lifter<CGAL::GMP_arithmetic_kernel>();
#endif
    return EXIT_SUCCESS;
}
// EOF
