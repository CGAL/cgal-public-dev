// Copyright (c) 2009, 2010, 2011, 2012 Max-Planck-Institut Saarbruecken (Germany).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Eric Berberich <eric.berberich@cgal.org> 
//                 Pavel Emeliyanenko <asm@mpi-sb.mpg.de>

#ifndef CGAL_ALGEBRAIC_KERNEL_2_BI_SOLVE_2_FLAGS_H
#define CGAL_ALGEBRAIC_KERNEL_2_BI_SOLVE_2_FLAGS_H 1

#include <CGAL/config.h>
#include <CGAL/Algebraic_kernel_d/flags.h>

// bisolve

#ifndef CGAL_BISOLVE_DISABLE_LOCALIZED_OPERATOR
#define CGAL_BISOLVE_DISABLE_LOCALIZED_OPERATOR 0 // 0 is default
#endif
#if CGAL_AK_D_SHOW_COMPILE_OPTIONS_AS_WARNING && CGAL_BISOLVE_DISABLE_LOCALIZED_OPERATOR
#warning Bi_solve_2: Disabled localized operator for x
#endif

#ifndef CGAL_BISOLVE_ENABLE_LENSE_FILTER_IN_T_TEST
#define CGAL_BISOLVE_ENABLE_LENSE_FILTER_IN_T_TEST 0 // 0 is default
#endif
#if CGAL_AK_D_SHOW_COMPILE_OPTIONS_AS_WARNING && CGAL_BISOLVE_ENABLE_LENSE_FILTER_IN_T_TEST
#warning Bi_solve_2: Lense filter in well-seperation is enabled
#endif
 
// certifier

#ifndef CGAL_BISOLVE_DISABLE_COMBINATORIAL_CERTIFICATION
#define CGAL_BISOLVE_DISABLE_COMBINATORIAL_CERTIFICATION 0 // 0 is default
#endif
#if CGAL_AK_D_SHOW_COMPILE_OPTIONS_AS_WARNING && CGAL_BISOLVE_DISABLE_COMBINATORIAL_CERTIFICATION
#warning Bi_solve_2: Combinatorial certification is disabled
#endif

#ifndef CGAL_BISOLVE_DISABLE_BIDIRECTIONAL_CERTIFICATION
#define CGAL_BISOLVE_DISABLE_BIDIRECTIONAL_CERTIFICATION 0 // 0 is default
#endif
#if CGAL_AK_D_SHOW_COMPILE_OPTIONS_AS_WARNING && CGAL_BISOLVE_DISABLE_BIDIRECTIONAL_CERTIFICATION
#warning Bi_solve_2: Certification in both directions is disabled
#endif

#ifndef CGAL_BISOLVE_ENABLE_ARCAVOID 
#define CGAL_BISOLVE_ENABLE_ARCAVOID 1 // default TODO?
#endif
#if CGAL_BISOLVE_ENABLE_ARCAVOID 
#undef CGAL_BISOLVE_DISABLE_BITSTREAM_FILTER
#define CGAL_BISOLVE_DISABLE_BITSTREAM_FILTER 1 // DISABLE BITSTREAM -- DO NOT CHANGE
#if CGAL_AK_D_SHOW_COMPILE_OPTIONS_AS_WARNING
#warning Bi_solve_2: Enabled Arcavoid numerical solver
#endif
#endif

#ifndef CGAL_BISOLVE_ENABLE_NTL_FACTORIZE
#define CGAL_BISOLVE_ENABLE_NTL_FACTORIZE 0 // 0 is default
#endif
#if CGAL_AK_D_SHOW_COMPILE_OPTIONS_AS_WARNING && CGAL_BISOLVE_ENABLE_NTL_FACTORIZE
#warning Bi_solve_2: factorization by NTL is enabled
#endif


#ifndef CGAL_BISOLVE_DISABLE_BITSTREAM_FILTER
#define CGAL_BISOLVE_DISABLE_BITSTREAM_FILTER 0 // 0 is default
#endif
#if CGAL_AK_D_SHOW_COMPILE_OPTIONS_AS_WARNING && CGAL_BISOLVE_DISABLE_BITSTREAM_FILTER && !CGAL_BISOLVE_ENABLE_ARCAVOID
#warning Bi_solve_2: Bitstreamfilter is disabled
#endif

#ifndef CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS
#define CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS (CGAL_BISOLVE_DISABLE_BITSTREAM_FILTER && !CGAL_BISOLVE_ENABLE_ARCAVOID)
#endif
#if CGAL_AK_D_SHOW_COMPILE_OPTIONS_AS_WARNING && CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS
#warning Bi_solve_2: Active intervals disabled
#endif

#ifndef CGAL_ACK_BITSTREAM_USES_E08_TREE
#define CGAL_ACK_BITSTREAM_USES_E08_TREE 1 // do not change
#endif

// set precision according to the results of the norm test
#ifndef CGAL_BISOLVE_USE_ADJUSTABLE_PRECISION
#define CGAL_BISOLVE_USE_ADJUSTABLE_PRECISION 1 // 1 is default
#endif

// verbosity

#ifndef CGAL_BISOLVE_DEBUG
#define CGAL_BISOLVE_DEBUG 0 // 0 is default
#endif

#ifndef CGAL_BISOLVE_VERBOSE
#define CGAL_BISOLVE_VERBOSE 0 // 0 is default
#endif

#if CGAL_BISOLVE_VERBOSE
#define Bisolve_out(x) std::cout << x;
#define dbl(x) CGAL::to_double(x)
#define bfi(x) CGAL::lower(CGAL::convert_to_bfi(x))
#define STILL_ALIVE std::clog << __LINE__ << "\n";
#else
#define Bisolve_out(x) static_cast< void >(0);
#define STILL_ALIVE static_cast< void >(0);
#endif

// telemetry

#ifndef CGAL_BISOLVE_TELEMETRY
#define CGAL_BISOLVE_TELEMETRY 0 // 0 is default
#endif

#if CGAL_BISOLVE_TELEMETRY
#define Bisolve_telemetry_code(x) x;
#else
#define Bisolve_telemetry_code(x) 
#endif


#endif // CGAL_ALGEBRAIC_KERNEL_2_BI_SOLVE_2_FLAGS_H
// EOF
