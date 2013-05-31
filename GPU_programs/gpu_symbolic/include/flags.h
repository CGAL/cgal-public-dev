// ============================================================================
//
// Copyright (c) 2001-2008 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// this file is not part of any library ;-)
//
// ----------------------------------------------------------------------------
//
// Library       : CUDA MP
//
// File          : 
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

//!NOTE: this file is supposed to overwrite the respective files
//! from all libraries included. It contains macro definitions to control the
//! compilation of algorithms in a single library

#ifndef _FLAGS_H_
#define _FLAGS_H_
 
#warning unified library flags!!

#define CUMP_PREFETCH_FROM_CUDA_ARRAY 0 // use 2D texturing
#define CUMP_USE_PAGELOCKED_MEM       1 // allocate page-locked mem
#define CUMP_USE_PTX_ASSEMBLY         1 // use PTX assembly instead of
                                        // intrinsics
#define CUMP_USE_ATOMICS              1 // use atomic intructions
#define CUMP_VERBOSE                  0 // verbose output

#ifndef CUMP_USE_32BIT_MODULI_SET
#warning default moduli set selector
#define CUMP_USE_32BIT_MODULI_SET     1 // use 24 or 32-bit moduli set
#endif

#define CUMP_BENCHMARK_CALLS          0 // benchmark separate kernels
#define CUMP_MEASURE_MEM              0 // timing for GPU-host memory transfer

#define CUMP_DEVICE_MODULI_SZ     16382 // maximal # of moduli supported
#define CUMP_MODULI_STRIDE            4 // # of words per modulu

#define CUMP_COMPILE_GCD_KERNEL    1

#define CUMP_USE_PGCD_KERNEL       1 // use new pgcd kernel
#define CUMP_USE_GCD_BLOCK_KERNEL  1 // compile gcd block or monolithic kernel

#define CUMP_USE_STREAMS     1  // whether to use several streams

#define CUMP_RUN_GCD_ONLY     0 // whether to run only gcd kernel (without MRC)
#define CUMP_MRC_HOST_RECOVER 1 // whether to run MRC reconstrunction (host)

#define CUMP_GCD_USE_COPRIME_CHECK 1 // whether to use fast coprimality check
#define CUMP_GCD_TERMINATE_FLAG -1u // terminate flag
#define CUMP_GCD_RAW_OFFSET  16 // offset used by the kernel to return
                                // the size of polynomial

#define CUMP_RUN_RESULTANTS_ONLY      0 // compute only resultants
#define CUMP_RUN_RES_AND_INTERPOLATE  0 // resultants + interpolation
#define CUMP_USE_CRA_KERNEL           1 // apply the CRA (debug only)

#if (CUMP_RUN_RESULTANTS_ONLY + CUMP_RUN_RES_AND_INTERPOLATE + CUMP_USE_CRA_KERNEL > 1)
#error Flags must be mutually exclusive
#endif

#define CUMP_COMPILE_RESULTANTS_KERNEL 1

#define CUMP_COMPILE_DEBUG_KERNEL  0 // compile with debug support

#endif // _FLAGS_H_
