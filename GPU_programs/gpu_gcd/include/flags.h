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

#ifndef _FLAGS_H_
#define _FLAGS_H_

#warning Local flags!!

#define CUMP_PREFETCH_FROM_CUDA_ARRAY 0 // use 2D texturing
#define CUMP_USE_PAGELOCKED_MEM       1 // allocate page-locked mem
#define CUMP_USE_PTX_ASSEMBLY         1 // use PTX assembly instead of
                                        // intrinsics
#define CUMP_USE_ATOMICS              1 // use atomic intructions
#define CUMP_VERBOSE                  1 // verbose output
 
#define CUMP_USE_32BIT_MODULI_SET     1 // use 24 or 32-bit moduli set

#define CUMP_BENCHMARK_CALLS          1 // benchmark separate kernels
#define CUMP_MEASURE_MEM              1 // timing for GPU-host memory transfer

#define CUMP_DEVICE_MODULI_SZ     16384 // maximal # of moduli supported
#define CUMP_MODULI_STRIDE            4 // # of words per modulu

#define CUMP_COMPILE_GCD_KERNEL       1

//! NOTE NOTE: the flag \c CUMP_USE_PGCD_KERNEL overwrites
//! \c CUMP_USE_GCD_BLOCK_KERNEL be careful !!
//! differences are: PGCD kernel saves the results back not in reversed order
//! in addition it requires input data padding (nu/nv_ofs4) 
#define CUMP_USE_PGCD_KERNEL       1 // use new pgcd kernel
#define CUMP_USE_GCD_BLOCK_KERNEL  1 // compile gcd block or monolithic kernel
#define CUMP_COMPILE_DEBUG_KERNEL  1 // compile with debug support

#define CUMP_USE_STREAMS     1  // whether to use several streams

#define CUMP_RUN_GCD_ONLY     0 // whether to run only gcd kernel (without MRC)
#define CUMP_MRC_HOST_RECOVER 1 // whether to run MRC reconstrunction (host)

#define CUMP_GCD_USE_COPRIME_CHECK 0 // whether to use fast coprimality check
#define CUMP_GCD_TERMINATE_FLAG -1u // terminate flag
#define CUMP_GCD_RAW_OFFSET  16 // offset used by the kernel to return 
                                // the size of polynomial

// NOTE: be careful when switching to full iterations !!
// NOTE also that the # of lite iterations depends on the selected chunk size 
#define DEBUG_ITER_IDX  26   // index of the last iteration to run

#define DEBUG_BLOCK_IDX  0   // # of iterations to run

#endif // _FLAGS_H_
