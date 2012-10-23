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

#warning Using local res flags !!

#define CUMP_PREFETCH_FROM_CUDA_ARRAY 0 // use 2D texturing
#define CUMP_USE_PAGELOCKED_MEM       1 // allocate page-locked mem
#define CUMP_USE_PTX_ASSEMBLY         1 // use PTX assembly instead of
                                        // intrinsics
#define CUMP_USE_ATOMICS              1 // use atomic intructions
#define CUMP_USE_32BIT_MODULI_SET     1 // use 24 or 32-bit moduli set

#define CUMP_VERBOSE  1

#define CUMP_RUN_RESULTANTS_ONLY      0 // compute only resultants
#define CUMP_RUN_RES_AND_INTERPOLATE  0 // resultants + interpolation
#define CUMP_USE_CRA_KERNEL      1 // apply the CRA (debug only)
#define CUMP_MRC_HOST_RECOVER    1 // recover integers from MR digits

//! NOTE remark that in \c CUMP_RUN_RES_AND_INTERPOLATE mode, GPU-host results
//! might differ because results are not multiplied by modular inverse
//! (which happens in the CRA kernel)

#if (CUMP_RUN_RESULTANTS_ONLY + CUMP_RUN_RES_AND_INTERPOLATE + CUMP_USE_CRA_KERNEL > 1)
#error Flags are mutually exclusive
#endif

#define CUMP_BENCHMARK_CALLS 1  // benchmark separate kernels
#define CUMP_MEASURE_MEM     1  // measure time for GPU-host memory transfer

#define CUMP_DEVICE_MODULI_SZ 16382
#define CUMP_MODULI_STRIDE    4     // 4 elements per modulus

#define CUMP_COMPILE_RESULTANTS_KERNEL 1
#define CUMP_COMPILE_DEBUG_KERNEL      1

#endif // _FLAGS_H_
