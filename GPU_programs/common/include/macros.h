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

#ifndef _MACROS_H_
#define _MACROS_H_

#include <include/flags.h>
#include <gmp.h>

typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int limb;

typedef unsigned long ulong;
typedef unsigned long long uint64;
typedef long long int64;

#if CUMP_USE_32BIT_MODULI_SET
typedef double fp_limb;
#else
typedef float fp_limb;
#endif

#if CUMP_VERBOSE
#define STILL_ALIVE printf("%d\n", __LINE__);
#define CUMP_out(x) std::cerr << x;
#define CUMP_out2(x...) printf(x);
#else
#define STILL_ALIVE static_cast< void >(0);
#define CUMP_out(x) static_cast< void >(0);
#define CUMP_out2(x, ...) static_cast< void >(0);
#endif

#define WS 32
#define HF 16
#define LOG_WARP_SIZE 5

#if !CUMP_USE_ATOMICS
#warning Compiling in 1.0 compatibility mode: no atomic instruction support
#endif

#define CU_SYNC __syncthreads();

#if CUDART_VERSION >= 2000
#define CU_MEMFENCE __threadfence();
#else
#define CU_MEMFENCE static_cast< void >(0);
#endif

// Fermi has native 32-bit integer mul 
#if defined __CUDA_ARCH__ && __CUDA_ARCH__ >= 200

#define UMUL_PAD     0U
#define UMUL(x, y) (x) * (y)

#define CUMP_LAUNCH_BOUNDS(t, b) __launch_bounds__(t, b)
#else
#define UMUL_PAD     0x1000000U // 2^24
#define UMUL(x, y) __umul24((x), (y)) // GT200 and lower use 24-bit mul
// #define UMUL_PAD     0
// #define UMUL(x, y) (x) * (y)
#define CUMP_LAUNCH_BOUNDS(t, b)
#endif

// computes r = a - b subop c unsigned using extended precision
#define VSUBx(r, a, b, c, subop) \
     asm volatile("vsub.u32.u32.u32." subop " %0, %1, %2, %3;" :  \
                "=r"(r) : "r"(a) , "r"(b), "r"(c));

// computes r = a + b subop c unsigned using extended precision
#define VADDx(r, a, b, c, subop) \
     asm volatile("vadd.u32.u32.u32." subop " %0, %1, %2, %3;" :  \
                "=r"(r) : "r"(a) , "r"(b), "r"(c));

#if CUMP_USE_PTX_ASSEMBLY

#define UMUL24HI(c, a, b) \
        asm volatile("mul24.hi.u32 %0, %1, %2;" : \
                "=r"(c) : "r"(a) , "r"(b));

// add & output carry flag
#define UADDO(c, a, b) \
    asm volatile("add.cc.u32 %0, %1, %2;" : \
                "=r"(c) : "r"(a) , "r"(b));

// add with carry & output carry flag
#define UADDC(c, a, b) \
    asm volatile("addc.cc.u32 %0, %1, %2;" : \
                "=r"(c) : "r"(a) , "r"(b));

#else // otherwise use compiler-defined intrinsic
#warning using compiler intrinsics

#define UMUL24HI(c, a, b) \
        c = __umul24hi(a, b);

    //((__umulhi(a, b) << 16) + (__umul24(a, b) >> 16))
    //(__umul24(a >> 16, b) + __umul24(a & 0xffff, b >> 16))

#define UADDC(s, c, a, b) s = __uaddc((c), 0);
#define USETC(s, a, b) s = __uaddc(0, 0);
#define UADDO(s, a, b) s = __uaddo((a), (b));

#endif // CUMP_USE_PTX_ASSEMBLY

#if CUMP_USE_32BIT_MODULI_SET
#warning using 32-bit moduli set
#endif
 
#if CUDART_VERSION < 2000   // CUDA 1.0 does not like exceptions
#define CUMP_THROW(x) exit(1);
#else
#define CUMP_THROW(x) throw x();
#endif

// #define ADDC(s, c, a, b) s = c + ((a) < (b));
// #define SETC(s, a, b) s = ((a) < (b));
// #define ADDO(s, a, b) s = (a) + (b);

// #define CUMP_ERROR_CHECK { \
//     cudaError err = cudaGetLastError();                                   \
//     if( cudaSuccess != err) {                                                \
//         fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
//                 __FILE__, __LINE__, cudaGetErrorString( err) );              \
//         exit(EXIT_FAILURE);                                                  \
//     } }

#define CUMP_SAFE_CALL( call) {                                    \
    cudaError err = call;                                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString( err) );              \
        exit(EXIT_FAILURE);                                                  \
    } }

#define CUMP_SAFE_CALL_drv( call) {                                    \
    CUresult err = call;                                                    \
    if( CUDA_SUCCESS != err) {                                                \
        fprintf(stderr, "Cuda drv error in file '%s' in line %i: %d.\n",    \
                __FILE__, __LINE__, err);              \
        exit(EXIT_FAILURE);                                                  \
    } }


#define CUMP_CHECK_ERROR(errorMessage) do {                                 \
    cudaError_t err = cudaGetLastError();                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
                errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
        exit(EXIT_FAILURE);                                                  \
    }                                                                        \
    err = cudaThreadSynchronize();                                           \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
                errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
        exit(EXIT_FAILURE);                                                  \
    } } while (0)

#define CU_SAVEn(addr, idx, ofs, a, _n) do { \
    for(unsigned _i = 0; _i < _n; _i++) { \
        (addr)[(idx) + _i * (ofs)] = (a)[_i]; \
    }  } while(0);

#define CU_LOADn(addr, idx, ofs, a, _n) do { \
    for(unsigned _i = 0; _i < _n; _i++) { \
        (a)[_i] = (addr)[(idx) + _i * (ofs)]; \
    }  } while(0);

#define CU_LOAD_GLOBALn_dummy(addr, idx, ofs, a, _n) do { \
    for(unsigned _i = 0; _i < _n; _i++) { \
        (a)[_i] = (limb&)(addr) + _i; \
    }  } while(0);

#define CU_LOAD_GLOBALn(addr, idx, ofs, a, _n) do { \
    for(unsigned _i = 0; _i < _n; _i++) { \
        (a)[_i] = (addr)[(idx) + _i * (ofs)]; \
    }  } while(0);

#define CU_BEGIN_TIMING() { \
    cudaThreadSynchronize(); \
    for(unsigned __iii = 0; __iii < n_iters + 1; __iii++) {

#if CUDART_VERSION < 2000 
#define CU_END_TIMING( ms ) \
        if(__iii == 0) { \
            cudaThreadSynchronize(); \
        } \
    } \
    ms /= n_iters; \
    }

#else
#define CU_END_TIMING( ms ) \
        if(__iii == 0) { \
            cudaThreadSynchronize(); \
            cudaEventRecord( e_start, 0 ); \
        } \
    } \
    cudaEventRecord( e_end, 0 ); \
    cudaEventSynchronize( e_end ); \
    cudaEventElapsedTime( &ms, e_start, e_end ); \
    ms /= n_iters; \
    }

#endif // CUDART_VERSION

#if (CUMP_BENCHMARK_CALLS && CUDART_VERSION >= 2000)

#define LOCAL_TM_BEGIN(evt1) \
    if(__iii == 1) { cudaEventRecord(evt1, 0); }
    

#define LOCAL_TM_END(evt2) \
    if(__iii == 1) { cudaEventRecord(evt2, 0); }

#define FETCH_TM(ms, evt1, evt2) do { \
    cudaEventSynchronize(evt2); \
    cudaEventElapsedTime(&ms, evt1, evt2); \
    } while(0);
#else // CUMP_BENCHMARK_CALLS

#define LOCAL_TM_BEGIN(...) \
    (void)0;
    
#define LOCAL_TM_END(...) \
    (void)0;

#define FETCH_TM(...) \
    (void)0;
#endif // CUMP_BENCHMARK_CALLS

#define FETCH_FROM_CUDA_ARRAY
#ifdef FETCH_FROM_CUDA_ARRAY
#define TEX1D(tex, index) tex1D((tex), (index))
#define TEX2D(tex, x, y)  tex2D((tex), (x), (y))
#else
#define TEX1D(tex, index) tex1Dfetch((tex), (index))
#define TEX2D(tex, x, y) (void)0
#endif

#if CUMP_PREFETCH_FROM_CUDA_ARRAY
#define TEX_ROW_FETCH(_a, _, _lin_ofs, _stride) do { \
        unsigned _ofs = _lin_ofs + _stride, \
            _row_ofs = _ofs & (TILE_SZ*TILE_SZ*N_TILES_PER_ROW - 1), \
            _lidx = _row_ofs & TILE_SZ*TILE_SZ - 1, \
            _lx = __umul24(_row_ofs >> 2*LOG_TILE_SZ, 0x1000000 + TILE_SZ) + \
                (_lidx & TILE_SZ - 1), _ly; \
        _row_ofs = _ofs / (TILE_SZ*TILE_SZ*N_TILES_PER_ROW); \
        _ly = __umul24(_row_ofs, 0x1000000 + TILE_SZ) + (_lidx / TILE_SZ); \
        (_a) = TEX2D(tiled_tex, _lx, _ly); \
    } while(0);

#define TEX_COLUMN_FETCH(_a, _, _lin_ofs, _stride) do { \
        unsigned _ofs = _lin_ofs + _stride, \
        _col_ofs = _ofs / (TILE_SZ*TILE_SZ*N_TILES_PER_COL), \
        _lx = __umul24(_col_ofs, 0x1000000 + TILE_SZ) + \
            (_ofs & TILE_SZ - 1), _ly; \
        _col_ofs = _ofs & (TILE_SZ*TILE_SZ*N_TILES_PER_COL - 1), \
        _ly = _col_ofs / TILE_SZ; \
        (_a) = TEX2D(tiled_tex, _lx, _ly); \
    } while(0);

#else
#define TEX_ROW_FETCH(_a, _g_in, _lin_ofs, _stride) do { \
        (_a) = (_g_in)[_lin_ofs + _stride]; \
    } while(0);

#define TEX_COLUMN_FETCH(_a, _g_in, _lin_ofs, _stride) \
        TEX_ROW_FETCH(_a, _g_in, _lin_ofs, _stride)

#endif 

#endif // _MACROS_H_
