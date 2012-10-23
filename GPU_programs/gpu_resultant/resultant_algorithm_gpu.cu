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
    
#ifndef _RESULTANT_ALGORITHM_GPU_CU_
#define _RESULTANT_ALGORITHM_GPU_CU_
 
#define LOG_TILE_SZ 6 // 2 ^ LOG_TILE_SZ texture dimension
#define TILE_SZ (1 << LOG_TILE_SZ)

#define N_TILES_PER_COL 16 // # of texture tiles in a single column

#if ((N_TILES_PER_COL - 1) & N_TILES_PER_COL)
#error N_TILES_PER_COL must be power-of-two!!
#endif

#if CUDART_VERSION < 2000 // NVCC 1.0 does not like exceptions
#undef __EXCEPTIONS
#endif

#include <gmp.h>
#include <fstream>
#include <iostream>

#include <include/macros.h>
#include <CGAL/GPU_algorithms/resultant_algorithm.h>

#include <resultant_kernel.cu>
#include <interpolate_CRA_kernel.cu>

namespace CGAL {

namespace internal {

// using namespace std; // workaround streampos bug..

#if CUMP_PREFETCH_FROM_CUDA_ARRAY
typedef texture<unsigned, 2, cudaReadModeElementType> Texture_descr;

Texture_descr g_tiled_tex;
cudaArray *g_tiled_tex_array = 0;

void alloc_tex_mem_col(cudaArray **tex_array, unsigned& n);

void copy_tex_mem_col(cudaArray *tex_array, Texture_descr *tex,
        const unsigned *data, unsigned n);
#endif

#if (CUDART_VERSION >= 2000)
extern cudaEvent_t e_start, e_end;
#endif

void GPU_resultant::device_static_setup() {

    DEV_bytes_allocated = 5*1024*1024;
    CUMP_SAFE_CALL(cudaMalloc(&DEV_mem_ptr, DEV_bytes_allocated));
}

//! \c CPU_mem_ptr collects all dynamic memory allocated either in ordinary
//! or page-locked memory triggered by the flag \c CUMP_USE_PAGELOCKED_MEM
//! \c DEV_mem_ptr collects all dynamic memory allocated on GPU, optionally
//! input operands can be allocated in texture space pointed to by
//! \c tiled_tex_array - must be bound to actual texture before use.
//! \c batch_size - # of inputs to allocate memory required for
//! parallel processing
//! \c word_size - bytes per one data element
//! \c nr_sz - the amount of input data to be transferred to GPU
//! \c prod_sz - size of product: allocated twice on CPU: 
//! \c aux_cpu_mem - auxiliary space allocated for CPU only
//! \c aux_dev_mem - auxiliary space allocated for device only
bool GPU_resultant::alloc_device_mem() {

    if(CUMP_N_MODULI > CUMP_DEVICE_MODULI_SZ / CUMP_MODULI_STRIDE) {
        printf("ERROR: insufficient device constant memory..\n");
        return false; 
    }

    mem_size_nr = nr_sz * word_size;
    mem_size_prod = prod_sz * batch_size * word_size;

    // don't use any sh*ty aux_batch_size, compute mem size without it
    mem_size_cpu_aux = aux_cpu_mem * word_size;
    mem_size_dev_aux = aux_dev_mem * word_size;

    unsigned new_dev_mem, new_host_mem;
#if CUMP_PREFETCH_FROM_CUDA_ARRAY
    // TODO: nr_sz must be aligned by TILE_SZ (holds by definition)
    if(g_tiled_tex_array != 0) {
        CUMP_SAFE_CALL(cudaFreeArray(g_tiled_tex_array));
    }

    alloc_tex_mem_col(&g_tiled_tex_array, nr_sz);
    mem_size_nr = nr_sz * word_size;

    new_dev_mem = mem_size_prod * 2 + mem_size_dev_aux;
#else
    // allocate write-to device-mapped memory
//     CUMP_SAFE_CALL(cudaHostAlloc(&DEV_mem_write_ptr,
//                 mem_size_nr + mem_size_dev_aux, cudaHostAllocMapped));

//!  mem_size_prod * 2 because resultants_kernel outputs
//! two values: modular resultant and denominator !!!
    new_dev_mem = mem_size_nr + mem_size_prod * 4 + mem_size_dev_aux;
#endif // CUMP_PREFETCH_FROM_CUDA_ARRAY

    if(new_dev_mem > DEV_bytes_allocated || DEV_mem_ptr == 0) {
        if(DEV_mem_ptr != 0)
            CUMP_SAFE_CALL(cudaFree(DEV_mem_ptr));
        DEV_bytes_allocated = new_dev_mem;
        CUMP_SAFE_CALL(cudaMalloc(&DEV_mem_ptr, new_dev_mem));
    }

    if(DEV_mem_ptr == 0) {
        printf("ERROR: unable to allocate device mem..\n");
        return false;
    }
    new_host_mem = mem_size_nr + mem_size_prod*2 + mem_size_cpu_aux;

    if(new_host_mem > CPU_bytes_alloced || CPU_mem_ptr == 0) {
#if CUMP_USE_PAGELOCKED_MEM
#warning using page-locked mem
        if(CPU_mem_ptr != 0)
            cudaFreeHost(CPU_mem_ptr);
        CUMP_SAFE_CALL(cudaMallocHost(&CPU_mem_ptr, new_host_mem));
#else
        if(CPU_mem_ptr != 0)
            free(CPU_mem_ptr);
        CPU_mem_ptr = malloc(new_host_mem);
#endif
        CPU_bytes_alloced = new_host_mem;
    }

    if(CPU_mem_ptr == 0) {
        printf("ERROR: unable to allocate CPU mem..\n");
        return false;
    }
    return true;
}

//! Input: \c const_mem: (CUMP_N_MODULI * CUMP_MODULI_STRIDE + 2)
//!        \c Mods: mem_size_dev_aux (Mods, Inks, Mus, InvMods)
//!        \c U: mem_size_nr (set of polynomial coefficients)
//! Output: \c R: mixed radix digits
void GPU_resultant::launch_kernel(const unsigned *const_mem,
        const unsigned *Mods, const unsigned *U, unsigned *R) {

    //! additional 2: for copying nu & nv
    CUMP_SAFE_CALL(cudaMemcpyToSymbol(dev_const_mem, const_mem,
            (CUMP_N_MODULI * CUMP_MODULI_STRIDE + 2) * sizeof(unsigned)));

    double flops; // estimated flop count per one block
    float ms;

#if (CUMP_BENCHMARK_CALLS && CUDART_VERSION >= 2000)
    const unsigned N_TIMERS = 4;
    cudaEvent_t tma[N_TIMERS], tmb[N_TIMERS];
    for(unsigned i = 0; i < N_TIMERS; i++) {
        cudaEventCreate(&tma[i]);
        cudaEventCreate(&tmb[i]);        
    }
#endif

#if CUMP_COMPILE_DEBUG_KERNEL
    unsigned n_iters = 10;
#else
    unsigned n_iters = 0;
#endif

#if CUMP_MEASURE_MEM
    CU_BEGIN_TIMING()
#endif

#if CUMP_PREFETCH_FROM_CUDA_ARRAY

    CUMP_out2("Fetching from texture mem: total of %.3f Kb\n",
        (double)mem_size_nr / 1024.0)
    //! + aux_dev_mem ??
    copy_tex_mem_col(g_tiled_tex_array, &g_tiled_tex, U, nr_sz); 
    unsigned *devU = 0, *devR = (unsigned *)DEV_mem_ptr;
#else
    //! \c devR - stores intermediate result, \c devR2 - final result
    //! (additional space for coalesced mem access)
    unsigned *devU = (unsigned *)DEV_mem_ptr, *devMods = devU + nr_sz,
        *devR = devMods + aux_dev_mem, *devR2 = devR + prod_sz * batch_size,
        *devR3 = devR2 + prod_sz * batch_size,
        *devR4 = devR3 + prod_sz * batch_size;
    cudaMemcpy(devU, U, mem_size_nr, cudaMemcpyHostToDevice);

    // NOTE: use async copy here ??
    cudaMemcpy(devMods, Mods, mem_size_dev_aux, cudaMemcpyHostToDevice);
    unsigned *devModIndex = 0;
    CUMP_SAFE_CALL(cudaGetSymbolAddress((void **)&devModIndex, dev_mod_index));
#endif // CUMP_PREFETCH_FROM_CUDA_ARRAY

#if !CUMP_MEASURE_MEM
    CU_BEGIN_TIMING()
#endif

    CUMP_SAFE_CALL(cudaMemset(devModIndex, 0,
             CUMP_N_MODULI * sizeof(unsigned)));

    flops = 0;
    LOCAL_TM_BEGIN(tma[0]);    
    resultant_kernel_dispatch(devR, devR2, devU);
    LOCAL_TM_END(tmb[0]);    

#if !CUMP_RUN_RESULTANTS_ONLY
    LOCAL_TM_BEGIN(tma[1]);
    // such a mess of address spaces - not desired
    mod_inverse_kernel1_dispatch(devR, devR2, devR3, devR4);
    LOCAL_TM_END(tmb[1]);
    
    LOCAL_TM_BEGIN(tma[2]);
    interpolate_kernel_dispatch(devR3, devR4, devU);
    LOCAL_TM_END(tmb[2]);

    LOCAL_TM_BEGIN(tma[1]);
    mod_inverse_kernel2_dispatch(devU, devU, devMods);
    LOCAL_TM_END(tmb[1]);
#endif
 
    LOCAL_TM_BEGIN(tma[3]);
#if CUMP_USE_CRA_KERNEL
    CRA_kernel_dispatch(devR4, devR3, devMods, devU);
#endif
    LOCAL_TM_END(tmb[3]);

#if !CUMP_MEASURE_MEM
    CU_END_TIMING(ms)
#endif

#if !CUMP_RUN_RESULTANTS_ONLY
   cudaMemcpy(R, devModIndex, CUMP_N_MODULI * sizeof(unsigned),
         cudaMemcpyDeviceToHost);

    bool failed = false;
    n_pts_failed = 0;
    for(unsigned i = 0; i < CUMP_N_MODULI; i++) {

        unsigned diff = n_real_pts - R[i]; // failed if R[i] < n_real_pts
        if((int)diff > 0) {
            failed = true;
            n_pts_failed = std::max(n_pts_failed, diff);
        }
    }
    if(failed) {
        printf("Failed: NSR: %d\n", n_pts_failed);
        CUMP_THROW(GPU_algorithm_exception);
    }
#endif

#if CUMP_USE_CRA_KERNEL
    cudaMemcpy(R, devR4, mem_size_prod, cudaMemcpyDeviceToHost);
//     cudaMemcpy(R, devR, mem_size_prod, cudaMemcpyDeviceToHost);
#elif CUMP_RUN_RESULTANTS_ONLY
    cudaMemcpy(R, devR, mem_size_prod, cudaMemcpyDeviceToHost);
#else // CUMP_RUN_RES_AND_INTERPOLATE
    cudaMemcpy(R, devR3, mem_size_prod, cudaMemcpyDeviceToHost);
//     cudaMemcpy(R, devR2, mem_size_prod, cudaMemcpyDeviceToHost);
#endif

#if CUMP_MEASURE_MEM
    CU_END_TIMING(ms)
#endif

#if CUMP_BENCHMARK_CALLS
    double s = ms / 1000.0;

    float res_ms, mod_inv_ms, int_ms, crt_ms;
    FETCH_TM(res_ms, tma[0], tmb[0]);
    FETCH_TM(mod_inv_ms, tma[1], tmb[1]);
    FETCH_TM(int_ms, tma[2], tmb[2]);
    FETCH_TM(crt_ms, tma[3], tmb[3]);
        
    CUMP_out2("batch_size: %d\n", batch_size);
    //NOTE: flops counter is corrupt..
    double Gflop = 1e-9 * flops * batch_size;
    double GB = 1e-9 * (double)(mem_size_nr + mem_size_prod);
        
    CUMP_out2("GPU time elapsed: %f ms; Gflop/s: %f; GB/s: %f\n\n", ms, Gflop/s,
            GB/s);
    CUMP_out2("res time: %f ms; mod_inv time: %f ms; interp time: %f ms; "
        "CRA time: %f ms\n", res_ms, mod_inv_ms, int_ms, crt_ms);

    CUMP_out2("\n\namount of data transferred to device: %.2f Kb\n"
            "amount of data transferred back to CPU: %.2f Kb\n",
            (double)(mem_size_nr + mem_size_dev_aux)/1024.0,
            (double)(mem_size_prod)/1024.0);

#if CUDART_VERSION >= 2000
    for(unsigned i = 0; i < N_TIMERS; i++) {
        cudaEventDestroy(tma[i]);
        cudaEventDestroy(tmb[i]);        
    }
#endif

    if(write_to_file) {
        std::ofstream out(benchmark_file, std::ios::app);
//         out << "# res: " << res_ms << " mod_inv: " << mod_inv_ms <<
//             " interp: " << int_ms << " CRA: " << crt_ms <<
//             " total_GPU: " << ms << "\n";
        out << " res: " << res_ms << " mod_inv: " << mod_inv_ms <<
            " interp: " << int_ms << " CRA: " << crt_ms <<
            " total_GPU: " << ms;
    }
#endif // CUMP_BENCHMARK_CALLS
}

void GPU_resultant::resultant_data_padding() {

#if 0
    unsigned min_nr;
    if(nr >= 32 && nr <= 31*2) {
        min_nr = 32, max_nr = 64;
        // due to saving offset, for this kernel the number of points must be
        // even !!!
        CUMP_N_BATCHES_PER_MOD += (CUMP_N_BATCHES_PER_MOD & 1);
        CUMP_out("*********** 32-thids kernel invoke\n");
    
    } else if(nr <= 63*2) {
        min_nr = 31*2+1, max_nr = 128;
        CUMP_out("*********** 64-thids kernel invoke\n");

    } else if(nr <= 95*2) {
        min_nr = 63*2+1, max_nr = 96*2;
        CUMP_out("*********** 96-thids kernel invoke\n");

    } else if(nr <= 127*2) {
        min_nr = 95*2+1, max_nr = 128*2;
        CUMP_out("*********** 128-thids kernel invoke\n");

    } else {
        CUMP_out2("unsupported poly degree: nr: %d\n", nr);
        CUMP_THROW(GPU_algorithm_exception);
    }
    
    prod_sz = 1; // each block outputs a single elem
    unsigned _nv = nr - nu;
    if(nr < min_nr || nr > max_nr || nu > max_nr/2 - 1 || _nv > nu
            || _nv > nr) {
        CUMP_out2("incorrect input: nr: %u; nu: %u\n", nr, nu);
        CUMP_THROW(GPU_algorithm_exception);
    }
#else
    // choose kernel size depending on maximal degree 'nu': to handle
    // unbalanced operands correctly
    // for nu < 16: use 4 resultants per kernel

    nv = nr - nu, max_uv = std::max(nv, nu);
    
    if(max_uv <= 31) {
        max_nr = 64;
        // due to saving offset, for this kernel the number of points must be
        // even !!!
        CUMP_N_BATCHES_PER_MOD += (CUMP_N_BATCHES_PER_MOD & 1);
        CUMP_out("*********** 32-thids kernel invoke\n");
    
    } else if(max_uv <= 63) {
        max_nr = 128;
        CUMP_out("*********** 64-thids kernel invoke\n");

    } else if(max_uv <= 95) {
        max_nr = 96*2;
        CUMP_out("*********** 96-thids kernel invoke\n");

    } else if(max_uv <= 127) {
        max_nr = 128*2;
        CUMP_out("*********** 128-thids kernel invoke\n");

    } else {
        CUMP_out2("unsupported poly degree: nr: %d\n", nr);
        CUMP_THROW(GPU_algorithm_exception);
    }
    
    prod_sz = 1; // each block outputs a single elem
    
    if(nr > max_nr /*|| _nv > nu*/ || nv > nr) {
        CUMP_out2("incorrect input: nr: %u; nu: %u\n", nr, nu);
        CUMP_THROW(GPU_algorithm_exception);
    }
#endif
} 

void GPU_resultant::resultant_kernel_dispatch(unsigned *devR, unsigned *devR2,
    unsigned *devU) {

    dim3 threads;
    // this is a full grid size (without streams)
    dim3 grid(CUMP_N_BATCHES_PER_MOD, CUMP_N_MODULI, 1);
    //n_blocks = grid.x * grid.y;
    unsigned shm_size;

#if 0
// lite kernel size is: (3 * 32 + 16) * 2 (for two warps)

#define DISPATCH_CALL(__n) \
    else if(nr <= __n*2 - 2) { \
        threads = dim3(__n, 1, 1); \
        shm_size = (3 * __n + 16) * sizeof(unsigned); \
        resultant_block_kernel< __n ><<<grid, threads, shm_size>>> \
                 (devR, devR2, devU, deg_x1, deg_x2, out_data_padding); \
    }

    if(nr <= 32*2 - 2) {
        threads = dim3(32, 2, 1);
        grid = dim3((CUMP_N_BATCHES_PER_MOD + 1)/2, CUMP_N_MODULI, 1);
        shm_size = (3 * 32 + 16) * 2 * sizeof(unsigned); 
        resultant_block_kernel< 32 ><<<grid, threads, shm_size>>> 
                 (devR, devR2, devU, deg_x1, deg_x2, out_data_padding); 
    }
    DISPATCH_CALL(64)
    DISPATCH_CALL(96)
    DISPATCH_CALL(128)
#else

    #define DISPATCH_CALL(__n) \
    else if(max_uv <= __n - 1) { \
        threads = dim3(__n, 1, 1); \
        shm_size = (3 * __n + 16) * sizeof(unsigned); \
        resultant_block_kernel< __n ><<<grid, threads, shm_size>>> \
                 (devR, devR2, devU, deg_x1, deg_x2, out_data_padding); \
    }

    if(max_uv <= 31) {
        threads = dim3(32, 2, 1);
        grid = dim3((CUMP_N_BATCHES_PER_MOD + 1)/2, CUMP_N_MODULI, 1);
        shm_size = (3 * 32 + 16) * 2 * sizeof(unsigned); 
        resultant_block_kernel< 32 ><<<grid, threads, shm_size>>> 
                 (devR, devR2, devU, deg_x1, deg_x2,
                     out_data_padding); 
    }
    DISPATCH_CALL(64)
    DISPATCH_CALL(96)
    DISPATCH_CALL(128)
#endif

#undef DISPATCH_CALL
}

// TODO TODO TODO: need doubles_kernel for interpolation 

//! determine data padding depending on the # of evaluation points (nbatches)
//! for resultant algorithm
//! this must reflect the # of interpolate kernels
void GPU_resultant::interpolate_data_padding() {

// it is assumed that this flag is used for sm20 architecture
#if CUMP_USE_32BIT_MODULI_SET 
    const unsigned MAX_THIDS = 1024;
#else
    const unsigned MAX_THIDS = 512;
#endif        

    unsigned nthids = ((CUMP_N_BATCHES_PER_MOD + 3) / 4 + 31) & ~31;
    if(nthids < 64)
        nthids = 64;
    else if(nthids > MAX_THIDS) {
        CUMP_out2("unsupported batch size: %d %d\n", CUMP_N_BATCHES_PER_MOD,
            MAX_THIDS);
        CUMP_THROW(GPU_algorithm_exception);
    }
    out_data_padding = nthids*4;

    CUMP_out2("nbatches %d; nthids: %d; padding: %d\n",
        CUMP_N_BATCHES_PER_MOD, nthids, out_data_padding);
}

//! \c devR - in: poly evaluations; out: coefficients
//! \c devInvDets - write out denominators
//! \c devXs - interpolation points
void GPU_resultant::interpolate_kernel_dispatch(unsigned *devR,
    const unsigned *devXs, unsigned *devInvDets) {

// TODO: use doubles_kernel when # of interpolation points <= 256

//! IMPORTANT: we use data padding according to CUMP_N_BATCHES_PER_MOD
//! but the # of threads (interpolation points) is given by \c n_real_pts !!

    unsigned shm_size;
    unsigned nthids = ((n_real_pts + 3) / 4 + 31) & ~31;
    if(nthids < 64)
        nthids = 64;

    dim3 threads(nthids, 1, 1);
    dim3 grid(CUMP_N_MODULI, 1, 1);
    shm_size = (3 * out_data_padding / 4) * sizeof(unsigned);

    CUMP_out2("interpolate kernel thids: %d; shmem_size: %f Kb\n",
            nthids, shm_size / 1024.0);

#define DISPATCH_CALL(__n__) else if(nmod4 == __n__) { \
        interpolate_quads_kernel< __n__ ><<<grid, threads, shm_size>>> \
            (devR, devXs, devInvDets, n_real_pts, out_data_padding); \
    }

    unsigned nmod4 = n_real_pts & 3;

    if(0);
    DISPATCH_CALL(0)
    DISPATCH_CALL(1)
    DISPATCH_CALL(2)
    DISPATCH_CALL(3)

#undef DISPATCH_CALL
}

void GPU_resultant::mod_inverse_kernel1_dispatch(unsigned *devR,
        unsigned *devR2, unsigned *devR3,
        unsigned *devR4) {

    dim3 threads(128, 1, 1);
    unsigned nx = (CUMP_N_BATCHES_PER_MOD + 127) / 128,
            shm_size = (256 + 16) * sizeof(unsigned);
    dim3 grid(nx, CUMP_N_MODULI, 1);

    printf("nx: %d; CUMP_N_MODULI: %d\n", nx, CUMP_N_MODULI);

    mod_inverse_kernel1<<< grid, threads, shm_size >>>(
        devR, devR2, devR3, devR4, CUMP_N_BATCHES_PER_MOD,
                 out_data_padding);
}

//! \c devR - write out inverse of denominators
//! \c devU - input denominators
//! \c devMods - set of moduli (as well as invk and mu), cannot be read
//! from constant mem because moduli are allocated on per-thread basis
void GPU_resultant::mod_inverse_kernel2_dispatch(unsigned *devR,
     const unsigned *devU, const unsigned *devMods) {

    dim3 threads(64, 1, 1);
    unsigned nx = (CUMP_N_MODULI + 63) / 64;
    dim3 grid(nx, 1, 1);
    
    // only problem is where to save modular inverses initially...?
    mod_inverse_kernel2<<< grid, threads >>>(devR, devU,
         devMods, mods_padding, CUMP_N_MODULI);
}

void GPU_resultant::CRA_kernel_dispatch(unsigned *devOut,
    const unsigned *devIn, const unsigned *devMods,
        const unsigned *devInvDets) {

    dim3 grid, threads;
//! IMPORTANT: we use data padding according to CUMP_N_BATCHES_PER_MOD
//! but the # of threads (interpolation points) is given by \c n_real_pts !!
    if(CUMP_N_MODULI <= 512) { // use doubles kernel up to 256 threads

         unsigned nx = ((CUMP_N_MODULI + 1) / 2 + 31) & ~31,
            shm_size = 16 * sizeof(unsigned);

        if(nx == 32) {
            threads = dim3(32, 2, 1);
            grid = dim3((n_real_pts + 1)/2, 1, 1);

            CUMP_out("**** invoking 32x2 CRA kernel ****\n");

            if(CUMP_N_MODULI & 1) {
                CRA_doubles_kernel< 32, false ><<< grid, threads, shm_size >>>
                    (devOut, devIn, devMods, devInvDets, mods_padding,
                         CUMP_N_MODULI, out_data_padding);
            } else {
                CRA_doubles_kernel< 32, true ><<< grid, threads, shm_size >>>
                    (devOut, devIn, devMods, devInvDets, mods_padding,
                         CUMP_N_MODULI, out_data_padding);
            }
        } else {
            if(nx > 512) {
                CUMP_out2("unsupported CRA size: %d\n", nx);
                CUMP_THROW(GPU_algorithm_exception);
            }

            CUMP_out("**** invoking CRA doubles kernel ****\n");

            threads = dim3(nx, 1, 1);
            grid = dim3(n_real_pts, 1, 1);

            // here BlockSz acts simply as a switch (not an actual block size)
            if(CUMP_N_MODULI & 1) {
                CRA_doubles_kernel< 64, false ><<< grid, threads, shm_size >>>
                    (devOut, devIn, devMods, devInvDets, mods_padding,
                         CUMP_N_MODULI, out_data_padding);
            } else {
                CRA_doubles_kernel< 64, true ><<< grid, threads, shm_size >>>
                (devOut, devIn, devMods, devInvDets, mods_padding,
                     CUMP_N_MODULI, out_data_padding);
            }
        }

    } else { // CUMP_N_MODULI > 512

        CUMP_out("**** invoking CRA quads kernel ****\n");

        unsigned nx = ((CUMP_N_MODULI + 3) / 4 + 31) & ~31,
            shm_size = 16 * sizeof(unsigned);

        if(nx > 512) {
            CUMP_out2("unsupported CRA size: %d\n", nx);
            CUMP_THROW(GPU_algorithm_exception);
        }

        threads = dim3(nx, 1, 1);
        grid = dim3(n_real_pts, 1, 1);
        CUMP_out("**** invoking quad CRA kernel ****\n");

        unsigned nmod4 = CUMP_N_MODULI % 4;
        // here BlockSz acts simply as a switch (not an actual block size)
        if(nmod4 == 0) {
            CRA_quads_kernel< 64, 0 ><<< grid, threads, shm_size >>>
                (devOut, devIn, devMods, devInvDets, mods_padding,
                     CUMP_N_MODULI, out_data_padding);
        } else if(nmod4 == 1) {
            CRA_quads_kernel< 64, 1 ><<< grid, threads, shm_size >>>
                (devOut, devIn, devMods, devInvDets, mods_padding,
                     CUMP_N_MODULI, out_data_padding);
        } else if(nmod4 == 2) {
            CRA_quads_kernel< 64, 2 ><<< grid, threads, shm_size >>>
                (devOut, devIn, devMods, devInvDets, mods_padding,
                     CUMP_N_MODULI, out_data_padding);
        } else {
            CRA_quads_kernel< 64, 3 ><<< grid, threads, shm_size >>>
                (devOut, devIn, devMods, devInvDets, mods_padding,
                     CUMP_N_MODULI, out_data_padding);
        }
    }
}

//! allocates \c n (in *words*) tiled texture memory \c tiled_tex_array in
//! column-major order
void alloc_tex_mem_col(cudaArray **tiled_tex_array, unsigned& n) {

    n = (n + TILE_SZ - 1) & ~(TILE_SZ - 1);

    unsigned n_tiles = (n + TILE_SZ*TILE_SZ-1) >> (2*LOG_TILE_SZ),
         w = (n_tiles + N_TILES_PER_COL-1) / N_TILES_PER_COL, h = n_tiles;

    if(h > N_TILES_PER_COL)
        h = N_TILES_PER_COL;

    CUMP_out2("tex_mem_size: %d words; n_tiles: %d: %d x %d\n", n, n_tiles,
            w, h);

    cudaChannelFormatDesc desc = cudaCreateChannelDesc< unsigned >();
        //(8, 0, 0, 0, cudaChannelFormatKindUnsigned);

    CUMP_SAFE_CALL(cudaMallocArray(tiled_tex_array, &desc,
         w * TILE_SZ * sizeof(unsigned), h * TILE_SZ));
}

#if CUMP_PREFETCH_FROM_CUDA_ARRAY
//! copies \c sz *words* of \c data to texture memory pointed to by
//! \c tex_array in column-major order and binds it to texture \c tex
void copy_tex_mem_col(cudaArray *tex_array,
        Texture_descr *tex, const unsigned *data, unsigned n) {

    const unsigned *src = data;
    for(int sz = n, ofsx = 0; sz > 0; ofsx += TILE_SZ * sizeof(unsigned)) {

        // process data in column-major order
        for(int ofsy = 0, i = 0; (sz > 0 && i < N_TILES_PER_COL);
            sz -= TILE_SZ*TILE_SZ, src += TILE_SZ*TILE_SZ, i++,
                    ofsy += TILE_SZ) {

            unsigned w = TILE_SZ * sizeof(unsigned),
                h = (sz >= TILE_SZ*TILE_SZ ? TILE_SZ : sz / TILE_SZ);
            if(h < TILE_SZ)
                CUMP_out2("last block height: %d\n", h);

        // TODO: you can copy full columns (ie in one loop)
            CUMP_SAFE_CALL(cudaMemcpy2DToArray(tex_array, ofsx, ofsy,
                 src, w, w, h, cudaMemcpyHostToDevice));
        }
    }
//     unsigned w = TILE_SZ * sizeof(unsigned), h = n / TILE_SZ;
//     CUMP_SAFE_CALL(cudaMemcpy2DToArray(tiled_tex_array, 0, 0,
//                   data, w, w, h, cudaMemcpyHostToDevice));

    cudaChannelFormatDesc desc = cudaCreateChannelDesc<unsigned>();
    CUMP_SAFE_CALL(cudaBindTextureToArray(tex, tex_array, desc));
}
#endif // CUMP_PREFETCH_FROM_CUDA_ARRAY

GPU_resultant::~GPU_resultant() {
    free_device_mem();
}

void GPU_resultant::free_device_mem() {

    CUMP_out2("GPU_res: free_device_mem\n")

    if(CPU_mem_ptr != 0) {
#if CUMP_USE_PAGELOCKED_MEM
        
        cudaFreeHost(CPU_mem_ptr);
#else
        free(CPU_mem_ptr);
#endif
        CPU_mem_ptr = 0;
    }

    if(DEV_mem_ptr != 0) {
        CUMP_SAFE_CALL(cudaFree(DEV_mem_ptr));
        DEV_mem_ptr = 0;
    }

#if CUMP_PREFETCH_FROM_CUDA_ARRAY
    if(g_tiled_tex_array != 0) {
        CUMP_SAFE_CALL(cudaFreeArray(g_tiled_tex_array));
        g_tiled_tex_array = 0;
    }
#endif

}

} // namespace internal 

} // namespace CGAL 

#endif // _RESULTANT_ALGORITHM_GPU_CU_
