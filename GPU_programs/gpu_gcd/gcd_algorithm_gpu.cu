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

// #include <include/cudart_version.h>

#if CUDART_VERSION < 2000
#undef __EXCEPTIONS // disable exceptions for CUDA 1.0 (bug workaround)
#endif

#include <iostream>
#include <fstream>

#include <include/macros.h>
#include <CGAL/GPU_algorithms/gcd_algorithm.h>
#include <include/const_mem_defs.h>

#ifdef CUDART_VERSION 
#if CUMP_USE_STREAMS && (CUDART_VERSION < 2000)
#error Streams are not supported on this CUDA platform
#endif
#endif

#if (CUMP_RUN_GCD_ONLY && CUMP_MRC_HOST_RECOVER)
#error Flags are mutually exclusive !!
#endif

#include <mod_reduce_kernel.cu>

#if (CUMP_USE_GCD_BLOCK_KERNEL && !CUMP_USE_PGCD_KERNEL)
#include <QR_GCD_block_kernel.cu>
#include <QR_GCD_lite_kernel.cu>
#endif

#if (!CUMP_USE_GCD_BLOCK_KERNEL && !CUMP_USE_PGCD_KERNEL)
#warning Using monolithic GCD kernel
#include <QR_GCD_monolithic_kernel.cu>
#endif

#if CUMP_USE_PGCD_KERNEL
#include <PGCD_kernel.cu>
#endif

#if !CUMP_RUN_GCD_ONLY
#include <MRC_block_kernel.cu>
#endif

namespace CGAL {

namespace internal {

#if CUDART_VERSION >= 2000
extern cudaEvent_t e_start, e_end;
#endif

cudaStream_t *g_streams = NULL;

void GPU_gcd::device_static_setup() {

    DEV_bytes_allocated = 5*1024*1024;
    CUMP_SAFE_CALL(cudaMalloc(&DEV_mem_ptr, DEV_bytes_allocated));
}

bool GPU_gcd::alloc_device_mem() {

    if(n_moduli > CUMP_DEVICE_MODULI_SZ / CUMP_MODULI_STRIDE) {
        printf("ERROR: insufficient device constant memory..\n");
        return false; 
    }

    mem_size_in = data_sz * word_size;
    mem_size_prod = prod_sz * word_size;
    mem_size_cpu_aux = aux_cpu_mem * word_size;
    mem_size_dev_aux = aux_dev_mem * word_size;

    unsigned new_dev_mem, new_host_mem;
    // allocate write-to device-mapped memory
//     CUMP_SAFE_CALL(cudaHostAlloc(&DEV_mem_write_ptr,
//                 mem_size_in + mem_size_dev_aux, cudaHostAllocMapped));

    new_dev_mem = mem_size_in + mem_size_prod + mem_size_dev_aux;

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

    new_host_mem = mem_size_in + mem_size_prod*2 + mem_size_cpu_aux;
#if CUMP_COMPILE_DEBUG_KERNEL
//     new_host_mem += mem_size_prod; // alloc space for reference solution only
                                   // in debug mode
#endif

    if(new_host_mem > CPU_bytes_alloced || CPU_mem_ptr == 0) {
#if CUMP_USE_PAGELOCKED_MEM
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

void GPU_gcd::launch_kernel(unsigned *const_mem,
    const unsigned *Mods, const unsigned *U, unsigned *R, bool& coprime) {

#if CUMP_COMPILE_GCD_KERNEL

#if CUMP_USE_GCD_BLOCK_KERNEL
    unsigned *devU = (unsigned *)DEV_mem_ptr,
        *devCPU_In = devU + data_padding * n_moduli,
        *tmpLimbs = devCPU_In + mods_padding,
        *devR = devCPU_In + data_sz,
        *devIn = devR + prod_sz, *devOut = devIn + block_io_sz,
        *devInvDets = devOut + block_io_sz,
        *devMods = devInvDets + mods_padding;
#else
    //! NOTE change the size of devR to 'prod_sz'
    unsigned *devU = (unsigned *)DEV_mem_ptr,
        *devR = devU + data_sz,
        *devMods = devR + data_sz, *devOuts = devMods + aux_cpu_mem;
#endif

    const_mem[NU] = nu, const_mem[NV] = nv;
    const_mem[NU_OFS4] = nu_ofs4, const_mem[NV_OFS4] = nv_ofs4;
    const_mem[MAX_NU] = max_nu;
    const_mem[DATA_PAD] = data_padding;
    const_mem[GCD_SZ] = gcd_sz;
    const_mem[N_BLOCKS] = n_blocks; // keep original # of blocks needed to
                                    // calculate proper batch offset
//     const_mem[DATA_IN] = (unsigned)devU;
//     const_mem[DATA_OUT] = (unsigned)devR;
    const_mem[MODS_PAD] = mods_padding;

    cudaMemcpyToSymbol(dev_const_mem, const_mem,
            (n_moduli * c_stride + ConstMemDataSz) * word_size);

//     printf("devU: %x; devR: %x\n", devU, devR);
            
//     g_streams = new cudaStream_t[n_streams];
//     for(unsigned i = 0; i < n_streams; i++) {
//         cudaStreamCreate(&g_streams[i]);
//     }

#if CUMP_COMPILE_DEBUG_KERNEL
    unsigned n_iters = 0;
#else
    unsigned n_iters = 0;
#endif
    float ms = 0;

#if !CUMP_USE_GCD_BLOCK_KERNEL
    unsigned *devOutsFlag = 0;
    cudaGetSymbolAddress((void **)&devOutsFlag, dev_outs_flag);
#else
/*    unsigned *devIters = 0;
    cudaGetSymbolAddress((void **)&devIters, dev_iters);
    cudaMemset(devIters, 0, CUMP_DEVICE_MODULI_SZ / CUMP_MODULI_STRIDE *
             sizeof(unsigned));*/
#endif

#if (CUMP_BENCHMARK_CALLS && CUDART_VERSION >= 2000)
    const unsigned N_TIMERS = 5;
    cudaEvent_t tma[N_TIMERS], tmb[N_TIMERS];
    for(unsigned i = 0; i < N_TIMERS; i++) {
        cudaEventCreate(&tma[i]);
        cudaEventCreate(&tmb[i]);        
    }
#endif

    unsigned *mem_out;
#if CUMP_MEASURE_MEM
    CU_BEGIN_TIMING()
#endif

    cudaMemset(devU, 0, data_padding * n_moduli * sizeof(unsigned));
    cudaMemcpy(devCPU_In, U, mem_size_in, cudaMemcpyHostToDevice);

#if !CUMP_USE_GCD_BLOCK_KERNEL
//     cudaMemset(devR, 0xcc, mem_size_in);
    cudaMemset(devOutsFlag, 0, CUMP_DEVICE_MODULI_SZ / CUMP_MODULI_STRIDE *
             sizeof(unsigned));
#else
//! zero out all IO memory to make sure the last block does not read garbage
//! during ``lite'' iterations
#if !CUMP_USE_PGCD_KERNEL
    cudaMemset(devIn, 0, block_io_sz * 2 * sizeof(unsigned));
#endif
#endif
    //! NOTE: use async copy here ?? this is not used now
    cudaMemcpy(devMods, Mods, aux_cpu_mem * word_size,
         cudaMemcpyHostToDevice);

#if !CUMP_MEASURE_MEM
    CU_BEGIN_TIMING()
#endif

    LOCAL_TM_BEGIN(tma[0]);
    mod_reduce_kernel_dispatch(tmpLimbs, devU, devMods);
    LOCAL_TM_END(tmb[0]);

    LOCAL_TM_BEGIN(tma[1]);
#if CUMP_USE_GCD_BLOCK_KERNEL
    // use 'R' as a page-locked buf during mem transfer
#if CUMP_USE_PGCD_KERNEL
    QR_gcd_kernel_dispatch(mem_out, devR, devOut, devU, R);
#else
    QR_gcd_kernel_dispatch(mem_out, devIn, devOut, devU, R);
#endif
    
#else
    QR_gcd_kernel_dispatch(mem_out, devR, devOuts, devU, R);
#endif
    LOCAL_TM_END(tmb[1]);

#if !CUMP_RUN_GCD_ONLY
    // get the size of gcd in the first word
    cudaMemcpy(R, devR, n_moduli*sizeof(unsigned), cudaMemcpyDeviceToHost);

    for(unsigned i = 1; i < n_moduli; i++) {
//         printf("%d gcd sz: %d\n", i, R[i]);
        if(R[i] != R[0]) {
            std::cerr << "Failed: gcd_real_sz\n";
            CUMP_THROW(GPU_algorithm_exception);
        }
    }
    gcd_real_sz = R[0];
    
    CUMP_out2("\n############ gcd_real_sz: %d\n", gcd_real_sz)
    LOCAL_TM_BEGIN(tma[2]);

    coprime = (gcd_real_sz == 1);
    if(coprime) {
        return;
    }

    unsigned *devLcfs = devR + mods_padding, *devGCDLcf = devCPU_In;
        
#if CUMP_USE_PGCD_KERNEL
    devLcfs += gcd_real_sz - 1;  //! PGCD kernel saves results *not* in reversed order
#endif 

    mod_inverse_kernel1_dispatch(devInvDets, devLcfs, devMods, devGCDLcf);
    LOCAL_TM_END(tmb[2]);

    //! NOTE NOTE: here we use devIn because it has larger size: block_io_sz*2
    LOCAL_TM_BEGIN(tma[3]);
    MRC_kernel_dispatch(devIn, devR + mods_padding, devMods, devInvDets);
    LOCAL_TM_END(tmb[3]);
#endif // !CUMP_RUN_GCD_ONLY

#if !CUMP_MEASURE_MEM
    CU_END_TIMING(ms)
#endif

#if CUMP_RUN_GCD_ONLY
//! NOTE remark that the first 'mods_padding' elements of 'devR'
//! keep the degrees of images
    cudaMemcpy(R, devR, mem_size_prod, cudaMemcpyDeviceToHost);
//    cudaMemcpy(R, mem_out, mem_size_prod, cudaMemcpyDeviceToHost);
#else
    //! try \c mem_size_prod if not working
    //! NOTE here we use devIn because it has larger size: block_io_sz*2
    cudaMemcpy(R, devIn, gcd_real_sz * mods_padding * sizeof(unsigned),
             cudaMemcpyDeviceToHost);
#endif // CUMP_RUN_GCD_ONLY

#if CUMP_MEASURE_MEM
    CU_END_TIMING(ms)
#endif

//     if(g_streams != NULL) {
//         for(unsigned i = 0; i < n_streams; i++) {
//             cudaStreamDestroy(g_streams[i]);
//         }
//         delete []g_streams;
//     }

#if CUMP_BENCHMARK_CALLS
    double s = ms / 1000.0;
    double GB = 1e-9 * (double)(mem_size_in + mem_size_prod);

    float mod_cvt_ms, gcd_ms, mod_inv_ms, mrc_ms;
    FETCH_TM(mod_cvt_ms, tma[0], tmb[0]);
    FETCH_TM(gcd_ms, tma[1], tmb[1]);
    FETCH_TM(mod_inv_ms, tma[2], tmb[2]);
    FETCH_TM(mrc_ms, tma[3], tmb[3]);

    CUMP_out2("GPU time elapsed: %f ms; GB/s: %f\n\n", ms, GB/s);
    CUMP_out2("mod_cvt: %f ms; gcd: %f ms; mod_inv: %f ms; mrc: %f ms\n\n",
            mod_cvt_ms, gcd_ms, mod_inv_ms, mrc_ms);
    CUMP_out2("\n\namount of data transferred to device: %.2f Kb\n"
            "amount of data transferred back to CPU: %.2f Kb\n",
            (double)(mem_size_in + aux_cpu_mem * word_size)/1024.0,
            (double)(mem_size_prod)/1024.0);

    if(write_to_file) {
        std::ofstream out(benchmark_file, std::ios::app);
        out << " gpu_total: " << ms << " gcd: " << gcd_ms << " mod_inv: " <<
             mod_inv_ms << " mrc: " << mrc_ms ;
    }

#if CUDART_VERSION >= 2000
    for(unsigned i = 0; i < N_TIMERS; i++) {
        cudaEventDestroy(tma[i]);
        cudaEventDestroy(tmb[i]);        
    }
#endif

#endif // CUMP_BENCHMARK_CALLS

#else
#warning sqfree_factorize_kernel: dummy compilation
#endif // CUMP_COMPILE_GCD_KERNEL
}

// choosing the chunk size and # of blocks depending on the input paramenters
void GPU_gcd::QR_gcd_kernel_data_padding() {

    chunk_sz = 128;
    // maximal number of blocks
    unsigned nmax = (nu+nv + chunk_sz-1) / (2*chunk_sz);
    // # of blocks is computed such that 'nu+1' elements are covered
    n_blocks = (nu + chunk_sz) / (2*chunk_sz) + 1;
    if(n_blocks > nmax)
        n_blocks = nmax;
//     n_last_block = nr - n_blocks * 2*chunk_sz + chunk_sz;
}

void GPU_gcd::mod_reduce_kernel_dispatch(const unsigned *devU,
        unsigned *devR, const unsigned *devMods) {

    unsigned shm_size = (48) * word_size, chunks_per_cf = 1, n_thids;

    if(n_moduli < 64) // at least 32 threads per block
        n_thids = max(n_moduli, 32);
    else { // otherwise divide in 64-thread chunks
        chunks_per_cf = (n_moduli + 63) / 64;
        n_thids = 64;
    }

    unsigned n_coeffs = nu + nv + 2;
    // for the first time assume #mods = #threads
    dim3 grid(n_coeffs, chunks_per_cf);
    dim3 threads(n_thids, 1);

    CUMP_out2("\n##### mod_reduce kernel: chunks_per_cf: %d; nthids: %d\n",
            chunks_per_cf, n_thids);

    if(n_thids == 64) {
        mod_reduce_kernel< 64 ><<< grid, threads, shm_size >>>(
                devR, devU, devMods, limbs_f, limbs_g, n_moduli);
    } else { // otherwise template parameter does not matter
        mod_reduce_kernel< 32 ><<< grid, threads, shm_size >>>(
                devR, devU, devMods, limbs_f, limbs_g, n_moduli);
    }
}

//! dispatching \c QR_gcd_block_kernel calls
//! \c devU - data copied from CPU
//! \c devIn & \c devOut - memory spaces for "ping-ponging"
#if CUMP_USE_GCD_BLOCK_KERNEL
void GPU_gcd::QR_gcd_kernel_dispatch(unsigned *& mem_out,
        unsigned *devIn, unsigned *devOut, unsigned *devU,
            unsigned *page_locked_buf) {

#if CUMP_USE_PGCD_KERNEL
/** NOTE NOTE NOTE: make sure nu is not very large !! **/
/// NOTE note also that quad kernel loads data with a particular offset !!

//     unsigned *devU = (unsigned)const_mem[_DATA_IN],
//             *devR = (unsigned)const_mem[_DATA_OUT];

    if(0) {//nu < 1024 && nv < 256) {
        dim3 threads(max_nv, 1);
        dim3 grid(n_moduli, 1);
        unsigned shm_size = (nu + 8) * sizeof(unsigned);

        CUMP_out2("\n********** PGCD_lite_kernel: #thids: %d shm_size: "
            "%.2f (kb)\n", threads.x, ((float)shm_size/1024.0f));

        // QR_gcd_kernel_dispatch(mem_out, devR, devOuts, devU, R);
        /** DATA_IN = devU = devU
            DATA_OUT = devR = devIn*/
        
        PGCD_lite_kernel< false ><<< grid, threads, shm_size >>>(
                devU, devIn);

    } else {
        // NOTE: if the # of threads is not divisible by warp, it might that
        // our thread split does not work properly because warp will be
        // split in the middle
        unsigned nthids = ((nv + 4) / 4 + 31) & ~31, x2thids = 1, gridx;
    // NOTE NOTE NOTE: if you use this you also have to split shared
    // mem arrays between threads !!!
//         if(n_moduli != 1) {
//             if(nthids == 64) {
//                 x2thids = 2;
//             } else if(nthids == 32) {
//                 x2thids = 4;
//             }
//         }
        gridx = (n_moduli + x2thids-1) / x2thids;

        dim3 threads(nthids, x2thids);
        dim3 grid(gridx, 1);

        if(nthids >= 1024) {
            CUMP_out2("PGCD_quad_kernel: unsupported data size\n");
            CUMP_THROW(GPU_algorithm_exception)
        }
        
        const unsigned CacheLn = 32;
        unsigned shm_size = (nthids + CacheLn + 8) *
                    x2thids * sizeof(unsigned);

        CUMP_out2("\n********** PGCD_quad_kernel: gridx: %d; #thids: %dx%d shm_size: "
            "%.2f (kb)\n", gridx, threads.x, threads.y, ((float)shm_size/1024.0f));

//         printf("devR: %x\n", devIn);
            
        if(x2thids == 1)
            PGCD_quad_kernel< false ><<< grid, threads, shm_size >>>(
                    devU, devIn);
        else
            PGCD_quad_kernel< true ><<< grid, threads, shm_size >>>(
                    devU, devIn);
    }
    return;
#else
    if(nu + nv <= 512) {
        // block size is multiple of 64
        const unsigned block_sz = (nu + nv + 63) & ~63;

        dim3 threads(block_sz, 1);
        unsigned shm_size = (block_sz + 16) * sizeof(unsigned);
        dim3 grid(n_moduli, 1);

        CUMP_out2("********** QR_gcd_lite_kernel: nr: %d block_sz: %d\n",
            nu+nv, block_sz);

        QR_gcd_lite_kernel<<< grid, threads, shm_size >>>();
        return;
    }

    const unsigned block_sz = chunk_sz * 2;
    dim3 threads(block_sz, 1);
    // NOTE: can we reduce the memory usage somehow ??
    unsigned shm_size = (4*chunk_sz + 16) * sizeof(unsigned),
        n_used_streams = n_streams, mods_per_stream, mods_aux,
        stream_shift = 0, i = 0;

    if(n_moduli < n_used_streams)
        n_used_streams = n_moduli;
    mods_per_stream = n_moduli / n_used_streams; // # of mods per one stream
    mods_aux = n_moduli % n_used_streams;   // additional mods left

    CUMP_out2("QR_gcd_kernel: ChunkSz: %d; n_blocks: %d; "
            "block_io_sz: %d (32-bit words); n_streams: %d \n",
            chunk_sz, n_blocks, block_io_sz, n_used_streams)

        // # of blocks during ``lite'' iterations is given by 'nu' and it
    // does not change with iterations
    unsigned start_idx = 0, *p_in = devIn, *p_out = devOut,
        nb = n_blocks;

//!     QR_gcd_block_kernel< x, false > <<< grid, threads, shm_size, \
//!              g_streams[i] >>>(p_out, p_in, start_idx, 0, stream_shift); \

#define GCD_BLOCK_KERNEL_DISPATCH1(x) \
    case x: \
        QR_gcd_block_kernel< x, false > <<< grid, threads, shm_size >>> \
            (p_out, p_in, start_idx, 0, stream_shift); \
        break;
#define GCD_BLOCK_KERNEL_DISPATCH2(x) \
    case x: \
        QR_gcd_block_kernel< x, true ><<< grid, threads, shm_size >>> \
            (p_out, p_in, start_idx, shift, stream_shift); \
        break;       
#define GCD_LAST_BLOCK_DISPATCH(x) \
    case x: \
        QR_gcd_last_block< x ><<< grid, threads, shm_size >>> \
            (p_out, p_in, start_idx, stream_shift); \
        break;

    while(1) {

        stream_shift = 0;
        for(i = 0; i < n_used_streams; i++) {
        unsigned mps = mods_per_stream;
        if(i < mods_aux)
            mps++;
//         CUMP_out2("stream: %d; mps: %d\n", i, mps);
        dim3 grid(mps, nb);

        switch(chunk_sz) {
//         GCD_BLOCK_KERNEL_DISPATCH1(32)
//         GCD_BLOCK_KERNEL_DISPATCH1(64)
        GCD_BLOCK_KERNEL_DISPATCH1(128)
//         GCD_BLOCK_KERNEL_DISPATCH1(256)
        default:
            printf("unsupported chunk size: %d\n", chunk_sz);
            CUMP_THROW(GPU_algorithm_exception);
        }

        stream_shift += mps;
        } // for(i)

        //for(i = 0; i < n_used_streams; i++) {

//         cudaError_t err = cudaThreadSynchronize();
//         printf("%d stream: %s\n", 0, cudaGetErrorString(err));

/*        cudaMemcpy((void *)page_locked_buf, p_out, 16 * sizeof(unsigned),
                    cudaMemcpyDeviceToHost);

        unsigned j=0;
        while(1) {
            i = 0;
            cudaError_t err = cudaStreamQuery(g_streams[i]);
            printf("%d stream: %s\n", j++, cudaGetErrorString(err));

            if(err == cudaSuccess)
                break;
        }
        printf("=========================== device ready =================\n");*/
        //}

        start_idx += chunk_sz;
        std::swap(p_in, p_out);
//         if(i == DEBUG_ITER_IDX+1)
//             break;
        if(start_idx >= nv) // indicates the end of ``lite'' iterations
            break;
    }
    // NOTE NOTE: observe that p_in & p_out are already swapped here !!
    unsigned shift = (nv - (start_idx-chunk_sz)) % chunk_sz, N = 0;
//     CUMP_out2("shift = %d\n", shift);

//! unless \c shift is exactly 0, we do not count the full iteration because
//! it is split in 2 parts
    if(shift == 0) {
        start_idx += chunk_sz;
    }
    shift++; // this is to distinguish the case where shift == 0

#if 1
    unsigned iters = 0;
    while(1) {

        stream_shift = 0;
        for(i = 0; i < n_used_streams; i++) {
        unsigned mps = mods_per_stream;
        if(i < mods_aux)
            mps++;

//         CUMP_out2("Full start_idx: %d: N: %d; n_blocks: %d; "
//             " stream: %d; mps: %d\n", start_idx, N, nb, i, mps);
        dim3 grid(mps, nb);

        switch(chunk_sz) {
//         GCD_BLOCK_KERNEL_DISPATCH2(32)
//         GCD_BLOCK_KERNEL_DISPATCH2(64)
        GCD_BLOCK_KERNEL_DISPATCH2(128)
//         GCD_BLOCK_KERNEL_DISPATCH2(256)
        default:
            printf("unsupported chunk size: %d\n", chunk_sz);
            CUMP_THROW(GPU_algorithm_exception);
        }

        stream_shift += mps;
        } // for(i)

        // shift == 0 indicates that we are running not the first fill iter
        N = (nu + nv - start_idx);
        start_idx += chunk_sz, shift = 0;

        iters++;
        if((iters % 3) == 0) {
            cudaMemcpy((void *)page_locked_buf, p_out, 16 * sizeof(unsigned),
                    cudaMemcpyDeviceToHost);
            if(page_locked_buf[0] == CUMP_GCD_TERMINATE_FLAG) {
                mem_out = p_out;
                printf("EVENT FIRED in step: %d\n", iters);
                return;
            }
        }

        if(N <= 3*chunk_sz)
            break;
        nb = (N + chunk_sz - 1) / (2*chunk_sz);
        std::swap(p_in, p_out);
    }

    CUMP_out2("last iteration: N = %d\n", N);
    std::swap(p_in, p_out);

    stream_shift = 0;
    for(i = 0; i < n_used_streams; i++) {
        unsigned mps = mods_per_stream;
        if(i < mods_aux)
            mps++;
        dim3 grid(mps, 1);

        switch(chunk_sz) {
//         GCD_LAST_BLOCK_DISPATCH(32)
//         GCD_LAST_BLOCK_DISPATCH(64)
        GCD_LAST_BLOCK_DISPATCH(128)
//         GCD_LAST_BLOCK_DISPATCH(256)
        default:
            printf("unsupported chunk size: %d\n", chunk_sz);
            CUMP_THROW(GPU_algorithm_exception);
        }
        stream_shift += mps;
    } // for(i)
#else
    std::swap(p_in, p_out);
#endif
    // NOTE: here we should be very carefull about the grid size because
    // memory is allocated accounting for the ``lite'' iterations..
    mem_out = p_out;
#undef GCD_BLOCK_KERNEL_DISPATCH1    
#undef GCD_BLOCK_KERNEL_DISPATCH2
#undef GCD_LAST_BLOCK_DISPATCH

#endif  // CUMP_USE_PGCD_KERNEL
}

#else // !CUMP_USE_GCD_BLOCK_KERNEL
void GPU_gcd::QR_gcd_kernel_dispatch(unsigned *& mem_out,
       unsigned *devR, unsigned *devOuts, unsigned *devU) {

    const unsigned block_sz = chunk_sz * 2;

    CUMP_out2("QR_gcd_kernel launch: ChunkSz: %d; n_blocks: %d; n_last: %d\n",
            chunk_sz, n_blocks, n_last_block)

    dim3 threads(block_sz, 1);
    dim3 grid(n_moduli, n_blocks);
    unsigned shm_size = (4*chunk_sz + 16) * sizeof(unsigned);

    if(chunk_sz == 32) {
        QR_gcd_monolithic_kernel< 32 ><<< grid, threads, shm_size >>>(
            devR, devOuts, devU);
    } else
        throw "NYI";
    
    mem_out = devR;
}
#endif // CUMP_USE_GCD_BLOCK_KERNEL

void GPU_gcd::MRC_kernel_dispatch(unsigned *devOut, const unsigned *devIn,
     const unsigned *devMods, const unsigned *devInvDets) {

#if !CUMP_RUN_GCD_ONLY
    dim3 grid, threads;
    unsigned shm_size;

    if(gcd_real_sz == 0 || gcd_real_sz > nv + 1) {
        printf("Incorrect gcd size: %d\n", gcd_real_sz);
        CUMP_THROW(GPU_algorithm_exception);
    }

    if(n_moduli <= 64) {
        shm_size = 16 * sizeof(unsigned);
        threads = dim3(32, 2, 1);
        grid = dim3((gcd_real_sz + 1) / 2, 1, 1);
        CUMP_out2("**** invoking 32x2 MRC kernel ****\n");

        if(n_moduli & 1) {
            MRC_doubles_kernel< 32, false ><<< grid, threads, shm_size >>>
                 (devOut, devIn, devMods, devInvDets, mods_padding,
                         n_moduli, gcd_sz);
        } else {
            MRC_doubles_kernel< 32, true ><<< grid, threads, shm_size >>>
                 (devOut, devIn, devMods, devInvDets, mods_padding,
                     n_moduli, gcd_sz);
        }
        return;
    }

#define MRC_KERNEL_DISPATCH(x) \
  MRC_block_kernel< x ><<< grid, threads, shm_size >>>(devOut, devIn, \
        devMods, devInvDets, mods_padding, n_moduli, gcd_sz);
    
    grid = dim3(gcd_real_sz, 1, 1);
    shm_size = 16 * sizeof(unsigned);

    if(n_moduli <= 128) {
        threads = dim3(64, 1, 1);
        shm_size += 64*2 * sizeof(unsigned);
        MRC_KERNEL_DISPATCH(2)
        return;
    }
    // process the main cases when n_moduli > 128:
    // compute optimal # of chunks, so that the last chunk has 
    // the minimal number of idle threads with the priority given to
    // 4 chunks (because of smaller block size)
    int CHZ3 = ((n_moduli + 2) / 3 + 63) & ~63,
        nlast3 = (int)n_moduli - (3 - 1) * CHZ3,
        CHZ4 = ((n_moduli + 3) / 4 + 63) & ~63,
        nlast4 = (int)n_moduli - (4 - 1) * CHZ4;

    if(nlast4 > CHZ4)
        nlast4 = -1;
    CUMP_out2("CHZ3: %d nlast3: %d; CHZ4: %d nlast4: %d\n",
            CHZ3, nlast3, CHZ4, nlast4);
            
    unsigned n_thids, chz;
    if(CHZ3 > 512 || nlast3 < 0 || nlast4 > nlast3) {
        if(CHZ4 > 512) { // in case # of moduli is too large -> quit
            int CHZ5 = ((n_moduli + 4) / 5 + 63) & ~63;
            if(CHZ5 > 512) {
                printf("FATAL: unsupported # of mods: %d\n", n_moduli);
                CUMP_THROW(GPU_algorithm_exception)
            }
            n_thids = CHZ5, chz = 5;
        } else {
            n_thids = CHZ4, chz = 4;
            if(nlast4 < 0 || nlast4 > CHZ4) {
                printf("FATAL: nlast4: %d\n", nlast4);
                CUMP_THROW(GPU_algorithm_exception)
            }
        }
    } else {
        n_thids = CHZ3, chz = 3;
        if(nlast3 < 0 || nlast3 > CHZ3) {
            printf("FATAL: nlast3: %d\n", nlast3);
            CUMP_THROW(GPU_algorithm_exception)
        }
    }

/*      [65..128] = 64 x 2 (CHUNK_SZ x N_CHUNKS)
        [129..192] = 64 x 3
        [193..256] = 64 x 4
        [257..320] = 64 x 5
        [321..384] = 128 x 3
        [385..512] = 128 x 4
        [513..768] = 256 x 3
        [769..1024] = 256 x 4
        [1025..1152] = 384 x 3
        [1153..1536] = 384 x 4
        [1537..1792] = 448 x 4
        [1792..2048] = 512 x 4*/

    CUMP_out2("########### MRC_kernel: n_moduli: %d; n_thids: %d; CHZ: %d\n",
            n_moduli, n_thids, chz);
   
    threads = dim3(n_thids, 1, 1);
    shm_size += n_thids*2 * sizeof(unsigned);
    switch(chz) {
    case 3: MRC_KERNEL_DISPATCH(3)
        break;
    case 4: MRC_KERNEL_DISPATCH(4)
        break;
    case 5: MRC_KERNEL_DISPATCH(5)
        break;
    default:
        printf("MRC_block_kernel: chz: %d\n", chz);
        CUMP_THROW(GPU_algorithm_exception);      
    }
//     CRA_quads_kernel< 64, 1 ><<< grid, threads, shm_size >>>(devOut, devIn,
//         devMods, devInvDets, mods_padding, n_moduli, gcd_sz);

#endif // !CUMP_RUN_GCD_ONLY
}

//! \c devGCDLcf - a set of residues of gcd(lc(f),lc(g))
void GPU_gcd::mod_inverse_kernel1_dispatch(unsigned *devOut,
    const unsigned *devIn, const unsigned *devMods,
            const unsigned *devGCDLcf) {

#if !CUMP_RUN_GCD_ONLY
    const unsigned ChunkSz = 64;
    unsigned n = (n_moduli + ChunkSz-1) / ChunkSz, 
        shm_size = 16 * sizeof(unsigned);

    dim3 threads(ChunkSz, 1, 1);
    dim3 grid(n, 1, 1);

    mod_inverse_kernel1< ChunkSz ><<< grid, threads, shm_size >>>(
        devOut, devIn, devMods, devGCDLcf, mods_padding, n_moduli, gcd_sz);
#endif
}

void GPU_gcd::free_device_mem() {

    CUMP_out("GPU_gcd: free_device_mem..\n");

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
