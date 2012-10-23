// ============================================================================
//
// Copyright (c) 2001-2009 Max-Planck-Institut Saarbruecken (Germany).
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

#ifndef _INTERPOLATE_CRA_KERNEL_CU_
#define _INTERPOLATE_CRA_KERNEL_CU_

#define RUN_INTERPOLATE_DOUBLES 0
#define RUN_INTERPOLATE_QUADS   1
#define RUN_CRA_KERNEL          1
#define RUN_TEST_RED_KERNEL     0

#include <include/macros.h>
#include <zmod_dev.cu>

//! computes the product of 2^LDN elements (min LDN = 5) stored in shared mem
//! \c data , number of working threads: 2^(LDN-1)
//! NOTE: \c data must be padded by ones
//! the result is returned by \c 0th thid !!!
//! \c ThidCheck : whether to constrain the number of working thids
//! (set to \c false if several prefix sums are run independently for
//! different warp groups)
//! \c postscan - memory storage for postscan
template < int N, bool ThidCheck >
__device__ __forceinline__ unsigned prefix_mul(unsigned *data, unsigned *postscan, volatile unsigned m,
            volatile fp_limb invk) {

    unsigned thid = threadIdx.x;
    volatile unsigned *t = (volatile unsigned *)data + UMUL(thid >> 5,
            UMUL_PAD + (1 << 6)) + (thid & WS-1);

    if(thid < 32) // maximum 32 postscan elements: 32*WS = 2048
        postscan[thid] = 1; // fill with identities

    unsigned a0 = 1; // this is for data padding
    if(ThidCheck && thid < N / 2) {
 
// TODO: the first reduction is only required for even # of warps
// the last (odd) warp does not perform it
        if(N / 2 >= WS && (((N / 2 & 1) == 0 || (int)thid < N / 2 - WS))) {
            a0 = mul_m(t[0], t[WS], m, invk);
            t[0] = a0;
        } else {
            a0 = t[0];
        }
        a0 = mul_m(a0, t[HF], m, invk);
        t[0] = a0;
        a0 = mul_m(a0, t[8], m, invk);
        t[0] = a0;
        a0 = mul_m(a0, t[4], m, invk);
        t[0] = a0;
        a0 = mul_m(a0, t[2], m, invk);
        t[0] = a0;
        a0 = mul_m(a0, t[1], m, invk);
        t[0] = a0;
    }

    // now we need to multiply N / 64 more elements
    int nElems = (N + 63) / 64, nElemsPadded = nElems;
    // N <= 64: nothing to be done
    // N <= 128: use 1 thread
    // N <= 256: use 2 threads
    // N <= 512: use 4 threads ..

    if((nElems & nElems - 1) != 0) { // compile-time evaluation
        if(nElems < 4)
            nElemsPadded = 4;
        else if(nElems < 8)
            nElemsPadded = 8;
        else if(nElems < 16)
            nElemsPadded = 16;
        else if(nElems < 32)
            nElemsPadded = 32;
    }

    volatile unsigned *t2 = postscan;
    if(nElems / 2 >= 1) {

    // relocate data to avoid bank conflicts during the second reduction
        if((thid & WS-1) == 0) { // every first thid in a warp
  // NOTE NOTE: read-write hazard.. dangerous ??
            unsigned id = thid >> 5;
            t2[id] = a0; // at most 16 thids
        }
        CU_SYNC // we cannot sync inside a data-dependent branch
    }

    if((int)thid < nElems / 2) {

        // each thread will handle 2 elements with WS*2 stride
        t2 += thid;

        nElemsPadded /= 2; // nElems is at least 1 here
        a0 = t2[0]; // nElems 64 -> first stride
        // here nElems = 4
        a0 = mul_m(a0, t2[nElemsPadded], m, invk);

        nElemsPadded /= 2;// == 2
        if(nElemsPadded >= 1) {
            t2[0] = a0; 
            a0 = mul_m(a0, t2[nElemsPadded], m, invk);
        }

        nElemsPadded /= 2; // == 1
        if(nElemsPadded >= 1) { // strides decrease by 2 each step
            t2[0] = a0; 
            a0 = mul_m(a0, t2[nElemsPadded], m, invk);
        }

        nElemsPadded /= 2;
        if(nElemsPadded >= 1) { // strides decrease by 2 each step
            t2[0] = a0; 
            a0 = mul_m(a0, t2[nElemsPadded], m, invk);
        }

        // up to nElems == 32, ie., N = 32*64 = 2048
        nElemsPadded /= 2;
        if(nElemsPadded >= 1) { // strides decrease by 2 each step
            t2[0] = a0; 
            a0 = mul_m(a0, t2[nElemsPadded], m, invk);
        }
    }

    return a0;
}

//! Vandermonde interpolation using 256 threads, i.e. 21 regs per thread
//! to run 3 blocks per SM
//! \c n - number of interpolation points: 256 < n <= 512
//! \c g_U layout: \c n points

//! NOTE NOTE NOTE: interpolation points g_U must be in reversed order now..
template < unsigned LDN >
__global__ void interpolate_doubles_kernel(unsigned *g_R, const unsigned *g_U,
        unsigned n) {

#if RUN_INTERPOLATE_DOUBLES
    extern __shared__ unsigned shared[];

    unsigned thid = threadIdx.x, bidx_x = blockIdx.x;//, bidx_y = blockIdx.y;

    unsigned *r = shared;
    const unsigned m_stride = 4; // constant memory stride
    const unsigned N = 1 << LDN;
    unsigned *mods = dev_const_mem;

    //! remember 2 is here: this is to match memory layout for resultant kernel
    mods += 2 + UMUL(bidx_x, UMUL_PAD + m_stride);
    volatile unsigned m = mods[0];
    volatile fp_limb invk, invm;

#if CUMP_USE_32BIT_MODULI_SET
    invk = __hiloint2double(mods[3], mods[2]);
#else
    invk = __int_as_float(mods[2]);
    invm = __int_as_float(mods[3]);
#endif

    volatile unsigned mx100 = UMUL(100, m);
    // fool the compiler: 24-bit mul masks out the highest bit
    //volatile unsigned _2e23 = 0x4b000000 + UMUL(thid, UMUL_PAD);

    // TODO: must be saved with padding even if n < N !!
    unsigned ofs = UMUL(bidx_x, UMUL_PAD + N);
    const unsigned *in_U = g_U + ofs;

    unsigned n_thids = (n + 1) / 2; // # of working thids (+1 to handle odd
                                    // numbers)
    uint2 a, b, w;
    if(thid < n_thids) // read out only relevant data
        b = ((const uint2 *)in_U)[thid];

    // NOTE As & Bs have per-thread indexing !!!!
    const unsigned BlockSz = N / 2; // 2 elements processed by single thread
    unsigned *Bs = r + thid, *As = Bs + BlockSz,
        *xs = r + BlockSz*2; // xs is not per-thread indexing !!!

    //!//////////////////////////////////////////////////////
    //! for the first time xs[i] = i + 1
    volatile int j = n - 1;
    unsigned a0 = 1, a_prev;

    // xs[thid] = thid + 1;
    Bs[0] = b.x; // shift down all elements, s.t. Bs[0] becomes b0

    CU_SYNC

    b.x = b.y, b.y = 0;

    a.x = thid + 1;
    if(a.x < n_thids) { // what to do in case of odd n ??
                        // force even number of evaluation points always ??
        b.y = Bs[1]; // except the last element
    }
    
    a.x = 1, a.y = a.x;
    w.x = a.x, w.y = a.x; // fill with identities to make sure reduction
            // works correctly even with less # of active threads

    while(1) {

    if(thid == n_thids - 1) {
        // so we have 0 * a0 - 1 * b0 = -b0
        if(n & 1) { // handle odd # of interpolation points
            a.x = 1, b.x = m;
        } else {
            a.y = 1, b.y = m;
        }
    }

    //! shift all index ranges by -(j+1)
    // r[0] == b0
    unsigned x0 = r[0];

    b.x = sub_mul_reduce_m(b.x, a0, a.x, x0, m, invk, invm, mx100);
    b.y = sub_mul_reduce_m(b.y, a0, a.y, x0, m, invk, invm, mx100);

    //! this just defines the order of evaluation points - can be changed later
    if(j == thid)
        w.x = a0;
    else if(j == thid + BlockSz)
        w.y = a0;

    // b0 not needed anymore
    unsigned b0 = UMUL(thid, UMUL_PAD + 2) + 1;
    a0 = j - b0; // == n_m_j - (thid * 2 + 1)
    //! this just defines the order of evaluation points - can be changed later
    b0 += n - j - 1;

    if(thid == (j >> 1)) {
        a_prev = m;
    }

    x0 = j + 1; // nn - j // xs[j]
    unsigned t1 = m, mul1 = x0, t2 = m, mul2 = x0;

    if(j == 0) // 0
        break;
    j--;

    // that is, if interpolation points are decreasing -> you can
    // benefit from that..=> use mul_small_m
    if((int)a0 >= 0) {
        mul1 -= (n - b0);//__usad(m, b0 + 1, x0);// m + x0 - (b0 + 1);
    }

    if((int)a0 > 0)
        mul2 -= (n - 1 - b0);//__usad(m, b0 + 2, x0);//m + x0 - (b0 + 2);

    if((int)a0 == 0) {
        t2 = a_prev;
    }

    if((int)a0 < 0) {
        t1 = a_prev;
        t2 = a.x;
    }

    a.x = submul_m(t1, mul1, a.x, m, invk);
//     a.x = mul_small_m(mul1, a.x, m, invk);
//     a.x = sub_m(t1, a.x, m);

    a.y = submul_m(t2, mul2, a.y, m, invk);
//     a.y = mul_small_m(mul2, a.y, m, invk);
//     a.y = sub_m(t2, a.y, m);

    //! imagine you have to perform all these instructions at *each* iteration
    //! does it cost one mul_m ?? (taking into account that you will have
    //! to execute reduction with 256 threads at the end..

    // NOTE: possible read-write hazard on the next iteration
    As[0] = a.x, Bs[0] = b.x, a_prev = a.x;

    CU_SYNC

    a.x = a.y, b.x = b.y;
    if(thid + 1 < n_thids) { // TODO: can merge these two conditions into one
        a.y = As[1]; // read out elements shifted by one (except the last
        b.y = Bs[1]; // thid)
    }

    a0 = r[BlockSz]; // == As[-thid]

    } // while(1)

    w.x = mul_m(w.x, w.y, m, invk); // multiply per-thid results
    xs[thid] = w.x;
    CU_SYNC

    // LDN-1 because we already premultiplied by elements pairwise
    // 
    a0 = prefix_mul< N / 2, true > (xs, r, m, invk);

//     if(thid == 0) {
//         g_R[thid] = a0;
//     }
//     return;

    if(thid < n_thids) {
        ((uint2 *)(g_R + ofs))[thid] = b;
    }

#else
#warning interpolate_doubles_kernel: dummy compilation
#endif // RUN_INTERPOLATE_DOUBLES
}

//! Vandermonde interpolation inner loop using 4 elements per thread
//! \c b - interpolation data (gets overwritten by the results)
//! \c w - returns denominator factors to be multiplied
//! \c n - number of interpolation points
//! \c n_thids - number of active thids
//! \c r - shared memory working space (3 * \c block_sz words)
//! \c Mod4 - parametrizes kernel by (n mod 4) where \c n is the number
//! of evaluation points
template < unsigned Mod4 >
__device__ __forceinline__ unsigned __interpolate_quads_internal(uint4& b, unsigned *r,
        const unsigned *g_Xs, unsigned n, unsigned n_thids,
        unsigned block_sz, volatile unsigned m,
            volatile fp_limb invk, volatile fp_limb invm) {

    unsigned thid = threadIdx.x;
    //! NOTE As, Bs and xs have per-thread indexin !!
    unsigned *Bs = r + thid, *As = Bs + block_sz, *xs = As + block_sz;

    volatile int j = n - 1;
    unsigned a0 = 1, a_prev;

    uint4 a, w;
    if(thid < n_thids) {
        w = ((uint4 *)g_Xs)[thid];
// better to wait with saving to shared mem, so that we can overlap mem
// access with ALUs..
    }

    // shift the vectors down
    xs[0] = w.x, w.x = w.y, w.y = w.z, w.z = w.w, w.w = 0;
    Bs[0] = b.x, b.x = b.y, b.y = b.z, b.z = b.w, b.w = 0;
    a.x = thid + 1;

// NOTE: it's better to place this sync as futher as possible
    CU_SYNC
    if(a.x < n_thids) { // what to do in case of odd n ??
                        // force even number of evaluation points always ??
        b.w = Bs[1]; // except the last element

        w.w = xs[1]; // read "the tail" of the next thread
    }

    a.x = 1, a.y = a.x, a.z = a.x, a.w = a.x;
    unsigned det = 1; // this will store the denominator

    while(1) {

    unsigned x0;

    if(thid == n_thids - 1) {
        b.w = det; // last thid computes denominator
        a.w = 0;
    }

    //! shift all index ranges by -(j+1)
    x0 = r[0]; // r[0] == b0

    unsigned mx100 = UMUL(100, m);
    b.x = sub_mul_reduce_m(b.x, a0, a.x, x0, m, invk, invm, mx100);
    b.y = sub_mul_reduce_m(b.y, a0, a.y, x0, m, invk, invm, mx100);
    b.z = sub_mul_reduce_m(b.z, a0, a.z, x0, m, invk, invm, mx100);
    b.w = sub_mul_reduce_m(b.w, a0, a.w, x0, m, invk, invm, mx100);

    if(thid == n_thids - 1) {
        det = b.w;

        if(Mod4 == 1) { // this is compile-time branch
// NOTE: it might happen that x0 == 0 then b == m: no explicit reduction is
// required because b will be reduced later on in subsequent operations
            a.x = 1, b.x = m - x0;
        } else if(Mod4 == 2) {
            a.y = 1, b.y = m - x0;
        } else if(Mod4 == 3) {
            a.z = 1, b.z = m - x0;
        } else {
            a.w = 1, b.w = m - x0; // so we have m * a0 - 1 * b0 = -b0
        }
    }

    //! this just defines the order of evaluation points - can be changed later
    x0 =  //j+1;
        r[block_sz*2]; // == xs[-thid]

    // b0 not needed anymore
    a0 = j - UMUL(thid, UMUL_PAD + 4) - 1;
//     unsigned b0 = n - 1 - a0;

    if(thid == (j >> 2)) {
        a_prev = m;
    }

    if(j == 0)
        break;
    j--;

    CU_SYNC

    Bs[0] = b.x; // b.x is free now

/**
 if((int)a0 > 2) { 4 first

    } else if((int)a0 == 2) { // 3 first, 1 second

    } else if((int)a0 == 1) { // 2 first, 2 second

    } else if((int)a0 == 0) { // 1 first, 3 second

    a0 < 0: 4 second
*/

//! do ping-ponging registers: mul1 <--> t1
    unsigned t1 = m, mul1 = x0;

    if((int)a0 >= 0) {
        // TODO: any way to replace this by __sad ??
        // mul[i] is either x0 or x0 - x[i]
        mul1 -= w.x; 
    }
    //! inlined "sub_m" to use predicate exec instead of branching
    mul1 = umin(mul1, mul1 + m);

    if((int)a0 < 0) {
        t1 = a_prev;
    }

    b.x = a.x;
    a.x = submul_m(t1, mul1, a.x, m, invk);

    mul1 = m; // mul1 is now t2

    if((int)a0 < 0) {
        mul1 = b.x;
    }

    if((int)a0 == 0) { // why it works without this as well ?
        mul1 = a_prev;
    }

    t1 = x0; // t1 is now mul2, b.x is free here
    if((int)a0 >= 1) {
        t1 -= w.y; // (n - 1 - b0);
    }
    t1 = umin(t1, t1 + m);

    b.x = a.y;
    a.y = submul_m(mul1, t1, a.y, m, invk);

    t1 = m;
    if((int)a0 < 1) {
        t1 = b.x; // t1 is t3 now
    }

    mul1 = x0; // mul1 is mul3 again
    if((int)a0 >= 2) {
        mul1 -= w.z;//(n - 2 - b0);
//         mul1 = sub_m(mul1, w.z, m);
    }
    mul1 = umin(mul1, mul1 + m);

    if((int)a0 == 1) {
        t1 = a_prev;
    }

    b.x = a.z;
    a.z = submul_m(t1, mul1, a.z, m, invk);

    mul1 = m; // mul1 is t4
    if((int)a0 == 2) {
        mul1 = a_prev;
    }

    if((int)a0 < 2) {
        mul1 = b.x;
    }    

    t1 = x0; // t1 is mul4
    if((int)a0 >= 3) {
        t1 -= w.w;
//         t1 = sub_m(t1, w.w, m);
    }
    t1 = umin(t1, t1 + m);
    a.w = submul_m(mul1, t1, a.w, m, invk);

    //! imagine you have to perform all these instructions at *each* iteration
    //! does it cost one mul_m ?? (taking into account that you will have
    //! to execute reduction with 256 threads at the end)..

    // NOTE: possible read-write hazard on the next iteration
    As[0] = a.x, a_prev = a.x;

// NOTE: take into account that the # of xs in use decreases with each
// iteration
    xs[0] = w.x;

    CU_SYNC

    // shift down the variables
    a.x = a.y, b.x = b.y, a.y = a.z;
    b.y = b.z; a.z = a.w, b.z = b.w;

    w.x = w.y, w.y = w.z, w.z = w.w;

        //TODO: can use a0 & b0 here for temporary variables..
    if(thid + 1 < n_thids) { // TODO: can merge these two conditions into one
        a.w = As[1]; // read out elements shifted by one (except the last
        b.w = Bs[1]; // thid)
        w.w = xs[1];
    }

//     CU_SYNC

    a0 = r[block_sz]; // == As[-thid]

    } // while(1)

    return det;
}

// interpolation kernel processing 4 elements per thread
//! \c Mod4 = n mod 4
template < unsigned Mod4 > 
__global__ void CUMP_LAUNCH_BOUNDS(1024, 1)
interpolate_quads_kernel(unsigned *g_R, const unsigned *g_Xs,
        unsigned *g_InvDets, unsigned n, const unsigned out_data_padding) {

#if RUN_INTERPOLATE_QUADS
    extern __shared__ unsigned shared[];

    unsigned thid = threadIdx.x, bidx_x = blockIdx.x;//, bidx_y = blockIdx.y;
//! bidx_x - indexes over moduli

    unsigned *r = shared;
    //const unsigned N = 1 << LDN;
    const unsigned m_stride = 4; // constant memory stride
    unsigned *mods = dev_const_mem;// + UMUL(bidx_x, UMUL_PAD + m_stride);

    //! remember 2 is here !!!!!
    mods += 2 + UMUL(bidx_x, UMUL_PAD + m_stride);
    volatile unsigned m = mods[0];
    volatile fp_limb invk, invm;

#if CUMP_USE_32BIT_MODULI_SET
    invk = __hiloint2double(mods[3], mods[2]);
#else
    invk = __int_as_float(mods[2]);
    invm = __int_as_float(mods[3]);
#endif

    unsigned ofs = UMUL(bidx_x, out_data_padding);
    unsigned n_thids = (n + 3) / 4; // # of working thids
    //! offsets added! do not forget
    g_R += ofs, g_Xs += ofs;

    uint4 b;
    if(thid < n_thids) // read out only relevant data
        b = ((uint4 *)g_R)[thid];

//     uint4 w;
//     if(thid < n_thids) {
//         w = ((uint4 *)g_Xs)[thid];
//     }
//     CU_SYNC

    // 4 elements per thread => N/4 max thids
    unsigned block_sz = out_data_padding / 4;
    unsigned det = __interpolate_quads_internal< Mod4 >(b, r, g_Xs, n, n_thids,
        block_sz, m, invk, invm);

    if(thid < n_thids) {
        ((uint4 *)g_R)[thid] = b;
    }

    // number of elements: BlockSz * 4
    // LDN-1 because we have performed one mul already
    if(thid == n_thids - 1) {
        g_InvDets[bidx_x] = det;
    }

#else
#warning interpolate_quads_kernel: dummy compilation
#endif // RUN_INTERPOLATE_QUADS
}

//! 64-thread kernel computing modular inverses
__global__ void mod_inverse_kernel2(unsigned *g_R, const unsigned *g_U,
        const unsigned *g_Mods, unsigned mods_padding, unsigned n_mods) {

     extern __shared__ unsigned shared[];
//! number of blocks: nmods / 64
    unsigned thid = threadIdx.x, bidx_x = blockIdx.x;
//! bidx_x: indexes over moduli
    //! again offset-dependent reading in kernel is bad
    unsigned ofs = bidx_x * 64;
    //! (m, mu, inv_mod, invk)
    const unsigned *g_Mus = g_Mods + mods_padding;

    if(ofs + thid >= n_mods) // to protect from "garbage" inputs
        return;

    unsigned m = g_Mods[thid + ofs], mu = g_Mus[thid + ofs],
            inva = g_U[thid + ofs];
    // NOTE: beware if inva == 0 the algorithm does not hang up
    inva = montgomery_inverse(inva, m, mu);

// save with the same offset: available for in-place mod
    g_R[thid + ofs] = m - inva; // store -inva^-1
}

template < bool IsEven >
__device__ __forceinline__ void __CRA_shmem_update2(uint2 w, uint2 m, unsigned *r,
            volatile int idx) {

    if(IsEven) { // compile-time decision
        if(idx == 0) {
            r[0] = w.y, r[1] = m.y; 
        }
    } else { 
        if(idx == 1) {
            r[0] = w.x, r[1] = m.x;
        }
    }
}

//! single iteration of the inner CRT loop parametrized by (n_mods mod 2)
//! \c r_in and \c r_out - shared mem entries for "ping-ponging"
template < bool IsEven >
__device__ __forceinline__ void __CRA_doubles_internal(uint2& w, uint2& M, uint2 m,
        fp_limb2 invk, const unsigned *r_in, unsigned *r_out,
                volatile int idx) {

    if(IsEven) { // compile-time decision
        if(idx <= 0) {
            w.x = submul_m(w.x, r_in[0], M.x, m.x, invk.x);
            M.x = mul_m(M.x, r_in[1], m.x, invk.x); 
            w.y = submul_m(w.y, r_in[0], M.y, m.y, invk.y);
            M.y = mul_m(M.y, r_in[1], m.y, invk.y); 
        }
    } else {
        if(idx <= 1) {
            w.x = submul_m(w.x, r_in[0], M.x, m.x, invk.x);
            M.x = mul_m(M.x, r_in[1], m.x, invk.x); 
        }
        if(idx <= 0) {
            w.y = submul_m(w.y, r_in[0], M.y, m.y, invk.y);
            M.y = mul_m(M.y, r_in[1], m.y, invk.y); 
        }   
    }
    // share updated elements between all working threads
    __CRA_shmem_update2< IsEven >(w, m, r_out, idx);
}

//! \c devR - out data (Mixed-Radix digits)
//! \c devU - in data (residues)
//! \c mods_padding - moduli padding
//! \c g_Mods layout: (m, invk, mu, inv_mod) stored with \c mods_padding
//! \c n_mods - # of moduli: determines # of working thids

//! \c BlockSz works as a switch: BlockSz == WS: process two batches
//! at a time (split up by the 2nd block dimension)
//! \c IsEven indicates whether \c n_mods is even
template < unsigned BlockSz, bool IsEven >
__global__ void CRA_doubles_kernel(unsigned *g_R, const unsigned *g_U,
    const unsigned *g_Mods, const unsigned *g_InvDets, unsigned mods_padding,
            unsigned n_mods, unsigned out_data_padding) {

#if RUN_CRA_KERNEL
    extern __shared__ unsigned shared[];
//! blockIdx.x - indexes over batches
    unsigned thid = threadIdx.x, bidx_x = blockIdx.x,
        nthids = (n_mods + 1) / 2; // # of working thids

    unsigned split = (BlockSz == WS ? threadIdx.y : 0);
    unsigned *r = shared + split * 8;
    //! possibly in texture mem ??
    //! moduli in descreasing order: m[0] > m[1] > ... > m[n-1]
    uint2 m;
    fp_limb2 invk;
    m = ((uint2 *)g_Mods)[thid];

    //! (m, mu, inv_mod, invk)
    //! bad that reading is offset-dependent (could be hard to debug this..)
    const unsigned *g_InvMods = g_Mods + mods_padding * 2,
                   *g_InvKs = g_InvMods + mods_padding;

    if(thid < nthids) {
#if CUMP_USE_32BIT_MODULI_SET
        uint2 lo = ((uint2 *)g_InvKs)[thid],
              hi = ((uint2 *)(g_InvKs + mods_padding))[thid];

        invk.x = __hiloint2double(hi.x, lo.x);
        invk.y = __hiloint2double(hi.y, lo.y);
#else
        invk = ((fp_limb2 *)g_InvKs)[thid];
#endif
    }

    //! read is uncoalesced, write back must be coalesced
    unsigned ofs;
    if(BlockSz == WS)
        ofs = UMUL(UMUL(thid, out_data_padding) + bidx_x,
                UMUL_PAD + 2) + split;
    else
        ofs = UMUL(thid, out_data_padding*2) + bidx_x;

    uint2 inv_mod, w, M;
    if(thid < nthids) {
        inv_mod = ((uint2 *)g_InvMods)[thid];
        w.x = g_U[ofs];     // *awful* uncoalesced mem access..grrr
        w.y = g_U[ofs + out_data_padding];
    }

    if(BlockSz == WS)
        ofs = UMUL(UMUL(bidx_x, UMUL_PAD + 2) + split, mods_padding) +
                thid*2;
    else
        ofs = UMUL(bidx_x, mods_padding) + thid*2;

    uint2 det = ((uint2 *)g_InvDets)[thid];
    // # of denominators is moduli-bound
    w.x = mul_m(w.x, det.x, m.x, invk.x);
    w.y = mul_m(w.y, det.y, m.y, invk.y);

    // inv_mod[n - 1] = 1
    // inv_mod[n - 2] = m[n - 1]^-1 mod m[n - 2]
    // inv_mod[n - 3] = m[n - 1]m[n - 2]^-1 mod m[n - 3] ..
// NOTE: in fact we can handle n + 1 residues with n threads
    volatile int idx = UMUL(thid, UMUL_PAD + 2) + 1 - (n_mods - 1);
    __CRA_shmem_update2< IsEven >(w, m, r, idx);

    if(BlockSz > WS)
        CU_SYNC 

    idx++;
    if(idx <= 1) {
        w.x = sub_m(w.x, r[0], m.x);
        w.x = mul_m(w.x, inv_mod.x, m.x, invk.x);
        M.x = mul_m(r[1], inv_mod.x, m.x, invk.x);
    }
    if(n_mods <= 2)
        goto Lexit;

    if(idx <= 0) { 
        w.y = sub_m(w.y, r[0], m.y);
        w.y = mul_m(w.y, inv_mod.y, m.y, invk.y);
        M.y = mul_m(r[1], inv_mod.y, m.y, invk.y);
    }

    __CRA_shmem_update2< !IsEven >(w, m, r + 2, idx);

    idx++;
    {
    volatile int i = n_mods - 3;
    while(1) { // main loop unrolled by the factor of two

        if(BlockSz > WS)
            CU_SYNC
    
        // ping-ponging shared mem access in each iteration to save on syncs
        __CRA_doubles_internal< IsEven >(w, M, m, invk, r + 2, r, idx);

        idx++;
        if(BlockSz > WS)
            CU_SYNC

        __CRA_doubles_internal< !IsEven >(w, M, m, invk, r, r + 2, idx);
        
         if(i <= 1)
             break;
        i -= 2; idx++;
    } // while
    }
    if(BlockSz > WS)
        CU_SYNC

    //w.x*m.y + w.y (w.x is the highest MR digit)
//    uint2 u;
//     mad_24_wide((unsigned&)u.y, (unsigned&)u.x, w.x, m.y, w.y);
//     u.y >>= 16, u.y += (u.x < w.y);

Lexit:
    if(thid < nthids) {
        ((uint2 *)(g_R + ofs))[0] = w;
    }
#endif // RUN_CRA_KERNEL
}

template < unsigned Mod4 >
__device__ __forceinline__ void __CRA_shmem_update4(uint4 w, uint4 m, unsigned *r,
            volatile int idx) {

    if(Mod4 == 0) { // compile-time decision
        if(idx == 0) {
            r[0] = w.w, r[1] = m.w;
        }
    } else if(Mod4 == 3) {
        if(idx == 1) {
            r[0] = w.z, r[1] = m.z;
        }
    } else if(Mod4 == 2) {
        if(idx == 2) {
            r[0] = w.y, r[1] = m.y;
        }
    } else if(Mod4 == 1) {
        if(idx == 3) {
            r[0] = w.x, r[1] = m.x;
        }
    }
}

//! single iteration of the inner CRT loop parametrized by (n_mods mod 2)
//! \c r_in and \c r_out - shared mem entries for "ping-ponging"
template < unsigned Mod4 >
__device__ __forceinline__ void __CRA_quads_internal(uint4& w, uint4& M, uint4 m,
        fp_limb2 invk1, fp_limb2 invk2, const unsigned *r_in, unsigned *r_out,
                volatile int idx) {

    switch(Mod4) { // compile-time decision
    case 0:
       if(idx <= 0) {
           w.x = submul_m(w.x, r_in[0], M.x, m.x, invk1.x);
           M.x = mul_m(M.x, r_in[1], m.x, invk1.x);
           w.y = submul_m(w.y, r_in[0], M.y, m.y, invk1.y);
           M.y = mul_m(M.y, r_in[1], m.y, invk1.y);
           w.z = submul_m(w.z, r_in[0], M.z, m.z, invk2.x);
           M.z = mul_m(M.z, r_in[1], m.z, invk2.x);
           w.w = submul_m(w.w, r_in[0], M.w, m.w, invk2.y);
           M.w = mul_m(M.w, r_in[1], m.w, invk2.y);
       }
       break;
    case 3:
       if(idx <= 0) {
           w.w = submul_m(w.w, r_in[0], M.w, m.w, invk2.y);
           M.w = mul_m(M.w, r_in[1], m.w, invk2.y);
       }
       if(idx <= 1) {
           w.x = submul_m(w.x, r_in[0], M.x, m.x, invk1.x);
           M.x = mul_m(M.x, r_in[1], m.x, invk1.x);
           w.y = submul_m(w.y, r_in[0], M.y, m.y, invk1.y);
           M.y = mul_m(M.y, r_in[1], m.y, invk1.y);
           w.z = submul_m(w.z, r_in[0], M.z, m.z, invk2.x);
           M.z = mul_m(M.z, r_in[1], m.z, invk2.x);
       }
       break;
    case 2:
       if(idx <= 0) {
           w.w = submul_m(w.w, r_in[0], M.w, m.w, invk2.y);
           M.w = mul_m(M.w, r_in[1], m.w, invk2.y);
           w.z = submul_m(w.z, r_in[0], M.z, m.z, invk2.x);
           M.z = mul_m(M.z, r_in[1], m.z, invk2.x);
       }
       if(idx <= 2) {
           w.x = submul_m(w.x, r_in[0], M.x, m.x, invk1.x);
           M.x = mul_m(M.x, r_in[1], m.x, invk1.x);
           w.y = submul_m(w.y, r_in[0], M.y, m.y, invk1.y);
           M.y = mul_m(M.y, r_in[1], m.y, invk1.y);
       }
       break;
    case 1:
       if(idx <= 0) {
           w.w = submul_m(w.w, r_in[0], M.w, m.w, invk2.y);
           M.w = mul_m(M.w, r_in[1], m.w, invk2.y);
           w.z = submul_m(w.z, r_in[0], M.z, m.z, invk2.x);
           M.z = mul_m(M.z, r_in[1], m.z, invk2.x);
           w.y = submul_m(w.y, r_in[0], M.y, m.y, invk1.y);
           M.y = mul_m(M.y, r_in[1], m.y, invk1.y);
       }
       if(idx <= 3) {
           w.x = submul_m(w.x, r_in[0], M.x, m.x, invk1.x);
           M.x = mul_m(M.x, r_in[1], m.x, invk1.x);
       }
       break;
    }
    // share updated elements between all working threads
    __CRA_shmem_update4< Mod4 >(w, m, r_out, idx);
}

// same as above but processes 4 moduli per thread
template < unsigned BlockSz, unsigned Mod4 >
__global__ void CRA_quads_kernel(unsigned *g_R, const unsigned *g_U,
    const unsigned *g_Mods, const unsigned *g_InvDets, unsigned mods_padding,
            unsigned n_mods, unsigned out_data_padding) {

#if RUN_CRA_KERNEL
    extern __shared__ unsigned shared[];
//! blockIdx.x - indexes over batches
    unsigned thid = threadIdx.x, bidx_x = blockIdx.x,
        nthids = (n_mods + 3) / 4; // # of working thids

    unsigned split = (BlockSz == WS ? threadIdx.y : 0);
    unsigned *r = shared + split * 8;
    //! possibly in texture mem ??
    //! moduli in descreasing order: m[0] > m[1] > ... > m[n-1]
    uint4 m;
    fp_limb2 invk1, invk2;
    if(thid < nthids)
        m = ((uint4 *)g_Mods)[thid];

    //! (m, mu, inv_mod, invk)
    //! bad that reading is offset-dependent (could be hard to debug this..)
    const unsigned *g_InvMods = g_Mods + mods_padding * 2,
                   *g_InvKs = g_InvMods + mods_padding;

    if(thid < nthids) {
#if CUMP_USE_32BIT_MODULI_SET
        uint4 lo = ((uint4 *)g_InvKs)[thid],
              hi = ((uint4 *)(g_InvKs + mods_padding))[thid];

        invk1.x = __hiloint2double(hi.x, lo.x);
        invk1.y = __hiloint2double(hi.y, lo.y);
        invk2.x = __hiloint2double(hi.z, lo.z);
        invk2.y = __hiloint2double(hi.w, lo.w);
#else
        uint4 d = ((uint4 *)g_InvKs)[thid];
        invk1.x = __int_as_float(d.x);
        invk1.y = __int_as_float(d.y);
        invk2.x = __int_as_float(d.z);
        invk2.y = __int_as_float(d.w);
#endif
    }

    //! read is uncoalesced, write back must be coalesced
    unsigned ofs;
    if(BlockSz == WS) {
        ofs = UMUL(thid, out_data_padding*4) + bidx_x * 2 + split;
    } else
        ofs = UMUL(thid, out_data_padding*4) + bidx_x;

    uint4 inv_mod, w, M;
    if(thid < nthids) {
        inv_mod = ((uint4 *)g_InvMods)[thid];
        w.x = g_U[ofs];     // *awful* uncoalesced mem access..grrr
        w.y = g_U[ofs + out_data_padding];
        w.z = g_U[ofs + UMUL(out_data_padding, 2)];
        w.w = g_U[ofs + UMUL(out_data_padding, 3)];
    }

    if(BlockSz == WS) {
        ofs = UMUL(UMUL(bidx_x, UMUL_PAD + 2) + split, mods_padding) +
                 thid*4;
    } else
        ofs = UMUL(bidx_x, mods_padding) + thid*4;

    uint4 det = ((uint4 *)g_InvDets)[thid];
    // # of denominators is moduli-bound
    w.x = mul_m(w.x, det.x, m.x, invk1.x);
    w.y = mul_m(w.y, det.y, m.y, invk1.y);
    w.z = mul_m(w.z, det.z, m.z, invk2.x);
    w.w = mul_m(w.w, det.w, m.w, invk2.y);

    volatile int idx = UMUL(thid, UMUL_PAD + 4) + 3 - (n_mods - 1);
    __CRA_shmem_update4< Mod4 >(w, m, r, idx);

    if(BlockSz > WS)
        CU_SYNC 

    idx++;
    if(idx <= 3) {
        w.x = sub_m(w.x, r[0], m.x);
        w.x = mul_m(w.x, inv_mod.x, m.x, invk1.x);
        M.x = mul_m(r[1], inv_mod.x, m.x, invk1.x);
    }

    if(idx <= 2) {
        w.y = sub_m(w.y, r[0], m.y);
        w.y = mul_m(w.y, inv_mod.y, m.y, invk1.y);
        M.y = mul_m(r[1], inv_mod.y, m.y, invk1.y);
    }

    if(idx <= 1) {
        w.z = sub_m(w.z, r[0], m.z);
        w.z = mul_m(w.z, inv_mod.z, m.z, invk2.x);
        M.z = mul_m(r[1], inv_mod.z, m.z, invk2.x);
    }

    if(idx <= 0) {
        w.w = sub_m(w.w, r[0], m.w);
        w.w = mul_m(w.w, inv_mod.w, m.w, invk2.y);
        M.w = mul_m(r[1], inv_mod.w, m.w, invk2.y);
    }

    __CRA_shmem_update4< (Mod4 + 3) % 4 >(w, m, r + 2, idx);

    idx++;
    {
    volatile int i = n_mods - 3;
    while(1) { // main loop unrolled by the factor of four

        if(BlockSz > WS)
            CU_SYNC

        // ping-ponging shared mem access in each iteration to save on syncs
        __CRA_quads_internal< (Mod4 + 2) % 4 >(w, M, m, invk1, invk2,
                r + 2, r, idx);

        idx++;
        if(BlockSz > WS)
            CU_SYNC

        __CRA_quads_internal< (Mod4 + 1) % 4 >(w, M, m, invk1, invk2,
                r, r + 2, idx);

        idx++;
        if(BlockSz > WS)
            CU_SYNC

        __CRA_quads_internal< Mod4 >(w, M, m, invk1, invk2, r + 2, r, idx);

        idx++;
        if(BlockSz > WS)
            CU_SYNC

        __CRA_quads_internal< (Mod4 + 3) % 4 >(w, M, m, invk1, invk2,
                r, r + 2, idx);

         if(i <= 3)
             break;
        i -= 4; idx++;
    } // while
    }
    if(BlockSz > WS)
        CU_SYNC

    //w.x*m.y + w.y (w.x is the highest MR digit)
//    uint2 u;
//     mad_24_wide((unsigned&)u.y, (unsigned&)u.x, w.x, m.y, w.y);
//     u.y >>= 16, u.y += (u.x < w.y);

Lexit:
    if(thid < nthids) {
        ((uint4 *)(g_R + ofs))[0] = w;
    }
#endif // RUN_CRA_KERNEL
}

#if 0
// use 128 threads per block 
__global__ void test_red_kernel(unsigned *g_R, const unsigned *g_U) {

    extern __shared__ unsigned shared[];

    unsigned thid = threadIdx.x, bidx_x = blockIdx.x;
    unsigned *r = shared;
    unsigned ofs = bidx_x << 7, thid_in_warp = thid & WS-1;

    unsigned a = (g_U + ofs)[thid];

    volatile unsigned *t = (volatile unsigned *)r + HF + UMUL(thid >> 5,
            WS + HF + 1) + thid_in_warp;

    t[-HF] = 0;
    t[0] = a;

    a = a + t[-HF], t[0] = a;

//     (g_R + ofs)[thid] = a;
//     return;
    a = a + t[-8], t[0] = a;
    a = a + t[-4], t[0] = a;
    a = a + t[-2], t[0] = a;
    a = a + t[-1], t[0] = a;
    
    CU_SYNC

    volatile unsigned *t2 = r + HF + UMUL(WS*4 >> 5, WS + HF + 1);

    if(thid < 4) {

        unsigned loc_ofs = HF + WS-1 + UMUL(thid, WS + HF + 1);
        unsigned a2;

        volatile unsigned *ps = t2 + thid;
        ps[-2] = 0;

        a2 = r[loc_ofs]; ps[0] = a2;
        a2 = a2 + ps[-2], ps[0] = a2;
        a2 = a2 + ps[-1], ps[0] = a2;
    }

    CU_SYNC

    a = a + t2[(thid >> 5) - 1];

//     r[0] = 0;
    
    asm volatile("mov.u32 %r11, shared;" : );
    asm volatile("red.shared.add.u32 [%r11], %0;" :
                "+r"(a) : );

//     CU_SYNC

    a = r[0];

    (g_R + ofs)[thid] = a;//r[thid];
}
#endif // RUN_TEST_RED_KERNEL

#endif // _INTERPOLATE_CRT_KERNEL_CU_
