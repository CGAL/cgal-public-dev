// ============================================================================
//
// Copyright (c) 2001-2010 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// this file is not yet part of any library
//
// ----------------------------------------------------------------------------
//
// Library       : 
//
// File          : 
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================
 
#ifndef _MRC_BLOCK_KERNEL_CU_
#define _MRC_BLOCK_KERNEL_CU_

#include <zmod_dev.cu>      // modular arithmetic includes

//! Mixed Radix Conversion (MRC) algorithms

//! loads \c n_load residues and corresponding moduli set from gmem
//! \c g_ofs: offset to read from the global mem
//! \c chunk_ofs: offset of this chunk 
//! residues from \c g_U are read with stride \c in_data_padding
__device__ __forceinline__ void __MRC_gmem_load(
    unsigned& w, unsigned& m, fp_limb& invk,
     unsigned& inv_mod, unsigned thid, unsigned chunk_ofs, unsigned chunk_sz,
     unsigned mods_padding, unsigned in_data_padding, const unsigned *g_U,
    const unsigned *Mods, const unsigned *g_InvMods, const unsigned *g_InvKs,
    const unsigned *g_InvDets) {

    if((int)chunk_ofs < 0) {
        chunk_sz += chunk_ofs;
        chunk_ofs = 0;
    }

    if(thid < chunk_sz) {
        chunk_ofs += thid;

        m = Mods[chunk_ofs];
#if CUMP_USE_32BIT_MODULI_SET
        unsigned lo = g_InvKs[chunk_ofs],
                 hi = (g_InvKs + mods_padding)[chunk_ofs];
        invk = __hiloint2double(hi, lo);
#else
        invk = ((fp_limb *)g_InvKs)[chunk_ofs];
#endif
        inv_mod = g_InvMods[chunk_ofs];

        unsigned det = g_InvDets[chunk_ofs], ofs;
        // total offset: g_ofs + (chunk_ofs + thid) * in_data_padding
        ofs = UMUL(chunk_ofs, in_data_padding)/* + CUMP_GCD_RAW_OFFSET*/;
        w = g_U[ofs];
        w = mul_m(w, det, m, invk); // not a good solution because you load
                             // data and immediately use it.. no ALU overlap
    } else {
        w = 0; m = 2;
    }
}

//! performs the first step of MRC algorithm
//! just multiplies w by inv_mods and sets M = inv_mod
__device__ __forceinline__ void __MRC_init(unsigned *w, unsigned *M,
     unsigned *m, fp_limb *invk, unsigned inv_mod, const unsigned k) {

// NOTE: this even can be used for 0th thid because inv_mod == 1 for it !!
    w[k] = mul_m(w[k], inv_mod, m[k], invk[k]);
    M[k] = inv_mod;
}

//! computes MRC digits in a "stair-like" fashion
//! \c w[0] digits are done, \c w[1] need additional step
//! performs exactly \c n_iters iterations
//! only those digits for which \c idx < 0 get updated !
//! digit for which \c idx == 0 is the first one to be shared and *not* updated
//! if \c LastChunk is set, only threads for which \c active is on 
//! update the values of \c wy
__device__ __forceinline__ void __MRC_stairs(unsigned *w, unsigned* M,
      unsigned *m, fp_limb *invk, unsigned *r1, unsigned *r2, 
        volatile int idx, int n_iters, bool is_active) {

    int i = n_iters;
    // declaring volatile increases register pressure ??
// NOTE: makes sense to unroll the loop by the factor of 2 ?
// no need to swap r1 and r2 then !
    while(i > 0) {
        if(idx == 0) {
            r1[0] = w[0], r1[1] = m[0];
        }
        idx++, i--;
        CU_SYNC
    
        if(idx <= 0) {
            w[0] = submul_m(w[0], r1[0], M[0], m[0], invk[0]);
            M[0] = mul_m(M[0], r1[1], m[0], invk[0]);
        }

        if(is_active) {
            w[1] = submul_m(w[1], r1[0], M[1], m[1], invk[1]);
            M[1] = mul_m(M[1], r1[1], m[1], invk[1]); 
        }
        unsigned *t = r1; r1 = r2, r2 = t;
    }
}

//! computes MRC digits in a "stair-like" fashion
__device__ __forceinline__ void __MRC_stairs_final(unsigned *w, unsigned *M,
      unsigned *m, fp_limb *invk, const unsigned k, unsigned *r1, 
        unsigned *r2, volatile int idx, int n_iters) {
// TODO: need to unroll by the factor of 2 and use half the threads
    int i = n_iters;
    while(i > 0) {
        if(idx == 0) {
            r1[0] = w[k], r1[1] = m[k];
        }
        idx++, i--;
        CU_SYNC
    
        if(idx <= 0) {
            w[k] = submul_m(w[k], r1[0], M[k], m[k], invk[k]);
            M[k] = mul_m(M[k], r1[1], m[k], invk[k]);
        }
        unsigned *t = r1; r1 = r2, r2 = t;
    }
}

//! evaluates a block of \c chunk_sz MR digits using the set of already
//! computed digits \c ys with resp. moduli \c ms stored in shared mem
//! if \c LastChunk is set, only threads for which \c active is on 
//! update the values of \c wy
__device__ __forceinline__ void __MRC_block(unsigned *w, unsigned *M,
      unsigned *m, fp_limb *invk, const unsigned k, unsigned *ys,
      unsigned *ms, unsigned n_iters) {

    int i = 0;
    while(1) { // unrolled by the factor of 2, can do more..
        w[k] = submul_m(w[k], ys[i], M[k], m[k], invk[k]);
        M[k] = mul_m(M[k], ms[i], m[k], invk[k]);
        
        w[k] = submul_m(w[k], ys[i+1], M[k], m[k], invk[k]);
        M[k] = mul_m(M[k], ms[i+1], m[k], invk[k]);

        i += 2;
        if(i == n_iters)
            break;
    }
}

/*!**************************************************************************

Mixed radix conversion (MRC) algorithm: input is split in
chunks[0..n_chunks-1] depending on the number of threads per block

MRC_load(chunk[0], chunk[1])

for(i = 0; i < n_chunks - 1; i++) {
    MRC_stairs();

    if(i > 0)
        MRC_load(chunk[i + 1]);
    for(j = 0; j < i; j++) {
        MRC_block();
    }
}
MRC_stair_last();

nchunks = 3
i = 0..1
when i == 1: load chunk[2]

nchunks = 4
i = 0..2
when i == 1: load chunk[2]
when i == 2: load chunk[3]: j=0..1 (two iters)

****************************************************************************/

template < unsigned N_CHUNKS >
__global__ void /*CUMP_LAUNCH_BOUNDS(512, 2)*/
MRC_block_kernel(unsigned *g_R, const unsigned *g_U,
    const unsigned *g_Mods, const unsigned *g_InvDets, 
    unsigned mods_padding, unsigned n_mods, unsigned in_data_padding) {

    extern __shared__ unsigned shared[];
//! blockIdx.x - indexes over batches
    unsigned thid = threadIdx.x;
    const unsigned chunk_sz = blockDim.x; // # of threads == chunk size
    volatile unsigned n_last = n_mods - UMUL(N_CHUNKS, chunk_sz) + chunk_sz;
   
//     unsigned split = 0; //(BlockSz == WS ? threadIdx.y : 0);
//     unsigned *r = shared, *r1 = r + 2, *r2 = r, *ys = r + 4, 
//         *ys2 = ys + chunk_sz, *ms = ys2 + chunk_sz, *ms2 = ms + chunk_sz;
    unsigned *r = shared, *r1 = r + 2, *r2 = r, *ys = r + 4, 
        *ms = ys + chunk_sz;
    bool active = true;

    //! moduli in descreasing order: m[0] > m[1] > ... > m[n-1]
    unsigned inv_mod, m[4], w[4], M[2];
    fp_limb invk[2];
    //! attention: \c g_R address is per-thread now !!
    g_R += UMUL(blockIdx.x, mods_padding) + thid; 
    g_U += blockIdx.x; // in-offset

    //! (m, mu, inv_mod, invk)
    const unsigned *Mods = g_Mods, *g_InvMods = Mods + mods_padding * 2,
               *g_InvKs = g_InvMods + mods_padding;
    
    //! NOTE: observe that \c w[0] are loaded from the high addresses !!
    __MRC_gmem_load(w[0], m[0], invk[0], inv_mod, thid, n_mods - chunk_sz,
         chunk_sz, mods_padding, in_data_padding, g_U, Mods, g_InvMods,
             g_InvKs, g_InvDets);

    // inv_mod[n - 1] = 1 inv_mod[n - 2] = m[n - 1]^-1 mod m[n - 2]
    // inv_mod[n - 3] = m[n - 1]m[n - 2]^-1 mod m[n - 3] ..

    if(thid != chunk_sz - 1) { // leading thid shares its associated residue
        __MRC_init(w, M, m, invk, inv_mod, 0);
    }

    // no need to wait for SYNC then
    __MRC_gmem_load(w[1], m[1], invk[1], inv_mod, thid, n_mods - chunk_sz*2,
         chunk_sz, mods_padding, in_data_padding, g_U, Mods, g_InvMods,
         g_InvKs, g_InvDets);

    __MRC_init(w, M, m, invk, inv_mod, 1);

    int i = 0;
    #pragma unroll
    while(1) {

        volatile int idx = thid - (chunk_sz - 1);
        active = (i < N_CHUNKS - 2 || thid < n_last);
        __MRC_stairs(w, M, m, invk, r1, r2, idx, chunk_sz, active);

        if(N_CHUNKS == 2 || i++ == N_CHUNKS - 2)
            break;  

        //! remember that MR digits are computed in reversed order !
        //! that is ys[chunk_sz-1] is the 1st (lowest)
        //! and ys[0] is the highest digit
        if(N_CHUNKS == 3 || i == 1) { // save to ys/ms on the first run
            ys[chunk_sz-1 - thid] = w[0], ms[chunk_sz-1 - thid] = m[0];

        } else if(N_CHUNKS >= 4 && i == 2) {
//             ys2[chunk_sz-1 - thid] = w[0], ms2[chunk_sz-1 - thid] = m[0];
            w[2] = w[0], m[2] = m[0];
        } else if(N_CHUNKS == 5 && i == 3) { // means i == 3 ??
            w[3] = w[0], m[3] = m[0];
        }
        CU_SYNC

        unsigned ofs = n_mods - UMUL(chunk_sz, i);
        g_R[ofs] = w[0];
        ofs -= chunk_sz; 
 
         // do not overwrite w[0]'s to overlap gmem load-store with ALUs
        w[0] = w[1], M[0] = M[1], m[0] = m[1], invk[0] = invk[1];
        __MRC_gmem_load(w[1], m[1], invk[1], inv_mod, thid, ofs - chunk_sz,
            chunk_sz, mods_padding, in_data_padding, g_U, Mods, g_InvMods,
                 g_InvKs, g_InvDets);
        __MRC_init(w, M, m, invk, inv_mod, 1);

// N_CHUNKS = 3: i = 0..1; j = [0..i-1]
// N_CHUNKS = 4: i = 0..2
// N_CHUNKS = 5: i = 0..3

        active = (i < N_CHUNKS - 2 || thid < n_last);
        unsigned n_iters = chunk_sz;//UMUL(i, chunk_sz);
        // here we use consecutive memory from ys and than from ys2
        if(active)
            __MRC_block(w, M, m, invk, 1, ys, ms, n_iters);

        unsigned x, y;
        if(N_CHUNKS == 5 && i == 2) {
            x = ys[chunk_sz-1 - thid], y = ms[chunk_sz-1 - thid];
            // keep ws/ms for the next iteration
        }

        if(N_CHUNKS >= 4 && i >= 2) {
            CU_SYNC // wait until previous MRC block is done otherwise 
                    // shared mem might be corrupt

            ys[chunk_sz-1 - thid] = w[2], ms[chunk_sz-1 - thid] = m[2];
            CU_SYNC
            if(active)
                __MRC_block(w, M, m, invk, 1, ys, ms, n_iters);
        } 

        if(N_CHUNKS == 5) {
            if(i == 2) {
                CU_SYNC
                ys[chunk_sz-1 - thid] = x, ms[chunk_sz-1 - thid] = y;
                CU_SYNC
            } else if(i == 3) {
                CU_SYNC
                ys[chunk_sz-1 - thid] = w[3], ms[chunk_sz-1 - thid] = m[3];
                CU_SYNC
                if(active)
                    __MRC_block(w, M, m, invk, 1, ys, ms, n_iters);
            }
        }
    }
    
    g_R[n_last] = w[0];
    volatile int idx = thid - (n_last - 1);
    //  chunk_sz-1 iterations because 1st digit wz[0] is correct (only shared)
    __MRC_stairs_final(w, M, m, invk, 1, r1, r2, idx, n_last-1);

    if(thid < n_last)
        g_R[0] = w[1];
}

/*! **************************************************************************
********* old version: single block MRC algorithm
******************************************************************************/

template < bool IsEven >
__device__ __forceinline__ void __MRC_shmem_update2(uint2 w, uint2 m, unsigned *r,
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
__device__ __forceinline__ void __MRC_doubles_internal(uint2& w, uint2& M, uint2 m,
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
    __MRC_shmem_update2< IsEven >(w, m, r_out, idx);
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
__global__ void MRC_doubles_kernel(unsigned *g_R, const unsigned *g_U,
    const unsigned *g_Mods, const unsigned *g_InvDets, unsigned mods_padding,
            unsigned n_mods, unsigned in_data_padding) {

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
        ofs = UMUL(UMUL(thid, in_data_padding) + bidx_x,
                UMUL_PAD + 2) + split;
    else
        ofs = UMUL(thid, in_data_padding*2) + bidx_x;

    uint2 inv_mod, w, M;
    if(thid < nthids) {
        inv_mod = ((uint2 *)g_InvMods)[thid];

//         ofs += CUMP_GCD_RAW_OFFSET;
        w.x = g_U[ofs];     // *awful* uncoalesced mem access..grrr
        w.y = g_U[ofs + in_data_padding];
    }

    if(BlockSz == WS)
        ofs = UMUL(UMUL(bidx_x, UMUL_PAD + 2) + split, mods_padding) +
                thid*2;
    else
        ofs = UMUL(bidx_x, mods_padding) + thid*2;

    if(thid < nthids) {
        uint2 det = ((uint2 *)g_InvDets)[thid];
//      # of denominators is moduli-bound
        w.x = mul_m(w.x, det.x, m.x, invk.x);
        w.y = mul_m(w.y, det.y, m.y, invk.y);
    }

    // inv_mod[n - 1] = 1
    // inv_mod[n - 2] = m[n - 1]^-1 mod m[n - 2]
    // inv_mod[n - 3] = m[n - 1]m[n - 2]^-1 mod m[n - 3] ..
// NOTE: in fact we can handle n + 1 residues with n threads
    volatile int idx = UMUL(thid, UMUL_PAD + 2) + 1 - (n_mods - 1);

    __MRC_shmem_update2< IsEven >(w, m, r, idx);

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

/*!********************************************************************
!********************************************************************
!********************************************************************/
//     goto Lexit;

    __MRC_shmem_update2< !IsEven >(w, m, r + 2, idx);

    idx++;
    {
    volatile int i = n_mods - 3;
    while(1) { // main loop unrolled by the factor of two

        if(BlockSz > WS)
            CU_SYNC
    
        // ping-ponging shared mem access in each iteration to save on syncs
        __MRC_doubles_internal< IsEven >(w, M, m, invk, r + 2, r, idx);

        idx++;
        if(BlockSz > WS)
            CU_SYNC

        __MRC_doubles_internal< !IsEven >(w, M, m, invk, r, r + 2, idx);
        
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
}

#if 0
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

//! single iteration of the inner CRT loop parametrized by (n_mods mod 4)
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

    if(thid < nthids) {
        ((uint4 *)(g_R + ofs))[0] = w;
    }
}
#endif

//! computes modular inverse of leading coefficients to obtain monic gcd
//! \c g_GCDLcf - set of residues of gcd(lc(f),lc(g))
template < unsigned ChunkSz >
__global__ void mod_inverse_kernel1(unsigned *g_Out, const unsigned *g_In,
        const unsigned *g_Mods, const unsigned *g_GCDLcf, 
        unsigned mods_padding, unsigned n_mods, unsigned in_data_padding) {

    extern __shared__ unsigned shared[];
    unsigned thid = threadIdx.x, bidx_x = blockIdx.x;

//! bidx_x: indexes over chunks of respective size
    unsigned thid_ofs = UMUL(bidx_x, ChunkSz + UMUL_PAD) + thid;
    unsigned ofs = UMUL(thid_ofs, in_data_padding), a, m, mu, lcf;
    fp_limb invk;

    //! (m, mu, inv_mod, invk)
    const unsigned *Mods = g_Mods, *Mus = Mods + mods_padding,
            *InvKs = Mods + UMUL(mods_padding, 3);

    if(thid_ofs >= n_mods)
        return;

#if CUMP_USE_32BIT_MODULI_SET
    unsigned lo = InvKs[thid_ofs],
             hi = (InvKs + mods_padding)[thid_ofs];
    invk = __hiloint2double(hi, lo);
#else
    invk = ((fp_limb *)InvKs)[thid_ofs];
#endif
    // read leading elements with stride 'in_data_padding'
    a = g_In[ofs/* + CUMP_GCD_RAW_OFFSET*/];
    m = Mods[thid_ofs], mu = Mus[thid_ofs];
    lcf = g_GCDLcf[thid_ofs];
    
    if(a != 0) 
        a = montgomery_inverse(a, m, mu);
    g_Out[thid_ofs] = mul_m(a, lcf, m, invk);
}

#endif // _MRC_BLOCK_KERNEL_CU_
