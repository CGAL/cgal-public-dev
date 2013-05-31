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
 
#ifndef _MOD_REDUCE_KERNEL_CU_
#define _MOD_REDUCE_KERNEL_CU_

#include <zmod_dev.cu>      // modular arithmetic includes


// divides [u1,u0] by d
__device__ unsigned __div2by1(unsigned u1, unsigned u0,
         volatile unsigned d, volatile unsigned v) {

//WRONG WRONG WRONG:
    unsigned q0 = u1 * v + u0, r;
    unsigned q1 = __umulhi(u1, v) + u1 + (q0 < u0); // addc
    q1++;

    r = u0 - q1*d; // 1 mad
    if(r > q0) {
        r = r + d;
    }
    if(r >= d) {
        r = r - d;
    }
    return r;
//     uint64 xx = (((uint64)u1)<<32)+u0;
//     uint64 yy = xx % d;
//     return (unsigned)yy;
}

__device__ void __gload_chunk(unsigned *L, const unsigned *g_In,
    unsigned thid, unsigned chunk_ofs, unsigned chunk_sz) {

    if((int)chunk_ofs < 0) {
        chunk_sz += chunk_ofs;
        chunk_ofs = 0;
    }

    if((int)thid < (int)chunk_sz) {
        L[thid] = g_In[chunk_ofs + thid];
    }
}

  //NOTE NOTE: need to adapt for 24-bit arithmetic !!
// devR, devU, devMods, limbs_f, limbs_g, n_moduli);

 template < unsigned NTHIDS > /*CUMP_LAUNCH_BOUNDS(512, 4)*/
__global__ void mod_reduce_kernel(unsigned *g_Out, const unsigned *g_In,
        const unsigned *g_Mods, const unsigned limbs_f,
        const unsigned limbs_g, const unsigned n_mods) {

    extern __shared__ unsigned shared[];
    unsigned thid = threadIdx.x, bidx_x = blockIdx.x,
            bidx_y = blockIdx.y;

//! bidx_x: indexes over different large integers
//! bidx_y: indexes over batches for the same integer

    const unsigned m_stride = 4; // const mem  stride
    unsigned nu = dev_const_mem[NU], nv = dev_const_mem[NV];

    if(bidx_x <= nu) { // this block works with f coefficients
        g_In += UMUL(limbs_f, bidx_x) + 1;
    } else { // otherwise with g coefficients
        g_In += UMUL(limbs_f, nu + 1) + UMUL(limbs_g, bidx_x - (nu + 1)) + 1;
    }

    const unsigned ChunkSz = 16;
    // L0 and L1 are two cache lanes of size ChunkSz each
    unsigned *r = shared, *L0 = r + 1, *L1 = L0 + ChunkSz, *L;

    if(thid == 0) {
        r[0] = g_In[-1]; // reads in the actual # of limbs
    }

    CU_SYNC
    unsigned n_limbs = abs((int)r[0]),
        chunk_ofs = n_limbs - ChunkSz, c = 0;

    if(n_limbs == 0) // indicates zero coefficient => quit immediately
        goto Lexit;
    
    {
    //! NOTE NOTE: it is assumed that n_limbs is at least 16 !!
    __gload_chunk(L0, g_In, thid, chunk_ofs, ChunkSz);

    CU_SYNC

    chunk_ofs -= ChunkSz;
    // load the second chunk as soon as possible
    __gload_chunk(L1, g_In, thid, chunk_ofs, ChunkSz);

    unsigned m = 0, d, v, ofs;
    //! ATTENTION: thid-based offset !!
    ofs = UMUL(NTHIDS + UMUL_PAD, bidx_y) + thid;
    if(ofs < n_mods)
        m = g_Mods[ofs];
    d = m*2, v = -d;

    uint64 xx = 0xffffffffffffffffULL;
    uint64 yy = (xx / d);// - (1ULL<<32);
    v = (unsigned)yy;

    int i, j;
#if 0
    c = 0;
    for(i = n_limbs - 1; i >= 0; i--) {
        if(ofs < n_mods) {
            uint64 xx = (((uint64)c) << 32) + g_In[i];
            c = (unsigned)(xx % m);
        }
    }
#else
    for(i = n_limbs - 1; i >= 0; i -= ChunkSz) {

        if(ofs < n_mods) // make sure we do not wait for "dummy" threads
                         // for underpopulated blocks
        for(j = min(ChunkSz-1, i); j >= 0; j--) {
            c = __div2by1(c, L0[j], d, v);
        }

        CU_SYNC // formally we need to sync here
        L = L0, L0 = L1, L1 = L;
        chunk_ofs -= ChunkSz;
        __gload_chunk(L1, g_In, thid, chunk_ofs, ChunkSz);
    }

    if(c >= m)
        c -= m;
#endif
    if((int)r[0] < 0)
        c = m - c;
    }
       
Lexit:
    const unsigned max_nu = dev_const_mem[MAX_NU],
        data_pad = dev_const_mem[DATA_PAD],
        nu_ofs4 = dev_const_mem[NU_OFS4],
        nv_ofs4 = dev_const_mem[NV_OFS4];

    if(bidx_x <= nu) {
        g_Out += bidx_x + nu_ofs4;
    } else {
        g_Out += (int)(max_nu + bidx_x - (nu + 1) + nv_ofs4);
    }

    unsigned mod_ofs = UMUL(NTHIDS + UMUL_PAD, bidx_y) + thid;
    g_Out += UMUL(mod_ofs, data_pad);

    if(mod_ofs < n_mods) // to prevent dummy threads write garbage to gmem
        g_Out[0] = c;
}

#endif // _MOD_REDUCE_KERNEL_CU_
