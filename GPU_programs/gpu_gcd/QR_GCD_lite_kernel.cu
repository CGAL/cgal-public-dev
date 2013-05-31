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
 
#ifndef _QR_GCD_LITE_KERNEL_CU_
#define _QR_GCD_LITE_KERNEL_CU_

#include <zmod_dev.cu>      // device modular arithmetic includes

template < bool FullUpdate >
__device__ void __GJG_update_lite(uint4& G, unsigned& t, int& di, int& j,
        int stop_i, unsigned *r, unsigned *g_row, 
        const unsigned ChunkSz, unsigned thid, unsigned m, fp_limb invk, 
        fp_limb invm, volatile unsigned mx100) {

    unsigned *L = r + thid; // space for shift-downs
    const unsigned BlockSz = ChunkSz * 2;
    
    int step = 1;
    unsigned a0, b0, c0, d0; 
    while(1) {

        // read only from the upper half of threads when !UpdateRun
        if(thid == BlockSz - 1 - j) {
            g_row[0] = G.x, g_row[1] = G.y;
            if(FullUpdate) {
                g_row[2] = G.z, g_row[3] = G.w;
 //! NOTE NOTE it seems that we need to look for GCD only when 
 //! start_idx + j >= nu ??
                if(G.x == G.z && G.y == G.w)
                    g_row[4] = 1; // indicates that we have found a gcd
            }
        }
        CU_SYNC

        if(FullUpdate) {
            if(g_row[4] == 1)
                break;
        }

        if(di >= 0) {
            a0 = g_row[0], b0 = g_row[1];
            t = add_mul_reduce_m(G.x, a0, G.y, b0, m, invk, invm, mx100);

            if(FullUpdate) { // compile-time decision
                unsigned s;
                c0 = g_row[2], d0 = g_row[3];
                s = add_mul_reduce_m(G.z, c0, G.w, d0, m, invk, invm, mx100);
                t = sub_m(t, s, m);
            }
            //! NOTE: here we cannot save 't' in L[j] for UpdateRun because
            //! all threads participate, hence L[j] goes out of boundaries
            L[j] = t;   // in fact 'dets' are in registers 't'
        }
        CU_SYNC

        if(di > 0) {
            unsigned l, s; 

            l =  L[j + 1];
            s = sub_m(l, t, m);
            l = r[BlockSz - 1];
                
            G.x = add_mul_reduce_m(G.x, l, s, a0, m, invk, invm, mx100);
            G.y = add_mul_reduce_m(G.y, l, s, b0, m, invk, invm, mx100);

            if(!FullUpdate && thid == 0) {
                G.z = mul_m(G.z, l, m, invk);
            }

            if(FullUpdate) {
                G.z = add_mul_reduce_m(G.z, l, s, c0, m, invk, invm, mx100);
                G.w = add_mul_reduce_m(G.w, l, s, d0, m, invk, invm, mx100);
            }
        }
        j++, di -= step;
        if(j == stop_i)   // after k iterations: j == k    
            break;
    }
}

//! single-kernel version
// template < unsigned ChunkSz >
__global__ void /*CUMP_LAUNCH_BOUNDS(512, 4)*/
QR_gcd_lite_kernel() {

    extern __shared__ unsigned shared[];
//! bidx_x - moduli index
//! bidx_y - index within a batch for single modulus
    unsigned thid = threadIdx.x, bidx_x = blockIdx.x;
    unsigned *r = shared, *mods = dev_const_mem;

    const unsigned BlockSz = blockDim.x, // # of threads == block size
            ChunkSz = BlockSz / 2, m_stride = 4; // const mem  stride

    // L needs BlockSz words + 16 for flags
    unsigned *g_row = r, *misc = g_row + 5,
        *L = r + 6 + thid; // == r + 6
    // g_row: represents the top generator row: [a, b, c, d] - 4 elts + 1 flag

    /*const*/ unsigned nu = dev_const_mem[NU], nv = dev_const_mem[NV];

    mods += UMUL(bidx_x, UMUL_PAD + m_stride) + ConstMemDataSz;
    /*volatile*/ unsigned m = mods[0];
    /*volatile*/ fp_limb invk, invm;

    volatile unsigned mx100;
#if CUMP_USE_32BIT_MODULI_SET
    mx100 = 0, invm = 0;
    invk = __hiloint2double(mods[3], mods[2]);
#else
    mx100 = UMUL(100, m);
    invk = __int_as_float(mods[2]);
    invm = __int_as_float(mods[3]);
#endif
    
    if(thid == 0) {
        g_row[4] = 0; // zero-out g_row flag before loading the main data
    }
    CU_SYNC // really necessary ??

    uint4 G; // represents four generator columns (a,b,c,d)
    uint2 ld; 
    unsigned t;

    /*const */ unsigned max_nu = dev_const_mem[MAX_NU],
            block_ofs = UMUL(bidx_x, dev_const_mem[DATA_PAD]);
    const unsigned *In0 = (unsigned *)dev_const_mem[DATA_IN],
            *f_in = In0 + block_ofs, *g_in = f_in + max_nu;

    G.x = 0, G.y = 0;
    unsigned ofs1 = thid + 1 - ChunkSz*2;
    t = ofs1 + nu;
    if((int)t >= 0) {
        G.x = f_in[t];
    }
    t = ofs1 + nv;
    if((int)t >= 0) {
        G.y = g_in[t];
    }

    ld.x = 0, ld.y = 0;
    t = ofs1 + nu + nv;
    if((int)t >= 0 && t <= nu) {
        ld.x = f_in[t];
    }
    if((int)t >= 0 && t <= nv) {
        ld.y = g_in[t];
    }

    if(thid == 0) {
        G.z = 1;
    }

    int di, j = 0, i_count = nv;
    di = BlockSz - 1 - thid; 

    __GJG_update_lite< false >(G, t, di, j, i_count, r + 6, g_row, ChunkSz,
         thid, m, invk, invm, mx100);

    i_count = nu + nv; // run the remaining full iterations
    if(thid == 0) {
        misc[0] = G.z; // share the product of 'dets' between all threads
    }
    CU_SYNC

    t = misc[0]; 
    G.z = mul_m(ld.x, t, m, invk);
    G.w = mul_m(ld.y, t, m, invk);
    
    // continue iterations from where we stopped
    // NOTE: in fact gcd detection is not required for all iterations: only
    // during the last 'nv' iters
    i_count = nv + nu;
//     i_count = j+3;
    int ii = j;
    __GJG_update_lite< true >(G, t, di, j, i_count, r + 6, g_row, ChunkSz,
         thid, m, invk, invm, mx100);

    block_ofs = UMUL(bidx_x, dev_const_mem[GCD_SZ]) +
                         dev_const_mem[MODS_PAD];
    unsigned *Out0 = (unsigned *)dev_const_mem[DATA_OUT];

    if(g_row[4] != 0) {
       
        if(thid == 0) 
            Out0[bidx_x] = i_count - j + 1;

        ofs1 = BlockSz - 1 - thid;
        if(ofs1 <= nv)
            Out0[block_ofs + ofs1] = L[0];

    } else {
        if(thid == 0) {
            Out0[bidx_x] = 1; // set the gcd size to 1*/
            Out0[block_ofs] = 1;
        }
    }
}

#endif // _QR_GCD_LITE_KERNEL_CU_

