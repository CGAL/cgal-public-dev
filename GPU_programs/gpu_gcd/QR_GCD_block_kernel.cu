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
 
#ifndef _QR_GCD_BLOCK_KERNEL_CU_
#define _QR_GCD_BLOCK_KERNEL_CU_

#include <zmod_dev.cu>      // device modular arithmetic includes

// memory to save the product of 'dets'
__device__ unsigned dev_prods[CUMP_DEVICE_MODULI_SZ / CUMP_MODULI_STRIDE];

// global iteration counter
// __device__ unsigned dev_iters[CUMP_DEVICE_MODULI_SZ / CUMP_MODULI_STRIDE];

//! inner loop of generator updates
//! \c FullUpdate indicates that all four generator columns need to be updated

//! 0th thread returns a product of dets in a register \c G.z (when neither
//! \c FullUpdate nor \c UpdateRun is specified)

//! \c t for threads >= ChunkSz is a set of 'dets' in *reversed* order, i.e.,
//! the first 'det' is given by thread \c BlockSz-1

//! \c L[0] for threads < \ChunkSz is a set of 'outs'

//! returns the number of iterations run
//! shared memory \c r must have enough space for ChunkSz*2 + 1 elements
template < unsigned ChunkSz, bool FullUpdate, bool UpdateRun >
__device__ void __GJG_update(uint4& G, unsigned& t, int& di, int& j,
        int stop_i, const uint4& G0, unsigned *r, const unsigned *Ldets,
        const unsigned *Lins, unsigned *g_row, unsigned thid,
        /*volatile*/ unsigned m, /*volatile*/ fp_limb invk, /*volatile*/ fp_limb invm,
        /*volatile*/ unsigned mx100) {

    unsigned *L = r + thid; // space for shift-downs
    unsigned BlockSz = ChunkSz * 2;
    
    int step = 1;
    unsigned a0, b0, c0, d0; 
    while(1) {

        // read only from the upper half of threads when !UpdateRun
        if(!UpdateRun && thid == BlockSz - 1 - j) {
            g_row[0] = G.x, g_row[1] = G.y;
            if(FullUpdate) {
                g_row[2] = G.z, g_row[3] = G.w;
 //! NOTE NOTE it seems that we need to look for GCD only when 
 //! start_idx + j >= nu ??
                if(G.x == G.z && G.y == G.w)
                    g_row[4] = 1; // indicates that we have found a gcd
            }
        }
        if(UpdateRun && thid == BlockSz - 1 - j) {
            g_row[0] = G0.x, g_row[1] = G0.y;
            if(FullUpdate) {
                g_row[2] = G0.z, g_row[3] = G0.w;
            }
        }
        CU_SYNC

        if(FullUpdate && !UpdateRun) {
            if(g_row[4] == 1)
                break;
/*            if(j == stop_i) // break the loop here because we need to detect
                break;      // the gcd at the very last iteration*/
        }

        if(UpdateRun || di >= 0) {
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
            if(!UpdateRun)
                L[j] = t;   // in fact 'dets' are in registers 't'
            else
                L[0] = t;
        }
        CU_SYNC

        if(UpdateRun || di > 0) {
            unsigned l, s; 

            if(!UpdateRun) {
                l =  L[j + 1];
            } else {
                if(thid == BlockSz - 1)
                    l = Lins[j];
                else
                    l = L[1];
            }
            s = sub_m(l, t, m);

            if(!UpdateRun)
                l = r[BlockSz - 1];
            else
                l = Ldets[j]; 
                
            G.x = add_mul_reduce_m(G.x, l, s, a0, m, invk, invm, mx100);
            G.y = add_mul_reduce_m(G.y, l, s, b0, m, invk, invm, mx100);

            if(!FullUpdate && !UpdateRun && thid == 0) {
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

    // mul the last 'det' only in case it hasn't been collected before
    if(!FullUpdate && !UpdateRun && thid == 0 && j == ChunkSz) {
       G.z = mul_m(G.z, r[BlockSz - 1], m, invk); 
    }
}

template < bool FullUpdate >
__device__ void __init_gmem_load(uint4& G0, uint4& G, uint2& ld,
        const unsigned *f_in, const unsigned *g_in,
        unsigned ofs1, unsigned ofs2, const unsigned nu, const unsigned nv,
        unsigned *shared_det, unsigned thid, /*volatile*/ unsigned m,
        /*volatile*/ fp_limb invk) { 

        G0.x = 0, G0.y = 0;
        // the leading block needs to load only for lower threads
        unsigned t = ofs1 + nu;
        if((int)t >= 0 && (!FullUpdate || t <= nu)) {
            G0.x = f_in[t];
        }

        t = (!FullUpdate ? ofs1 + nv : ofs1 + nu);
        if((int)t >= 0 && (!FullUpdate || t <= nv)) {
            G0.y = g_in[t];
        }

        ld.x = 0, ld.y = 0;
        t = ofs2 + nu;
        if((int)t >= 0 && (!FullUpdate || t <= nu)) {
            ld.x = f_in[t];
        }

        t = (!FullUpdate ? ofs2 + nv : ofs2 + nu);
        if((int)t >= 0 && (!FullUpdate || t <= nv)) {
            ld.y = g_in[t];
        }
        if(!FullUpdate && thid == 0) {
            G0.z = 1;
        }

        if(FullUpdate) {
            t = shared_det[0]; 
            G0.z = mul_m(G0.x, t, m, invk);
            G0.w = mul_m(G0.y, t, m, invk);
            G.z = mul_m(ld.x, t, m, invk);
            G.w = mul_m(ld.y, t, m, invk);
        }
}

// __launch_bounds__(maxThreadsPerBlock, minBlocksPerMultiprocessor)

//! GCD by QR-factorization of Sylvester matrix
//! block gcd kernel executes in several launches to provide the data flow
//! between threads blocks without the need for explicit synchronization

//! \c g_U memory layout (\c DATA_PAD elements per one lane):
//! f[0] f[1] .. f[nu] 0 .. 0 (padding) g[0] g[1] .. g[nv] 0 .. 0 (padding)
//! position of g[0] is marked by \c max_nu
//! number of lanes == # of moduli
//! \c g_BlockIO is used for inter-block communication

//! each block processes two chunks of size \c ChunkSz
//! the first one is the current ``head'' of matrix generator, the second is
//! the working set for this block. # of threads per block = 2*ChunkSz
//! NOTE: \c ChunkSz must be aligned by the half-warp size
template < unsigned ChunkSz, bool FullUpdate >
__global__ void CUMP_LAUNCH_BOUNDS(256, 4)
QR_gcd_block_kernel(unsigned *g_Out, const unsigned *g_In,
                const unsigned start_idx, const unsigned n_ofs,
                const unsigned stream_shift = 0) {

    extern __shared__ unsigned shared[];
//! bidx_x - moduli index
//! bidx_y - index within a batch for single modulus
    unsigned thid = threadIdx.x, bidx_x = blockIdx.x + stream_shift,
            bidx_y = blockIdx.y;
    unsigned *r = shared, *mods = dev_const_mem;

    // memory consumption: BlockSz + 1 + 2*ChunkSz + 11 = 4*ChunkSz + 12
    unsigned BlockSz = ChunkSz * 2, m_stride = 4; // const mem  stride
    // L needs BlockSz words, Lins --> Ldets needs ChunkSz words
    unsigned *L = r + thid, *Lins = r + BlockSz, *Ldets = Lins + ChunkSz;
    
    unsigned *g_row = Ldets + ChunkSz, *misc = g_row + 5;
    // g_row: represents the top generator row: [a, b, c, d] - 4 elts + 1 flag

    // start_idx == 0 indicates that this is a first run and the data has to be
    // read from f_in and g_in, otherwise it is read from g_R in a format:
    // 4 columns (a,b,c,d) of size (nu + nv) each

    // account for reserved space for formal parameters
    const unsigned nu = dev_const_mem[NU], nv = dev_const_mem[NV],
            n_blocks = dev_const_mem[N_BLOCKS]; // batch size per one modulus

    // (n_blocks + 1) * ChunkSz*4 per modulus: same for ``lite'' & full iters
    unsigned io_ofs = UMUL(bidx_x, UMUL(n_blocks,
                4*BlockSz + UMUL_PAD) + ChunkSz*4);

    mods += UMUL(bidx_x, UMUL_PAD + m_stride) + ConstMemDataSz;
    /*volatile*/ unsigned m = mods[0];
    /*volatile*/ fp_limb invk, invm;

    /*volatile*/ unsigned mx100;
#if CUMP_USE_32BIT_MODULI_SET
    mx100 = 0, invm = 0;
    invk = __hiloint2double(mods[3], mods[2]);
#else
    mx100 = UMUL(100, m);
    invk = __int_as_float(mods[2]);
    invm = __int_as_float(mods[3]);
#endif
     
    if(FullUpdate && thid == 0) {
        g_row[4] = 0; // zero-out g_row flag before loading the main data
    }
    CU_SYNC // really necessary ??

    uint4 G, G0; // represents four generator columns (a,b,c,d)
    uint2 ld; // accessing the register can cause memory stall => do not
              // touch them until the right time comes
    unsigned t;

    if((!FullUpdate && start_idx == 0) || (FullUpdate && n_ofs != 0)) {

        /*const */ unsigned max_nu = dev_const_mem[MAX_NU],
                block_ofs = UMUL(bidx_x, dev_const_mem[DATA_PAD]);
        const unsigned *In0 = (unsigned *)dev_const_mem[DATA_IN],
                 *f_in = In0 + block_ofs, *g_in = f_in + max_nu;
        //! NOTE NOTE: here we need to add offsets to f_in & g_in

        unsigned ofs1, ofs2 = thid + 1 - UMUL(bidx_y, BlockSz + UMUL_PAD) -
                 ChunkSz*3;

        if(thid >= ChunkSz) { // load the current update sequence
            ofs1 = thid + 1 - ChunkSz*2; // upper threads: ChunkSz*2
        } else {
            ofs1 = ofs2 + ChunkSz*2; // load one block above
        }

        if(FullUpdate) { // 'sync' before loading data from global mem
            if(thid == 0) { 
                misc[0] = dev_prods[bidx_x];
            }
            CU_SYNC
            ofs1 += n_ofs-1; ofs2 += n_ofs-1; // advance pointers accordingly
        }

    // NOTE: preload data to different register sets, so that there is a
    // a sufficient overlap between memory operations and ALU
        __init_gmem_load< FullUpdate >(G0, G, ld, f_in, g_in, ofs1, ofs2, 
                nu, nv, misc, thid, m, invk);
    }

    if((!FullUpdate && start_idx != 0) || (FullUpdate && n_ofs != 0)) {
        // load data from g_R: here we need to distinguish FullUpdate cases
        g_In += io_ofs;
        unsigned ofs1, ofs2;

        ofs2 = UMUL(bidx_y, 2*BlockSz + UMUL_PAD);
        if(!FullUpdate && thid == 0) { //! NOTE: possible racing condition
            G0.z = dev_prods[bidx_x];
        }
        if(thid >= ChunkSz) {
            ofs1 = -2u*ChunkSz;
        } else {
            ofs1 = ofs2;
            ofs2 += 4*ChunkSz;
        }

        uint2 tmp; 
        tmp = ((uint2 *)(g_In + ofs1))[thid];
        G0.x = tmp.x, G0.y = tmp.y;
        ld = ((uint2 *)(g_In + ofs2))[thid];

    } else if(FullUpdate) { // this branch executes when start_idx != 0
        g_In += io_ofs;
        unsigned ofs1, ofs2;

        if(thid == 0) {
            misc[0] = g_In[0];
        }

        CU_SYNC
        // NOTE: shall we also write g_Out[0] in this case?
        if(misc[0] == CUMP_GCD_TERMINATE_FLAG) {
            if(bidx_y == 0 && thid == 0)
                g_Out[io_ofs] = CUMP_GCD_TERMINATE_FLAG;
            return;
        }

        ofs2 = UMUL(bidx_y, 4*BlockSz + UMUL_PAD);
        if(thid >= ChunkSz) {
            ofs1 = -4u*ChunkSz;
        } else {
            ofs1 = ofs2;
            ofs2 += 8*ChunkSz;
        }

        G0 = ((uint4 *)(g_In + ofs1))[thid];
        G = ((uint4 *)(g_In + ofs2))[thid];
    }
    // G0^up: the current updating sequence
    // G0^low: working set to generate 'Louts'
    // G: the main working set for this block (temporarily loaded to 'ld')

    int di, j, i_count = ChunkSz;

    if(!FullUpdate) { // indicates that we run the last ``lite'' iteration
        if(start_idx + ChunkSz > nv) 
            i_count = nv - start_idx;
    } 

    // NOTE: observe that the leading block does not load G^low
    if(thid < ChunkSz) 
        di = ChunkSz - 1 - thid; // current working set
    else
        di = BlockSz - 1 - thid; // precompute updating sequence

    j = 0;
    if(FullUpdate && n_ofs != 0) { 
        j = n_ofs-1;
        di -= j;
    }
    
    __GJG_update< ChunkSz, FullUpdate, false >(G0, t, di, j, i_count, G0, r,
             Ldets, Lins, g_row, thid, m, invk, invm, mx100);

    // Louts are generated by lower threads while Ldets are generated by
    // upper threads
    if(thid < ChunkSz) {
        Lins[thid] = L[0];  // and 'Louts' of the lower half threads
    } else {
        Ldets[BlockSz - 1 - thid] = t; // save 'Ldets' in reversed order
    }

    if(!FullUpdate || n_ofs != 0) {
        G.x = ld.x, G.y = ld.y;
    }

    if(FullUpdate && g_row[4] == 1) {
        // run the # of iterations given by 'j' parameter
        i_count = j; 
    }

    j = 0; // no need for 'di' here
    if(FullUpdate && n_ofs != 0) {
        j = n_ofs-1;
    }
    CU_SYNC

    __GJG_update< ChunkSz, FullUpdate, true >(G, t, di, j, i_count, G0, r,
             Ldets, Lins, g_row, thid, m, invk, invm, mx100);

    // if we have run all iterations and gcd is not yet detected: check the
    // last possibility: namely that gcd is detected after all i_count
    // iterations have been run: in that case the columns of G are no longer
    // linearly-independent

    //! FIXME: there is still a possibility that 
    if(FullUpdate && thid >= BlockSz - 32 && g_row[4] == 0) {
        bool what = __any(G.x != 0 || G.y != 0);// || G.z != 0_;
        if(what) {
            what = __all(G.x == G.z && G.y == G.w);
            if(thid == BlockSz - 1 && what)
                g_row[4] = 1;
        }
    }  
    CU_SYNC

    if(!FullUpdate) {

#if 0
    if(0) { 
        unsigned block_ofs = UMUL(bidx_x, BlockSz);
        g_Out += block_ofs;

        if(bidx_y == DEBUG_BLOCK_IDX)
            g_Out[BlockSz-1 - thid] = G.x; // for lower threads
    } else  
#endif
    {
        unsigned ofs1 = UMUL(bidx_y, 2*BlockSz + UMUL_PAD);

        if(bidx_y == 0 && thid == 0) { //! carefull !! race condition ? 
            dev_prods[bidx_x] = G0.z;  
        }

        g_Out += io_ofs;
        if(thid < ChunkSz)
            ofs1 += 2*ChunkSz;
        else
            ofs1 -= 2*ChunkSz;

        if(i_count != ChunkSz) { // if i_count == ChunkSz even at the last
                    // ``lite'' iteration => no need to save anything
            if(bidx_y == 0 && thid >= ChunkSz) {
                //! a thread offset: g_Out + ofs + thid*2
                ((uint2 *)(g_Out - 2*ChunkSz))[thid] = make_uint2(G0.x, G0.y);
            }
            ofs1 += 2*ChunkSz; // make enough space for G0's
        } 
        ((uint2 *)(g_Out + ofs1))[thid] = make_uint2(G.x, G.y);
    } 
                 
    } else { // FullUpdate

        if(g_row[4] != 0 /*n_iter == DEBUG_ITER_IDX*/ ) {

            unsigned block_ofs = UMUL(bidx_x, dev_const_mem[GCD_SZ]) +
                             dev_const_mem[MODS_PAD],
                     *Out0 = (unsigned *)dev_const_mem[DATA_OUT];
            unsigned nu = dev_const_mem[NU], nv = dev_const_mem[NV],
                    ofs = UMUL(bidx_y, BlockSz + UMUL_PAD) + ChunkSz;
// leading gcd coeffs are given by L[j-1..ChunkSz-1] elements during generating
// run, the remaining gcd elements are given by L[0] from respective blocks;
// hence the leading block must save two parts: leading L as well as L
// from its block. For the very last block only the lower threads [0..n_last-1]
// provide the valid data
// here j - is the value of the loop counter during the break

            if(bidx_y == 0) { // NOTE: do we need to check the first block ??
                // can just use upper halfwarp to save the data
                if(thid < ChunkSz && thid >= i_count - 1) {
                    Out0[block_ofs + ofs-1u - thid] = Lins[thid];
                }
                // size of gcd polynomial (not the degree!!)
                j = nu + nv + ChunkSz + 1 - (start_idx + i_count); 
                if(thid == ChunkSz) 
                    Out0[bidx_x] = j;

                if(thid == BlockSz-1) {  // set the flag indicating that a gcd
                    g_Out[io_ofs] = CUMP_GCD_TERMINATE_FLAG; // is computed
                }
            }

            ofs += BlockSz - thid - i_count, j = L[0];
            if(ofs <= nv /*&& j != 0*/)
                Out0[block_ofs + ofs] = j;
            
        } else {
            unsigned ofs1 = UMUL(bidx_y, 4*BlockSz + UMUL_PAD);
            g_Out += io_ofs;
            if(thid < ChunkSz)
                ofs1 += 4*ChunkSz;
            else
                ofs1 -= 4*ChunkSz;

            ((uint4 *)(g_Out + ofs1))[thid] = G; // save the whole bunch here
        }
    }

/** algorithm outline:
    start_i = 0;
    while(1) {

        // NOTE: use shmem entry for loop counter ?? instead of registers
        stop_i = ChunkSz;
        if(start_i + ChunkSz > nv)
            stop_i = nv - start_i;

        j = 0;
        if(start_i < nv) { // or stop_i > 0
            GJG_up_lite(UseOuts, stop_i)
        }
        // this is a place when to switch to full updates
        if(j == stop_i || start_i == nv) {
            update_generators()
        }

        //! stop_i = ChunkSz if we are doing only lite iterations now
        //! 0 < stop_i < ChunkSz if we in the middle
        // stop_i <= 0 if we already switched to full iterations
        if((int)stop_i < ChunkSz)
            GJG_up_full((UseOuts, ChunkSz)

        load_outs()

        if(start_i < nv)
            GJG_up_lite(UseIns, stop_i)
        if((int)stop_i < ChunkSz)
            GJG_up_full((UseIns, ChunkSz)
        start_i += N
    }
**/
}

/*!***************************************************************************
******************************************************************************/

//! performs the last iteration when only a single block left,
//! that is, N <= 3*ChunkSz; always assumes \c FullUpdate = true
template < unsigned ChunkSz >
__global__ void CUMP_LAUNCH_BOUNDS(256, 4)
QR_gcd_last_block(unsigned *g_Out, const unsigned *g_In,
            const unsigned start_idx, const unsigned stream_shift = 0) {

    extern __shared__ unsigned shared[];
//! bidx_x - moduli index
//! bidx_y == 0: **last block**
    unsigned thid = threadIdx.x, bidx_x = blockIdx.x + stream_shift;  
    unsigned *r = shared, *mods = dev_const_mem;

    // memory consumption: BlockSz + 1 + 2*ChunkSz + 11 = 4*ChunkSz + 12
    unsigned BlockSz = ChunkSz * 2, m_stride = 4; // const mem  stride
    // L needs BlockSz words, Lins --> Ldets needs ChunkSz words
    unsigned *L = r + thid, *Lins = r + BlockSz, *Ldets = Lins + ChunkSz,
         *g_row = Ldets + ChunkSz, *misc = g_row + 5;
    // g_row: represents the top generator row: [a, b, c, d] - 4 elts + 1 flag

    // account for reserved space for formal parameters
    /*const */unsigned n_blocks = dev_const_mem[N_BLOCKS]; 

    // (n_blocks + 1) * ChunkSz*4 per modulus: same for ``lite'' & full iters
    unsigned io_ofs = UMUL(bidx_x, UMUL(n_blocks,
                4*BlockSz + UMUL_PAD) + ChunkSz*4);

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

    uint4 G, G0; // represents four generator columns (a,b,c,d)
    unsigned t;

    g_In += io_ofs;
    unsigned ofs1, ofs2;

    if(thid == 0) {
        misc[0] = g_In[0];
    }

    CU_SYNC
    if(misc[0] == -1u) {
        if(thid == 0) // do we need to do this ?? nope
            g_Out[io_ofs] = -1u;
        return;
    }

    // here we have only single block left, hence bidx_y == 0
    ofs2 = 0;
    if(thid >= ChunkSz) {
        ofs1 = -4u*ChunkSz;
    } else {
        ofs1 = ofs2;
        ofs2 += 8*ChunkSz;
    }
    G0 = ((uint4 *)(g_In + ofs1))[thid];
    G = ((uint4 *)(g_In + ofs2))[thid];

    int di, j, i_count = ChunkSz;

    // NOTE: observe that the leading block does not load G^low
    if(thid < ChunkSz) 
        di = ChunkSz - 1 - thid; // current working set
    else
        di = BlockSz - 1 - thid; // precompute updating sequence

    j = 0;
// G0^up: generating sequence == G0^low because this is the last block
// hence we can also generate the Louts by the last block
    __GJG_update< ChunkSz, true, false >(G0, t, di, j, i_count, G0, r,
             Ldets, Lins, g_row, thid, m, invk, invm, mx100);

    // Louts are generated by lower threads while Ldets are generated by
    // upper threads
    if(thid < ChunkSz) {
        Lins[thid] = L[0];  // and 'Louts' of the lower half threads
    } else {
        Ldets[BlockSz - 1 - thid] = t; // save 'Ldets' in reversed order
    }

    if(g_row[4] == 1) {
        // run # of iterations given by 'j' parameter
        i_count = j; 
    }

    j = 0;
    CU_SYNC

//!********** if i_count == 0 (causing launch failure) ***********************
//!************* this indicates that algorithm failed to detect **************
//! a GCD during the main iterations !!!
    __GJG_update< ChunkSz, true, true >(G, t, di, j, i_count, G0, r,
             Ldets, Lins, g_row, thid, m, invk, invm, mx100);

    // if we have run all iterations and gcd is not yet detected: check the
    // last possibility: namely that gcd is detected after all i_count
    // iterations have been run: in that case the columns of G are no longer
    // linearly-independent
    if(thid == BlockSz - 1 && g_row[4] == 0) {
    // the upper thread is guaranteed to contain the correct data
        if(G.x == G.z && G.y == G.w) 
            g_row[4] = 1;
    }
    CU_SYNC

    if(g_row[4] == 0) {
        // in this case we have to run the remaining updates,
        // otherwise just output the gcd..
        di = BlockSz - 1 - thid; // precompute updating sequence
        i_count = BlockSz; // NOTE: should be according to the size
        j = 0;             // of the last block
    
        CU_SYNC

        __GJG_update< ChunkSz, true, false >(G, t, di, j, i_count, G, r,
             Ldets, Lins, g_row, thid, m, invk, invm, mx100);

        if(thid == 0 && g_row[4] != 0)
            g_row[4] = 2; // stupid way to distinguish the cases..
        CU_SYNC 
    }

    unsigned block_ofs = UMUL(bidx_x, dev_const_mem[GCD_SZ]) +
                             dev_const_mem[MODS_PAD],
            *Out0 = (unsigned *)dev_const_mem[DATA_OUT];

    if(g_row[4] != 0) {
        unsigned ofs = ChunkSz, nu = dev_const_mem[NU], nv = dev_const_mem[NV];

        if(g_row[4] == 1) { // bidx_y == 0 here
            if(thid < ChunkSz && thid >= i_count - 1)
                Out0[block_ofs + ofs-1u - thid] = Lins[thid];

            j = nu + nv + ChunkSz + 1 - (start_idx + i_count); // size of gcd
        } else {
            i_count = ChunkSz + 1; // to compensate for Out0 += ChunkSz ofs
            j = nu + nv + 1 - (start_idx + j); // size of gcd (not degree !!!)
        }

        if(thid == ChunkSz)
            Out0[bidx_x] = j; // save the size of gcd ***/

        ofs += BlockSz - thid - i_count; // shift all data by 'i_count - 1'
        if(ofs <= nv)
            Out0[block_ofs + ofs] = L[0];
            
   } else {
        // NOTE: this branch means that polynomials are coprime:
        if(thid == 0) {
            Out0[bidx_x] = 1; // set the gcd size to 1*/
            Out0[block_ofs] = 1;
        }
   }
}

#endif // _QR_GCD_BLOCK_KERNEL_CU_

