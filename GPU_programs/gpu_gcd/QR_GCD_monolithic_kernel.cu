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
 
#ifndef _QR_GCD_MONOLITHIC_KERNEL_CU_
#define _QR_GCD_MONOLITHIC_KERNEL_CU_

#include <zmod_dev.cu>      // device modular arithmetic includes

__device__ unsigned dev_outs_flag[CUMP_DEVICE_MODULI_SZ / CUMP_MODULI_STRIDE];

template < unsigned ChunkSz >
__device__ void __shmem_setup(unsigned *r, unsigned *& Lins,
        unsigned *& Ldets, unsigned *& g_row, unsigned *& break_i,
        unsigned *& stop_i, unsigned *& count_i, unsigned *& shift_f,
        unsigned *& misc) {

    const unsigned BlockSz = ChunkSz*2;
    Lins = r + BlockSz, Ldets = Lins + ChunkSz, g_row = Ldets + ChunkSz;
    break_i = g_row + 8, stop_i = break_i + 1, count_i = stop_i + 1,
    shift_f = count_i + 1, misc = shift_f + 1;
}

//! load and premultiply the remaining generator columns c & d by 'det'
__device__ void __setup_generators(uint4& G,
        const unsigned *f_in, const unsigned *g_in, unsigned det,
        const unsigned read_ofs, const unsigned nu, const unsigned nv,
        unsigned m, volatile fp_limb invk) {

    G.z = 0, G.w = 0;  // read the remaining columns c & d

    unsigned ofs = read_ofs + nu + nv;
    // for the most blocks this will be just dummy instructions without
    // running the real 'load' 
    if((int)ofs > 0 && ofs <= nu) { // we do not load 0th element
        G.z = f_in[ofs];
    }
    if((int)ofs > 0 && ofs <= nv) {
        G.w = g_in[ofs];
    }

    G.z = mul_m(G.z, det, m, invk);  // mul columns by the denominator
    G.w = mul_m(G.w, det, m, invk);
}

//! exchanges data between thread blocks: block in [0..n_blocks-2] writes out
//! its data. Then blocks [1..n_blocks-1] read data from respective previous
//! blocks; data is loaded from \c thread-indexed shared mem pointed to by
//! \c Louts and saved to shared mem \c Lins
//! only lower \c ChunkSz threads participate

//! \c break_i flag indicates that we do not load new update sequence
template < unsigned ChunkSz, bool FullUpdate >
__device__ void __load_outs(uint4& G, unsigned *Lins,
    const unsigned *Louts, unsigned *g_Outs, unsigned *g_Head,
    unsigned bidx_x, unsigned bidx_y, unsigned *count_i, unsigned *break_i,
        unsigned thid) {

    volatile unsigned *pmutex = (volatile unsigned *)dev_outs_flag + bidx_x;
    // all except the last block
    if(bidx_y < gridDim.y-1) {
        if(thid < ChunkSz) { // write outs
        //! NOTE: better add bidx_y*ChunkSz to g_Outs offset directly ??
            g_Outs[thid + UMUL(bidx_y, ChunkSz + UMUL_PAD)] = Louts[0];
        }
        if(bidx_y == count_i[0] / 2) {
            if(count_i[0] & 1) { // odd counter: save the lower part of G
            if(thid < ChunkSz) {
                if(FullUpdate)
                    ((uint4 *)g_Head)[thid] = G;
                else if(!break_i[0])
                    ((uint2 *)g_Head)[thid] = make_uint2(G.x, G.y);
            } } else { // even counter: save the upper part of G
            if(thid >= ChunkSz) {
                if(FullUpdate)
                    ((uint4 *)g_Head)[thid - ChunkSz] = G;
                else if(!break_i[0])
                    ((uint2 *)g_Head)[thid - ChunkSz] = make_uint2(G.x, G.y);
            } }
        }

        CU_MEMFENCE
        CU_SYNC

        // all blocks in a group use a single mutex
        if(thid == 0) {
#ifdef CUMP_USE_ATOMICS
            __uAtomicAdd((unsigned *)pmutex, 1);
#else
            *((unsigned *)pmutex) += 1;
#endif
        }
        // leading block does not participate on 'even' iterations
        if(bidx_y == count_i[0] / 2 && (count_i[0] & 1) == 0)
            return;
    }
    if(thid == 0) {
        while((int)pmutex[0] < bidx_y);
    }
    CU_SYNC

// probably leading block does not have to wait on a mutex to read G
// but in general does not matter

    if(thid < ChunkSz) { // load the next leading rows to *lower threads*
        if(FullUpdate)
            G = ((uint4 *)g_Head)[thid];
        else if(!break_i[0]) {
            uint2 in = ((uint2 *)g_Head)[thid];
            G.x = in.x, G.y = in.y; 
        }
    } else { 
        Lins[thid - ChunkSz] = g_Outs[thid +
                    UMUL(bidx_y-2, ChunkSz + UMUL_PAD)];
    }
}

//! inner loop of generator updates
//! \c FullUpdate indicates that all four generator columns need to be updated

//! 0th thread returns a product of dets in a register \c G.z (when neither
//! \c FullUpdate nor \c UpdateRun is specified)

//! \c t for threads >= ChunkSz is a set of 'dets' in *reversed* order, i.e.,
//! the first 'det' is given by thread \c BlockSz-1

//! \c L[0] for threads < \ChunkSz is a set of 'outs'

//! returns the number of iterations run
//! shared memory \c r must have enough space for ChunkSz*2 + 1 elements

//! \c ID == 0: not update run; \c ID == 1: update run;
//! \c ID == 2: full updates for lower half; lite updates for upper half
template < unsigned ChunkSz, bool FullUpdate, bool UpdateRun >
__device__ void __GJG_update(uint4& G, unsigned& t, int& di, int& j,
        int stop_i, const uint4& G0, unsigned *r, const unsigned *Ldets,
        const unsigned *Lins, unsigned *g_row, unsigned *shift_f,
        unsigned thid, unsigned m, volatile fp_limb invk,
        volatile fp_limb invm, volatile unsigned mx100) {

    unsigned *L = r + thid; // space for shift-downs
    const unsigned BlockSz = ChunkSz * 2;

// updates of 'di' in !UpdateRun mode (split by 2 thread chunks)
//         1st  2nd  <-- update location (di += step)
// lower:   0   dec
// upper:  inc   0

    int step = 1;
    if(!UpdateRun) {
        if(thid < ChunkSz)
            g_row += 4;
        else
            step = 0;
    }

    unsigned a0, b0, c0, d0, s; /*j = start_j[0];*/
    while(1) {

        if(di == 0) { // not sure if this is a good solution, although
                      // we split conditions on the warp boundary
            if(!UpdateRun && thid < ChunkSz) {
                g_row[0] = G.x, g_row[1] = G.y;
                if(FullUpdate) {
                    g_row[2] = G.z, g_row[3] = G.w;
//                 if(a == c && b == d)
//                     g_row[4] = 1; // flag indicating exit: gcd computed
                }
            } else {
                g_row[0] = G0.x, g_row[1] = G0.y;
                if(FullUpdate) {
                    g_row[2] = G0.z, g_row[3] = G0.w;
                }
            }
        }
        CU_SYNC

//         if(FullUpdate && g_row[4] == 1)
//             break;

        if((UpdateRun && shift_f[0]) || di >= 0) {
            a0 = g_row[0], b0 = g_row[1];
            t = add_mul_reduce_m(G.x, a0, G.y, b0, m, invk, invm, mx100);

            if(FullUpdate) { // compile-time decision
                c0 = g_row[2], d0 = g_row[3];
                s = add_mul_reduce_m(G.z, c0, G.w, d0, m, invk, invm, mx100);
                t = sub_m(t, s, m);
            }
            //! FIXME: when shift_f[0] && UpdateRun is set, saving 't' in L[j]
            //! might be dangerous because thread index goes out of bounds
            L[j] = t;   // in fact 'dets' are in registers 't'
        }
        CU_SYNC

        if(!UpdateRun) { // 'di' needs to be changed for upper threads in 
            di += (step ^ 1);  // !UpdateRun mode because we use 'Lins' here
        }
    
        if((UpdateRun && shift_f[0]) || di > 0) {
            unsigned l; 

            if((!UpdateRun || shift_f[0]) && thid == BlockSz - 1) {
                l = Lins[j]; // for leading thread
            } else
                l = L[j + 1]; 

            s = sub_m(l, t, m);

            if(!UpdateRun && thid < ChunkSz) { // checked for !UpdateRun
                l = r[ChunkSz - 1];
            } else
                l = Ldets[j]; // r[BlockSz - 1 - j]; 
                
            G.x = add_mul_reduce_m(G.x, l, s, a0, m, invk, invm, mx100);
            G.y = add_mul_reduce_m(G.y, l, s, b0, m, invk, invm, mx100);

            // NOTE: G.z is now computed twice..
            if(!FullUpdate && !UpdateRun && thid == 0) {
                G.z = mul_m(G.z, l, m, invk);
            }

            if(FullUpdate) {
                G.z = add_mul_reduce_m(G.z, l, s, c0, m, invk, invm, mx100);
                G.w = add_mul_reduce_m(G.w, l, s, d0, m, invk, invm, mx100);
            }
        }
        j++, di -= step;

        if(j == stop_i) // after k iterations: j == k
            break;
    }

    if(!FullUpdate && !UpdateRun && thid == 0) {
        G.z = mul_m(G.z, r[ChunkSz - 1], m, invk); // mul the last one
    }

//     if(thid == 0) {
//         start_j[0] = j; // update the counter
//     }
    //! if we break the loop earlier: gcd is given by in register 't' by
    //! threads [0 .. i+1]
}

template < unsigned ChunkSz, bool FullUpdate >
__device__ void __main_loop_internal(uint4& G, uint4& G0, int& di, int& j,
        unsigned nv, unsigned *g_Head, unsigned *g_Outs,
        unsigned *r, unsigned thid, unsigned m, volatile fp_limb invk,
        volatile fp_limb invm, volatile unsigned mx100) {

    const unsigned BlockSz = ChunkSz*2;
    unsigned *L = r + thid, *Lins, *Ldets, *g_row,
            *break_i, *stop_i, *count_i, *shift_f, *misc;
    // g_row: represents the top generator row: [a, b, c, d] - 4 elements

    __shmem_setup< ChunkSz >(r, Lins, Ldets, g_row, break_i, stop_i, count_i,
        shift_f, misc);

    const unsigned bidx_x = blockIdx.x, bidx_y = blockIdx.y;
    unsigned t;

    while(1) {

    if(thid == 0 && (!break_i[0] || FullUpdate)) {
        j = UMUL(++count_i[0], ChunkSz + UMUL_PAD);
        t = ChunkSz;
        if(j + ChunkSz > nv)
            t = nv - j;
        stop_i[0] = t;
    }
    CU_SYNC  // sync to broadcast start_i[0] for all thids

    j = 0;
    if(stop_i[0] == 0) // handling the boundary case: no need for lite iters.
        break;

// G0: upper: leading rows of the current generator (needed for update)
// G0: lower: leading rows of the (next) generator
// G - current data with upper threads need to be updated using Lins

    __GJG_update< ChunkSz, FullUpdate, false >(G, t, di, j, stop_i[0], G0, r,
             Ldets, Lins, g_row, shift_f, thid, m, invk, invm, mx100);
//! \c legend: G.z for thid == 0 contains the product of all dets
//! only for upper half of threads
//! t - set of 'dets' in reversed order, that is, dets[0] is for the last thid
//! L[0] for lower threads: set of Louts to be used for shift-through

/** ****************************************************************/

// here we should have G^low - new update sequence prepared
//     if(count_i[0] == 1)
//         break;
/** ****************************************************************/

    if(break_i[0]) // indicates that we are about to switch to the
        break;     // full iterations

    if(FullUpdate) {
        if(thid >= ChunkSz)
            j = L[-ChunkSz]; // protect Louts from modification *upper half*
    }
// moving data: G^low -> G0^up (bulge update sequence)
// G0^low -> G^low (working set of this block)

    CU_SYNC
    if(thid < ChunkSz) {
        Ldets[thid] = G.x;
        Lins[thid] = G.y;
        if(FullUpdate) {
            L[0] = G.z, L[ChunkSz] = G.w;
        }
    }
    CU_SYNC

    if(thid >= ChunkSz) {
        G0.x = Ldets[thid - ChunkSz];
        G0.y = Lins[thid - ChunkSz];
        if(FullUpdate) {
            G0.z = L[-ChunkSz], G0.w = L[0];
        }
    } else {
        // what if we do not execute this in case of break_i
        G.x = G0.x, G.y = G0.y;
        if(FullUpdate) {
            G.z = G0.z, G.w = G0.w;
        }
    }

    CU_SYNC // NOTE: use memory in such a way to avoid CU_SYNC ?

    if(thid < ChunkSz) {
        Ldets[ChunkSz - 1 - thid] = t; // save in reversed order
    } else { // ChunkSz ofs because Louts are for upper threads only
        if(!FullUpdate)
            Lins[thid - ChunkSz] = L[-ChunkSz];
        else
            Lins[thid - ChunkSz] = j;
    }

    di = BlockSz - 1 - thid, j = 0;

    unsigned lead_i = (bidx_y == count_i[0] / 2);
    if(thid == 0) {
        // meaning that shift-through is only used for 'even' interations
        shift_f[0] = lead_i & !(count_i[0] & 1);
    }
    CU_SYNC

    // at this point count_i[0] is already 1
    __GJG_update< ChunkSz, FullUpdate, true >(G, t, di, j, stop_i[0], G0, r,
             Ldets, Lins, g_row, shift_f, thid, m, invk, invm, mx100);

/** ****************************************************************/
     if(count_i[0] == 1)
         break;
/** ****************************************************************/

    if(thid < ChunkSz)
        di = ChunkSz - 1 - thid;
    else 
        di = thid - (BlockSz - 1); // UseIns for upper threads

    if(bidx_y == count_i[0] / 2 && !(count_i[0] & 1)) {
        if(thid >= ChunkSz)
            di = -BlockSz; // disable upper threads
    }
   
    if(stop_i[0] < ChunkSz) { // indicating the end of lite iterations

        if(thid < ChunkSz) // disable lower threads
            di = -BlockSz; // skip generating the update sequence
        if(thid == 0)  // set the flag indicating the last iteration
            break_i[0] = 1;
    }

    if(thid < ChunkSz) { // save G to protect from modifications
        G0.x = G.x, G0.y = G.y; 
        if(FullUpdate) {
            G0.z = G.z, G0.w = G.w;
        }
    }

    CU_SYNC // sync because 'count_i' is to be changed ?? need this ?

//! legend: \c count_i is even: lead G^up -> G^low
//!         \c count_i is even: lead G^low -> G^low
    __load_outs< ChunkSz, FullUpdate >(G, Lins, L, g_Outs, g_Head, bidx_x,
             bidx_y, count_i, break_i, thid);
    // G0^up: current updating sequence
    // G^up: working set to be updates
    // G^low: the next updating sequence

//     if(bidx_y == count_i[0]) { // leading thread block uploads its results
//                                // to global memory
//         break;
//     }
    }
}

//! GCD by QR-factorization of Sylvester matrix
//! monolithic GCD kernel uses thread memory fence for interblock
//! synchronization hence only single kernel launch is necessary (experimental)

//! \c g_U memory layout (\c data_padding elements per one lane):
//! f[0] f[1] .. f[nu] 0 .. 0 (padding) g[0] g[1] .. g[nv] 0 .. 0 (padding)
//! position of g[0] is marked by \c nv_start
//! number of lanes == # of moduli
//! \c g_BlockIO is used for inter-block communication

//! each block processes two chunks of size \c ChunkSz
//! the first one is the current ``head'' of matrix generator, the second is
//! the working set for this block. # of threads per block = 2*ChunkSz
//! NOTE: \c ChunkSz must be aligned by the half-warp size
template < unsigned ChunkSz >
__global__ void QR_gcd_monolithic_kernel(unsigned *g_R, unsigned *g_BlockIO,
      const unsigned *g_U) {
// const parameters are likely to be placed to constant mem by the compiler

    extern __shared__ unsigned shared[];
//! bidx_x indexes over moduli
//! bidx_y - batch index 
    unsigned thid = threadIdx.x, bidx_x = blockIdx.x, bidx_y = blockIdx.y;
    unsigned *r = shared, *mods = dev_const_mem;
    // memory consumption: BlockSz + 1 + 2*ChunkSz + 11 = 4*ChunkSz + 12
    const unsigned BlockSz = ChunkSz * 2;
    // L needs BlockSz words, Lins --> Ldets needs ChunkSz words
    unsigned *L = r + thid, *Lins, *Ldets, *g_row,
            *break_i, *stop_i, *count_i, *shift_f, *misc;
    // g_row: represents the top generator row: [a, b, c, d] - 4 elements

    __shmem_setup< ChunkSz >(r, Lins, Ldets, g_row, break_i, stop_i, count_i,
        shift_f, misc);

    const unsigned nu = mods[NU], nv = mods[NV],
            max_nu = mods[MAX_NU], max_nv = mods[MAX_NV],
            n_blocks = mods[N_BLOCKS]; // batch size per one modulus

    unsigned block_ofs = UMUL(bidx_x, max_nu + max_nv),
             outs_ofs = UMUL(bidx_x, UMUL(gridDim.y + 3, ChunkSz + UMUL_PAD));
    const unsigned *f_in = g_U + block_ofs, *g_in = f_in + max_nu,
        m_stride = 4; // const memory stride

    g_BlockIO += outs_ofs;
    unsigned *g_Head = g_BlockIO, *g_Outs = g_Head + ChunkSz*4;

    mods += UMUL(bidx_x, UMUL_PAD + m_stride) + ConstMemDataSz;
    volatile unsigned m = mods[0];
    volatile fp_limb invk, invm;

#if CUMP_USE_32BIT_MODULI_SET
    invk = __hiloint2double(mods[3], mods[2]);
#else
    invk = __int_as_float(mods[2]);
    invm = __int_as_float(mods[3]);
#endif
    // fool the compiler: 24-bit mul masks out the highest bit
    volatile unsigned _2e23 = 0x4b000000 + UMUL(thid, UMUL_PAD),
#if CUMP_USE_32BIT_MODULI_SET
        mx100 = 0;
#else
        mx100 = UMUL(100, m);
#endif
    uint4 G, G0; // represents four generator columns (a,b,c,d)
    G.x = 0, G.y = 0;

    unsigned read_ofs = thid + 1 - ChunkSz, t;

    if(thid < ChunkSz) { // only lower threads load data
        t = read_ofs + nu;
        if((int)t >= 0) {
            G.x = f_in[t];
        }
        t = read_ofs + nv;
        if((int)t >= 0) {
            G.y = g_in[t];
        }
    }

    read_ofs = thid + 1 - UMUL(bidx_y, BlockSz + UMUL_PAD) - ChunkSz*3;

    G0.x = 0, G0.y = 0;
    t = read_ofs + nu;
    if((int)t >= 0) {
        G0.x = f_in[t];
    }
    t = read_ofs + nv;
    if((int)t >= 0) {
        G0.y = g_in[t];
    }

    if(thid >= ChunkSz) { // load upper register set
        G.x = G0.x, G.y = G0.y;
    }
    G.z = 1; // set G.z to 1 to collect the product of 'dets' by thread 0

    if(thid == 0) {
        count_i[0] = -1u;
        break_i[0] = 0; // flag indicating that we finish the iteration
    }

    int di = -BlockSz, j;
    if(thid < ChunkSz) {
        di = ChunkSz - 1 - thid; // lower threads work
    }

    // TODO: before going on with the full iterations,
    // test throughly the inner loop for ``lite'' iterations, to make
    // sure it terminates correctly in all cases..

    __main_loop_internal< ChunkSz, false >(G, G0, di, j, nv,
        g_Head, g_Outs, r, thid, m, invk, invm, mx100);

//     if(thid == 0)
//         misc[0] = G.z;
// 
//     CU_SYNC

    if(bidx_y == 0) {
//        g_R[block_ofs + thid] = Lins[thid];
        g_R[block_ofs + BlockSz-1 - thid] = G.x; // for lower threads
    }

    /** generating bulge sequence & iteration type are completely orthogonal */
    //! legend: G0^up - current update sequence
    //! G0^low (G^low) + G^up - current working set

    //! NOTE in fact you can load the next sequence to G^low (no need for
    //! special flags): if you can exchange them properly

#if 0
    CU_SYNC // ??
    // have to transfer from upper register set to lower ones
    if(thid >= ChunkSz) {
        Ldets[thid - ChunkSz] = G0.x;
        Lins[thid - ChunkSz] = G0.y;
    }
    CU_SYNC

    if(thid < ChunkSz) {
        G.x = Ldets[thid]; // in case we have loaded next full updates
        G.y = Lins[thid];  // need to exchange registers here, ie.:
                           // G0.x = G.x, G0.y = G.y;
    }

    //! NOTE: do not fogret about \c shift_f[0] !!
    {
    // at this point 'j == stop[i]' marks the start of ''full'' iterations
    if(thid == 0) { // share 'det' for all threads
        misc[0] = G.z;
        r[ChunkSz] = 1;
    }

    CU_SYNC
// now things are completely orthogonal: generate new full updates & apply them

    if(thid < ChunkSz) { // lower threads only: load the remaining
                         // update sequence
       
        read_ofs = thid + 1 - ChunkSz;
        __setup_generators(G, f_in, g_in, misc[0], read_ofs, nu, nv, m, invk);
    }
 
    read_ofs = thid + 1 - UMUL(bidx_y, BlockSz + UMUL_PAD) - ChunkSz*3;
    // load data to G0 (which should be free, then copy the upper half to G)
    __setup_generators(G0, f_in, g_in, misc[0], read_ofs, nu, nv, m, invk);
    if(thid >= ChunkSz) {
        G.z = G0.z, G.w = G0.w;
    }

// here: G^low: current update sequence
// G0^low - working set first part
// G^up - working set second part

    j = stop_i[0]; 
    if(thid < ChunkSz)
        di = ChunkSz - 1 - thid - j; 
    else
        di = -BlockSz; // skip upper threads

    // run the remaining full iterations
    __GJG_update< ChunkSz, true, false >(G, t, di, j, ChunkSz, G0, r,
             Ldets, Lins, g_row, shift_f, thid, m, invk, invm, mx100);

    if(thid >= ChunkSz)
        j = L[-ChunkSz]; // protect Louts from modification *upper half*

    CU_SYNC
    if(thid < ChunkSz) {
        Ldets[thid] = G.x;
        Lins[thid] = G.y;
        L[0] = G.z, L[ChunkSz] = G.w;
    } 
    CU_SYNC

    if(thid >= ChunkSz) {
        G0.x = Ldets[thid - ChunkSz];
        G0.y = Lins[thid - ChunkSz];
        G0.z = L[-ChunkSz], G0.w = L[0];
    } else
        G = G0;

    CU_SYNC // NOTE: use memory in such a way to avoid CU_SYNC ?

    if(thid < ChunkSz) {
        Ldets[ChunkSz - 1 - thid] = t; // save in reversed order
    } else { // ChunkSz ofs because Louts are for upper threads only
        Lins[thid - ChunkSz] = j;
    }

    if(thid == 0) // setup shift-through flag
        shift_f[0] = (bidx_y == count_i[0]);

    CU_SYNC

    j = stop_i[0];
    di = BlockSz - 1 - thid - j;
    __GJG_update< ChunkSz, true, true >(G, t, di, j, ChunkSz, G0, r,
             Ldets, Lins, g_row, shift_f, thid, m, invk, invm, mx100);
    }

    // G0^up holds the current updating sequence
    if(bidx_y == 1) {
//        g_R[block_ofs + thid] = Lins[thid];
        g_R[block_ofs + BlockSz-1 - thid] = G.x; // for lower threads
    }
#endif

#if 0

    if(start_i[0] < nv) {
 
//! NOTE: supposedly it multiplies by 'res' 2 times, hence
//! either way you cannot get the correct output..
//! TASK: probably a good solution to detach ``lite'' iterations from the
//! full ones completely - this way we can reduce the shared memory transfer
//! and register usage

        // switch to full iterations
        if(j == stop_i[0] || start_i[0] == nv) {
            if(thid == 0) // share 'det' for all threads
                misc[0] = G.z;

            CU_SYNC
            if(thid < ChunkSz) {
                read_ofs = thid + 1 - ChunkSz;
                __setup_generators(G, f_in, g_in, misc[0], read_ofs,
                     nu, nv, m, invk);
            }

            read_ofs = thid + 1 - UMUL(bidx_y, BlockSz + UMUL_PAD) - ChunkSz*3;
            // load data to G0 (which should be free, then copy the upper half
            // to G)
            __setup_generators(G0, f_in, g_in, misc[0], read_ofs, nu,
                     nv, m, invk);
            if(thid >= ChunkSz) {
                G.z = G0.z, G.w = G0.w;
            }
        }
    }

    if((int)stop_i[0] < ChunkSz) // run the remaining ChunkSz - j iters.
        __GJG_update< ChunkSz, true, false >(G, t, di, j, ChunkSz, G0, r,
             Ldets, Lins, g_row, shift_f, thid, m, invk, invm, mx100);

    if(count_i[0] == 1)
        break;

    if(thid < ChunkSz)
        g_row -= 4;

    if(thid >= ChunkSz)
        j = L[-ChunkSz]; // protect Louts from modifications *upper half*

// moving data: G^low -> G0^up (bulge update sequence)
// G0^low -> G^low (working set of this block)
// hence, originally, data should be loaded to G0 ??

    CU_SYNC
    if(thid < ChunkSz) {
        Ldets[thid] = G.x;
        Lins[thid] = G.y;
        L[0] = G.z, L[ChunkSz] = G.w;
    } 
    CU_SYNC

    if(thid >= ChunkSz) {
        G0.x = Ldets[thid - ChunkSz];
        G0.y = Lins[thid - ChunkSz];
        G0.z = L[-ChunkSz], G0.w = L[0];
    } else
        G = G0;

    CU_SYNC // NOTE: use memory in such a way to avoid CU_SYNC ?

    if(thid < ChunkSz) {
        Ldets[ChunkSz - 1 - thid] = t; // save in reversed order
    } else { // ChunkSz ofs because Louts are for upper threads only
        Lins[thid - ChunkSz] = j;
    }

    if(thid == 0) // setup shift-through flag
        shift_f[0] = (bidx_y == count_i[0]);

    CU_SYNC

    di = BlockSz - 1 - thid, j = 0;
    if(start_i[0] < nv)
        __GJG_update< ChunkSz, false, true >(G, t, di, j, stop_i[0], G0, r,
             Ldets, Lins, g_row, shift_f, thid, m, invk, invm, mx100);
    
    if((int)stop_i[0] < ChunkSz) // run the remaining ChunkSz - j iters.
        __GJG_update< ChunkSz, true, true >(G, t, di, j, ChunkSz, G0, r,
             Ldets, Lins, g_row, shift_f, thid, m, invk, invm, mx100);

//     if(bidx_y == 0) {
//         g_R[block_ofs + BlockSz-1 - thid] = G.x;
//     }
//     return;

    if(thid < ChunkSz) {
        G0 = G; // save G to protect from modifications
    }
    // copy L to Lins from the previous block
    // loads lower part of G of the first block to the same register of this
    // block
    __load_outs< ChunkSz >(G, Lins, L, g_Outs, g_Head, bidx_x, bidx_y,
                count_i, thid);
    // G0 upper half: current updating sequence
    // G0 lower half: the next updating sequence to be generated

    if(bidx_y == count_i[0]) { // leading thread block uploads its results
                               // to global memory
        return;
    }

    if(thid == 0) {  // NOTE: check if we need this her..
        r[ChunkSz] = 1; //! set the ``magic'' symbol to avoid zeroing
                    //! out the leading element
    }
    CU_SYNC

// G0: upper: leading rows of the current generator (needed for update)
// G0: lower: leading rows of the (next) generator
// G - current data with upper threads need to be updated using Lins

    if(thid < ChunkSz) {
        di = ChunkSz - 1 - thid;
    } else {
        di = thid - (BlockSz - 1); // UseIns for upper threads
    }

    }  /** ************* end of the main loop *********************/
#endif


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

#endif // _QR_GCD_MONOLITHIC_KERNEL_CU_

