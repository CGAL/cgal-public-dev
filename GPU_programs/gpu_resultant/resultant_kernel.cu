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

#ifndef _RESULTANTS_KERNEL_CU_
#define _RESULTANTS_KERNEL_CU_

#define RUN_RESULTANT_BLOCK     1

#include <include/macros.h>
#include <zmod_dev.cu>

// one array is for global offsets per each modulus
// another array is to count # of failed points
__device__ unsigned dev_mod_index[CUMP_DEVICE_MODULI_SZ / CUMP_MODULI_STRIDE];

//! type B iterations inner loop
//! result is returned in shared mem
template < unsigned BlockSz >
__device__ __forceinline__
void __resultant_typeB_internal(unsigned a, unsigned c2,
    unsigned d2, unsigned *ACs_base, unsigned *As, unsigned *Cs, unsigned *pmul, unsigned *pb0,
    unsigned *pd0, unsigned nu, volatile unsigned m, volatile fp_limb invk,
    volatile fp_limb invm, volatile unsigned mx100) {

    unsigned thid = threadIdx.x;
    //! second type iterations (nu times)
    unsigned lc = 1, t1, a0 = 1, b0 = m - 1, c0 = 1, d0;
    unsigned b;

    if(BlockSz > WS)
        CU_SYNC

    // nu elements should be remapped
    if(thid < nu) {
        As[0] = a, Cs[0] = c2,
        pmul[thid] = d2;
    }

    if(BlockSz > WS)
        CU_SYNC

    volatile unsigned i = nu - 1;

    if(thid < nu) { // reverse element order
        b = i - thid * 2;
        a = As[b], c2 = Cs[b], d2 = pmul[i - thid];
    } else
        pmul[thid] = 1; // pad with identities for reduction

    b = 0;
    if(thid == i) { // nu - 1?
        a = a0, b = b0, pd0[0] = d2;
    }

    if(BlockSz > WS)
        CU_SYNC

    // here BlockSz == # of threads
    // iterate while i < BlockSz/2

/*    unsigned niters = nu;*/
    // after nu - BlockSz/2 iterations we are left to perform
    // nu - (nu - BlockSz/2) = BlockSz/2 iterations => all BlockSz threads
    // are occupied

// NOTE: if nu == WS (initially) -> you can skip first type iters completely
    if(BlockSz == WS || nu > WS)
    while(1) {

    d0 = pd0[0]; // possible hazard on the first iteration (but very unlikely)

    if(BlockSz > WS)
        CU_SYNC

    // valid index range: i+1 .. nu + nv - 1
    if(thid <= i) {

        t1 = sub_mul_reduce_m(a, c0, b, d0, m, invk, invm, mx100);
        As[0] = mul_m(t1, a0, m, invk);

        t1 = sub_mul_reduce_m(b, a0, a, b0, m, invk, invm, mx100);
        b = mul_m(t1, lc, m, invk);

        t1 = sub_mul_reduce_m(c2, a0, d2, b0, m, invk, invm, mx100);
        Cs[0] = mul_m(t1, lc, m, invk);

        t1 = sub_mul_reduce_m(d2, c0, c2, d0, m, invk, invm, mx100);
        d2 = mul_m(t1, a0, m, invk);
    }

    if(thid == 0) {
// TODO you can multiply and even store in shared mem inside
        t1 = mul_m(lc, lc, m, invk);
        lc = mul_m(a0, t1, m, invk); // la == a0
        pmul[i] = lc;
    }
    i--;

    if(thid == i) { // next thread saves b & d because previous values are 0s
        pb0[0] = b, pd0[0] = d2;
    }

    if(BlockSz > WS)
        CU_SYNC

// TODO: can quit earlier because element remapping occurs anyway
    if(thid <= i) { // read previous elements
        a = As[1];
        c2 = Cs[1]; // beware of the read-write hazard on the next iteration

        unsigned *pac0 = ACs_base + i + 1; //ACs_base + i - 1;
        a0 = pac0[0], c0 = pac0[BlockSz];
        b0 = pb0[0], lc = pmul[i + 1]; // read out lc directly from pmul
    } /*else if(thid == 0)
        a0 = ACs_base[i - 1];*/

    if(BlockSz == WS) {
        if(i == -1u)
            break;
    } else if(i == WS - 1) // switch to iters without sync
        break;

    } // while(1)

    if(BlockSz > WS)
        CU_SYNC // need this sync ??
/**
    unsigned *ACs_base = denom + 1, *As = ACs_base + thid, *Cs = As + BlockSz,
        *pmul = ACs_base + BlockSz * 2;
*/

//! iterations without sync (1st warp only)
    if(thid < WS && BlockSz > WS)
    while(1) {

    d0 = pd0[0]; // possible hazard on the first iteration (but very unlikely)
    // valid index range: i+1 .. nu + nv - 1
    if(thid <= i) {

        t1 = sub_mul_reduce_m(a, c0, b, d0, m, invk, invm, mx100);
        As[0] = mul_m(t1, a0, m, invk);

        t1 = sub_mul_reduce_m(b, a0, a, b0, m, invk, invm, mx100);
        b = mul_m(t1, lc, m, invk);

        t1 = sub_mul_reduce_m(c2, a0, d2, b0, m, invk, invm, mx100);
        Cs[0] = mul_m(t1, lc, m, invk);

        t1 = sub_mul_reduce_m(d2, c0, c2, d0, m, invk, invm, mx100);
        d2 = mul_m(t1, a0, m, invk);
    }

    if(thid == 0) {
        t1 = mul_m(lc, lc, m, invk);
        lc = mul_m(a0, t1, m, invk); // la == a0
        pmul[i] = lc;
    }
    i--;

    if(thid == i) { // next thread saves b & d because previous values are 0s
        pb0[0] = b, pd0[0] = d2;
    }

    if(i == -1u)
        break;

    if(thid <= i) { // read previous elements
        a = As[1];
        c2 = Cs[1]; // beware of the read-write hazard on the next iteration

        unsigned *pac0 = ACs_base + i + 1; //ACs_base + i - 1;
        a0 = pac0[0], c0 = pac0[BlockSz];
        b0 = pb0[0], lc = pmul[i + 1]; // read out lc directly from pmul
    }

    } // while(1)

#if 0
    // here we have BlockSz/2 upper threads working
    if(thid >= BlockSz/2)
    // NOTE: make sure As does not get overwritten
        ACs_base[thid - BlockSz/2] = d;

    // first half work with a/b, second half work with d/c

    // hence d & c should be passed to the second group of threads
    // c is already in shared mem

    // As is partly not used
    // nope: as is not used completely as long as no iterations made
    // NOTE: BlockSz/2 <= nu (by default)
    if(thid >= niters) {   }
#endif
}

//! \c x_deg1/2 - maximal degree of a polynomial in X
//! TODO: pass \c nu/nv also as parameters ?
//! \c g_U layout:
//! f[0][x_deg1], f[1][x_deg1],.. f[nu][x_deg1] - first lane
//! f[0][x_deg1-1], f[1][x_deg1-1],.. f[nu][x_deg1-1] - second lane
//! ...
//! f[0][0], f[1][0],.. f[nu][0] - last lane
//! total of (x_deg1 + 1) * 64 values (with padding)

//! g[0][x_deg2], g[1][x_deg2],.. g[nv][x_deg2]
//! g[0][x_deg2-1], g[1][x_deg2-1],.. g[nv][x_deg2-1]
//! ...
//! g[0][0], g[1][0],.. g[nv][0] 
//! total of (x_deg2 + 1) * 64 values (with padding)

//! memory consumption: BlockSz * 3 + 16 (BlockSz = WS*2, WS*3, WS*4)
//! (BlockSz * 3 + 16) * 2 for BlockSz = WS
//! supported block sizes: 32, 64, 96, 128
//! when BlockSz = 32 two resultants by two warps are computed
//! separately
//! \c out_data_padding - resultants for each modulus are stored with
//! padding for better memory coalescing
//! \c g_R - saves resultants, \c g_R2 - saves denominators
template < unsigned BlockSz >
__global__ void CUMP_LAUNCH_BOUNDS(128, 8)
resultant_block_kernel(unsigned *g_R, unsigned *g_R2, const unsigned *g_U,
        unsigned x_deg1, unsigned x_deg2, unsigned out_data_padding) {

#if RUN_RESULTANT_BLOCK
    extern __shared__ unsigned shared[];

// TODO TODO: need more efficient thread work distribution:
// when nu and nv are unbalanced
// currently f[0][x]..f[nu][x] are stored with padding max_nr/2
// which is the same as BlockSz; as well as g[0][x]..g[nv][x] is stored with
// padding max_nr/2

// change this to:
// + single lane of size BlockSz*2 storing
// f[0][x]..f[nu][x] g[0][x]..g[nv][x]: better for unbalanced 'nu' and 'nv'
// - requires padding for unbalanced 'degx1' and 'degx2': that is choosing
// degx = max(degx1, degx2)
// + more effiective for small nu and nv (boundary case) because during the
// first run we evaluate BlockSz points and during the second run can
// restrict # of threads to evaluate the rest points (while in the current
// implementation excess threads work with "garbage")

    unsigned thid = threadIdx.x, bidx_x = blockIdx.x, bidx_y = blockIdx.y;
    unsigned split = (BlockSz == WS ? threadIdx.y : 0);

//! bidx_x: indexes over evaluation points
//! bidx_y: indexes over moduli

    unsigned *r = shared;
    if(BlockSz == WS) // separate shared mem space transparently
        r += UMUL(split, BlockSz * 3 + 16);

    const unsigned m_stride = 4; // constant memory stride
    unsigned *mods = dev_const_mem;

    volatile unsigned nu = mods[0], nv = mods[1];
    mods += 2 + UMUL(bidx_y, UMUL_PAD + m_stride);
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
    volatile unsigned _2e23 = 0x4b000000 + UMUL(thid, UMUL_PAD);

    // generator G(a, b); generator B(c, d)
    // a = (x_deg1 + 1) * max_nr/2; b = (x_deg2 + 1) * max_nr/2
    unsigned a = UMUL(x_deg1, UMUL_PAD + BlockSz) + BlockSz,
         b = UMUL(x_deg2, UMUL_PAD + BlockSz) + BlockSz;
    unsigned c, d, c2, d2;

    unsigned ofs = UMUL(bidx_y, a + b) + thid;
    //! BlockSz * (x_deg1 + xdeg2 + 2) elements per modulus with zero padding
    //! BlockSz == max_nr / 2

    volatile unsigned i;
#if 1
    if(BlockSz > WS)
        c2 = bidx_x + 1; //gridDim.x - bidx_x;
    else
        c2 = UMUL(bidx_x, UMUL_PAD + 2) + split + 1;
            //UMUL(gridDim.x - bidx_x, UMUL_PAD + 2) - split;

    TEX_COLUMN_FETCH(c, g_U, ofs, 0)
    for(i = BlockSz; i < a; i += BlockSz) { // also serves as an offset
        c = mul_small_m(c, c2, m, invk);
        TEX_COLUMN_FETCH(d, g_U, ofs, i)
        c = add_m(c, d, m);
    }

    ofs += a;
    TEX_COLUMN_FETCH(d, g_U, ofs, 0)
    for(i = BlockSz; i < b; i += BlockSz) { // also serves as an offset
        d = mul_small_m(d, c2, m, invk);
        TEX_COLUMN_FETCH(d2, g_U, ofs, i)
        d = add_m(d, d2, m);
    }
#else
    d2 = umin(a, b)*2;

    c = g_U[ofs], d = g_U[ofs + BlockSz/2]; // prefetch 2 at a time
    for(i = BlockSz; i < d2; i += BlockSz) {

        c2 = g_U[ofs + i];
        a = g_U[ofs + i + BlockSz/2];

        c = mul_small_m(c, bidx_x, m, invk);
        d = mul_small_m(d, bidx_x, m, invk);

        c = add_m(c, c2, m);
        d = add_m(d, a, m);
    }
#endif
    
    if(BlockSz > WS)
        ofs = UMUL(bidx_y, out_data_padding) + bidx_x;
    else    
        ofs = UMUL(bidx_y, out_data_padding) + UMUL(bidx_x, UMUL_PAD + 2)
             + split;
       
    //! r layout: pu[nu], pv[nv], pv[nv]^-1, pu[nu]^-1 - different order !!
    if(thid == nv) {
        r[1] = d; // pv[nv]
    }
    if(thid == nu) { // save leading elements
        r[0] = c; // pu[nu]
    }

    if(thid == nu) {
        volatile unsigned mu = mods[1];
        r[3] = c;
        if(c != 0)
            r[3] = montgomery_inverse(c, m, mu);
    }

    if(BlockSz > WS)
        CU_SYNC

//     if(thid == nu) {
//         g_R[ofs] = mul_m(c, r[3], m, invk);
//     }
//     return;

    //! if leading coefficients are zero: return with zero-denominator
    //! correct this ...
    if(r[1] == 0 || r[3] == 0) {
        if(thid == 0) {
//             g_R[ofs] = 0x77776;
            g_R2[ofs] = 0; 
        }
        return;
    }

    i = BlockSz - nv, b = nu - thid;
    // values pu/pv are now in registers - we multiply directly by
    // inverses before saving to shared mem
    
    unsigned *pu = r + 4, *pv = pu + nu;
    if(thid < nu) {
        pu[thid] = mul_m(c, r[3], m, invk);
    }

//     if(thid - 1 < nv) { // valid thids: 1..nv inclusive
//         // d's are stored in reversed order
//         pv[nv - thid] = d;//mul_m(d, r[2], m, invk);
//     }
    if(thid < nv) { // except nv-th this
        pv[thid] = d;  //mul_m(d, r[2], m, invk);
    }

    if(BlockSz > WS)
        CU_SYNC

    c = (thid == i); // 0th thid if nv == WS*2 (max)
    d = c, a = 0, c2 = 0;

    if(d == 1)
        d = r[1];

// NOTE: here for thid > i pu[i + b] can go out of index range if nu < nv..
    if(thid > i) { // the first thid == nv has 1 assigned
        if((int)(i + b) >= 0) // BlockSz - nv + nu - thid (thid > BlockSz - nv)
            c = pu[i + b];
        d = pv[BlockSz - thid]; // it's easy for d as we load exactly nv values
    }
    if((int)(b - nv) >= 0) { // load the rest nu elements
        c2 = pu[b - nv];  //! b = nu - thid here !!
    }

    if(thid == 0) // process it additionally
        d2 = sub_m(pv[0], 1, m);

    else if((int)b > 0)  // except zero thid
        d2 = sub_m(0, pu[b], m);

    // r shmem layout (two first elements: for leading coeffs)
    // 0  1  2  3  4
    // b0 d0 lc det denom (la == a0)
    unsigned *pb0 = r + 2, *pd0 = pb0 + 1, *plc = pd0 + 1, *det = plc + 1,
            *denom = det + 1;

    //! do not forget that all variables already indexed by \c thid
    unsigned *ACs_base = denom + 1, *As = ACs_base + thid, *Cs = As + BlockSz,
        *pmul = ACs_base + BlockSz * 2;

    if(thid == i) {
        pd0[0] = d; // keep d0, although we now that d0 == 1 at the beginning
    }

    //! first type iterations (nv times)
    while(1) {

    if(BlockSz > WS) // no need to sync if single warp works
        CU_SYNC

    b = pd0[0]; // read out d0 from shared mem
    if(thid >= i) {
        As[0] = c;
        d = submul_m(d, b, c, m, invk);
    }
    // all thids work with the second part (any exceptions here ?)
    As[BlockSz] = c2;
    d2 = submul_m(d2, b, c2, m, invk);

    if(thid == BlockSz - i)
        a = b; // keep d0 in a

    i++;
    if(BlockSz > WS)
        CU_SYNC

    if(thid >= i) { // read previous elements
        c = As[-1]; 
    }
    c2 = As[BlockSz - 1];

    if(thid == i)   // no read-write hazard here because in the last iter
        pd0[0] = d; // i == BlockSz => nothing is written
    if(i == BlockSz)   
        break;
    } // while(1)

    __resultant_typeB_internal< BlockSz >(a, c2, d2, ACs_base, As, Cs,
         pmul, pb0,  pd0, nu, m, invk, invm, mx100);

    // need to break on the first i >= nu-32 which is divisible by 32 (warp)

    if(BlockSz > WS)
        CU_SYNC

    // lower warp: Cs = ACs_base + BlockSz (computes det)
    // upper warp: pmul = ACs_base + BlockSz*2 - WS = ACs_base + WS*3
    // (computes denom)

    a = Cs[0];
    if(thid >= nu) // nu elements to be multiplied => pad remaining with 1s
        a = 1;      
    if(thid < nv) { // r[0]: lu; computes lu^nv
        a = mul_m(a, r[0], m, invk);
    }
    Cs[0] = a;

    if(BlockSz > WS)    
        CU_SYNC

    unsigned *data = ACs_base + BlockSz; // == Cs
    volatile unsigned *t = data;

#if 0
    if(thid == 0) {
        a = 1;
        for(int ii = 0; ii < BlockSz; ii++) {
            a = mul_m(a, data[ii], m, invk);
        }
    }
#else

    a = thid & WS-1;
    if(BlockSz == WS*3) {

        if(thid < WS*2) {
            t += UMUL(thid >> 5, UMUL_PAD + BlockSz) + a;
        } else {
            // third warp is split on half-warp basis
            t += WS*2 + UMUL(a >> 4, UMUL_PAD + BlockSz) + (a & HF-1);
        }
    } else if(BlockSz == WS) {
        // use 2 half-warps of each warp to evaluate 4 prefix muls
        t += UMUL(a >> 4, UMUL_PAD + BlockSz) + (a & HF-1);

    } else {
        t += UMUL(thid >> 5, UMUL_PAD + (1 << 6)) + a;
    }

    a = t[0];
    // if BlockSz == WS*3 then thid < WS*2
    // otherwise no conditions
    if(BlockSz > WS)
    if(BlockSz != WS*3 || thid < WS*2) {
        a = mul_m(a, t[WS], m, invk);
        t[0] = a;
    } 

    a = mul_m(a, t[HF], m, invk);
    t[0] = a;
    a = mul_m(a, t[8], m, invk);
    t[0] = a;
    a = mul_m(a, t[4], m, invk);
    t[0] = a;
    a = mul_m(a, t[2], m, invk);
    t[0] = a;
    a = mul_m(a, t[1], m, invk);
    t[0] = a;

    if(BlockSz > WS*2) {

        CU_SYNC

        if((int)thid < 2) {
            // each thread will handle 2 elements with WS*2 stride
            t = data + UMUL(thid, UMUL_PAD + BlockSz);
            a = t[0]; // nElems 64 -> first stride
            a = mul_m(a, t[WS*2], m, invk);
        }
    }
#endif
    //! results in:
    //! BlockSz = WS: det = a (thid: 0); denom = a (thid: HF)
    //! BlockSz = WS*2: det = a (thid: 0); denom = a (thid: WS)
    //! BlockSz = WS*3 or WS*4: det = a (thid: 0); denom = a (thid: 1) !!!

    if(thid == 0) {
        g_R[ofs] = a; // this is resultant
    }
    if(BlockSz == WS) {
        if(thid == HF)
            g_R2[ofs] = a; // this is a denominator

    } else if(BlockSz == WS*2) { // compile-time decision
        if(thid == WS)
            g_R2[ofs] = a; // this is a denominator

    } else if(BlockSz == WS*3 || BlockSz == WS*4) {
        if(thid == 1)
            g_R2[ofs] = a; // this is a denominator
    }
#else
#warning resultant_block_kernel: dummy compilation
#endif // RUN_RESULTANT_BLOCK
}

//! locates "holes" (zeros) in the input sequence \c x
//! memory: 49 words per 2 warps: 49*2 words (for 128 threads) + HF for postscan
//! returns thread offsets of a new (compact) sequence
//! \c scanbase[sum_ofs] stores the # of valid entries
__device__ __forceinline__ unsigned __stream_compact128_internal(
    unsigned *scanbase, unsigned x, unsigned& sum_ofs) {

    unsigned thid = threadIdx.x, thid_in_warp = thid & WS - 1;
    //! using "alternating" scans layout, ie., layouts of two consecutive
    //! warps are placed to the same mem region
    volatile ushort *scan = (volatile ushort *)(scanbase + HF + thid_in_warp +
            UMUL(WS + HF + 1, thid >> 6)) + ((thid >> 5) & 1);

    ushort s = (x != 0);

    scan[-HF*2] = 0; // identity symbol
    scan[0] = s;

    // beware of sign propagation for signed shorts
    s = s + scan[-2], scan[0] = s;
    s = s + scan[-4], scan[0] = s;
    s = s + scan[-8], scan[0] = s;
    s = s + scan[-16], scan[0] = s;
    s = s + scan[-32], scan[0] = s;

    sum_ofs = HF + UMUL(WS + HF + 1, (WS*4) >> 6);
    unsigned *postscan = scanbase + sum_ofs;
    sum_ofs += 3;

    if((thid & WS-1) == WS-1)
        postscan[thid >> 5] = s;

    CU_SYNC

    if(thid < 4) {
        volatile int *temp = (volatile int *)postscan + thid; 
        temp[-4] = 0; // load identity symbol
        unsigned t = temp[0];
        t = t + temp[-1], temp[0] = t;
        t = t + temp[-2], temp[0] = t;
//         t = t + temp[-4], temp[0] = t;
    }

    CU_SYNC

    s = scan[0];
    // postscan[-1] equals to identity
    s = s + postscan[(thid >> 5) - 1];

// NOTE: make sure 'scanbase' does not overlap with the main data storage,
// otherwise sync is required
    return (unsigned)s - 1; // -1 because of inclusive prefix sum
}

//! computes modular inverses, multiplies with data elements and writes back
//! modular resultants
//! # of thids: 128
//! memory: 49*2 + 16 + 4 + 128 words
//! \c g_R - reads resultants, writes out divided by denominators
//! \c g_R2 - reads denominators, writes out evaluation points
//! \c npoints - # of evaluation points (for indexing)
__global__ void mod_inverse_kernel1(unsigned *g_R, unsigned *g_R2,
        unsigned *g_R3, unsigned *g_R4,
        unsigned npoints, unsigned out_data_padding) {

    extern __shared__ unsigned shared[];
//! number of blocks: nmods x padding / 128
    unsigned thid = threadIdx.x, bidx_x = blockIdx.x, bidx_y = blockIdx.y;

//! bidx_x: indexes over evaluation points
//! bidx_y: indexes over moduli

    const unsigned m_stride = 4, block_sz = 128; 
    unsigned *r = shared, *scanbase = r + block_sz, *s_ofs = scanbase + block_sz;

    unsigned *mods = dev_const_mem + 2 + UMUL(bidx_y, UMUL_PAD + m_stride);
    volatile unsigned m = mods[0], mu = mods[1];
    volatile fp_limb invk, invm;

#if CUMP_USE_32BIT_MODULI_SET
    invk = __hiloint2double(mods[3], mods[2]);
#else
    invk = __int_as_float(mods[2]);
    invm = __int_as_float(mods[3]);
#endif

    unsigned block_ofs = UMUL(bidx_x, UMUL_PAD + block_sz);
    unsigned ofs = UMUL(bidx_y, out_data_padding) + thid;

    unsigned a = 0, inva = 0;
    // do not read in garbage data for the last block
    if(block_ofs + thid < npoints) { 
        inva = g_R2[ofs + block_ofs], a = g_R[ofs + block_ofs];
    }

    inva = montgomery_inverse(inva, m, mu);
    a = mul_m(a, inva, m, invk);

    // call stream compaction after montgomery inverse for better parallelism
//! + 128 words for \c r itself (to avoid overlap)
    unsigned sum_ofs, s;
    s = __stream_compact128_internal(scanbase, inva, sum_ofs);
    unsigned data_sz = scanbase[sum_ofs];

    CU_SYNC

    //a = a | (thid << 24); // use 7 MSB for local thread offset
    if(inva != 0) { // save only non-zero elements
        r[s] = a;
        scanbase[s] = thid;
    }

    CU_SYNC

    a = r[thid];
    unsigned xs = scanbase[thid] + 1 + block_ofs;
//     unsigned xs = (a >> 24) + 1 + block_ofs;
//     a &= 0xffffff;

    // # of points to be saved (can be < block_sz for the last block)
    if(thid == 0) {
#if CUMP_USE_ATOMICS
        // read in a block write offset from the global variable
        s_ofs[0] = __uAtomicAdd((unsigned *)dev_mod_index + bidx_y, data_sz);
#endif
    }

    CU_SYNC 
#if CUMP_USE_ATOMICS
    ofs += s_ofs[0];

//     ofs += block_ofs;
//     a = scanbase[0];
//     unsigned bofs = bidx_x * block_sz;
//     if(bidx_x == gridDim.x - 1)
//         bofs = 0;
//     else
//         bofs = npoints - (bidx_x+1) * block_sz;
// 
//     ofs += bofs;
#else
    ofs += block_ofs;
#endif

//!NOTE NOTE: if you change the order of evaluation points you *must*
//! divide by denominator otherwise GPU-host results will differ !!
    if(thid < data_sz) {
        g_R3[ofs] = a;
        g_R4[ofs] = xs;
    }
}

#endif // _RESULTANTS_KERNEL_CU_
