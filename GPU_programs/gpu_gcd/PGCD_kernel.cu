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
 
#ifndef _PGCD_KERNEL_CU_
#define _PGCD_KERNEL_CU_

#include <zmod_dev.cu>      // device modular arithmetic includes

//! dividend is stored in shared mem \c L 
//! divisor in registers \c g , leading coeff of the divisor: \c lcG[0]
//! thread 'nv' is responsible for computing the powers of \c lcG[0]
//! which are placed in \c L[i_exp]
//! at each iteration \c nu is decreased by the amount specified by \c skip[0]
//! \c L[nu] is the leading coeff of the dividend
//! \c nu = -1 indicates that gcd is computed (in registers \c g )
template < bool UseSync >
__device__ __forceinline__ void __pgcd_lite_internal(unsigned *L, unsigned *lcG,
        unsigned *skip, unsigned g, const unsigned i_exp, unsigned& nu,
        unsigned nv, unsigned thid, unsigned m, fp_limb invk,
        fp_limb invm, volatile unsigned mx100) {

    unsigned t = nu - nv + thid;
    if(thid == nv) {
        t = i_exp; // points to special entry
        lcG[0] = g; // share the leading element btw all threads
        L[t] = g;  // L[-1] stores d^e
        g = 0;
        skip[0] = 1;
    }
    CU_SYNC

    unsigned lcg, lcf;
    unsigned j = 0;
    while(1) {

        lcg = lcG[0], lcf = L[nu];
        if(thid < skip[0]) {
           lcg = L[i_exp];
        }

        if(UseSync)
            CU_SYNC

        if(thid <= nv) {
            unsigned f = L[t];
            L[t] = sub_mul_reduce_mX(f, lcg, lcf, g, m,
                invk, invm, mx100);
        }

        if(UseSync)
            CU_SYNC

        if(thid == 0) {
            // search for the 1st non-zero leading element
            unsigned s = nu - 1;
            while((int)s >= 0 && L[s] == 0) {
                s--;
            }
            skip[0] = nu - s; // count # of elements skipped
        }
        if(UseSync)
            CU_SYNC

        // NOTE NOTE: you can do iterations without zero-checks
        // but then you need to check after each pseudo-div
        // if the dividend vanishes completely: can be done using
        // vote intrinsics
        nu -= skip[0];
//         j++;
//         if(j == 47)
//             return;

        if(thid < nv) {
            t -= skip[0];
        }
        if((int)nu < (int)nv)
            break;

    }
    if((int)t >= 0 && thid < skip[0] && thid < nv) {
        L[t] = mul_m(L[t], lcg, m, invk);
    }
}

__device__ __forceinline__ void __pgcd_lite_loop(unsigned *L, unsigned *lcG, unsigned *skip,
        unsigned& g, unsigned& nu, unsigned& nv,
        unsigned thid, unsigned m, fp_limb invk,
        fp_limb invm, volatile unsigned mx100) {

    unsigned k = 0;

    while(1) { 
        const unsigned i_exp = -1u;

        //! the dividend is contained in L[0..nu] where nu < nv
        //! NOTE: no need to sync after the call
//         if(k == 5) mx100 = 1;

        __pgcd_lite_internal< true >(L, lcG, skip, g, i_exp, nu, nv, thid, m,
                invk, invm, mx100);

        if((int)nu <= 0) // nu == -1 indicates that gcd is computed
            break;

//         if(++k == 1)
//             break;

        // swap the dividend and divisor to continue iterations:
        unsigned t = g;
        if(thid <= nu) {
            g = L[thid];        
        }

        CU_SYNC
        if(thid <= nv) {
            L[thid] = t;
        }
        t = nv, nv = nu, nu = t;
    }
    if(nu == 0) {
        nv = 0;
        g = 1;
    }
}

//! default implementation of PGCD algorithm
template < bool SplitThreads >
__global__ void /*CUMP_LAUNCH_BOUNDS(512, 4)*/
PGCD_lite_kernel(const unsigned *In0, unsigned *Out0) {

    extern __shared__ unsigned shared[];
//! bidx_x - moduli index
//! bidx_y - index within a batch for single modulus
    unsigned thid = threadIdx.x, bidx_x = blockIdx.x;
    unsigned *r = shared, *mods = dev_const_mem;
   // leading coeffs of the dividend and divisor
    // L[-1] is used to store d^e
    unsigned *lcG = r, *skip = lcG + 1, *L = skip + 2;

    if(SplitThreads)
        bidx_x = UMUL(bidx_x, blockDim.y) + threadIdx.y;

    const unsigned m_stride = 4; // const mem  stride
    /*const*/ unsigned nu = dev_const_mem[NU], nv = dev_const_mem[NV];
    mods += UMUL(bidx_x, UMUL_PAD + m_stride) + ConstMemDataSz;
    /*volatile*/ unsigned m = mods[0];
    volatile fp_limb invk, invm;
    volatile unsigned mx100;

    MODARITHM_INIT(thid)

    unsigned max_nu = dev_const_mem[MAX_NU],
            block_ofs = UMUL(bidx_x, dev_const_mem[DATA_PAD]);
//     const unsigned *In0 = (unsigned *)dev_const_mem[DATA_IN],
    const unsigned *f_in = In0 + block_ofs, *g_in = f_in + max_nu;

    const unsigned nu_ofs4 = dev_const_mem[NU_OFS4],
                   nv_ofs4 = dev_const_mem[NV_OFS4];
    f_in += nu_ofs4, g_in += nv_ofs4;

    // chunk_size = max_nv (that is nv+1 aligned by 16 boundary)
    // nu - degree of the dividend, nv - degree of the divisor
    unsigned t, g, max_nv = blockDim.x;

    // load a whole dividend to shared mem by chunks of size max_nv
    t = 0;
    while(t <= nu) {
        if(t + thid <= nu) {
            L[t + thid] = f_in[t + thid];
        }
        t += max_nv;
    }

    if(thid <= nv) {
        g = g_in[thid]; // # of threads == nv
    }

    __pgcd_lite_loop(L, lcG, skip, g, nu, nv, thid, m, invk, invm, mx100);

    block_ofs = UMUL(bidx_x, dev_const_mem[GCD_SZ]) +
                         dev_const_mem[MODS_PAD];
//     unsigned *Out0 = (unsigned *)dev_const_mem[DATA_OUT];

    if(thid == 0)
        Out0[bidx_x] = nv + 1; // size of a gcd (not the degree!!)

    if(thid <= nv)
        (Out0 + block_ofs)[thid] = g;//L[thid];

}

//! ``lite'' version to run the remaining iterations using using all threads
//! F - new dividend of size 'nu', G - new divisor of size 'nv'
__device__ __forceinline__ void __pgcd_quad_lite(unsigned *L, unsigned *lcF, unsigned *lcG, unsigned *Out0,
    uint4 F, uint4 G, unsigned nu, unsigned nv, unsigned last_thid,
     unsigned thid, unsigned FIRST_THID, unsigned bidx_x, unsigned m,
     volatile fp_limb invk, volatile fp_limb invm, volatile unsigned mx100) {

    unsigned ofs = thid * 4, t, g;
    // NOTE NOTE: you should be careful here because thid and last_thid
    // can be shifted (need to shift them back)
    if(thid <= last_thid) { //! save the divisor in registers
        L[ofs] = G.x;       // 4-way bank conflicts.. grr
        L[ofs + 1] = G.y;
        L[ofs + 2] = G.z;
        L[ofs + 3] = G.w;
    }
    if(thid == last_thid) {
        L[ofs + 4] = lcG[0];
    }
    CU_SYNC
    t = (last_thid + 1) * 4 - nv;
    if(thid + FIRST_THID <= nv)
        g = L[t + thid + FIRST_THID];

    CU_SYNC     //! save the dividend in shared mem
    if(thid <= last_thid) { // NOTE NOTE: how are you sure that shmem offsets
        L[ofs] = F.x;       // are correct ??
        L[ofs + 1] = F.y;
        L[ofs + 2] = F.z;
        L[ofs + 3] = F.w;
    }
    if(thid == last_thid) {
        L[ofs + 4] = lcF[0];
    }
    t = (last_thid + 1) * 4 - nu;
    L += t;
    // [F, lcF] needs to be put in shared mem while [G, lcG] goes into
    // register space
    thid += FIRST_THID;
    
    unsigned *skip = lcF;
    CU_SYNC

    __pgcd_lite_loop(L, lcG, skip, g, nu, nv, thid, m, invk, invm, mx100);

    unsigned block_ofs = UMUL(bidx_x, dev_const_mem[GCD_SZ]) +
                         dev_const_mem[MODS_PAD];
//     unsigned *Out0 = (unsigned *)dev_const_mem[DATA_OUT];

//     ofs = (last_thid + 1) * 4 - nu;

    if(thid == 0)
        Out0[bidx_x] = nv + 1; // size of a gcd (not the degree!!)

    if(thid <= nv)
        (Out0 + block_ofs)[thid] = g;//nv+1;//L[thid];
}

template < bool FirstRun >
__device__ __forceinline__ void __pgcd_quad_internal(unsigned *L, unsigned *lcF, unsigned *lcG,
    uint4& F, uint4 G, unsigned& nu, unsigned nv, unsigned last_thid,
    unsigned *cache, const unsigned CacheLn, unsigned ofs,
    const unsigned *f_in, unsigned thid, unsigned m, volatile fp_limb invk,
        volatile fp_limb invm, volatile unsigned mx100) {

    unsigned j = 0;
        
    while(1) {
        

    CU_SYNC

    //! G.w for thid = nthids-1 is g[nv - 1] - next after the lcoeff
    //! F.w for thid = nthids-1 is f[nu - 1] where 'nu' is the current degree
    //! of the dividend
    if(thid <= last_thid) {
        
        if(lcF[0] != 0) {

        unsigned lcg = lcG[0], lcf = lcF[0]; // leading element

        if(FirstRun) {
        L[thid] = sub_mul_reduce_mX(F.w, lcg, lcf, G.w, m, invk, invm, mx100);
        F.w = sub_mul_reduce_mX(F.z, lcg, lcf, G.z, m, invk, invm, mx100);

        // 0th thid uses d^e instead of d
        F.x = sub_mul_reduce_mX(F.x, lcg, lcf, G.x, m, invk, invm, mx100);

        if(thid == 0) { // register selection: no branch statement
            lcg = F.x;
        }
        F.z = sub_mul_reduce_mX(F.y, lcg, lcf, G.y, m, invk, invm, mx100);
        F.y = F.x;
        } else { // !FirstRun

        L[thid] = sub_mul_reduce_mX(F.w, lcg, lcf, G.w, m, invk, invm, mx100);
        F.w = sub_mul_reduce_mX(F.z, lcg, lcf, G.z, m, invk, invm, mx100);
        F.z = sub_mul_reduce_mX(F.y, lcg, lcf, G.y, m, invk, invm, mx100);
        F.y = sub_mul_reduce_mX(F.x, lcg, lcf, G.x, m, invk, invm, mx100);
        }

        } else { //! lcF[0] == 0
            if(FirstRun && thid == 0)
                F.y = mul_m(F.y, F.x, m, invk);
            L[thid] = F.w, F.w = F.z, F.z = F.y, F.y = F.x;
        }
    }
//     if(lcF[0] != 0) //! NOTE NOTE debug only
//         ii++;

    CU_SYNC

    if(thid != 0 && thid <= last_thid) {
        F.x = L[thid - 1];
    }

    nu--;
    if((int)nu < (int)nv) {// NOTE beware if nv == 0 !!
        break;
    }
    
    if(thid == last_thid) { // in fact lcr[0] is saved in L[thid]
        lcF[0] = L[thid];
    }
    if(FirstRun && thid == 0) {
        F.y = cache[j];
    }
    j++;

//     if(ii == 1) {// here F.y == F.x for 0th thread
//         break;
//     }

    if(FirstRun && j == CacheLn) { // load another cache line and reset 'j'
        ofs -= CacheLn;
        if(thid < CacheLn) {
            j = 0;
            if((int)ofs >= 0)
                j = f_in[ofs];
            cache[(int)(CacheLn - 1 - thid)] = j;
        }
        j = 0;
    }

    } // while(1)
}

//! scans \c F to find the first non-vanishing leading coeff, modifies \c nu
//! appropriately. Sets \c nu=-1 if \c F vanishes indicating that gcd is ready
__device__ __forceinline__ void __lcf_scan(unsigned *L, unsigned *lcF, uint4& F,
        unsigned& nu, unsigned nv, unsigned last_thid,
            unsigned thid, unsigned FIRST_THID) {

//     unsigned nu_safe = nu;
    if(L[last_thid] != 0)  // leading element is non-zero => return
        goto Lexit;

    CU_SYNC // sync because the last thid also writes to 'L'

//     while(1) { // need to shift data at least once since L[last_thid] == 0
//         if(thid <= last_thid) {
//             L[thid] = F.w;
//             F.w = F.z, F.z = F.y, F.y = F.x, F.x = 0;
//         }
//         CU_SYNC
// 
//         // equivalent to: thid != 0 && thid <= last_thid
//         if(thid-1 < last_thid) {
//             F.x = L[thid - 1];
//         }
// 
//         nu--;
//         if((int)nu < 0) {
//             return;
//         }
// 
//         if(L[last_thid] != 0) {
//             goto Lexit;
//         }
// 
//         CU_SYNC
// 
//         if(nu_safe - nu >= 4) // run several iterations before warp voting
//             break;
//     }
#if 0
    //! ofs == (nv + 1) 4-aligned - nu
    //! ofs = ((nv + 4) & ~3) - (nv - 1);
   if(0)
    { //! FIXME FIXME pay attention that for shifted 'thid' and 'last_thid'
      //! there could be problems with voting..
        unsigned x = 0, t2, nwarps;
        t2 = thid + FIRST_THID;

        if(t2 < WS)
            L[t2] = 0;
        CU_SYNC

        if(thid <= last_thid) {
            // fast check whether all elements of F vanished
            x = (F.x + F.y != 0) || (F.z + F.w != 0);
        }
        x = __any(x);

        CU_SYNC

        if((t2 & WS-1) == 0) // warp-based index
            L[t2 / WS] = x;
        // note that 'thid' might be negative for some threads => need to
        // use real thid id's 
        
        CU_SYNC

//         if(thid == 0) {
//             int i = 0;
//             for(i = 0; i < WS; i++) {
//                 if(L[i] != 0)
//                     break;
//             }
//             L[0] = (i != WS);
//         }

        //! NOTE NOTE total number of warps <= 32 (1024 threads max)
        nwarps = (blockDim.x + WS - 1) >> 5;
        if(t2 < WS) {
            x = 0;
            if(t2 < nwarps) {
                x = L[t2];
                // can just check for zero using this and do the remaining
                // in a loop - very unlikely that it will influence the perofrmance
                // point is you need to scan the results of ballot anyway
            }
            L[t2] = __any(x);
        }
// 
        CU_SYNC
        if(L[0] == 0) { // indicates that GCD is found
            nu = -1u;
            return;
        }
    }
#endif

    // NOTE you can optionally use the results of 'any' operation
    // to quickly find the warp from which to start scanning
    // scan all elements until we find a non-zero
    while(1) { // need to shift data at least once since L[last_thid] == 0
        if(thid <= last_thid) {
            L[thid] = F.w;
            F.w = F.z, F.z = F.y, F.y = F.x, F.x = 0;
        }
        CU_SYNC

        if(thid-1 < last_thid) {
            F.x = L[(int)(thid - 1)];
        }

        nu--;
        if((int)nu < 0) {
            return;
        }

        if(L[last_thid] != 0)
            goto Lexit;       

        CU_SYNC
    }

Lexit:
    // here we have that L[last_thid] != 0
    if(thid == 0) {
        lcF[0] = L[last_thid];
    }
}

//! PGCD algorithm with 4 elements processed per thread
template < bool SplitThreads >
__global__ void CUMP_LAUNCH_BOUNDS(1024, 1)
PGCD_quad_kernel(const unsigned *In0, unsigned *Out0) {

    extern __shared__ unsigned shared[];
//! bidx_x - moduli index
    unsigned *r = shared, thid = threadIdx.x, bidx_x = blockIdx.x;
    const unsigned *mods = dev_const_mem;

    if(SplitThreads) {
        bidx_x = UMUL(bidx_x, blockDim.y) + threadIdx.y;
        // need to split 'r' as well
    }

    // const mem  stride & size of cache line
    const unsigned m_stride = 4, CacheLn = 32;
    // two cache lines here
    unsigned *lcF = r, *lcG = lcF + 1, *cache = lcG + 1,
            *L = cache + CacheLn; //! beware that we use offset L[thid+1] !!

    unsigned nu = dev_const_mem[NU], nv = dev_const_mem[NV];
    mods += UMUL(bidx_x, UMUL_PAD + m_stride) + ConstMemDataSz;
    /*volatile*/ unsigned m = mods[0];
    volatile fp_limb invk, invm;
    volatile unsigned mx100;

    MODARITHM_INIT(thid)

    unsigned max_nu = dev_const_mem[MAX_NU],
        block_ofs = UMUL(bidx_x, dev_const_mem[DATA_PAD]);

//     const unsigned *In0 = (unsigned *)dev_const_mem[DATA_IN];
    const unsigned *f_in = In0 + block_ofs, *g_in = f_in + max_nu;

    // chunk_size = max_nv (that is nv+1 aligned by 16 boundary)
    // nu - degree of the dividend, nv - degree of the divisor
    // index of the last busy thread:
    unsigned last_thid = nv / 4; // (nv + 4) / 4 - 1; 
    uint4 F, G; // register sets for dividend and divisor

//     if(last_thid == 0)
//         last_thid = 1;

// NOTE NOTE: later we can make sure that data is originally stored with
// padding, i.e.: nv = 9
// thid:      0 0 0 0 1 1 1 1 2 2 2 2
//            e e 0 1 2 3 4 5 6 7 8 9
//            x y z w x y z w x y z w

    unsigned ofs, t = ((nu + 4) & ~3) - ((nv + 4) & ~3);

    if(thid <= last_thid) {
        // data is aligned in such a way that g[nv] is loaded in G.w
        G = ((uint4 *)g_in)[thid];
        // data offset to make sure f[nu] is loaded in F.w last thread
        // provided that both arrays are properly aligned
        F = ((uint4 *)(f_in + (int)t))[thid];
    }
    
    // there is no any global loop counter:
    // the parameters 'nu' and 'nv' change independently
    // of the local cache index position
    if(thid == last_thid) {
        lcF[0] = F.w;   // this is lcf(f)
        lcG[0] = G.w;   // this is lcf(g)
    }
//! NOTE NOTE NOTE: you can incorporate data shift directly to the algorithm

    L[thid + 1] = G.w, G.w = G.z, G.z = G.y, G.y = G.x;
    CU_SYNC

    G.x = L[thid];
    L[thid] = F.w, F.w = F.z, F.z = F.y, F.y = F.x;

    if(thid == 0) {
        F.x = 1, G.x = 0; // 0th thread collects powers of 'd'
    }

    CU_SYNC

    if(thid != 0 && thid <= last_thid) {
        F.x = L[thid - 1];
    }

    CU_SYNC
    // NOTE NOTE: make sure 't' is not modified before
    // t is location where you started loading data from
    // load umin(t, CacheLn) elements from interval [0..t-1]
    ofs = t - CacheLn + thid;
//     unsigned ofs3 = t;
    if(thid < CacheLn) {
        // load in reversed order to facilitate indexing
        // cache[0] = f[t-1]; cache[1] = f[t-2]; .. cache[15] = f[t-16]
        t = 0;
        if((int)ofs >= 0)
            t = f_in[ofs];
        cache[(int)(CacheLn - 1 - thid)] = t;
    }
    
    __pgcd_quad_internal< true >(L, lcF, lcG, F, G, nu, nv, last_thid,
            cache, CacheLn, ofs, f_in, thid, m, invk, invm, mx100);

    // NOTE NOTE: lcF[0] is not set to the next element after the call
    // hence we should use L[last_thid] instead !!!

    if(thid == 0) { // zero out aux elements before scanning
        F.x = 0, F.y = 0;
    }
    unsigned FIRST_THID = 0/*, CNT = 0*/;

//     int STOP = 2000;
//     int ccnt = 0;
    // <-- loop
    while(1) {

    CU_SYNC // sync because F.x is loaded from shared mem

    //! perform zero check for the dividend:
    //! at this point we have nu = nv - 1 => total of 'nv' elements in \c F
    __lcf_scan(L, lcF, F, nu, nv, last_thid, thid, FIRST_THID);

//     ccnt++;
    if(/*ccnt == STOP ||*/ (int)nu < 0) { // G contains gcd
//         G = F;
        goto Lexit;
    }

    //! FIXME FIXME FIXME: if smth does not work try to disable thid offset
    //! it is necessary that 'nu' and 'nv' are zero padded to have the same
    //! size
    //! there is a subtle problem: at the beginning of the next pgcd call
    //! it always holds that nu = nv + 1 (or vice a versa): hence
    //! the divisor vector must be padded with zeros to prevent a garbage
    //! in the dividend
    t = last_thid - ((nv + 4) / 4 - 1); // compute new 'last_thid'
    thid -= t, last_thid -= t;
    FIRST_THID += t;

#if 0
    //! FIXME FIXME FIXME: quad_lite iterations are still not 100% working
    //! disable this in case of problems
    if(nv < blockDim.x) {  //! repartition data and run the ``lite'' version

        CU_SYNC
        // F is a new divisor, G is a new dividend
        __pgcd_quad_lite(L, lcG, lcF, Out0, G, F, nv, nu, last_thid,
              thid, FIRST_THID, bidx_x, m, invk, invm, mx100);
        return;
    }
#endif

    //! __pgcd_quad_internal starts with CU_SYNC => no need to sync before
    __pgcd_quad_internal< false >(L, lcG, lcF, G, F, nv, nu, last_thid,
            cache, CacheLn, ofs, f_in, thid, m, invk, invm, mx100);

    CU_SYNC // sync because F.x is loaded from shared mem
    __lcf_scan(L, lcG, G, nv, nu, last_thid, thid, FIRST_THID);

//     ccnt++;
    if(/*ccnt == STOP ||*/ (int)nv < 0) { // F contains gcd
        G = F; nv = nu;
        lcG = lcF;
        goto Lexit;
    }

    // use here 'nu + 1' threads so that 0th thread can "sweep" the garbage
    // in the dividend
    t = last_thid - ((nu + 4) / 4 - 1); // compute new 'last_thid'
    thid -= t, last_thid -= t;
    FIRST_THID += t;

    __pgcd_quad_internal< false >(L, lcF, lcG, F, G, nu, nv, last_thid,
            cache, CacheLn, ofs, f_in, thid, m, invk, invm, mx100);

    }
    // loop -->

Lexit:

// //     CU_SYNC // need this sync ??

//     shift the data in G to the left by the amount 't'
    t = (last_thid + 1) * 4 - nv;
    //! FIXME NOTE NOTE: it seems that 't' is never zero hence we can use
    //! this to ``shift-in'' the missing leading coefficient
    if(thid == last_thid) {
        L[thid + 1] = lcG[0];
    }

    
    while(t != 0) {
        CU_SYNC

        if(thid <= last_thid) {
            L[thid] = G.x;
            G.x = G.y, G.y = G.z, G.z = G.w;
        }
        CU_SYNC

        if(thid <= last_thid)
            G.w = L[thid + 1];
        t--;
    }

    block_ofs = UMUL(bidx_x, dev_const_mem[GCD_SZ]) +
                         dev_const_mem[MODS_PAD];
//     unsigned *Out0 = (unsigned *)dev_const_mem[DATA_OUT];

    if(thid == 0)
        Out0[bidx_x] = // (unsigned)(Out0 - Out0_);
                nv + 1; // size of a gcd (not the degree!!)

    if(thid <= last_thid)
        ((uint4 *)(Out0 + block_ofs))[thid] = G;

//     G.x = last_thid;  //G.y = nv;
//     if(thid <= last_thid)
//         ((uint4 *)(Out0 + block_ofs))[thid] = G;
}

#endif // _PGCD_KERNEL_CU_
