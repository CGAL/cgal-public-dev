
//! runs \c stop_j - \c j iterations of pseudo-division
//! dividend is stored in shared mem \c L , addressed as \c L[thid-j]
//! divisor in registers \c g leading element of divisor: \c d[0]
//! thread 'nv' is responsible for computing the powers of \c d[0]
//! which are placed in \c L[i_exp]
//! at iteration \c j , \c L[nv-j] is the leading element of the dividend
//! \c skip[0] should point to the number of elements skipped since the
//! previous iterations (at least 1)
// template < bool FullUpdate >
__device__ void __prim_div(unsigned *L, unsigned *d, unsigned *skip,
        unsigned& j, unsigned g, unsigned stop_j,
        const unsigned i_exp, const unsigned nv,
        unsigned thid, unsigned m, fp_limb invk,
        fp_limb invm, volatile unsigned mx100) {

    unsigned t = thid - j;
    if(thid == nv) {
        t = i_exp; // points to special entry
    }

    while(1) {

        unsigned lcd = d[0];
        if(thid < skip[0]) {
           lcd = L[i_exp];
        }

        unsigned f = L[t], lcr = L[nv - j];
        L[t] = sub_mul_reduce_m(f, lcd, lcr, g, m,
            invk, invm, mx100);

        CU_SYNC
   
        if(thid == 0) {
//              L[nv - j] = 0;
            // search for the 1st non-zero leading element
            unsigned s = j + 1; // skipping the 1st one
            while(L[nv - s] == 0 && s < stop_j) {
                s++;
            }
            skip[0] = s - j; // count # of elements skipped
        }
        CU_SYNC

        j += skip[0];
        if(j >= stop_j)
            return;

        if(thid < nv) {
            t -= skip[0];
        }
    }
}

__global__ void /*CUMP_LAUNCH_BOUNDS(512, 4)*/
prim_gcd_lite_kernel() {

    extern __shared__ unsigned shared[];
//! bidx_x - moduli index
//! bidx_y - index within a batch for single modulus
    unsigned thid = threadIdx.x, bidx_x = blockIdx.x;
    unsigned *r = shared, *mods = dev_const_mem;
   // leading coeffs of the dividend and divisor
    // L[-1] is used to store d^e
    unsigned *d = r, *skip = d + 1, *L = skip + 2; 

    const unsigned m_stride = 4; // const mem  stride

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

    /*const */ unsigned max_nu = dev_const_mem[MAX_NU],
            block_ofs = UMUL(bidx_x, dev_const_mem[DATA_PAD]);
    const unsigned *In0 = (unsigned *)dev_const_mem[DATA_IN],
            *f_in = In0 + block_ofs, *g_in = f_in + max_nu;

//! NOTE: the actual number of threads we need is 'nv'
//! however operands can be unbalanced - in this case we need
//! many iterations to get the data loaded !!

// chunk_size = max_nv (that is nv+1 aligned by 16 boundary)
// since in any case we need to align block size with that
// and also need to preload 'nv + 1' elements all at once

    const unsigned CacheLn = 32;
    // nu - degree of the dividend, nv - degree of the divisor
    unsigned t, g;

#if 1
    // load first part of the dividend to shared mem
    t = CacheLn + thid;
    // load umin(nv+1, nu+1) elements using nv+1 threads
    // starting from the top: f_in[nu], f_in[nu-1], f_in[nu-2] ...
    unsigned ofs = thid + nu - nv;
    if((int)ofs >= 0) {
        L[t] = f_in[ofs];
    }
#endif

    g = g_in[thid]; // # of threads == nv
    if(thid == nv) {
        d[0] = g; // share the leading element btw all threads
    }
    
    CU_SYNC

#if 1
    n_left = nu - nv;
    // load umin(n_left, CacheLn) elements using CacheLn thids
    if(thid < CacheLn) {
        ofs = thid + n_left - CacheLn;
        if((int)ofs >= 0) {
            L[thid] = f_in[ofs];
        }
    }
    n_left -= CacheLn;
#endif

    // # of iterations that can be run is:
    // unsigned n_max_iters = umin(nu - nv, CacheLn);
    if(thid == nv) {
        L[-1u] = g;  // L[-1] stores d^e
        g = 0;
        skip[0] = 1;
    }
    CU_SYNC

    unsigned j = 0, stop_j = umin(nu - nv + 1, CacheLn);
    const unsigned i_exp = -CacheLn-1;

    unsigned *L2 = L + CacheLn;
//     stop_j = 10;

    // j always equals stop_j at the end
    __prim_div(L2, d, skip, j, g, stop_j, i_exp, nv,
        thid, m, invk, invm, mx100);

    if(j < nu - nv + 1) {
        CU_SYNC
        // relocate L[] back before the cache line
        t = L[thid];
        CU_SYNC
        L[thid + CacheLn] = t;
        // no need to sync here

        if(thid < CacheLn) {
            ofs = thid + n_left - CacheLn;
            if((int)ofs >= 0) {
                L[thid] = f_in[ofs];
            }
            n_left -= CacheLn;
        }
        j += CacheLn; // move j back  
    }

    CU_SYNC

    // we either need to find the next g which is not zero
    // or just look over all of them..

    block_ofs = UMUL(bidx_x, dev_const_mem[GCD_SZ]) +
                         dev_const_mem[MODS_PAD];
    unsigned *Out0 = (unsigned *)dev_const_mem[DATA_OUT];

    if(thid <= nv)
        Out0[block_ofs + thid] = j;//L2[thid - j + 1];
}