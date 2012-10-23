
template < unsigned Mod4 >
__device__ unsigned __interpolate_quads_internal(uint4& b, unsigned *r,
        const unsigned *g_Xs, unsigned n, unsigned n_thids,
        unsigned block_sz, volatile unsigned m,
            volatile fp_limb invk, volatile fp_limb invm) {

    extern __shared__ unsigned shared[];

    unsigned thid = threadIdx.x, bidx_x = blockIdx.x;//, bidx_y = blockIdx.y;
//! bidx_x - indexes over moduli

    unsigned *r = shared;
    //const unsigned N = 1 << LDN;
    const unsigned m_stride = 4; // constant memory stride
    unsigned *mods = dev_const_mem;// + UMUL(bidx_x, UMUL_PAD + m_stride);

    unsigned ofs = UMUL(bidx_x, out_data_padding);
    unsigned last_thid = (n + 3) / 4 - 1; // index of the last working thid
    //! offsets added! do not forget
    g_R += ofs, g_Xs += ofs;

    uint4 b;
    if(thid <= last_thid) // read out only relevant data
        b = ((uint4 *)g_R)[thid];

    // 4 elements per thread => N/4 max thids
    unsigned block_sz = ((n + 3) / 4 + 31) & ~31;

    unsigned thid = threadIdx.x;
    //! NOTE As, Bs and xs have per-thread indexing !!
    unsigned *Bs = r + thid, *As = Bs + block_sz, *xs = As + block_sz;

    volatile int j = n - 1;
    unsigned a0 = 1, a_prev;

    uint4 a, w;
    if(thid <= last_thid) {
        w = ((uint4 *)g_Xs)[thid];
// better to wait with saving to shared mem, so that we can overlap mem
// access with ALUs..
    }

    // shift the vectors down
    xs[0] = w.x, w.x = w.y, w.y = w.z, w.z = w.w, w.w = 0;
    Bs[0] = b.x, b.x = b.y, b.y = b.z, b.z = b.w, b.w = 0;

    CU_SYNC
    if(thid  < last_thid) { // what to do in case of odd n ??
                        // force even number of evaluation points always ??
        b.w = Bs[1]; // except the last element
        w.w = xs[1]; // read "the tail" of the next thread
    }

    a.x = 1, a.y = a.x, a.z = a.x, a.w = a.x;
    unsigned det = 1; // this will store the denominator

    while(1) {

    unsigned x0;

    if(thid == last_thid) {
        b.w = det; // last thid computes denominator
        a.w = 0;
    }

    //! shift all index ranges by -(j+1)
    x0 = r[0]; // r[0] == b0

    b.x = sub_mul_reduce_m(b.x, a0, a.x, x0, m, invk, invm, mx100);
    b.y = sub_mul_reduce_m(b.y, a0, a.y, x0, m, invk, invm, mx100);
    b.z = sub_mul_reduce_m(b.z, a0, a.z, x0, m, invk, invm, mx100);
    b.w = sub_mul_reduce_m(b.w, a0, a.w, x0, m, invk, invm, mx100);

    if(thid == last_thid) {
        det = b.w;
        if(Mod4 == 1) { // this is compile-time branch
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
    x0 =  r[block_sz*2]; // == xs[-thid]

    // b0 not needed anymore
    k = j - 4*thid - 1;
    if(thid == j / 4) {
        a_prev = m;
    }

    if(j == 0)
        break;
    j--;

    CU_SYNC

    Bs[0] = b.x; // b.x is free now

/**
 if((int)k > 2) { 4 first

    } else if((int)k == 2) { // 3 first, 1 second

    } else if((int)k == 1) { // 2 first, 2 second

    } else if((int)k == 0) { // 1 first, 3 second

    k < 0: 4 second
*/

//! do ping-ponging registers: mul1 <--> t1
    unsigned t1 = 0, s1 = x0;

    if((int)k >= 0) {
        s1 -= w.x;
    } else {
        t1 = a_prev;
    }

    b.x = a.x;
    a.x = t1 - s1 * a.x;

    s1 = m; // s1 is now t2

    if((int)k < 0) {
        s1 = b.x;
    }

    if((int)k == 0) { // why it works without this as well ?
        s1 = a_prev;
    }

    t1 = x0; // t1 is now mul2, b.x is free here
    if((int)k >= 1) {
        t1 -= w.y; // (n - 1 - b0);
    } 

    b.x = a.y;
    a.y = s1 - t1 * a.y;

    t1 = 0;
    if((int)k < 1) {
        t1 = b.x; // t1 is t3 now
    }

    s1 = x0; // s1 is mul3 again
    if((int)k >= 2) {
        s1 -= w.z;
    }

    if((int)k == 1) {
        t1 = a_prev;
    }

    b.x = a.z;
    a.z = t1 - s1 * a.z;

    s1 = 0; // s1 is t4
    if((int)k == 2) {
        s1 = a_prev;
    }

    if((int)k < 2) {
        s1 = b.x;
    }    

    t1 = x0; // t1 is mul4
    if((int)k >= 3) {
        t1 -= w.w;
    }
    a.w = s1 - t1 * a.w;

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

    if(thid < last_thid) { // TODO: can merge these two conditions into one
        a.w = As[1]; // read out elements shifted by one (except the last
        b.w = Bs[1]; // thid)
        w.w = xs[1];
    }

    a0 = r[block_sz]; // == As[-thid]

    } // while(1)

    return det;
}


__global__ void interpolate_defaul_alg() {

    extern __shared__ unsigned shared[];

    unsigned thid = threadIdx.x, bidx_x = blockIdx.x;//, bidx_y = blockIdx.y;

    unsigned *r = shared;
    const unsigned *in_U = g_U + ofs;

    unsigned a, b, w;
    if(thid < n) {// read out only relevant data
        b = g_Y[thid]; // read in y-values
        w = g_X[thid]; // read in x-values
    }

    const unsigned block_sz = n; // data size
    unsigned *Bs = r + thid, *As = Bs + block_sz,
        *xs = r + block_sz*2, *det = xs + block_sz;

    // j -> (n - 1) - j
    volatile unsigned j = 0; // loop counter
    int k;

    Bs[0] = b; // keep the current leading elements
    As[0] = 1; // a's are initially ones
    xs[thid] = w;

    if(thid == 0)
        det[0] = 1;

    unsigned last_thid = n - 1;
    
    CU_SYNC

    while(1) {

        if(thid < last_thid) {
            a = As[1]; b = Bs[1];
        } else if(thid == last_thid) {
            a = 0, b = det[0];  // the last thread computes denominator
        }

        unsigned b0 = r[0], a0 = r[n];
// TODO: need to use spare thread to compute l_int
// we do not have more threads for update
        b = b * a0 - a * b0;

        if(thid == last_thid) {
            a = 1, det[0] = b, b = -b0; // b[j + n] = -b0
        }

        s = 0, t = 0;
        k = thid + 1 + j; // thread indicator
        if(k < n) {
            s = a, t = xs[k];
        } else if(k > n) {
            s = As[0], t = 1;
        }
        a = s * t - a * xs[j];

        if(++j == n)
            break;

        CU_SYNC // in order to avoid racing conditions

        As[0] = a; Bs[0] = b; // shift-down the vectors
        CU_SYNC

    } // while(1)

    return b / det[0]; // return the polynomial coefficients 
}

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