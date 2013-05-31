
#ifndef _BLOCK_UPDATES_GJG_C_
#define _BLOCK_UPDATES_GJG_C_

#include "include/generator_defs.h"
#include <include/modular_arithm.h>

namespace CGAL {

namespace internal {

// bool g_save = false;
// std::vector< zmod > g_debug;

//! first phase of generator updates: compute L factors of \c gs and \c ms and save them to \c gs_up and \c ms_up resp.
//! \c [g_start;g_beyond) - index range for \c gs
//! \c [i_start;i_beyond) - index range for \c ms
template < class NT >
void _gjg_update_1(const Generator< NT >& gs, const Generator< NT >& ms,
       GJG_update< NT >& gs_up, GJG_update< NT >& ms_up, NT& det,
            unsigned g_start, unsigned g_beyond,
            unsigned i_start, unsigned i_beyond, bool full_update = true) {

    unsigned j = g_start, i; // j is a leading element index
    NT a0 = gs.a[j], b0 = gs.b[j], c0 = gs.c[j], d0 = gs.d[j];

    for(i = g_start; (int)i < (int)g_beyond; i++) {
        NT t = a0 * gs.a[i] + b0 * gs.b[i];
        if(full_update)
            t = t - c0 * gs.c[i] - d0 * gs.d[i];
        gs_up[i] = t;
    }
    det = gs_up[j];

    for(i = i_start; (int)i < (int)i_beyond; i++) {
        NT t = a0 * ms.a[i] + b0 * ms.b[i];
        if(full_update)
            t = t - c0 * ms.c[i] - d0 * ms.d[i];
        ms_up[i] = t;
    }
}

//! \c full_update - indicates that full iterations are performed
template < class NT >
void _gjg_update_2(Generator< NT >& gs, Generator< NT >& ms,
       const GJG_update< NT >& gs_up, const GJG_update< NT >& ms_up, 
       const NT& det, unsigned g_start, unsigned g_beyond,
            unsigned i_start, unsigned i_beyond, bool full_update = true) {

    unsigned j = g_start, i, k; // j is a leading element index
    NT a0 = gs.a[j], b0 = gs.b[j], c0 = gs.c[j], d0 = gs.d[j];

    Generator< NT > *v = &gs;
    const GJG_update< NT > *v_up = &gs_up;
    unsigned i0 = g_start, i1 = g_beyond;

    for(k = 0; k < 2; k++) {
        // in fact det can be extracted from gs_up ??
        for(i = i0; (int)i < (int)i1; i++) {
            v->a[i] = v->a[i] * det - (*v_up)[i] * a0;
            v->b[i] = v->b[i] * det - (*v_up)[i] * b0;
            if(full_update) {
                v->c[i] = v->c[i] * det - (*v_up)[i] * c0;
                v->d[i] = v->d[i] * det - (*v_up)[i] * d0;
            }
        }
        v = &ms, v_up = &ms_up, i0 = i_start, i1 = i_beyond;
    }
}

template < class NT >
void _gjg_setup_full(Generator< NT >& gs, Generator< NT >& ms,
    const NT& res, unsigned g_start, unsigned g_beyond,
                   unsigned i_start, unsigned i_beyond) {

    Generator< NT > *v = &gs;
    unsigned i0 = g_start, i1 = g_beyond, i, k;

    for(k = 0; k < 2; k++) {
//         in fact det can be extracted from gs_up ??
        for(i = i0; (int)i < (int)i1; i++) {
            v->c[i] = v->c[i] * res;
            v->d[i] = v->d[i] * res;
        }
        v = &ms, i0 = i_start, i1 = i_beyond;
    }
}

//! shifts down the elements of \c v and optionally inserts elements \c a_in
//! at the top of \c v  
//! only elements in the range \c [i_start;i_beyond) affected 
template < class NT >
void _gjg_shiftdown(GJG_update< NT >& v, const NT& a_in, 
        unsigned i_start, unsigned i_beyond, bool insert) {

    unsigned i;
    for(i = i_beyond - 1; (int)i > (int)i_start; i--) { 
        v[i] = v[i] - v[i - 1];
    }
    if(insert) {
        v[i_start] = v[i_start] - a_in;
//         printf("inserting element: %d - %d\n:", v[i_start].x, a_in.x);
    }
}
        
//! runs \c n steps of schur algorithm starting with iteration \c iter_idx
//! updates are generated \c gs and applied to \c ms
//! full updates start with iteration \c switch_idx
//! if \c ins <> 0, it is used to update elements from \c ms otherwise \c ms
//! is updated in a ``stair-like'' fashion
//! \c shift_through : elements shifted out from \c gs are inserted at the top
//! of \c ms (in other words, gs == ins reversed)
template < class NT >
void schur_block(Generator< NT >& gs, Generator< NT >& ms,
       NT& res, NT& denom, unsigned iter_idx, unsigned switch_idx,
        unsigned n, bool shift_through,
        GJG_update< NT > *ins = 0, GJG_update< NT > *outs = 0) {

    unsigned g_sz = gs.a.size(), m_sz = ms.a.size(), i, j,
        i_start = 0, i_beyond = m_sz;

    bool use_ins = (ins != 0), use_outs = (outs != 0);
    if(use_ins) { // in case ins is used we update in reverse ``stair-like''
                  // fashion
        i_beyond = 1;
    }

    for(j = 0, i = iter_idx; j < n; j++, i++) {

        bool full_iter = (i >= switch_idx);
        GJG_update< NT > gs_up(g_sz, NT(0)), ms_up(m_sz, NT(0));

        if(i == switch_idx) {
        printf("########## setting up full iters: det = %d\n", res.x);
        //! NOTE NOTE: when switching to full iterations we always update
        //! the full range of \c ms if \c use_ins == false
            _gjg_setup_full(gs, ms, res, j, g_sz, 0, (use_ins ? 0 : m_sz));
        }

        NT det;
        _gjg_update_1(gs, ms, gs_up, ms_up, det, j, g_sz, i_start, i_beyond,
                full_iter);

        NT a_in(0), _(0);
        denom = denom * res * res, res = res * det;

        if(shift_through) {
            a_in = gs_up[g_sz - 1];
        } else if(use_ins) {
            a_in = (*ins)[j];
        } else {
            i_start++; 
//             printf("%d: stair-like: [%d; %d)\n", j, i_start, i_beyond);
        }

        if(use_outs) {
            (*outs)[j] = ms_up[m_sz - 1];
//             printf("outs[%d] = %d\n", j, outs->a[j].x);
        }
        _gjg_shiftdown(gs_up, _, j, g_sz, false);
        _gjg_shiftdown(ms_up, a_in, 0, m_sz, use_ins | shift_through);

        // should this be j+1.. g_sz ??
        _gjg_update_2(gs, ms, gs_up, ms_up, det, j, g_sz, i_start, i_beyond,
                full_iter);

        if(use_ins) {
           i_beyond++;
            // NOTE guard is needed if n > m_sz which can occur for the last
            if(i_beyond > m_sz) // block
                i_beyond = m_sz;
        }
    }
}
    
template < class NT, unsigned N >
void block_schur_QR_test(const std::vector< NT >& v1,
         const std::vector< NT >& v2, unsigned *out = 0,
            unsigned o_stride = 1) {

    unsigned n = v1.size()-1, m = v2.size()-1, r = n + m, i, j;

    if(n < m) {
        printf("sylvester_QR_gcd: incorrect parameters !!\n");
        throw 1;
    }
    unsigned n_blocks = (r + N - 1) / N - 1,
        n_last = r - (n_blocks * N), // # of elements for the last block
        rN = (n_blocks + 1) * N;   // allocate padded generators

    printf("block_schur_QR N: %d; r: %d; n_blocks: %d; n_last: %d\n",
            N, r, n_blocks, n_last);

    zmod v1lc = v1[n];
    zmod inv1(mod_inverse(v1lc.x, zmod::modulus)), inv2(1);

    std::vector< NT > a(r, NT(0)), b(r, NT(0)), c(r, NT(0)), d(r, NT(0));
    for(i = 0; i <= n; i++) {
        a[i] = (i <= n ? v1[n - i] : NT(0));
        b[i] = (i <= m ? v2[m - i] : NT(0));
        
        if(i < n)
            c[i + m] = v1[n - i];
        if(i < m)
            d[i + n] = v2[m - i];
    }

    typedef Generator< NT > GJG;
    typedef GJG_update< NT > GJG_up;
    
    typename std::vector< NT >::const_iterator ai = a.begin(),
            bi = b.begin(), ci = c.begin(), di = d.begin();
    
    GJG dummy, gs0(N, ai, bi, ci, di), gs;
    std::vector< GJG > msv(n_blocks);
    std::vector< GJG_up > outsv(n_blocks);
    GJG *ms = (GJG *)msv.data();
    GJG_up *outs = (GJG_up *)outsv.data();

    // distribute data btw blocks
    unsigned ofs;
    for(i = 0, ofs = N; i < n_blocks - 1; i++, ofs += N) {
        ms[i] = GJG(N, ai + ofs, bi + ofs, ci + ofs, di + ofs);
        outs[i] = GJG_up(N);
    }
    ms[n_blocks - 1] = GJG(n_last, ai + ofs, bi + ofs, ci + ofs, di + ofs);
    outs[n_blocks - 1] = GJG_up(n_last); // last block does not need outs

    /** **************** STEP 1 *******************/
    /** each block performs G_11..G_14 rotations **/
    /** the first block uses shift-through **/

    NT res0(1), denom0(1), res1, denom1;
    unsigned iter_idx = 0, switch_idx = m; // algorithm iteration counter

    g_debug = std::vector< NT >(N);
    // 2 blocks, 3*N data elements
#if 0
    // now we need to make a loop out of it and realize switch to the full
    // iterations

    j = 0;
    for(i = j; i < n_blocks; i++) {
        // limit the number of iterations for the last block
        unsigned n_iters = (i < n_blocks - 1 ? N : n_last);

        gs = gs0, res1 = res0, denom1 = denom0;
        GJG_up *pout = (i < n_blocks - 1 ? &outs[i] : 0);
        schur_block(gs, ms[i], res1, denom1, iter_idx, switch_idx,
              N, (i == j), (GJG_up *)0, pout);
    }

    for(i = 1; i < n_blocks; i++) {
        gs = gs0;
        res1 = res0, denom1 = denom0;
        g_save = (i == 1);
        schur_block(gs, ms[i], res1, denom1, iter_idx, switch_idx, N,
                 false, &outs[i - 1]);
    }

    return;
#else
    for(j = 0; j < n_blocks-1; j++) {

//         if(j == 5)
//             break;
 // NOTE: can we merge the two stages to avoid ``stair-like'' updates ??
        for(i = j; i < n_blocks; i++) {
            // limit the number of iterations for the last block
            unsigned n_iters = (i < n_blocks - 1 ? N : n_last);

            gs = gs0, res1 = res0, denom1 = denom0;
//! NOTE NOTE: we have to run N full iterations even for the last block
//! in order not to miss the transition to full iterations
        
            GJG_up *pout = (i < n_blocks - 1 ? &outs[i] : 0);
            schur_block(gs, ms[i], res1, denom1, iter_idx, switch_idx,
                N, (i == j), (GJG_up *)0, pout);
        }

        /** **************** STEP 2 *******************/
        /** each block except the first updates its ms[i] using 'outs' of the
            previous block **/
        for(i = j + 1; i < n_blocks; i++) {
            gs = gs0;
            res1 = res0, denom1 = denom0;
        // NOTE: here we need to perform exactly N iters even if the last
        // block has smaller size because leading elements must be updated
        // exactly N-1 times
            schur_block(gs, ms[i], res1, denom1, iter_idx, switch_idx, N,
                 false, &outs[i - 1]);
//             if(j == 2) {
//                 printf("\n########### block %d after %d full iterations\n",
//                     i, (j+1)*N);
//                 ms[i].print();
//             }
        }
        // each loop iteration equivalent to N steps of Schur algorithm
        iter_idx += N;
        res0 = res1, denom0 = denom1; // and determinant as well
        gs0 = ms[j]; // update leading block

    } // for(j)

//     schur_block(ms[n_blocks-2], ms[n_blocks-1], res1, denom1, iter_idx,
//              switch_idx, N, true);

    printf("############## truth: \n");
    ms[n_blocks-2].print();
    ms[n_blocks-1].print();

    // it remains to run the last 2*N iterations by the last block with
    // shift-through from the previous one (should be merged in a single call)
    schur_block(ms[n_blocks-2], ms[n_blocks-1], res1, denom1, iter_idx,
             switch_idx, N, true);

printf("TRUTH res = %d; denom = %d\n", res1.x, denom1.x);
    
    iter_idx += N;
    printf("\n########### ms2 last iterations: \n");
    ms[n_blocks-1].print(); // this is after performing 12 iterations

    // and finalize
    schur_block(ms[n_blocks-1], dummy, res1, denom1, iter_idx, switch_idx,
         n_last, false);
//     NT inv(mod_inverse(denom1.x, zmod::modulus));
//     det1 = det1 * v12lc_exp * inv;
#endif

    printf("TRUTH res = %d; denom = %d\n", res1.x, denom1.x);
}

} // namespace internal

} // namespace CGAL

#endif // _BLOCK_UPDATES_GJG_C_
