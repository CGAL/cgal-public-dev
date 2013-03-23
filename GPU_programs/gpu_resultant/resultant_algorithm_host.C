// ============================================================================
//
// Copyright (c) 2001-2008 Max-Planck-Institut Saarbruecken (Germany).
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

#include <fstream>
#include <getopt.h>

#define CGAL_MODULAR_FILTER_OFF

#include <include/macros.h>
#include <include/device_manager.h>

#include <CGAL/GPU_algorithms/resultant_algorithm.h>
#include <CGAL/GPU_algorithms/GPU_algorithm_templates.h>

#if CUMP_COMPILE_DEBUG_KERNEL
#include <CGAL/Algebraic_kernel_d/Polynomial_parser_d.h>
#endif

#include <include/modular_arithm.h>
#include <include/misc.h>
#include <matrix_algebra.C>

#include <CGAL/Timer.h>
 

namespace CGAL {

namespace internal {

void resultant_abort_handler(int sig)
{
    printf("resultant abort handler..\n");

    GPU_resultant& obj = GPU_resultant::instance();
    obj.~GPU_resultant();
}

GPU_resultant::GPU_resultant() {

    host_static_setup();
    GPU_device_manager& obj = GPU_device_manager::instance();
    obj.push_CUDA_context();
    device_static_setup();
    obj.pop_CUDA_context();
}

void GPU_resultant::host_static_setup() {

    CPU_mem_ptr = 0, DEV_mem_ptr = 0;
    batch_size = 1, word_size = 4;
    aux_cpu_mem = 0, aux_dev_mem = 0;
    CPU_bytes_alloced = 0, DEV_bytes_allocated = 0;
    pts_overrun = 8;

    CUMP_N_MODULI = -1u;
    CUMP_N_BATCHES_PER_MOD = -1u;
    no_cpu_run = true; // quit immediately after kernel launch

    write_to_file = false;
    fg_idx = 1;
    nu = -1u, nv = -1u;
    deg_x1 = 10, deg_x2 = 10;
    bits = 32;
    strcpy(benchmark_file, "benchmarks_out");

    // initializes device and loads moduli set
    GPU_device_manager& obj = GPU_device_manager::instance();
    obj.register_callback((GPU_DEVICE_CALLBACK)resultant_abort_handler);
}

bool GPU_resultant::setup(unsigned deg_y1_, unsigned deg_y2_,
        unsigned deg_x1_, unsigned deg_x2_, unsigned low_deg_,
        unsigned high_deg_, unsigned bits_) {

// TODO: additional handling required if f or g has a trivial factor y^k
// then we compute lower degree resultant and multiply it by appropriate
// power of y

    GPU_device_manager& obj = GPU_device_manager::instance();

    unsigned nmods, npts;
    high_deg = high_deg_, low_deg = low_deg_;
    nu = deg_y1_, nv = deg_y2_;
    deg_x1 = deg_x1_, deg_x2 = deg_x2_;
    bits = bits_;

    CUMP_out2("low_deg: %d; high_deg: %d; bits: %d\n", low_deg, high_deg,
             bits)
    double points_enlarge = 1.005;

    nmods = bits_to_moduli(bits, obj.mod_bits);
    if(nmods == 0) {
        std::cerr << "FATAL: moduli exhausted !!\n";
        return false;
    }
//     printf("nmoduli: %d; upper bound: %f\n", nmods, obj.mod_bits[nmods-1]);

    low_deg = 0; //! for the time being..
    n_real_pts = high_deg - low_deg + 1; // resultant degree to be recovered
    npts = (unsigned)std::ceil((double)n_real_pts * points_enlarge);

    
    // at least 8 additional points    
    npts = std::max(npts, n_real_pts + pts_overrun); 

    CUMP_out2("nmods: %d; n_real_pts: %d; npts: %d\n", nmods, n_real_pts, npts)
    nr = nu + nv;

    CUMP_N_MODULI = nmods;
    CUMP_N_BATCHES_PER_MOD = npts;

//     CUMP_N_MODULI += (CUMP_N_MODULI & 1); // in fact we can force the
                // number of moduli to be even (facilitates computation)

    mods_padding = (CUMP_N_MODULI + 15) & ~15; // next divisible by 32
    // invk can occupy 2 words (for double precision)
    aux_cpu_mem = mods_padding * 5; // (m, mu, inv_mod, invk) for CRA
    aux_dev_mem = aux_cpu_mem; 

    try {
        resultant_data_padding();   // compute max_nr (input data padding)
        interpolate_data_padding(); // compute output data padding
    } catch(GPU_algorithm_exception) {
        return false;
    }

    if(write_to_file) {
     //   strcpy(benchmark_file, "benchmarks_out");
        std::ofstream out(benchmark_file, std::ios::app);
        out << "deg_x1: " << deg_x1 << " deg_x2: " << deg_x2 <<
            " deg_y1: " << nu << " deg_y2: " << nv <<
            " nmods: " << nmods << " npts: " << npts;
//         out << "#fg_idx: " << fg_idx << "\n";
//         out << "#deg_x1: " << deg_x1 << " deg_x2: " << deg_x2 << "\n";
//         out << "#deg_y1: " << nu << " deg_y2: " << nv << "\n";
//         out << "#nmods: " << nmods << " npts: " << npts << "\n";
    }

    // padding both for moduli and eval points because data is transposed
    // in the CRT kernel
    batch_size = mods_padding * out_data_padding;
    prod_sz = 1; // each block outputs a single elem

    nr_sz = max_nr * (deg_x1 + deg_x2 + 2) * CUMP_N_MODULI / 2;

    if(!alloc_device_mem())
        return false; 

    return true;
}

bool GPU_resultant::run_gpu_part(const MPZ_vector_2& fv,
            const MPZ_vector_2& gv, MPZ_vector_1& r) {

    //! (m, mu, inv_mod, invk)
    unsigned *U = (unsigned *)CPU_mem_ptr,
        *R = U + nr_sz, *reference = R + prod_sz * batch_size, 
        *Mods = reference + prod_sz * batch_size,
        *Mus = Mods + mods_padding, *InvMods = Mus + mods_padding,
        *InvKs = InvMods + mods_padding;

    static unsigned const_mem[CUMP_DEVICE_MODULI_SZ];
    unsigned *pconst = const_mem;

    pconst[0] = nu, pconst[1] = nv;
    pconst += 2;

#if CUMP_BENCHMARK_CALLS
    CGAL::Timer rolex;
    float mod_reduce_tm, crt_cpu_tm;
    rolex.start();
#endif

    memset(U, 0, mem_size_nr);
    //! reduce polynomials and build up the moduli set without ``bad'' primes
    if(!RNS_conversion(fv, gv, U, (void *)pconst,
            Mods, InvKs, Mus)) {
        return false;
    }

#if CUMP_BENCHMARK_CALLS
    rolex.stop();
    mod_reduce_tm = rolex.time() * 1000.0f;
#endif

    CUMP_out2("nu: %d; nv: %d; deg_x1: %d; deg_x2: %d\n", nu, nv, deg_x1,
             deg_x2);

//! //////////////////////////////////////////////////////////////////////
//!TODO TODO: mod_inverse_seq must be precomputed
//! //////////////////////////////////////////////////////////////////////
    mod_inverse_seq(InvMods, Mods, CUMP_N_MODULI, 1);

    try {
        launch_kernel(const_mem, Mods, U, R);
    }
    catch(CGAL::internal::GPU_algorithm_exception) {
        return false;
    }

#if CUMP_MRC_HOST_RECOVER
#if CUMP_BENCHMARK_CALLS
    rolex.reset();
    rolex.start();
#endif

    CUMP_out2("RNS_recover\n");

    RNS_recover(R, Mods, r);

    CUMP_out2("RNS_recover done\n");

#if CUMP_BENCHMARK_CALLS
    rolex.stop();
    crt_cpu_tm = rolex.time() * 1000.0f;
#endif

#endif // CUMP_MRC_HOST_RECOVER

#if CUMP_BENCHMARK_CALLS
    if(write_to_file) {
        std::ofstream out(benchmark_file, std::ios::app);
//         out << "# mod_reduce: " << mod_reduce_tm
//                 << " CPU_CRA: " << crt_cpu_tm << "\n";
//     out << "#=========================================================\n";
        out << " mod_reduce: " << mod_reduce_tm
                << " CPU_CRA: " << crt_cpu_tm << "\n";
    }
    CUMP_out2("mod_reduce time: %f ms; CPU CRT time: %f ms\n",
            mod_reduce_tm, crt_cpu_tm);
#endif // CUMP_BENCHMARK_CALLS

    if(no_cpu_run)
        return true;

    //       memset(reference, 0, prod_sz*batch_size*sizeof(unsigned));
//       return checkme(R, reference, CUMP_N_BATCHES_PER_MOD,
//             out_data_padding, CUMP_N_MODULI, -1u);
//         return checkme(R, reference, n_real_pts,
//             out_data_padding, CUMP_N_MODULI, -1u);

    return reference_solution(const_mem, U, Mods, InvMods,
        reference, R);
}

//! verifies the result on the host: \c const_mem : constant memory pointer,
//! \c U - input set (polynomial coeffs), \c reference - memory set for
//! reference solution, \c GPU_res - GPU result
bool GPU_resultant::reference_solution(const unsigned *const_mem,
        const unsigned *U, const unsigned *Mods, const unsigned *InvMods,
        unsigned *reference, const unsigned *GPU_res) {

    const unsigned *pconst = const_mem + 2, *pU;
    unsigned *pref = reference;
    unsigned i, j, m;

    //! NOTE NOTE NOTE don't forget out_data_padding
    for(i = 0, pU = U; i < CUMP_N_MODULI; i++, pconst += CUMP_MODULI_STRIDE,
            pU += nr_sz / CUMP_N_MODULI) {

        m = pconst[0];
        zmod::set_mod(m);

        std::vector< unsigned > xs(CUMP_N_BATCHES_PER_MOD),
                 ys(CUMP_N_BATCHES_PER_MOD);

        CUMP_out2("m: %x\n", m);

        unsigned jidx = 0;
        for(j = 0; j < CUMP_N_BATCHES_PER_MOD; j++) {
            std::vector< zmod > v1m(nu + 1), v2m(nv + 1);
            zmod pt(j + 1);

            poly_eval2_pure(pU, pt, deg_x1, nu, max_nr/2,
                 v1m.begin());
            poly_eval2_pure(pU + max_nr/2 * (deg_x1 + 1), pt,
                     deg_x2, nv, max_nr/2, v2m.begin());

            bool failed;
            zmod det = resultant_modular_sylvester(v1m, v2m, failed);

            if(failed) {
#if CUMP_RUN_RESULTANTS_ONLY
                pref[j] = 0;    // save zero in case of ill-conditioned matrix
#endif
                printf("failed: modulus (%d): %x; point: %d\n", i, m, pt.x);
                continue; // we do not collect this point
            }
//             pref[jidx] = det.x;  pref[jidx] = pt.x;
            xs[jidx] = pt.x, ys[jidx] = det.x; jidx++;
#if CUMP_RUN_RESULTANTS_ONLY
            pref[j] = det.x;
#endif
        }

        if(jidx < n_real_pts) {
            printf("unable to recover resultant polynomial: %d %d\n",
                jidx, n_real_pts);
            exit(1);
        }

//         pref += out_data_padding;
//         std::swap(xs[31], xs[14]); std::swap(ys[31], ys[14]);

#if CUMP_USE_CRA_KERNEL // if using CRT kernel: save with mods_padding stride
        vandermonde_interpolate(xs.data(), ys.data(), pref, n_real_pts,
                 mods_padding);
        pref++;
#else

#if !CUMP_RUN_RESULTANTS_ONLY
        // we interpolate only with the required # of points
        vandermonde_interpolate(xs.data(), ys.data(), pref, n_real_pts, 1);
#endif

        if(write_to_file/* && i==0*/) {
//             printf("m = %d\n", m);
//             for(j = 0; j < CUMP_N_BATCHES_PER_MOD; j++) {
//                 printf("%d: %d\n", j, pref[j]);
//             }
//             printf("\n=======================\n");
        }
        pref += out_data_padding;

#endif // CUMP_USE_CRA_KERNEL
    }

    // now apply CRT to coefficients
#if CUMP_USE_CRA_KERNEL
    for(j = 0, pref = reference; j < n_real_pts; j++,
                pref += mods_padding) {

        compute_MR_digits_rev(pref, pref, 1, Mods, 1, InvMods, CUMP_N_MODULI);
    }

    return checkme(GPU_res, reference, CUMP_N_MODULI,
            mods_padding, n_real_pts, 0u, -1u);

#else
    return checkme(GPU_res, reference, n_real_pts,
           out_data_padding, CUMP_N_MODULI, 0u, -1u, true);

#endif // CUMP_USE_CRA_KERNEL

}

#if CUMP_COMPILE_RESULTANTS_KERNEL

//! reduces polynomial coefficients modulo set of primes, eliminates those
//! primes which divide leading coefficients of either polynomial
//! copies valid moduli to \c pconst
//! (must have enough space for \c CUMP_N_MODULI)
//! and duplicates to \c Mods, \c InvKs and \c Mus
bool GPU_resultant::RNS_conversion(const MPZ_vector_2& fv,
        const MPZ_vector_2& gv, unsigned *U, void *pconst_,
            unsigned *Mods, unsigned *InvKs, unsigned *Mus) {

    MOD_entry *pconst = (MOD_entry *)pconst_;
    unsigned *pU, *ppU;
    // first make sure that leading coeffs do not vanish
    const MPZ_vector_1& lcf = fv[nu], lcg = gv[nv];
    // addresses to save leading coeffs
    unsigned *p_lcf = U + nu, *p_lcg = U + (deg_x1 + 1) * max_nr / 2 + nv,
        n_mods(0);
    
    mpz_t rz;
    mpz_init(rz);

    GPU_device_manager& obj = GPU_device_manager::instance();
    unsigned r_stride = nr_sz / CUMP_N_MODULI, i, j;

    for(i = 0; i < obj.mod_table.size() && n_mods < CUMP_N_MODULI; i++) {

        const MOD_entry& me = obj.mod_table[i];
        unsigned m = me.m;
        bool zero_lc = true;
        for(j = deg_x1, pU = p_lcf; (int)j >= 0; j--, pU += max_nr/2) {

            if(lcf.size() <= j)
                continue;
            mpz_mod_ui(rz, (mpz_srcptr)&lcf[j], m);
            pU[0] = mpz_get_ui(rz);
            zero_lc &= (pU[0] == 0);   // count non-zero residues
        }
        if(zero_lc) {
            CUMP_out(m << ": vanishing lcoeff(f), skipping..\n");
            continue;
        }

        for(j = deg_x2, pU = p_lcg, zero_lc = true; (int)j >= 0; j--,
                pU += max_nr/2) {

            if(lcg.size() <= j)
                continue;
            mpz_mod_ui(rz, (mpz_srcptr)&lcg[j], m);
            pU[0] = mpz_get_ui(rz);
            zero_lc &= (pU[0] == 0);   // count non-zero residues
        }
        if(zero_lc) {
            CUMP_out(m << ": vanishing lcoeff(g), skipping..\n");
            continue;
        }
        // if we succeeded => advance p_lcf & p_lcg pointers
        p_lcf += r_stride, p_lcg += r_stride;
        pconst[n_mods] = me;
        Mods[n_mods] = m, Mus[n_mods] = me.mu;
        InvKs[n_mods] = (unsigned&)me.invk;

#if CUMP_USE_32BIT_MODULI_SET // save another half of invk here
        InvKs[n_mods + mods_padding] = (unsigned&)me.invm;
#endif
        n_mods++;
    }
    mpz_clear(rz);

    if(n_mods < CUMP_N_MODULI) {
        CUMP_out("RNS_conversion: moduli set exhausted..\n");
        return false;
    }

    // NOTE: leading coeffs are already computed: no need to recompute them
    for(j = deg_x1, pU = U; (int)j >= 0; j--, pU += max_nr/2) {
        for(i = 0, ppU = pU; i < nu; i++, ppU++) {
            if(fv[i].size() > j) {
                convert_to_RNS(ppU, r_stride, (mpz_srcptr)&fv[i][j],
                        Mods, CUMP_N_MODULI, 1);
            }
        }
    }

    for(j = deg_x2; (int)j >= 0; j--, pU += max_nr/2) {
        for(i = 0, ppU = pU; i < nv; i++, ppU++)
            if(gv[i].size() > j) {
                convert_to_RNS(ppU, r_stride, (mpz_srcptr)&gv[i][j],
                        Mods, CUMP_N_MODULI, 1);
            }
    }
    return true;
}

//! recovers large integer coefficients of a polynomial given a set of
//! Mixed-radix (MR) digits fro each coefficient
//! NOTE: MPZ_vector_1 must be disposed of after copying data back
bool GPU_resultant::RNS_recover(const unsigned *R, const unsigned *Mods,
            MPZ_vector_1& out) {

    const unsigned *pU, *pref;
    unsigned i, j;
//     typedef CORE::BigInt Integer;

    out = MPZ_vector_1(n_real_pts);
#if 0
    for(i = 0, pU = R; i < n_real_pts; i++, pU += mods_padding) {

        mpz_t pt;
        mpz_init_set_ui(pt, pU[0]);

        if(pU[0] > Mods[0]/2)
            mpz_sub_ui(pt, pt, Mods[0]);

        pref = pU + 1;
        for(j = 1; j < CUMP_N_MODULI; j++, pref++) {
            mpz_mul_ui(pt, pt, Mods[j]);

            mpz_t zz;
            mpz_init_set_ui(zz, pref[0]);
//             if(pref[0] > Mods[j]/2)
//                 mpz_sub_ui(zz, zz, Mods[j]);
            mpz_add(pt,pt,zz);
            mpz_clear(zz);
//            mpz_add_ui(pt, pt, pref[0]);
        }

//         mpz_init2(pt, 2 * BITS_PER_LIMB);
//         pt->_mp_d[0] = pU[0];
//         pt->_mp_d[1] = pU[1];
//         pt->_mp_size = 2;
//         pref = (limb *)pU + 2;
//         for(j = 2; j < CUMP_N_MODULI; j += 2, pref += 2) {
// 
//             mpz_t zz;
//             mpz_init_set_ui(zz, Mods[j]); // can be precomputed
//             mpz_mul_ui(zz, zz, Mods[j + 1]);
// 
//             mpz_mul(pt, pt, zz);
//             zz->_mp_d[0] = pref[0];
//             zz->_mp_d[1] = pref[1];
//             mpz_add(pt, pt, zz);
// 
//             mpz_clear(zz);
//         }

        std::cout << i << ": " << Integer(pt) << "\n";
        if(Integer(pt) != Integer(0)) {
//             std::cout << "wrong: " << i << ": " << Integer(pt) << "\n";
        }
        mpz_clear(pt);
    }
#else
 
    mpz_t *Zs = new mpz_t[CUMP_N_MODULI + 2],
            *zzs = Zs, *mzs = zzs + (CUMP_N_MODULI + 1) / 2;
    
    for(j = 0; j < (CUMP_N_MODULI + 1) / 2; j++) { // initially allocate space
                                             // for two words
        mpz_init2(zzs[j], sizeof(uint64) * 8);
        mpz_init2(mzs[j], sizeof(uint64) * 8);
    }

    const unsigned N_LIMBS = sizeof(uint64) * 8 / GMP_NUMB_BITS;
    for(i = 0, pU = R; i < n_real_pts; i++, pU += mods_padding) {

        unsigned n = CUMP_N_MODULI, n_iters = (n + 1) / 2;

        const unsigned *ppU, *pMods;
        for(j = 0, ppU = pU, pMods = Mods; j < n_iters; j++,
                    ppU += 2, pMods += 2) {

            uint64 x = ppU[0], y = pMods[0];
            if(2*j + 1 < n) {
                x = x * pMods[1] + (uint64)ppU[1];
                y = y * pMods[1];
            }
            if(x <= 0xffffffff) {
                mpz_set_ui(zzs[j], (unsigned)x);
            } else {
                memcpy(zzs[j]->_mp_d, &x, sizeof(uint64));
                zzs[j]->_mp_size = N_LIMBS;
            }
            memcpy(mzs[j]->_mp_d, &y, sizeof(uint64));
            mzs[j]->_mp_size = N_LIMBS;
        }
        if(pU[0] > Mods[0]/2) {
         // reverse zzs[0] in case of negative sign
            mpz_sub(zzs[0], zzs[0], mzs[0]);
        }

        //! process in a tree-like fashion; first we have:
        //! mzs[j] = m[2*j] * m[2*j + 1]
        //! zzs[j] = x[2*j + 1] + m[2*j + 1] * x[2*j]
        //! then:
        //! mzs[j] = m[4*j] * m[4*j + 1] * m[4*j + 2] * m[4*j + 3]
        //! zzs[j] = ... so on until we have only on a single "root node"
        while(n_iters != 1) {
            n = n_iters, n_iters = (n_iters + 1) / 2;

            mpz_t *pzzs = zzs, *pmzs = mzs;
            for(j = 0; j < n_iters; j++, pzzs += 2, pmzs += 2) {
                // # valid indices: 0..n_iters-1
                if(j*2 + 1 < n) {
                    mpz_mul(pzzs[0], pzzs[0], pmzs[1]);
                    mpz_add(zzs[j], pzzs[1], pzzs[0]);
                    mpz_mul(mzs[j], pmzs[1], pmzs[0]);
                } else { // otherwise just copy from previous level
                    mpz_set(zzs[j], pzzs[0]);
                    mpz_set(mzs[j], pmzs[0]);
                }
            }
        } // while

//         CUMP_out(i << ": " << Integer(zzs[0]) << "\n");
        mpz_init_set((mpz_ptr)&out[i], zzs[0]);
    }
    for(j = 0; j < CUMP_N_MODULI / 2; j++) {
        mpz_clear(zzs[j]);
        mpz_clear(mzs[j]);
    }

    delete []Zs;
#endif
    return true;
}
#endif // CUMP_COMPILE_RESULTANTS_KERNEL

bool GPU_resultant::internal_compute(const MPZ_vector_2& fv,
        const MPZ_vector_2& gv, MPZ_vector_1& r, unsigned deg_y1_,
        unsigned deg_y2_, unsigned deg_x1_, unsigned deg_x2_,
        unsigned low_deg_, unsigned high_deg_, unsigned bits_) {

    unsigned run = 0; // run = 1: tried swapped polynomials
                      // run = 2: tried increased # of excess points
    pts_overrun = 8;
    const MPZ_vector_2 *pf = &fv, *pg = &gv;

    GPU_device_manager& obj = GPU_device_manager::instance();
    bool ret = false;
    obj.push_CUDA_context(); // attach to the current context

    if(!setup(deg_y1_, deg_y2_, deg_x1_, deg_x2_, low_deg_, high_deg_, bits_))
         goto Lexit;

    ret = run_gpu_part(*pf, *pg, r);

// at this point we can either restart the algorithm it with increased
// # of excess points or swap polynomials..
//     CUMP_out("restarting with increased # of points..\n");
//     pts_overrun *= 2;
Lexit:
    obj.pop_CUDA_context();
    return ret;
}

#if CUMP_COMPILE_DEBUG_KERNEL

struct Args {
    unsigned CUMP_N_MODULI;
    unsigned CUMP_N_BATCHES_PER_MOD;
    bool no_cpu_run;
    unsigned nu;
    unsigned nv;
    char in_file[64];
    bool write_to_file;
    unsigned nentries;
    unsigned nbits;
    unsigned fg_idx;
    unsigned deg_x1;
    unsigned deg_x2;
    unsigned niters;
};

void print_args(const Args& a);

void parse(int argc, char **argv, Args& a);

template < class Poly_ >
bool read_from_file(const char *filename,  Poly_& f, Poly_& g) {

    CGAL::Polynomial_parser_d< Poly_ > parser;
    std::ifstream in(filename);
    int nread = 0;

    if(!in) {
        std::cerr << "Unable to open: " << filename << "\n";
        return false;
    }
    while(!in.eof()) {
        Poly_ p;
        std::string s;
        std::getline(in, s);

        if(s.length() == 0)
            continue;
        if(!parser(s, p)) {
          //std::cerr << "Parser error, trying another format..\n";
            try {
                std::stringstream ss(s);
                ss >> p;
            } catch(...) {
                std::cerr << "Invalid format of polynomial..skipping..\n";
                continue; 
            }
        }
        if(nread == 0) {
            f = p, nread = 1;
        } else if(nread == 1) {
            g = p, nread = 2;
            break;
        }
    }
    if(nread == 1) {
        g = CGAL::differentiate(f);
        nread++;
    }

    return (nread == 2);
}


#endif // CUMP_COMPILE_DEBUG_KERNEL

//! ./gpu_resultant --nu=31 --nv=31 --degx1=20 --degx2=15 --nbits=16 --nentries=1111 --usefile --idx=1
//! ./gpu_resultant --nu=27 --nv=20 --degx1=7 --degx2=7 --nbits=16
bool GPU_resultant::debug_run(int argc, char **argv) {

#if CUMP_COMPILE_DEBUG_KERNEL

    srand48(time(NULL));
    gmp_randinit_mt(rands);
    gmp_randseed_ui(rands, time(NULL));

    typedef CORE::BigInt Integer;
    typedef CGAL::Polynomial< Integer > Poly_int_1;
    typedef CGAL::Polynomial< Poly_int_1 > Poly_int_2;

//     return quick_run();
    Args a;
    a.CUMP_N_MODULI = -1u;
    a.CUMP_N_BATCHES_PER_MOD = -1u;
    a.no_cpu_run = false;
    a.nentries = 100;
    a.write_to_file = false;
    a.fg_idx = 1;
    a.nu = -1u, a.nv = -1u;
    a.deg_x1 = 12, a.deg_x2 = 11;
    a.nbits = 10;
    a.in_file[0] = 0;
    a.niters = 1;

    CGAL::set_pretty_mode(std::cout);
    CGAL::set_pretty_mode(std::cerr);

    parse(argc, argv, a);
    print_args(a);

    unsigned deg_y1 = 34, deg_y2 = 31;
    if(a.nu == -1u || a.nv == -1u) { 
        a.nu = deg_y1; a.nv = deg_y2;
    }
    CUMP_N_MODULI = a.CUMP_N_MODULI;
    CUMP_N_BATCHES_PER_MOD = a.CUMP_N_BATCHES_PER_MOD;
    no_cpu_run = a.no_cpu_run;
    write_to_file = a.write_to_file;
    fg_idx = a.fg_idx;

    CGAL::set_pretty_mode(std::cerr);
    CGAL::set_pretty_mode(std::cout);

    Poly_int_2 f, g;

#if 0
    CGAL::Polynomial_parser_d< Poly_int_2 > parser;    
    const char *fstr, *gstr;
    fstr = "x^10*y^11-xy+x-1";
    gstr = "x^14y^43+y^29+2";
/*    fstr = "x^35*y^35-1"; gstr = "x^33*y^33-1";*/
    parser(fstr, f);
    parser(gstr, g);
#else
    
    if(a.in_file == "" || !read_from_file(a.in_file, f, g)) {
        typedef Rand_coeff< Integer > Rand_coeff;
        f = CGAL::internal::generate_sparse_random_poly2< Integer, Rand_coeff >
            (a.nu, a.deg_x1, a.nentries, a.nbits);
        g = CGAL::internal::generate_sparse_random_poly2< Integer, Rand_coeff >
            (a.nv, a.deg_x2, a.nentries, a.nbits);
    }
#endif

    if(write_to_file) {
/*        std::ofstream out("maple_polys_in", std::ios::app);
        CGAL::set_pretty_mode(out);
        out << "\n#=====================================================\n";
        out << "f" << fg_idx << " := " << f << ":\n\n";
        out << "g" << fg_idx << " := " << g << ":\n\n";*/
        std::ofstream out(benchmark_file, std::ios::app);
        out << "bits: " << a.nbits << " ";
    }

    Poly_int_2 _f = f, _g = g;
    if(f.degree() > g.degree())
         std::swap(f, g);

    CUMP_out("f := " << f << ":\n\n")
    CUMP_out("g := " << g << ":\n\n")

    unsigned low_deg_, high_deg_, deg_x1_, deg_x2_, bits_;
    compute_resultant_bounds(f, g, low_deg_, high_deg_, bits_);

    MPZ_vector_2 fv, gv;
    construct_mpz_vector_2(f, fv, deg_x1_);
    construct_mpz_vector_2(g, gv, deg_x2_);

//!*******************************************************************
//      bits_ = 30;
//!*******************************************************************

    bool succeeded = true;
    MPZ_vector_1 r;
    while(a.niters--) {
        if(!internal_compute(fv, gv, r, f.degree(), g.degree(), deg_x1_,
                 deg_x2_, low_deg_, high_deg_, bits_)) {
            succeeded = false;
        }
    }

    if(succeeded) {
#if CUMP_MRC_HOST_RECOVER
        Poly_int_1 res =
               CGAL::internal::construct_polynomial_from_mpz< Integer >(r);

        Poly_int_1 truth = CGAL::resultant(_f, _g);

        if(truth != res) {
            CUMP_out("\nWrong resultant!!\n" << res << "\n" <<
            truth << "\n");
            CUMP_out("\nDIFF: " << (res - truth) << "\n");
        } else
            CUMP_out("CORRECT!!\n");
        CGAL::internal::dispose_mpz_vector(r);
#endif // CUMP_MRC_HOST_RECOVER
    }

    return succeeded;

#endif // CUMP_COMPILE_DEBUG_KERNEL
}

bool GPU_resultant::quick_run() {

#if 0

    srand48(time(NULL));
    gmp_randinit_mt(rands);
    gmp_randseed_ui(rands, time(NULL));

    unsigned timerCPU;
//     cutCreateTimer(&timerCPU);

    unsigned i, j;
    
//     unsigned deg_x1 = 5, deg_x2 = 7, nr = nu + nv,
//         nentries = 11, nbits = 23;

   CGAL::Polynomial_parser_d< Poly_int_2 > parser;

    const char *fstr = "x^35*y^35-1";
    const char *gstr = "x^33*y^33-1";

    CGAL::set_pretty_mode(std::cerr);
    CGAL::set_pretty_mode(std::cout);

    Poly_int_2 f, g;
    parser(fstr, f);
    parser(gstr, g);

    CUMP_out("f := " << f << ":\n\n")
    CUMP_out("g := " << g << ":\n\n")

    MPZ_vector_2 fv, gv;
    unsigned deg_x1_, deg_x2_;
    construct_mpz_vector(f, fv, deg_x1_);
    construct_mpz_vector(g, gv, deg_x2_);

    deg_x1_ = 35; deg_x2_ = 33;
    nu = 35; nv = 33;

 setup(f.degree(), g.degree(), deg_x1_, deg_x2_, 0, 100,
                 48);

    unsigned *pU, *ppU, *U = (unsigned *)CPU_mem_ptr,
        *R = U + nr_sz, *reference = R + prod_sz * batch_size, 
        *Mods = reference + prod_sz * batch_size,
        *InvKs = Mods + mods_padding,
        *Mus = InvKs + mods_padding,
        *InvMods = Mus + mods_padding;

    unsigned n_skips = 0;
    umod_t m((1<<24));

//     CUMP_N_MODULI = 3;
//     CUMP_N_BATCHES_PER_MOD = 10;
//     Mods = new unsigned[CUMP_N_MODULI];

    for(i = 0; i < CUMP_N_MODULI; i++) {
        get_next_mod(m, n_skips);
        zmod::set_mod(m);
        printf("%d: modulus: %d\n", i+1, m);
        Mods[i] = m;
        m--;
    }
    
  unsigned r_stride = nr_sz / CUMP_N_MODULI;
// 
memset(U, 0, mem_size_nr);

    for(j = deg_x1, pU = U; (int)j >= 0; j--, pU += max_nr/2) {
        for(i = 0, ppU = pU; i <= nu; i++, ppU++) {
            if(f[i].degree()+1 > j) {
                convert_to_RNS(ppU, r_stride, f[i][j].get_mp(),
                        Mods, CUMP_N_MODULI, 1);
            }
        }
    }

    for(j = deg_x2; (int)j >= 0; j--, pU += max_nr/2) {
        for(i = 0, ppU = pU; i <= nv; i++, ppU++) 
            if(g[nv - i].degree()+1 > j) { // reversed order !!
                convert_to_RNS(ppU, r_stride, g[nv - i][j].get_mp(),
                        Mods, CUMP_N_MODULI, 1);
            }
    }

  for(i = 0, pU = U; i < CUMP_N_MODULI; i++,
            pU += nr_sz / CUMP_N_MODULI) {

        m = Mods[i];

        zmod::set_mod(m);

    Vector_2 fp(nu + 1, Vector_1(1, zmod(0))),
            gp(nv + 1, Vector_1(1, zmod(0)));

//     fp[nu] =  Vector_1(nu + 1, zmod(0));
//     fp[nu][nu] = zmod(1);
//     fp[0][0] = -zmod(1);
//     gp[nv] =  Vector_1(nv + 1, zmod(0));
//     gp[nv][nv] = zmod(1);
//     gp[0][0] = -zmod(1);

        std::vector< unsigned > xs(CUMP_N_BATCHES_PER_MOD),
                 ys(CUMP_N_BATCHES_PER_MOD);

 printf("m: %d\n", m);

        unsigned jidx = 0;
        for(j = 0; j < CUMP_N_BATCHES_PER_MOD; j++) {
            std::vector< zmod > v1m(nu + 1), v2m(nv + 1);
            zmod pt(j + 1);//CUMP_N_BATCHES_PER_MOD - j);

            poly_eval2_pure((const limb *)pU, pt, deg_x1, nu, max_nr/2,
                 v1m.begin());
//             !NOTE: here order is reversed !!!
            poly_eval2_pure((const limb *)pU + max_nr/2 * (deg_x1 + 1), pt,
                     deg_x2, nv, max_nr/2, v2m.rbegin());    
//         poly_eval2_pure(fp, pt, v1m.begin());
//         poly_eval2_pure(gp, pt, v2m.begin());

        
//     printf("pt: %d\n", pt.x);
// // print_vector(fp[0]);
// 
//     print_vector(v1m);
//     print_vector(v2m);
    
//     return 1;
// exit(1);

//TODO: simplify resultant algorithm w.r.t. to the new developments: only
// one modular inverse 

            bool failed;
            zmod det = resultant_modular_sylvester(v1m, v2m, failed);
        printf("res: %d\n", det.x);
// exit(1);
            if(failed) {
                printf("failed: modulus (%d): %x; point: %d\n", i, m, pt.x);
                continue; // we do not collect this point
            }
            xs[jidx] = pt.x, ys[jidx] = det.x; jidx++;
//              pref[j] = det.x;
        }

    }
    return true;

#else
    return true;
#endif // CUMP_COMPILE_TEST_KERNEL
}


#if CUMP_COMPILE_DEBUG_KERNEL

void print_args(const Args& a) {

    printf("CUMP_N_MODULI: %d; CUMP_N_BATCHES_PER_MOD: %d;\n"
      "no_cpu_run: %d nentries: %d; nbits: %d\n",
       a.CUMP_N_MODULI, a.CUMP_N_BATCHES_PER_MOD, a.no_cpu_run,
                a.nentries, a.nbits);
    if(a.nu != -1u && a.nv != -1u)
        printf("nu: %d; nv: %d\n", a.nu, a.nv);
}

void parse(int argc, char **argv, Args& a) {

    static option parser_opts[] = {
        {"help", 0, 0, 0},          // 0
        {"nmods", 1, 0, 0},         // 1
        {"nbatches", 1, 0, 0},      // 2
        {"nocpurun", 0, 0, 0},      // 3
        {"nu", 1, 0, 0}, {"nv", 1, 0, 0}, // 4, 5
        {"usefile", 0, 0, 0},       // 6
        {"nentries", 1, 0, 0},      // 7
        {"nbits", 1, 0, 0},         // 8
        {"idx", 1, 0, 0},           // 9
        {"degx1", 1, 0, 0}, {"degx2", 1, 0, 0}, // 10, 11
        {"in", 1, 0, 0}, // 12
        {"niters", 1, 0, 0}, // 13
        {0, 0, 0, 0}
    };
    // NOTE: also take into account that 3exec script takes up to 9 options !!

    while(1) {

        int c, idx;
        c = getopt_long(argc, argv, "h", parser_opts, &idx);
        if(c == -1)
            break;

        if(c == '?') {
            exit(1);
        } else if(c == 'h' || (c == 0 && idx == 0)) {
            printf("\nUSAGE: %s", argv[0]);
            int i = 0;
            while(parser_opts[i].name != 0)
                printf(" [%s=]", parser_opts[i++].name);
            printf("\n\n");

            exit(1);
        }
        if(c != 0) {
            printf ("?? getopt returned character code 0%o ??\n", c);
            exit(1);
        }

        int argi;
        if(parser_opts[idx].has_arg && idx != 12)
            argi = atoi(optarg);
        switch(idx) {
        case 1:
            a.CUMP_N_MODULI = argi;
            break;
        case 2:
            a.CUMP_N_BATCHES_PER_MOD = argi;
            break;
        case 3:
            a.no_cpu_run = true;
            break;
        case 4:
            a.nu = argi;
            break;
        case 5:
            a.nv = argi;
            break;
        case 6:
            a.write_to_file = true;
            break;
        case 7:
            a.nentries = argi;
            break; 
        case 8:
            a.nbits = argi;
            break;
        case 9:
            a.fg_idx = argi;
            break;
        case 10:
            a.deg_x1 = argi;
            break;
        case 11:
            a.deg_x2 = argi;
            break;
        case 12:
            strncpy(a.in_file, optarg, 63);
            a.in_file[63] = 0;
            break;
        case 13:
            a.niters = argi;
            break;
        default:
            printf("unrecognized option..\n");
            exit(1);
        }
    }
}
#endif // CUMP_COMPILE_DEBUG_KERNEL

} // namespace internal 

} // namespace CGAL 
