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
    
#include <iostream>
#include <getopt.h>

#define CGAL_MODULAR_FILTER_OFF

#include <include/macros.h>
#include <include/const_mem_defs.h>
#include <include/device_manager.h>

#include <CGAL/GPU_algorithms/gcd_algorithm.h>
#include <CGAL/GPU_algorithms/GPU_algorithm_templates.h>

#include <CGAL/Timer.h>
#if CUMP_COMPILE_DEBUG_KERNEL
#include <CGAL/Algebraic_kernel_d/Polynomial_parser_d.h>
#endif 

#include <include/modular_arithm.h>
#include <include/misc.h>
#include <matrix_algebra.C>

#if CUMP_COMPILE_DEBUG_KERNEL
// #include <common/block_updates_GB.C>
#include <block_updates_GJG_new.C>
#endif    

namespace CGAL {

namespace internal {

void gcd_abort_handler(int sig)
{
    GPU_gcd& obj = GPU_gcd::instance();
    obj.~GPU_gcd();
}

GPU_gcd::GPU_gcd() {

    host_static_setup();
    GPU_device_manager& obj = GPU_device_manager::instance();
    obj.push_CUDA_context();
    device_static_setup();
    obj.pop_CUDA_context();
}

#if CUMP_COMPILE_DEBUG_KERNEL
CGAL::Timer tm_rns_cvt, tm_rns_rcv;
#endif

GPU_gcd::~GPU_gcd() {

#if CUMP_COMPILE_DEBUG_KERNEL
    std::cout << "\n########### tm_rns_cvt: " << tm_rns_cvt.time()*1000.f <<
       " ms; tm_rns_rcv: " << tm_rns_rcv.time()*1000.f << " ms\n";
#endif

    free_device_mem();
}


void GPU_gcd::host_static_setup() {

    CPU_mem_ptr = 0, DEV_mem_ptr = 0;
    word_size = sizeof(unsigned), c_stride = 4;   // constant memory stride
    aux_cpu_mem = 0, aux_dev_mem = 0;
    CPU_bytes_alloced = 0, DEV_bytes_allocated = 0;
    n_moduli = -1u, no_cpu_run = false;
    chunk_sz = 0, n_blocks = 0, gcd_real_sz = 0;
    n_streams = 1; // defaults to one stream

    strcpy(benchmark_file, "test_out");

    // initializes device and loads moduli set
    GPU_device_manager& obj = GPU_device_manager::instance();
    obj.register_callback((GPU_DEVICE_CALLBACK)gcd_abort_handler);
}

//! allocates data for univariate polynomials of degree \c deg_f_ and \c deg_g_
//! with estimated # of bits for gcd given by \c bits_
bool GPU_gcd::setup(unsigned deg_f_, unsigned deg_g_,
        unsigned bits_) {

// observe that # of coefficients is nu + 1 !!!
    nu = deg_f_, nv = deg_g_;

    //! requirenment: BlockSz/2 <= nu + nv <= BlockSz and nu >= nv
    //! and: all coefficients of f[nu..0] can be read by a single
    //! block read using BlockSz threads, i.e., max(nu) = BlockSz - 1

    if(/*nu + nv < BlockSz/2 || nu >= BlockSz ||*/ nu < nv) {
        CUMP_out2("Incorrect poly degrees: %d %d\n", nu, nv);
        return false;
    }

    GPU_device_manager& obj = GPU_device_manager::instance();

    n_moduli = bits_to_moduli(bits_, obj.mod_bits);
    if(n_moduli == 0) {
        std::cerr << "FATAL: moduli set exhausted !!\n";
        return false;
    }
#if CUMP_COMPILE_DEBUG_KERNEL
//     n_moduli = 1;
#endif
    CUMP_out2("============ nu: %d; nv: %d; bits: %d; n_moduli: %d\n",
            nu, nv, bits_, n_moduli);

    if(write_to_file) {
        std::ofstream out(benchmark_file, std::ios::app);
        out << "nu: " << nu << " nv: " << nv << " bits: " << bits_ <<
            " n_mods: " << n_moduli;
    }

    mods_padding = (n_moduli + 15) & ~15; // next divisible by 16
    // gcd_data_padding(); this supposedly will compute the estimated bounds
    //! NOTE: we have to allocate memory for nu+1 & nv+1 coeffs, resp.
    max_nu = (nu + 1 + 15) & ~15, max_nv = (nv + 1 + 15) & ~15;

#if CUMP_USE_PGCD_KERNEL
    // compute offsets
    nu_ofs4 = ((nu + 4) & ~3) - (nu + 1), nv_ofs4 = ((nv + 4) & ~3) - (nv + 1);
#else
    nu_ofs4 = 0, nv_ofs4 = 0; // no zero padding is needed for block kernel
#endif

    data_padding = max_nu + max_nv;     // data padding per each modulus

    QR_gcd_kernel_data_padding(); // get the chunk size, so that we can
                                  // compute the amount of memory for in-outs

    aux_dev_mem = 0;
    // required for MRC algorithm (not for the main algorithm)
    aux_cpu_mem = mods_padding * 5; // (m, mu, inv_mod, invk)
    // total amount of input data: mods_padding is used to save the residues
    // of gcd(lc_f, lc_g)
//     data_sz = data_padding * n_moduli + mods_padding;
    //! NOTE: later we can change data_sz appropriately 
//     data_sz += limbs_f * (nu + 1) + limbs_g * (nv + 1);
    data_sz = limbs_f * (nu + 1) + limbs_g * (nv + 1) + mods_padding;
    aux_dev_mem += data_padding * n_moduli;

    // ofs is needed to save the size of gcd computed; do we need chunk_sz ??
    gcd_sz = max_nv;
    //! gcd_sz + 1 because we reserve first 'mods_padding' elements to collect
    //! degrees of gcd images to sort out unlucky primes
    prod_sz = (gcd_sz + 1) * mods_padding;

#if CUMP_USE_GCD_BLOCK_KERNEL
    block_io_sz = (n_blocks * chunk_sz*2*4 + chunk_sz*4) * n_moduli;
    unsigned _ = mods_padding * (nv + 1);
    if(block_io_sz * 2 < _) { // enlarge io_sz if not enough memory for MRC
        block_io_sz = (_ / 2 + 15) & ~15;
    }
    // we allocate 2 batches of size block_io_sz: for devIn and devOut
    // mods_padding to save inverse of leading coefficients
    aux_dev_mem += mods_padding * 5 + mods_padding + block_io_sz * 2;
#else
    block_io_sz = 0;
    // the last one is for data sharing
    aux_dev_mem += mods_padding * 5 + ((n_blocks - 1) * chunk_sz +
        chunk_sz * 4) * n_moduli + data_sz;
    //!prod_sz = max_nv * n_moduli;     // gcd's degree is not larger than the
    // NOTE: in this case the memory for MRC algorithm is not set..
#endif

    CUMP_out2("##### gcd_sz: %d; data_sz: %d; prod_sz: %d; n_blocks: %d; "
        "block_io_sz: %d\n", gcd_sz, data_sz, prod_sz, n_blocks, block_io_sz);
//     CUMP_out2("### nu_ofs4: %d; nv_ofs4: %d\n", nu_ofs4, nv_ofs4);

    if(!alloc_device_mem())
        return false; 

    return true;
}

bool GPU_gcd::internal_compute(const MPZ_vector_1& fv, const MPZ_vector_1& gv,
             MPZ_vector_1& r, unsigned bits_) {

//     if(n_moduli != -1u)
//         bits_ = n_moduli * 31 - 1;
    GPU_device_manager& obj = GPU_device_manager::instance();
    bool ret = false;

    obj.push_CUDA_context(); // attach to the current context

    limbs_f = mpz_vector_bitlength(fv);
    limbs_g = mpz_vector_bitlength(gv);

//     CUMP_out2("limbs_f: %d; limbs_g: %d\n", limbs_f, limbs_g);
    unsigned N_BITS = GMP_NUMB_BITS, 
        FACTOR = N_BITS / (sizeof(unsigned) * 8);
    
    limbs_f = ((limbs_f + N_BITS - 1) / N_BITS * FACTOR + 16) & ~15;
    limbs_g = ((limbs_g + N_BITS - 1) / N_BITS * FACTOR + 16) & ~15;

    // padding to next divisible by 16 (+1 word to keep the actual size of cf)
//     limbs_f = ((limbs_f + N_BITS - 1) / N_BITS + 16) & ~15;
//     limbs_g = ((limbs_g + N_BITS - 1) / N_BITS + 16) & ~15;
    CUMP_out2("limbs_f: %d; limbs_g: %d\n", limbs_f, limbs_g);

#if CUMP_GCD_USE_COPRIME_CHECK
    if(limbs_f >= 7 || limbs_g >= 7) {

    if(!setup(fv.size()-1, gv.size()-1, 24))
        goto Lexit; 

    if(!run_gpu_part(fv, gv, r))
        goto Lexit;

    if(r.size() == 1) { // found coprime polynomials
        CUMP_out2("========= COPRIME CHECK succeeded ========\n")
        ret = true;
        goto Lexit;
    }

    dispose_mpz_vector(r); // otherwise continue as usual
    }
#endif

    if(!setup(fv.size()-1, gv.size()-1, bits_))
        goto Lexit;

    ret = run_gpu_part(fv, gv, r);
Lexit:
    obj.pop_CUDA_context(); // attach to the current context
    return ret;
}

bool GPU_gcd::run_gpu_part(const MPZ_vector_1& fv,
            const MPZ_vector_1& gv, MPZ_vector_1& r) {
//! U - host input data
//! R - GPU output
//! reference - host reference solution
//     new_host_mem = mem_size_in + mem_size_prod*2 + mem_size_cpu_aux;
    unsigned *U = (unsigned *)CPU_mem_ptr, *R = U + data_sz,
        *reference = R + prod_sz, *Mods = reference + prod_sz,
        *Mus = Mods + mods_padding, *InvMods = Mus + mods_padding,
        *InvKs = InvMods + mods_padding;

// NOTE NOTE: make sure that 'reference' and 'R' use only allowed mem-size
// that is, 'prod_sz' elements

    static unsigned const_mem[CUMP_DEVICE_MODULI_SZ];
    //! reserve 4 entries for formal parameters
    unsigned *pconst = const_mem + ConstMemDataSz;

// TODO: better use macros for this:
// #define CU_BENCHMARK(x) x
#if CUMP_BENCHMARK_CALLS
    CGAL::Timer rolex;
    float mod_reduce_tm, mrc_cpu_tm;
    rolex.start();
#endif

//     tm_rns_cvt.start();
    memset(U, 0, mem_size_in);
    //! reduce polynomials modulo set of primes and eliminate ``bad'' primes
    if(!RNS_conversion(fv, gv, U, (void *)pconst, Mods, InvKs, Mus)) {
        return false;
    }
#if CUMP_BENCHMARK_CALLS
    rolex.stop();
    mod_reduce_tm = rolex.time() * 1000.0f;
#endif

//! mod_inverse_seq can be precomputed
// new bottleneck: mod_inverse_seq must be somehow precomputed !!
// because it has quadratic complexity !!!
    mod_inverse_seq(InvMods, Mods, n_moduli, 1);
//     tm_rns_cvt.stop();

    bool coprime = false;
    try {
        launch_kernel(const_mem, Mods, U, R, coprime);
    }
    catch(GPU_algorithm_exception) {
        CUMP_out("gpu algorithm failed..\n");
        return false;
    }

#if CUMP_MRC_HOST_RECOVER
#if CUMP_BENCHMARK_CALLS
    rolex.reset();
    rolex.start();
#endif

//     tm_rns_rcv.start();
#if CUMP_MRC_HOST_RECOVER
    if(coprime) {
        r = MPZ_vector_1(1);
        mpz_init_set_ui((mpz_ptr)&r[0], 1);
    } else {
        RNS_recover(R, Mods, r);
    }    
#endif
//     tm_rns_rcv.stop();

#if CUMP_BENCHMARK_CALLS
    rolex.stop();
    mrc_cpu_tm = rolex.time() * 1000.0f;
#endif

#endif // CUMP_MRC_HOST_RECOVER

#if CUMP_BENCHMARK_CALLS
    CUMP_out2("mod_reduce time: %f ms; CPU MRC time: %f ms\n",
            mod_reduce_tm, mrc_cpu_tm);

    if(write_to_file) {
        std::ofstream out(benchmark_file, std::ios::app);
        out << " mod_reduce: " << mod_reduce_tm << " mrc_cpu: " <<
            mrc_cpu_tm << "\n";
    }
#endif // CUMP_BENCHMARK_CALLS

    if(no_cpu_run)
        return true;

    return reference_solution(fv, gv, pconst, U, Mods, InvMods, reference, R);
}

//! verifies the result on the host: \c const_mem : constant memory pointer,
//! \c U - input set (polynomial coeffs), \c reference - memory set for
//! reference solution, \c GPU_res - GPU result
bool GPU_gcd::reference_solution(const MPZ_vector_1& fv,
        const MPZ_vector_1& gv, const unsigned *const_mem,
        const unsigned *GCDlc, const unsigned *Mods, const unsigned *InvMods,
            unsigned *reference, const unsigned *GPU_res) {

#if CUMP_COMPILE_DEBUG_KERNEL
    typedef std::vector< zmod > Vector_mod_1;
    
    const unsigned *pconst = const_mem, *p_GCDlc = GCDlc;
    unsigned *pref = reference, *ppref;
    unsigned i, j, m;

    memset(pref, 0xff, prod_sz * sizeof(unsigned));

    Vector_mod_1 mod_f(nu + 1), mod_g(nv + 1);
    std::vector< unsigned > out(nu + nv);

    mpz_t mp_r;
    mpz_init(mp_r);
    unsigned gcd_host_sz = 0; // size of gcd
    for(i = 0; i < n_moduli; i++, pconst += c_stride, p_GCDlc++) {
        unsigned m = pconst[0];

//         CUMP_out2("%d: modulus: %d\n", i, m);
        zmod::set_mod(m);

        for(j = 0; j <= nu; j++) {
            mpz_mod_ui(mp_r, &fv[j], m);
            mod_f[j] = mpz_get_ui(mp_r);
        }
        for(j = 0; j <= nv; j++) {
            mpz_mod_ui(mp_r, &gv[j], m);
            mod_g[j] = mpz_get_ui(mp_r);
        }

        unsigned *ptr = pref, stride = 1;
#if !CUMP_RUN_GCD_ONLY          // save with mods_padding stride to allow
        stride = mods_padding;  // in-place data modification
        pref++;
#else
//! NOTE NOTE NOTE: need to fix this as now prod_sz is different !!!!!!!!!!

        pref += gcd_sz; /*chunk_sz*2*/
#endif

    /** **************************************************************/
    gcd_host_sz = pgcd_check< zmod >(mod_f, mod_g, ptr, stride);
    /** **************************************************************/

#if 0
        switch(chunk_sz) {
        case 32:
            gcd_host_sz = block_schur_QR_new< zmod, 32 >(mod_f, mod_g, ptr, stride);
            break;
        case 64:
            gcd_host_sz = block_schur_QR_new< zmod, 64 >(mod_f, mod_g, ptr, stride);
             break;
        case 128:
            gcd_host_sz = block_schur_QR_new< zmod, 128 >(mod_f, mod_g, ptr, stride);
            break;
        default:
            throw "WTF ???";
        }
#endif

#if !CUMP_RUN_GCD_ONLY
        unsigned div(mod_inverse(ptr[(gcd_host_sz-1)*stride], m));
        div = mul_mod(div, p_GCDlc[0], m); // mul with gcd of leading coeffs
        for(j = 0; j <= nv; j++, ptr += stride)
            ptr[0] = mul_mod(ptr[0], div, m);
#endif
    }
    mpz_clear(mp_r);

#if !CUMP_RUN_GCD_ONLY
    if(gcd_host_sz != gcd_real_sz) {
        CUMP_out2("******** gcd size differs: host: %d; gpu: %d\n", 
                gcd_host_sz, gcd_real_sz);
        return false;
    }
    for(j = 0, pref = reference; j < gcd_real_sz; j++, pref += mods_padding) {

        compute_MR_digits_rev(pref, pref, 1, Mods, 1, InvMods, n_moduli);
    }

    return checkme(GPU_res, reference, n_moduli,
            mods_padding, gcd_real_sz, 0u, -1u);

#else
    unsigned Mm = 0;//const_mem[0]; //! get the first modulus
//! NOTE NOTE: gcd_real_sz is set defined here !!!!
//! hence we use 'gcd_host_sz' instead
//! \c n elements per row stored with \c padding , # of rows: \c n_batches

    // need to compare only 'nu' elements
    return checkme(GPU_res + mods_padding, reference,
            gcd_sz, gcd_sz, n_moduli, Mm, -1u, false);

    unsigned useofs = 0; // CUMP_GCD_RAW_OFFSET
//     return checkme(GPU_res + mods_padding, reference,
//             gcd_host_sz, gcd_sz, n_moduli, Mm, -1u, false);

#endif // CUMP_RUN_GCD_ONLY
#else
  
    return true;
#endif
}

void print_limbs(unsigned *ptr, int n) {

    printf("[ ");
    if(n < 0)
        printf("-");
    for(int i = 0; i < std::abs(n); i++) {
        printf("%x ", ptr[i]);
    }
    printf(" ]\n");
}

//! reduces polynomial coefficients modulo set of primes, eliminates those
//! primes which divide leading coefficients of either polynomial
//! copies valid moduli to \c pconst (must have enough space for \c n_moduli)
//! and replicates to \c Mods, \c InvKs and \c Mus for MRC algorithm
bool GPU_gcd::RNS_conversion(const MPZ_vector_1& fv,
        const MPZ_vector_1& gv, unsigned *U, void *pconst_,
            unsigned *Mods, unsigned *InvKs, unsigned *Mus) {

    MOD_entry *pconst = (MOD_entry *)pconst_;
    // place to save residues of either polynomial
    unsigned *pgcd_lc = U, // residues of gcd(lcf,lcg)
        *plimbs = pgcd_lc + mods_padding; // raw data for coefficients
    mpz_srcptr lcf = &fv[nu], lcg = &gv[nv];
    unsigned n_m(0), i;
    
    mpz_t rz, mp_gcd;
    mpz_init(rz);
    mpz_init(mp_gcd);
    mpz_gcd(mp_gcd, lcf, lcg);

    GPU_device_manager& obj = GPU_device_manager::instance();

    for(i = 0; i < obj.mod_table.size() && n_m < n_moduli; i++) {

        const MOD_entry& me = obj.mod_table[i];
        unsigned m = me.m;
//         CUMP_out2("%d: modulus: %d %f\n", i, m, log((double)m)/log(2.0));

        if(mpz_divisible_ui_p((mpz_srcptr)lcf, m) != 0) {
            CUMP_out(m << ": vanishing lcoeff(f), skipping..\n");    
            continue;
        }

        if(mpz_divisible_ui_p((mpz_srcptr)lcg, m) != 0) {
            CUMP_out(m << ": vanishing lcoeff(g), skipping..\n");
            continue;
        }
        mpz_mod_ui(rz, mp_gcd, m);
        pgcd_lc[0] = mpz_get_ui(rz);
        pgcd_lc++;

        pconst[n_m] = me, Mods[n_m] = m, Mus[n_m] = me.mu;
        InvKs[n_m] = (unsigned&)me.invk;
#if CUMP_USE_32BIT_MODULI_SET // save another half of invk 
        InvKs[n_m + mods_padding] = (unsigned&)me.invm;
#endif 
        n_m++;
    }
    mpz_clear(rz);
    mpz_clear(mp_gcd);   

    if(n_m < n_moduli) {
        CUMP_out("RNS_conversion: moduli set exhausted..\n");
        return false;
    }

    const int FACTOR = (int)GMP_NUMB_BITS / (sizeof(unsigned) * 8);
    for(i = 0; i <= nu; i++, plimbs += limbs_f) {
        int sz = fv[i]._mp_size;
        
        // NOTE NOTE NOTE: alternatively we can also remove zero here
        // in case of 64 bits
        //!!  if(fv[i]._mp_d[sz-1] != 0)
        
        // NOTE NOTE NOTE: sz == 0 means the coefficient is 0 !! 
        memcpy(plimbs + 1, fv[i]._mp_d, std::abs(sz) * GMP_NUMB_BITS / 8);
        plimbs[0] = sz * FACTOR; // first word contains the data size (signed)

//         printf("size f %d: %d\n", i, (int)plimbs[0]);
//         print_limbs(plimbs + 1, plimbs[0]);
    }
    for(i = 0; i <= nv; i++, plimbs += limbs_g) {
        int sz = gv[i]._mp_size;
        memcpy(plimbs + 1, gv[i]._mp_d, std::abs(sz) * GMP_NUMB_BITS / 8);
        plimbs[0] = sz * FACTOR; // first word contains the data size (signed)

 /*       printf("size g %d: %d\n", i, (int)plimbs[0]);
        print_limbs(plimbs + 1, plimbs[0]);
 */   }

//     pf = U, pg = pf + max_nu; // leading elements already processed
//     pf += nu_ofs4, pg += nv_ofs4;
//     for(i = 0; i < nu; i++, pf++) {
//         convert_to_RNS(pf, data_padding, &fv[i], Mods, n_moduli, 1);
//     }
//     for(i = 0; i < nv; i++, pg++) {
//         convert_to_RNS(pg, data_padding, &gv[i], Mods, n_moduli, 1);
//     }

    return true;
}

//! recovers large integer coefficients of a polynomial given a set of
//! Mixed-radix (MR) digits fro each coefficient
//! NOTE: MPZ_vector_1 must be disposed of after copying data back
bool GPU_gcd::RNS_recover(const unsigned *R, const unsigned *Mods,
            MPZ_vector_1& out) {
#if CUMP_MRC_HOST_RECOVER

    const unsigned *pU, *pref;
    unsigned i, j;
    typedef CORE::BigInt Integer;

    out = MPZ_vector_1(gcd_real_sz);
    mpz_t *Zs = new mpz_t[n_moduli + 2],
            *zzs = Zs, *mzs = zzs + (n_moduli + 1) / 2;

    for(j = 0; j < (n_moduli + 1) / 2; j++) { // initially allocate space
                                             // for two words
        mpz_init2(zzs[j], sizeof(uint64) * 8);
        mpz_init2(mzs[j], sizeof(uint64) * 8);
    }

    const unsigned N_LIMBS = sizeof(uint64) * 8 / GMP_NUMB_BITS;
    for(i = 0, pU = R; i < gcd_real_sz; i++, pU += mods_padding) {

        unsigned n = n_moduli, n_iters = (n + 1) / 2;

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
#if CUMP_USE_PGCD_KERNEL
        mpz_init_set((mpz_ptr)&out[i], zzs[0]);
#else
        mpz_init_set((mpz_ptr)&out[gcd_real_sz - i - 1], zzs[0]);
#endif
    }
    for(j = 0; j < n_moduli / 2; j++) {
        mpz_clear(zzs[j]);
        mpz_clear(mzs[j]);
    }

    delete []Zs;
#endif // CUMP_MRC_HOST_RECOVER

    return true;
}

#if CUMP_COMPILE_DEBUG_KERNEL

struct Args {

    Args() : n_moduli(-1u), no_cpu_run(false),
        nu(-1u), nv(-1u), nentries(100), nbits(0),
        nstreams(-1u), ngcd(10), write_to_file(false) {
    }
    unsigned n_moduli;
    bool no_cpu_run;
    unsigned nu;
    unsigned nv;
    unsigned nentries;
    unsigned nbits;
    unsigned nstreams;
    unsigned ngcd;
    char in_file[512];
    bool write_to_file;
};

void print_args(const Args& a);

void parse(int argc, char **argv, Args& a);

template < class Poly_ >
bool read_from_file(const char *filename,  Poly_& f, Poly_& g) {

#define ONE_POLY 0

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

        if(s.length() == 0) {
            std::cerr << "zero length\n";
            continue;
        }

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
#if !ONE_POLY
        if(nread == 0) {
            f = p, nread = 1;
        } else if(nread == 1) {
            g = p, nread = 2;
            break;
        }
#else
        nread++, f = p;
        break; // read just one poly
#endif
    }
#if ONE_POLY
    g = f; g.diff();
    return (nread > 0);
#else
    return (nread == 2);
#endif
}

#endif // CUMP_COMPILE_DEBUG_KERNEL

bool GPU_gcd::debug_run(int argc, char **argv) {

#if CUMP_COMPILE_DEBUG_KERNEL

    typedef CORE::BigInt Integer;
    typedef CORE::BigRat Rational;
    typedef CGAL::Polynomial< Integer > Poly_1;

    srand48(time(NULL));
    gmp_randinit_mt(rands);
    gmp_randseed_ui(rands, time(NULL));

    CGAL::set_pretty_mode(std::cout);
    CGAL::set_pretty_mode(std::cerr);

    Args a;
    parse(argc, argv, a);
    print_args(a);

    unsigned deg_f = 10, deg_g = 11;
    if(a.nu == -1u || a.nv == -1u) { 
        a.nu = deg_f; a.nv = deg_g; // these are only used for random polys
    }
    if(a.nstreams != -1u)
        n_streams = a.nstreams;

    n_moduli = a.n_moduli;
    no_cpu_run = a.no_cpu_run;
    write_to_file = a.write_to_file;

    Poly_1 f, g;
#if 0
//! NOTE: these polynomials result in non-strongly regular case !!
//      const char *sgcd = "(18*x^5-27652*x^2+112)",
//         *sf="(x^11+x^10+x^9+x^8+345*x^2+4)", *sg = "(22*x^14+77*x^3+1231x^2-x+11)";

//     const char *sgcd = "(18*x^6-27652*x^4-x^3+112)",
//         *sf="(x^17-77*x^13+345*x^2-4)", *sg = "1";//222*x^4-(x^2+1)^3-023";

 /*   const char *sgcd = "344x^1000+(x^2-3x+1)^100-111-x+1",
        *sf = "2234523452345234523454*x^5000-(1*x^10-2*x-1)^400-111",
        *sg = "11234523452345*x^1000+(1*x^3-3*x+2)^500-(2*x-2)^200+123";
 */   // 4 blocks: nu + nv = 5*32 = 160

// dense case:
    const char *sgcd = "1(x+1*x^10+1)^18-23",
        *sf = "117x^101-(3311x^5-123x^2+1)^80-12",
        *sg = "888x^60 + 1123(x^4-891273x^2+1232)^60+888x-8977";
    
// sparse case:
//     const char *sgcd = "27897x^38-1(87878x-111111111111111111111111111111123123121*x^10-234x^3+1)^2-23",
//         *sf = "117x^101-(3311x^5-123x^2+1)^2-12",
//         *sg = "1122x^20 + 1123(x^4-891273x^2+1232)^2+333x-8977";

    CGAL::Polynomial_parser_d< Poly_1 > parser;
    Poly_1 gcd_fg;

    if(!parser(sf, f) || !parser(sg, g) || !parser(sgcd, gcd_fg)) {
        std::cerr << "Unable to parse poly..\n";
        return false;
    }
    f = f * gcd_fg; g *= gcd_fg;
// HACK HACK HACK
//     g = gcd_fg;

#else
    
    if(a.in_file == "" || !read_from_file(a.in_file, f, g)) {
        throw "NYI";
    }

    Poly_1 gcd_fg(1);
//     typedef Rand_coeff< Integer > RandCoeff;
//     f = generate_sparse_random_poly1< Integer, RandCoeff >
//             (a.nu, a.nentries, a.nbits);
//     g = generate_sparse_random_poly1< Integer, RandCoeff >
//             (a.nv, a.nentries, a.nbits);
//     gcd_fg = generate_sparse_random_poly1< Integer, RandCoeff >
//                 (a.ngcd, a.nentries, 1000);
// 
//     f = f * gcd_fg; g *= gcd_fg;

#endif
    if(f.degree() < g.degree()) 
         std::swap(f, g);

    f = CGAL::canonicalize(f);
    g = CGAL::canonicalize(g);

//     CUMP_out("sf: " << f << "\n")
//     CUMP_out("sg: " << g << "\n")
//     
    unsigned fi = 0, gi = 0;;
    while(f[fi] == Integer(0))
        fi++;
    while(g[gi] == Integer(0))
        gi++;
//      std::cout << "-------- " << fi << "; " << gi << "============\n";
    

//     CUMP_out("f:= " << f << ":\n")
//     CUMP_out("g:= " << g << ":\n")

    unsigned bits_;
    compute_gcd_bitlength(f, g, bits_);

    if(a.nbits != 0)
        bits_ = a.nbits;

    MPZ_vector_1 fv, gv, rv;
    CGAL::internal::construct_mpz_vector_1(f, fv);
    CGAL::internal::construct_mpz_vector_1(g, gv);

    bool success = internal_compute(fv, gv, rv, bits_);

//     if(a.in_file != "") {
//         std::stringstream ss;
//         CGAL::set_pretty_mode(ss);
//         ss << "[\"" << basename(a.in_file) <<  "\", [" << f << "," << g << "]],";
//         std::ofstream o("maple_out", std::ios::app);
//         o << ss.str();
//         o.close();
//     }

    if(success) {
#if CUMP_MRC_HOST_RECOVER    
    Poly_1 ggcd =
           CGAL::internal::construct_polynomial_from_mpz< Integer >(rv);

    ggcd = CGAL::canonicalize(ggcd);
    CUMP_out("gcd_check: " << ggcd << "\n\n")
    gcd_fg = CGAL::gcd(f,g);
    CUMP_out("gcd_truth: " << gcd_fg << "\n\n")

    if(ggcd == gcd_fg)
        CUMP_out("GCD correct!!\n")
    else {
        CUMP_out("GCD failed!!\n")
        CUMP_out("gcd_diff: " << (ggcd - gcd_fg) << "\n\n")
    }
/*
        std::stringstream ss;
        CGAL::set_pretty_mode(ss);
        ss << "FF:=[" << ggcd  << "]:";
        std::ofstream o("maple_FF", std::ios::app);
        o << ss.str();
        o.close();
    */

    CGAL::internal::dispose_mpz_vector(rv);
#endif
    } else
        CUMP_out("GPU algorithm FAILED\n");
    
#endif // CUMP_COMPILE_DEBUG_KERNEL

    return true;
}

#if CUMP_COMPILE_DEBUG_KERNEL

void print_args(const Args& a) {

    printf("\nPARAMS: n_moduli: %d; no_cpu_run: %d; nentries: %d; nbits: %d; "
        " nstreams: %d\n\n", a.n_moduli, a.no_cpu_run, a.nentries, a.nbits,
            a.nstreams);

    if(a.nu != -1u && a.nv != -1u)
        printf("nu: %d; nv: %d; ngcd: %d\n", a.nu, a.nv, a.ngcd);
}

void parse(int argc, char **argv, Args& a) {

    static option parser_opts[] = {
        {"help", 0, 0, 0},          // 0
        {"nmods", 1, 0, 0},         // 1
        {"nocpurun", 0, 0, 0},      // 2
        {"nu", 1, 0, 0}, {"nv", 1, 0, 0}, // 3, 4
        {"nentries", 1, 0, 0},      // 5
        {"nbits", 1, 0, 0},         // 6
        {"in", 1, 0, 0},            // 7
        {"usefile", 0, 0, 0},       // 8
        {"nstreams", 1, 0, 0},      // 9
        {"ngcd", 1, 0, 0},          // 10
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
            a.n_moduli = argi;
            break;
        case 2:
            a.no_cpu_run = true;
            break;
        case 3:
            a.nu = argi;
            break;
        case 4:
            a.nv = argi;
            break;
        case 5:
            a.nentries = argi;
            break; 
        case 6:
            a.nbits = argi;
            break;
        case 7:
            strncpy(a.in_file, optarg, 511);
            a.in_file[511] = 0;
            break;
        case 8:
            a.write_to_file = true;
            break;
        case 9:
            a.nstreams = argi;
            break;
        case 10:
            a.ngcd = argi;
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
