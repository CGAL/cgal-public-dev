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
    
#ifndef _TEST_RES_CU_
#define _TEST_RES_CU_
 
#define LOG_TILE_SZ 6 // 2 ^ LOG_TILE_SZ texture dimension
#define TILE_SZ (1 << LOG_TILE_SZ)

#define N_TILES_PER_COL 16 // # of texture tiles in a single column

#if ((N_TILES_PER_COL - 1) & N_TILES_PER_COL)
#error N_TILES_PER_COL must be power-of-two!!
#endif

#define CUMP_PREFETCH_FROM_CUDA_ARRAY 0 // whether to use 2D texturing
#define CUMP_USE_PAGELOCKED_MEM  1      // whether to allocate page-locked mem

#define CUMP_COMPILE_RESULTANT_KERNEL   0
#define CUMP_COMPILE_INTERPOLATE_KERNEL 0
#define CUMP_COMPILE_TEST_KERNEL        0

#define MEASURE_MEM 1   // whether to measure time for memory transfer
#define RUN_CRT 0       // whether to run CRT reconstruction

#if (RUN_CRT & CUMP_PREFETCH_FROM_CUDA_ARRAY)
#error Texturing is not yet implemented for the CRT!!
#endif

#define CRT_4_MODULI 1 // use CRT with 4 moduli

#include <include/macros.h>
#include <include/test_res.h>
#include <misc.cu>

#define CUMP_DEVICE_MODULI_SZ 16382

#if CUMP_COMPILE_RESULTANT_KERNEL | CUMP_COMPILE_INTERPOLATE_KERNEL
#include <resultants_kernel.cu>
#endif

void sigint_handler(int sig);

Test_res *test_obj = NULL;

Test_res::Test_res(int argc, char **argv) {

    signal(SIGINT, sigint_handler);

    CPU_mem_ptr = 0, DEV_mem_ptr = 0;
#if __CUDA_VERSION__ > 100
    streams = 0;
#endif

    batch_size = 1, word_size = 4;
    aux_cpu_mem = 0, aux_dev_mem = 0;

    CUMP_N_MODULI = 1;
    CUMP_N_BATCHES_PER_MOD = 1;

    CUMP_N_MAX_STREAMS = 1; // suggested # of streams
    CUMP_USE_STREAMS = 0;
    CUMP_N_STREAMS = 1;     // real # of streams
    NO_CPU_RUN = false; // quit immediately after kernel launch

    nentries = 1111;
    write_to_file = false;
    nu = -1u, nr = -1u;

    parse(argc, argv);
    print_args();

#if CUMP_COMPILE_RESULTANT_KERNEL
    setup_resultants();
#elif CUMP_COMPILE_INTERPOLATE_KERNEL
    setup_interpolate();
#endif
}

bool Test_res::run()
{
#if CUMP_COMPILE_RESULTANT_KERNEL
    return run_resultant_kernel();
#elif CUMP_COMPILE_INTERPOLATE_KERNEL
    return run_interpolate_kernel();
#elif CUMP_COMPILE_TEST_KERNEL
    return quick_run();
#else
    return true;
#endif
}

using namespace std; // workaround streampos bug..

//#define USE_26BIT_MODULI
#ifndef USE_26BIT_MODULI

#if CUMP_COMPILE_TEST_KERNEL
#define USE_DEFAULT_MUL_MOD // using 32-bit moduli for testing purposes
#endif

#include <modular_arithm.cu>
#else
#include <include/residue_type.h>
typedef Residue zmod;
#endif // USE_26BIT_MODULI

#include <matrix_algebra.cu>

void Test_res::setup_resultants() {

    m_stride = 4; // constant memory stride: # of elements per modulus

    // obtain data padding for interpolation kernels
    data_padding_dispatch(CUMP_N_BATCHES_PER_MOD);

    batch_size = CUMP_N_MODULI * padding;

    deg_x1 = 10, deg_x2 = 10;
    unsigned min_nr, _nu, _nr;

#if 1
    _nr = 126, _nu = 63;
    min_nr = 64, max_nr = 128; // needed for data padding 
#endif
#if 0
    _nr = 95*2, _nu = 95;
    min_nr = 128, max_nr = 96*2; // needed for data padding
#endif
#if 0
    _nr = 127*2, _nu = 127;
    min_nr = 96*2, max_nr = 256; // needed for data padding
#endif

    prod_sz = 1; // each block outputs a single elem

    if(nu == -1u || nr == -1u) { // take command line arguments first
        nu = _nu, nr = _nr;
    }

    unsigned nv = nr - nu;
    if(nr < min_nr || nr > max_nr || nu > max_nr/2 - 1 || nv > nu || nv > nr) {
        printf("incorrect input: nr: %u; nu: %u\n", nr, nu);
        exit(1);
    }

    nr_sz = max_nr * (deg_x1 + deg_x2 + 2) * CUMP_N_MODULI;
    setup();
}

void Test_res::setup_interpolate() {

    m_stride = 4;
    batch_size = CUMP_N_MODULI;

    max_nr = 256; // maximal # of interpolation points
    nr = max_nr - 11; // actual # of points

    // use max_nr to ensure data is properly padded in case of odd 'nr'
    nr_sz = max_nr * batch_size;
    prod_sz = max_nr;
    setup();
}

//! \c CPU_mem_ptr collects all dynamic memory allocated either in ordinary
//! or page-locked memory triggered by the flag \c CUMP_USE_PAGELOCKED_MEM
//! \c DEV_mem_ptr collects all dynamic memory allocated on GPU, optionally
//! input operands can be allocated in texture space pointed to by
//! \c tiled_tex_array - must be bound to actual texture before use.
//! \c batch_size - # of inputs to allocate memory required for
//! parallel processing
//! \c word_size - bytes per one data element
//! \c nr_sz - the amount of input data to be transferred to GPU
//! \c prod_sz - size of product: allocated twice on CPU: 
//! \c aux_cpu_mem - auxiliary space allocated for CPU only
//! \c aux_dev_mem - auxiliary space allocated for device only
void Test_res::setup() {

    printf("building prime table..\n");
    make_oddprime_bitarray(small_prime_limit, oddprime_array);

    if(CUMP_N_MODULI > CUMP_DEVICE_MODULI_SZ / m_stride) {
        printf("insufficient constant memory..\n");
        exit(1); //HACK HACK HACK
    }

    mem_size_nr = nr_sz * word_size;
    mem_size_prod = prod_sz * batch_size * word_size;

    // don't use any sh*ty aux_batch_size, compute mem size without it
    mem_size_cpu_aux = aux_cpu_mem * word_size;
    mem_size_dev_aux = aux_dev_mem * word_size;

#if CUMP_PREFETCH_FROM_CUDA_ARRAY
    // TODO: nr_sz must be aligned by TILE_SZ (holds by definition)
    alloc_tex_mem_col(&tiled_tex_array, nr_sz);
    mem_size_nr = nr_sz * word_size;
    CUMP_SAFE_CALL(cudaMalloc(&DEV_mem_ptr, mem_size_prod + mem_size_dev_aux));
#else
    tiled_tex_array = 0;

    // allocate write-to device-mapped memory
//     CUMP_SAFE_CALL(cudaHostAlloc(&DEV_mem_write_ptr,
//                 mem_size_nr + mem_size_dev_aux, cudaHostAllocMapped));
// 
//     // allocate read-from device memory
//     CUMP_SAFE_CALL(cudaMalloc(&DEV_mem_read_ptr, mem_size_prod));

    CUMP_SAFE_CALL(cudaMalloc(&DEV_mem_ptr, mem_size_nr +
             mem_size_prod + mem_size_dev_aux));
#endif // CUMP_PREFETCH_FROM_CUDA_ARRAY

    if(DEV_mem_ptr == 0) {
        printf("ERROR: unable to allocate device mem..\n");
        exit(1);
    }

#if CUMP_USE_PAGELOCKED_MEM
#warning using page-locked mem
    cudaMallocHost(&CPU_mem_ptr, mem_size_nr + mem_size_prod*2
            + mem_size_cpu_aux);
#else
    // one more mem_size_prod block for reference solution
    CPU_mem_ptr = malloc(mem_size_nr + mem_size_prod*2 + mem_size_cpu_aux);
#endif

    if(CPU_mem_ptr == 0) {
        printf("ERROR: unable to allocate CPU mem..\n");
        exit(1);
    }

#if __CUDA_VERSION__ == 200
    cudaEventCreate(&e_start);
    cudaEventCreate(&e_end);
#endif
}

//! determine data padding depending on the # of evaluation points (nbatches)
//! for resultant algorithm
//! this must reflect the # of interpolate kernels
void Test_res::data_padding_dispatch(unsigned nbatches) {

    unsigned nthids = ((nbatches + 3) / 4 + 31) & ~31;
    if(nthids < 64)
        nthids = 64;

    if(nthids <= 64) {
        padding = 64*4;
    } else if(nthids <= 96) {
        padding = 96*4;
    } else if(nthids <= 128) {
        padding = 128*4;
    } else if(nthids <= 192) {
        padding = 192*4;
    } else if(nthids <= 256) {
        padding = 256*4;
    } else if(nthids <= 384) {
        padding = 384*4;
    } else if(nthids <= 512) {
        padding = 512*4;
    } else {
        printf("unsupported batch size: %d\n", nbatches);
        exit(1);
    }

    printf("nbacthes %d; nthids: %d; padding: %d\n",
        nbatches, nthids, padding);
}

void Test_res::interpolate_kernel_dispatch(limb *devR, limb *devU) {

#if CUMP_COMPILE_RESULTANT_KERNEL | CUMP_COMPILE_INTERPOLATE_KERNEL
    dim3 threads;
    // this is a full grid size (without streams)
    dim3 grid(CUMP_N_MODULI, 1, 1);

    unsigned nn = CUMP_N_BATCHES_PER_MOD, shm_size;
    unsigned nthids = ((nn + 3) / 4 + 31) & ~31;
    if(nthids < 64)
        nthids = 64;
    threads = dim3(nthids, 1, 1);

#define DISPATCH_CALL(__n) \
    else if(nthids <= __n) { \
        const unsigned MaxN = __n * 4; \
        shm_size = (3 * MaxN / 4) * sizeof(limb); \
        interpolate_quads_kernel< MaxN ><<<grid, threads, shm_size>>> \
                 (devR, devU, nn); \
    }

    if(0);
    DISPATCH_CALL(64)
    DISPATCH_CALL(96)
    DISPATCH_CALL(128)
    DISPATCH_CALL(192)
    DISPATCH_CALL(256)
    DISPATCH_CALL(384)
    DISPATCH_CALL(512)
    else {
        printf("DISPATCH_CALL: unsupported nr: %d\n", nr);
        exit(1);
    }
#endif // CUMP_COMPILE_RESULTANT_KERNEL | CUMP_COMPILE_INTERPOLATE_KERNEL
}


bool Test_res::run_resultant_kernel() {

#if CUMP_COMPILE_RESULTANT_KERNEL

//modifyFPUStateX86(__FPU_CW_ROUND_MASK__, __FPU_CW_ROUND_CHOP__);
    typedef CORE::BigInt Integer;
    typedef Simple_matrix< Integer > Matrix;
    typedef Simple_vector< Integer > Vector;

    unsigned nv = nr - nu;
    unsigned *pU, *U = (unsigned *)CPU_mem_ptr,
        *R = U + nr_sz, *reference = R + prod_sz * batch_size;

    static unsigned const_mem[CUMP_DEVICE_MODULI_SZ];
    unsigned *pconst = const_mem;

    pconst[0] = nu, pconst[1] = nv;
    pconst += 2;
    srand48(time(NULL));

   // outer var: y; inner var: x
    std::vector< Poly_mod_2 > poly_mods1(CUMP_N_MODULI),
        poly_mods2(CUMP_N_MODULI);

    unsigned coeff_mod = 1000;
    unsigned i, j, k, u = (1u << 24), m = u - 2, nskips = lrand48() % 4;

    printf("nskips = %d\n", nskips);

    memset(U, 0, mem_size_nr);
    for(i = 0, pU = U; i < CUMP_N_MODULI; i++, pconst += m_stride, m--) {

        get_next_mod(m, nskips);
        printf("%d: modulus: %x\n", i, m);

        float invk = 65536.0f / m, invm = 1.0f / m;

        unsigned beta = 1 << 24;
        unsigned u1, u2;
        egcd(beta, m, u1, u2);
        if((int)u2 < 0) // compute m^-1 mod beta needed for Montgomery inv
            u2 += beta;

        pconst[0] = m, pconst[1] = (unsigned&)invk;
        pconst[2] = (unsigned&)invm; pconst[3] = beta - u2; 
        zmod::set_mod(m);

        Poly_mod_2 poly1, poly2;
        //TODO: you have to keep polynomials as well
        poly1 = generate_sparse_random_poly2(nu, deg_x1, nentries, coeff_mod);
        poly2 = generate_sparse_random_poly2(nv, deg_x2, nentries, coeff_mod);

        poly_mods1[i] = poly1, poly_mods2[i] = poly2;

        for(j = deg_x1; (int)j >= 0; j--, pU += max_nr/2) {
            for(unsigned k = 0; k <= nu; k++)
                if(poly1[k].size() > j) {
                    pU[k] = poly1[k][j].x;
                }
//             for(unsigned k = 0; k <= nv; k++)
//                 if(poly2[nv - k].size() > j) {
//                     pU[k + max_nr/2] = poly2[nv - k][j].x;
//                 }
        }
        for(j = deg_x2; (int)j >= 0; j--, pU += max_nr/2) {
            for(unsigned k = 0; k <= nv; k++)
                if(poly2[nv - k].size() > j) { // reversed order !!
                    pU[k] = poly2[nv - k][j].x;
                }
        }
//         std::cout << poly1 << "\n\n";
//          std::cout << poly2 << "\n\n";
    }

    //save_testcase("mytestcase", U, nr_sz);
    //load_testcase("mytestcase", U, nr_sz);
    //! additional 2: for copying nu & nv
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(dev_const_mem, const_mem,
            (CUMP_N_MODULI * m_stride + 2) * sizeof(limb)));

    unsigned batch_sz = batch_size, n_blocks, shm_size;
    double flops; // estimated flop count per one block
    float ms;

    FILE *fout;
    if(write_to_file) {
        fout = fopen("benchmarks_out", "aw");
    }

#if MEASURE_MEM
    unsigned n_iters = 10;
    BEGIN_TIMING()
#else
    unsigned n_iters = 2;
#endif

#if CUMP_PREFETCH_FROM_CUDA_ARRAY

    printf("Prefetching from texture mem: total of %.3f Kb\n",
        (double)mem_size_nr / 1024.0);

    copy_tex_mem_col(tiled_tex_array, U, nr_sz);

    limb *devU = 0, *devR = (limb *)DEV_mem_ptr;
#else

    limb *devU = (limb *)DEV_mem_ptr, *devR = devU + nr_sz;
    cudaMemcpy(devU, U, mem_size_nr, cudaMemcpyHostToDevice);

#endif // CUMP_PREFETCH_FROM_CUDA_ARRAY

#if !MEASURE_MEM
    BEGIN_TIMING()
#endif

    dim3 threads;
    // this is a full grid size (without streams)
    dim3 grid(CUMP_N_BATCHES_PER_MOD, CUMP_N_MODULI, 1);
    n_blocks = grid.x * grid.y; 

#if 1
    threads = dim3(64, 1, 1);
    shm_size = (3 * 64 + 16) * sizeof(limb);
    resultant_block_kernel< 64 ><<<grid, threads, shm_size>>>
                 (devR, devU, deg_x1, deg_x2, padding);

#endif
#if 0
    threads = dim3(96, 1, 1);
    shm_size = (3 * 96 + 16) * sizeof(limb);
    resultant_block_kernel< 96 ><<<grid, threads, shm_size>>>
                 (devR, devU, deg_x1, deg_x2);
#endif
#if 0
    threads = dim3(128, 1, 1);
    shm_size = (3 * 128 + 16) * sizeof(limb);
    resultant_block_kernel< 128 ><<<grid, threads, shm_size>>>
                 (devR, devU, deg_x1, deg_x2);
#endif
    flops = 0;

#if !MEASURE_MEM
    END_TIMING(ms)
#endif

    interpolate_kernel_dispatch(devR, devR);

    cudaMemcpy(R, devR, mem_size_prod, cudaMemcpyDeviceToHost);

#if MEASURE_MEM
    END_TIMING(ms)
#endif

    // 3*K term represents time spent for element-wise mul, modular inverse
    // and CRT reconstruct
    double s = ms / 1000.0;

    printf("batch_size: %d\n", batch_size);
    // 2 forward + 1 backward NTT + element-wise mul
    //NOTE: flops counter is corrupt..
    double Gflop = 1e-9 * flops * batch_size;
    double GB = 1e-9 * (double)(mem_size_nr + mem_size_prod);
        
    printf("GPU time elapsed: %f ms; Gflop/s: %f; GB/s: %f\n\n", ms, Gflop/s,
            GB/s);
    printf("\n\n amount of data transferred to device: %.2f Kb\n"
            "amount of data transferred back to CPU: %.2f Kb\n"
            "# of blocks executed: %d; shared mem size: %.2f Kb\n",
            (double)(mem_size_nr + mem_size_dev_aux)/1024.0,
            (double)(mem_size_prod)/1024.0,
            n_blocks, (double)shm_size/1024.0);

    if(write_to_file) {
        fprintf(fout, "max_nr: %d    nmods: %d     n_batches: %d     nu: %d    ms: %.2f\n", max_nr, CUMP_N_MODULI, CUMP_N_BATCHES_PER_MOD, nu, ms);
        fclose(fout);
    }

//     for(i = 0; i < batch_size; i++) {
//         if(R[i] == 0)
//             printf("zero resultant !! %d\n", (i / CUMP_N_BATCHES_PER_MOD));
//     }

    if(NO_CPU_RUN)
        return true;

    Vector_mod v1m(nu + 1), v2m(nv + 1); // NOTE: nu & nv - poly degrees

    limb *pref = (limb *)reference;
    pconst = const_mem + 2; // skip first two entries
    for(i = 0, pU = U; i < CUMP_N_MODULI; i++, pconst += m_stride,
            pref += padding) {

        m = pconst[0];
        zmod::set_mod(m);

        Vector_mod xs(CUMP_N_BATCHES_PER_MOD),
            ys(CUMP_N_BATCHES_PER_MOD);

        for(j = 0; j < CUMP_N_BATCHES_PER_MOD; j++, pU += max_nr) {

            zmod pt(CUMP_N_BATCHES_PER_MOD - j);
            poly_eval2(poly_mods1[i], pt, v1m.begin());
            poly_eval2(poly_mods2[i], pt, v2m.begin());

            //print_vector(v1m);
            //print_vector(v2m);
            zmod det = resultant_modular_sylvester(v1m, v2m, (limb *)R);
            if(det.x == 0) {
                printf("modulus (%d): %x; point: %x\n", i, m, j+1);
                //! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                det.x = 1; //! replace by one for the time being !!!!
                //! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            }
            xs[j] = pt, ys[j] = det;

//             zmod det2 = resultant_modular_sylvester_new(v1m, v2m, (limb *)R);
//             printf("results: (truth, test): %x %x\n",
//                     det.x, det2.x);
//             if(det2 != det)
//                 printf(" WRONG WRONG WRONG !!\n");
//             else
//                 printf(" CORRECT !!\n");

            //pref[0] = det.x;
#if 0
        mpz_t t;
        mpz_init(t);

        ppU = pU;
        Vector v1(v1m.size()), v2(v2m.size());
        v1[nu] = ppU[0], v2[nv] = ppU[1];
        ppU += 16;
        for(k = 0; k < nu; k++, ppU++)
            v1[k] = ppU[0];

        for(k = 0; k < nv; k++, ppU++)
            v2[k] = ppU[0];
        
        Matrix bz = hybrid_bezout_matrix(v1, v2);
        Integer truth = det_berkowitz(bz);

        mpz_mod_ui(t, truth.get_mp(), m);
        zmod truth_m(t->_mp_d[0]);

        if(det != truth_m)
            std::cout << "wrong determinant: " << det << "; truth: " << 
                truth_m << "\n";
        mpz_clear(t);
#endif
            if(det == zmod(0))
                std::cout << "zero resultant!\n";
        }
        Poly_mod_1 verify(CUMP_N_BATCHES_PER_MOD);
        vandermonde_interpolate(xs, ys, pref);
    }

    /*!********************************************************************/

    limb *ptruth = (limb *)reference, *pR = (limb *)R;

//     for(i = 0, pU = U; i < CUMP_N_MODULI; i++, pR += padding,
//             ptruth += padding) {
// 
//         limb *ppR = pR, *pptruth = ptruth;
//         Vector_mod poly = poly_res[i];
//         for(j = 0; j < CUMP_N_BATCHES_PER_MOD; j++, ppR++, pptruth++) {
// 
// //             printf("%d (GPU, truth): %x %x ", j, ppR[0],
// //                     poly[j].x);            
//             if(pptruth[0] != ppR[0])
//                 printf("DIFFERS!\n");
//             else
//                 ;//printf("\n");
//         }
//     }

    return checkme(R, reference, CUMP_N_BATCHES_PER_MOD,
            padding, CUMP_N_MODULI);
#else
#warning run_resultant_kernel: dummy compilation
    return true;
#endif // CUMP_COMPILE_RESULTANT_KERNEL
}

/*!*************************************************************************
 ***************** Vandermonde interpolate *********************************
 ***************************************************************************/
bool Test_res::run_interpolate_kernel() {

#if CUMP_COMPILE_INTERPOLATE_KERNEL

    unsigned *pU, *U = (unsigned *)CPU_mem_ptr,
        *R = U + nr_sz, *reference = R + prod_sz * batch_size;

    static unsigned const_mem[CUMP_DEVICE_MODULI_SZ];
    unsigned *pconst = const_mem + 2; //! this is to ensure we have
            //! the same memory layout as in case of resultant kernel

    srand48(time(NULL));

    unsigned n_entries = 111, coeff_mod = 1000,
        total_range = 5000, hole_range = 100; // probability of holes: 0.5
    unsigned i, j, k, u = (1u << 24), m = u - 2, nskips = lrand48() % 15;

    printf("nskips = %d\n", nskips);

    unsigned n_valid = 0;
    memset(U, 0, mem_size_nr);
//    for(i = 0, pU = U; i < nr_sz; i++, pU++) {
// 
//         int x = (lrand48() % total_range) - hole_range;
//         if(x < 0) // this indicates a hole
//             x = 0;
// 
//         n_valid += (x != 0);
//         printf("%d: %d\n", i, x);
//         pU[0] = x;
//     }

    std::vector< Poly_mod_1 > poly_mods(CUMP_N_MODULI),
            poly_res(CUMP_N_MODULI);

    unsigned deg_x = nr - 1, n_ex = 1170, coeff_modx = 2000;

#if 1
    for(i = 0, pU = U; i < CUMP_N_MODULI; i++, pconst += m_stride, m--,
            pU += max_nr) {

        get_next_mod(m, nskips);
        //m = 16776989;
        printf("%d: modulus: %x\n", i, m);
        //std::cerr << "modulus: " << m << "\n";

        float invk = 65536.0f / m, invm = 1.0f / m;

        unsigned beta = 1 << 24;
        unsigned u1, u2;
        egcd(beta, m, u1, u2);
        if((int)u2 < 0) // compute m^-1 mod beta needed for Montgomery inv
            u2 += beta;

        pconst[0] = m, pconst[1] = (unsigned&)invk;
        pconst[2] = (unsigned&)invm; pconst[3] = beta - u2;
        zmod::set_mod(m);

        Poly_mod_1 pv = generate_sparse_random_poly1(deg_x, n_ex, coeff_modx);
        poly_mods[i] = pv;

        Vector_mod xs(deg_x + 1), ys(deg_x + 1);
        zmod cur(1);

        unsigned *ppU = pU + nr - 1;
        for(j = 0; j <= deg_x; j++, ppU--) {

            //std::cout << cur << " ";
            xs[deg_x - j] = cur; // enforce interpolate points in decreasing
                                // order
            ys[deg_x - j] = poly_eval1(pv, xs[deg_x - j]);
            // strange but we have to save pU in reversed order
            ppU[0] = ys[deg_x - j].x;
    
            if(ppU[0] == 0)
                printf("ZERO !!\n");
            cur = cur + zmod((rand() % 1) + 1); // number of point skips
        }

//         for(j = 0; j <= deg_x; j++, pU++) {
//             printf("%x\n", xs[j].x);
//         }
//         return true;

        if(!NO_CPU_RUN) {
            Poly_mod_1 verify(deg_x + 1);
            vandermonde_interpolate(xs, ys, verify.data());
            poly_res[i] = verify;
        }
    }

    CUDA_SAFE_CALL(cudaMemcpyToSymbol(dev_const_mem, const_mem,
             (CUMP_N_MODULI * m_stride + 2) * sizeof(limb)));
#else

    for(i = 0, pU = U; i < CUMP_N_MODULI; i++, pconst += m_stride, m--,
            pU += max_nr) {

        unsigned *ppU = pU;
        for(j = 0; j <= deg_x; j++, ppU++) {

            ppU[0] = lrand48();
   
        }

    }
#endif

    unsigned batch_sz = batch_size, n_blocks, shm_size, padding;
    double flops; // estimated flop count per one block
    float ms;

#if MEASURE_MEM
    unsigned n_iters = 20;
    BEGIN_TIMING()
#else
    unsigned n_iters = 2;
#endif

#if CUMP_PREFETCH_FROM_CUDA_ARRAY

    printf("Prefetching from texture mem: total of %.3f Kb\n",
        (double)mem_size_nr / 1024.0);

    copy_tex_mem_col(tiled_tex_array, U, nr_sz);

    limb *devU = 0, *devR = (limb *)DEV_mem_ptr;
#else

    limb *devU = (limb *)DEV_mem_ptr, *devR = devU + nr_sz;
    cudaMemcpy(devU, U, mem_size_nr, cudaMemcpyHostToDevice);

#endif // CUMP_PREFETCH_FROM_CUDA_ARRAY

#if !MEASURE_MEM
    BEGIN_TIMING()
#endif

    dim3 threads;
    // this is a full grid size (without streams)
    dim3 grid(CUMP_N_MODULI, 1, 1);
    n_blocks = grid.x * grid.y;

    unsigned nthids = ((nr + 3) / 4 + 31) & ~31;
    if(nthids < 64)
        nthids = 64;
    threads = dim3(nthids, 1, 1);

    if(nthids <= 64) {
        const unsigned MaxN = 64*4; // maximal number of elements that can be
                                    // handles with a proper stride
        shm_size = (3 * MaxN / 4) * sizeof(limb);
        interpolate_quads_kernel< MaxN ><<<grid, threads, shm_size>>>
                 (devR, devU, nr);
        padding = MaxN;
    } else if(nthids <= 96) {
        const unsigned MaxN = 96*4; // maximal number of elements that can be
                                    // handles with a proper stride
        shm_size = (3 * MaxN / 4) * sizeof(limb);
        interpolate_quads_kernel< MaxN ><<<grid, threads, shm_size>>>
                 (devR, devU, nr);
        padding = MaxN;
    } else if(nthids <= 128) {

        const unsigned MaxN = 128*4; // maximal number of elements that can be
                                    // handles with a proper stride
        shm_size = (3 * MaxN / 4) * sizeof(limb);
        interpolate_quads_kernel< MaxN ><<<grid, threads, shm_size>>>
                 (devR, devU, nr);
        padding = MaxN;
    } else if(nthids <= 192) {
        const unsigned MaxN = 192*4; // maximal number of elements that can be
                                    // handles with a proper stride
        shm_size = (3 * MaxN / 4) * sizeof(limb);
        interpolate_quads_kernel< MaxN ><<<grid, threads, shm_size>>>
                 (devR, devU, nr);
        padding = MaxN;
    } else if(nthids <= 256) {

        const unsigned MaxN = 256*4; // maximal number of elements that can be
                                    // handles with a proper stride
        shm_size = (3 * MaxN / 4) * sizeof(limb);
        interpolate_quads_kernel< MaxN ><<<grid, threads, shm_size>>>
                 (devR, devU, nr);
        padding = MaxN;
    } else if(nthids <= 512) {
        const unsigned MaxN = 256*4; // maximal number of elements that can be
                                    // handles with a proper stride
        shm_size = (3 * MaxN / 4) * sizeof(limb);
        interpolate_quads_kernel< MaxN ><<<grid, threads, shm_size>>>
                 (devR, devU, nr);
        padding = MaxN;
    }

    // align by the warp boundary (well, it happens automatically anyway)
    printf("nr = %d nthids = %d padding: %d\n", nr, nthids, padding);

    //! 128 threads, 3*128 mem: DO NOT FORGET !!!
    threads = dim3(128, 1, 1);
    shm_size = (3*128) * sizeof(limb); 
//     interpolate_quads_kernel< 9 ><<<grid, threads, shm_size>>>
//                   (devR, devU, nr);

    threads = dim3(512, 1, 1);
    shm_size = (3*512) * sizeof(limb);
//     interpolate_quads_kernel< 11 ><<<grid, threads, shm_size>>>
//                   (devR, devU, nr);

//     interpolate_doubles_kernel< 10 ><<<grid, threads, shm_size>>>
//                   (devR, devU, nr);

    threads = dim3(128, 1, 1);
    shm_size = (3*128) * sizeof(limb);
//     test_red_kernel<<<grid, threads, shm_size>>> (devR, devU);
    
    flops = 0;

#if !MEASURE_MEM
    END_TIMING(ms)
#endif

    cudaMemcpy(R, devR, mem_size_prod, cudaMemcpyDeviceToHost);

#if MEASURE_MEM
    END_TIMING(ms)
#endif

    // 3*K term represents time spent for element-wise mul, modular inverse
    // and CRT reconstruct
    double s = ms / 1000.0;

    printf("batch_size: %d\n", batch_size);
    // 2 forward + 1 backward NTT + element-wise mul
    //NOTE: flops counter is corrupt..
    double Gflop = 1e-9 * flops * batch_size;
    double GB = 1e-9 * (double)(mem_size_nr + mem_size_prod);
        
    printf("GPU time elapsed: %f ms; Gflop/s: %f; GB/s: %f\n\n", ms, Gflop/s,
            GB/s);
    printf("\n\n amount of data transferred to device: %.2f Kb\n"
            "amount of data transferred back to CPU: %.2f Kb\n"
            "# of blocks executed: %d; shared mem size: %.2f Kb\n",
            (double)(mem_size_nr + mem_size_dev_aux)/1024.0,
            (double)(mem_size_prod)/1024.0,
            n_blocks, (double)shm_size/1024.0);

    if(NO_CPU_RUN)
        return true;

    limb *pref = (limb *)reference;

    pconst = const_mem + 2; 
    for(i = 0, pU = U; i < CUMP_N_MODULI; i++, pconst += m_stride,
            pref += padding) { // account for padding as well

        Poly_mod_1 pv = poly_res[i];
        for(j = 0; j < nr; j++) {
            pref[j] = pv[j].x;
        }
    }

//     for(i = 0, pU = U; i < CUMP_N_MODULI; i++,
//             pref += max_nr, pU += max_nr) { // account for padding as well
// 
//         limb sum = 0;
//         for(j = 0; j < nr; j++) {
//             sum += pU[j];
//             pref[j] = sum;
//         }
//     }

    limb *ptruth = (limb *)reference, *pR = (limb *)R;
    return //false;
        checkme(R, reference, nr, max_nr, batch_size);
#else
#warning run_interpolate_kernel: dummy compilation
    return true;
#endif // CUMP_COMPILE_INTERPOLATE_KERNEL
}

#ifdef USE_26BIT_MODULI
#include <include/ieee754_fp.h>
#endif

//NOTE NOTE TODO: common.mk: change to NVCCFLAGS +=
// otherwise optimizations does not work:
bool Test_res::quick_run() {

#if CUMP_COMPILE_TEST_KERNEL

//NOTE: can use _controlfp() to force to_nearest fp rounding mode

#ifdef USE_26BIT_MODULI
    modifyFPUStateX86(__FPU_CW_ROUND_MASK__, __FPU_CW_ROUND_CHOP__);
#endif

     srand48(time(NULL));

// TODO: generate some matrix with large entries
// evaluate its resultant using exact approach (gmp)
// and modular approach (CRT) as well as try representing numbers
// as polynomials and apply NTT for point evaluation

    typedef CORE::BigInt Integer;
    //typedef int Integer;
    typedef Simple_matrix< Integer > Matrix;
    typedef Simple_vector< Integer > Vector;

//     unsigned bitlen_v1 = 32, bitlen_v2 = 28;
    bool res = true;

    unsigned timerCPU;
    cutCreateTimer(&timerCPU);

    gmp_randstate_t rands;
    gmp_randinit_mt(rands);

#ifndef USE_26BIT_MODULI
    unsigned n_skips = 13, ii = 0;
    umod_t m((1<<31));
    while(1) {
        if(is_prime(m) && (m&3)==3 && n_skips-- == 0) {
            printf("%d modulus: %u %d\n", ++ii, m, (m&3));
            break;
        }
        m--;    
    }
    zmod::set_mod(m);
#endif

#if 0
    // degree           0   1   2   3   4  5    
//     int sample_v1[] = {-11, -2, 6, 8, 656, 12},
//         sample_v2[] = {-1, 5, 3, 4, -2, 0};
// 
//     Vector_mod v1mod(deg_v1 + 1), v2mod(deg_v2 + 1);
//     for(unsigned i = 0; i <= deg_v1; i++) {
//         //v1[i] = Integer(sample_v1[i]);
//         v1mod[i] = zmod(sample_v1[i]);
//         //mpz_urandomb(v1[i].get_mp(), rands, bitlen_v1);
//     }
//  
//     for(unsigned i = 0; i <= deg_v2; i++) {
//         //v2[i] = Integer(sample_v2[i]);
//         v2mod[i] = zmod(sample_v2[i]);
// //         mpz_urandomb(v2[i].get_mp(), rands, bitlen_v2);
//     }

    unsigned deg_x = 7, n_ex = 7, coeff_modx = 2000;
    Poly_mod_1 pv1 = generate_sparse_random_poly1(deg_x, n_ex, coeff_modx);

    Vector_mod xs(deg_x + 1), ys(deg_x + 1);
    zmod cur(1);
    for(unsigned i = 0; i <= deg_x; i++) {

        //std::cout << cur << " ";
        xs[i] = cur;
        
        ys[i] = poly_eval1(pv1, xs[i]);
        cur = cur + zmod((rand() % 40) + 1); // number of point skips
    }

    //std::cout << "\npoly = " << pv1 << "\n";

    //Matrix bezout = hybrid_bezout_matrix(v1, v2);
    //print_vector(v1mod); print_vector(v2mod);
    //print_matrix(bezout);

//     Integer det = det_berkowitz(bezout);
//     printf("++++ determinant: %s\n", det.get_str().c_str());
//     printf("hadamards bound: %.f\n", hadamards_bound(bezout));

    Vector_mod verify(deg_x + 1);

    cutStartTimer(timerCPU);
    vandermonde_interpolate(xs, ys, verify);

    cutStopTimer(timerCPU);

    std::cout << "Elapsed time: " << 
        (cutGetTimerValue(timerCPU)/1000.0f) << "\n";

    if(verify == pv1)
        std::cout << "CORRECT\n";
    else
        std::cout << "FAILED\n";    

    std::cout << verify[0] << " and " << pv1[0] << "\n";
    cutDeleteTimer(timerCPU); 
    return true;
#endif

    mpz_t t;
    mpz_init(t);

    //  TODO: store bivariate polynomials as an array of arrays
    // array in y with coefficients as polynomials in x
    // substitute x with some value [0,m-1] to obtain a univariate polynomial
    // in y. Compute their resultant and count how many failures are there
    // for random inputs..

    unsigned deg_x1 = 10, deg_x2 = 12, deg_y1 = 8, deg_y2 = 14;
    
    // outer var: y; inner var: x
    Poly_mod_2 poly1, poly2; 

    unsigned n_entries = 11, coeff_mod = 1000;
    poly1 = generate_sparse_random_poly2(deg_y1, deg_x1,
         n_entries, coeff_mod);
    poly2 = generate_sparse_random_poly2(deg_y2, deg_x2,
         n_entries, coeff_mod);   

     std::cout << "======= poly1: " << poly1 << "\n\n";
     std::cout << "======= poly2: " << poly2 << "\n\n";


    unsigned low_deg, high_deg;
    double bits;
//     compute_resultant_bounds< double >(poly1, poly2, low_deg, high_deg, bits);

    cutStartTimer(timerCPU);

#if 0
    Vector_mod v1m(deg_y1 + 1), v2m(deg_y2 + 1);
    for(unsigned evalp = 1; evalp <= 1; evalp++) {

         std::cout << "\n" << evalp << ": ";
        zmod x(evalp);
        poly_eval2(poly1, x, v1m.begin());
        poly_eval2(poly2, x, v2m.begin());

//         print_vector(v1m);
//         print_vector(v2m);

//         for(unsigned i = 0; i <= deg_v2; i++) {
//             mpz_mod_ui(t, v2[i].get_mp(), m);
//             v2m[i] = zmod(t->_mp_d[0]);
//         }
//         mpz_clear(t);
        
        zmod det = resultant_modular_sylvester(v1m, v2m);

        if(det == zmod(0))
            std::cout << "zero resultant!\n";
#if 1
        Vector v1(deg_y1 + 1), v2(deg_y2 + 1);
        for(unsigned i = 0; i <= deg_y1; i++)
            v1[i] = v1m[i].x;

        for(unsigned i = 0; i <= deg_y2; i++)
            v2[i] = v2m[i].x;

        Matrix bz = hybrid_bezout_matrix(v1, v2);
        Integer truth = det_berkowitz(bz);

//         std::cout << "truth: " << truth << "\n";
        
        mpz_mod_ui(t, truth.get_mp(), m);
        zmod truth_m(t->_mp_d[0]);

        if(det != truth_m)
            std::cout << "wrong determinant: " << det << "; truth: " << 
                truth_m << "\n";
#endif
    }
    cutStopTimer(timerCPU);

    std::cout << "Elapsed time: " <<
        (cutGetTimerValue(timerCPU)/1000.0f) << "\n";

    cutDeleteTimer(timerCPU);
#endif
    return res;

#else
#warning quick_run: dummy compilation
    return true;
#endif // CUMP_COMPILE_TEST_KERNEL
}

Test_res::~Test_res() {

    if(CPU_mem_ptr != 0) {
#if CUMP_USE_PAGELOCKED_MEM
        printf("freeing page-locked mem..\n");
        cudaFreeHost(CPU_mem_ptr);
#else
        free(CPU_mem_ptr);
#endif
    }

    if(DEV_mem_ptr != 0) {
        printf("freeing device mem..\n");
        CUDA_SAFE_CALL(cudaFree(DEV_mem_ptr));
    }

#if CUMP_PREFETCH_FROM_CUDA_ARRAY
    CUDA_SAFE_CALL(cudaFreeArray(tiled_tex_array));
#endif

#if __CUDA_VERSION__ == 200
    if(streams != NULL) {
        for(unsigned i = 0; i < CUMP_N_STREAMS; i++) {
            cudaStreamDestroy(streams[i]);
        }
        delete []streams;
    }
    cudaEventDestroy(e_start);
    cudaEventDestroy(e_end);
#endif
}

void Test_res::test_case(limb *A, limb *B, unsigned n, unsigned test_index,
         bool show_dump) {

    //printf("n = %d; test_index = %d\n", n, test_index);
}

void Test_res::print_args() {

    printf("CUMP_N_MODULI: %d; CUMP_N_BATCHES_PER_MOD: %d;\n"
      "CUMP_N_MAX_STREAMS: %d; NO_CPU_RUN: %d MEASURE_MEM: %d; RUN_CRT: %d; "
        "nentries: %d\n",
       CUMP_N_MODULI, CUMP_N_BATCHES_PER_MOD, CUMP_N_MAX_STREAMS, NO_CPU_RUN,
                MEASURE_MEM, RUN_CRT, nentries);
    if(nu != -1u && nr != -1u)
        printf("nr: %d; nu: %d\n", nr, nu);
}

#include <getopt.h>

void Test_res::parse(int argc, char **argv) {

    static option parser_opts[] = {
        {"help", 0, 0, 0},          // 0
        {"nmods", 1, 0, 0},         // 1
        {"nbatches", 1, 0, 0},      // 2
        {"nstreams", 1, 0, 0},      // 3
        {"usestreams", 1, 0, 0},    // 4
        {"nocpurun", 0, 0, 0},      // 5
        {"nu", 1, 0, 0}, {"nr", 1, 0, 0}, // 6, 7
        {"usefile", 0, 0, 0},       // 8
        {"nentries", 1, 0, 0},      // 9
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
        if(parser_opts[idx].has_arg)
            argi = atoi(optarg);
        switch(idx) {
        case 1:
            CUMP_N_MODULI = argi;
            break;
        case 2:
            CUMP_N_BATCHES_PER_MOD = argi;
            break;
        case 3:
            CUMP_N_MAX_STREAMS = argi;
            break;
        case 4:
            CUMP_USE_STREAMS = argi;
            break;
        case 5:
            NO_CPU_RUN = true;
            break;
        case 6:
            nu = argi;
            break;
        case 7:
            nr = argi;
            break;
        case 8:
            write_to_file = true;
            break;
        case 9:
            nentries = argi;
            break;
        default:
            printf("unrecognized option..\n");
            exit(1);
        }
    }
}

#if 0
    printf("modulus: %x\n", m);

    unsigned n_min = 10000, n_max = 0, xmin, xmax;

    unsigned rolex1, rolex2;
    cutCreateTimer(&rolex1);
    cutCreateTimer(&rolex2);

    for(i = 0; i < 50000; i++) {

        unsigned x = lrand48() % m;
        zmod zx(x);

    cutStartTimer(rolex1);
        unsigned r1 = mod_inverse_montgomery2(x, m, n1);
    cutStopTimer(rolex1);

    cutStartTimer(rolex2);
        unsigned trust = mod_inverse(x, m);
            //mod_inverse_montgomery(x, m, n2);
    cutStopTimer(rolex2);

        int n_iters = n1;
        if(n_iters < n_min) {
            n_min = n_iters, xmin = x; }
        if(n_iters > n_max) {
            n_max = n_iters, xmax = x; }

        if(n1 != n2)
            printf("wrong inverse: %u - %u %u - %u\n", n1, n2,
                r1, trust);
        //printf("n_iters = %d\n", n_iters);
    }
    printf("n_min: %d - %x; n_max: %d - %x\n", n_min, xmin, n_max, xmax);

    std::cout << "Elapsed time: " << 
        (cutGetTimerValue(rolex1)/1000.0f) << " and " <<
             (cutGetTimerValue(rolex2)/1000.0f) << "\n";

    return true;

#endif

void sigint_handler(int sig)
{
    printf("SIGINT received, bailing out..\n");

    if(test_obj != NULL) {
        delete test_obj;
        test_obj = NULL;
    }
    signal(SIGINT, SIG_DFL);
    raise(SIGINT);
}

#endif // _TEST_MUL_CU_
