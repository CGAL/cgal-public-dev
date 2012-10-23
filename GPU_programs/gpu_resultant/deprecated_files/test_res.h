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

#ifndef _TEST_RES_H_
#define _TEST_RES_H_


struct Test_res {

    Test_res(int argc, char **argv);

    ~Test_res();
    
    void test_case(limb *A, limb *B, unsigned n, unsigned test_index,
        bool show_dump = false);

    void setup();
    void setup_resultants();
    void setup_interpolate();

    bool run();
    bool run_resultant_kernel();
    bool run_interpolate_kernel();

    bool quick_run();

    void parse(int argc, char **argv);
    void print_args();

    void data_padding_dispatch(unsigned nbatches);
    void interpolate_kernel_dispatch(limb *devR, limb *devU);

protected:
    unsigned CUMP_N_MODULI, CUMP_N_BATCHES_PER_MOD, CUMP_USE_STREAMS,
        CUMP_N_MAX_STREAMS, CUMP_N_STREAMS, CUMP_N_BATCHES_PER_STREAM,
         CUMP_N_CRT_BATCHES, NO_CPU_RUN;

    // nv = nr - nv
    unsigned nentries; // maximal # of non-zero entries in polynomials
    unsigned nr, nu, deg_x1, deg_x2; // polynomial degrees
    unsigned max_nr, nr_sz; // maximal nr value and
                    // max # of poly coefficients for all moduli
    unsigned prod_sz; // size of the result (product)

    unsigned batch_size, m_stride, word_size; // OBSOLETE ?
    unsigned padding;

    unsigned mem_size_nr, mem_size_prod;
    unsigned mem_size_cpu_aux, mem_size_dev_aux;
    unsigned aux_cpu_mem, aux_dev_mem;

    void *CPU_mem_ptr, *DEV_mem_ptr;

    bool write_to_file; // whether to write benchmark results to file

#if __CUDA_VERSION__ > 100
    cudaStream_t *streams;
#endif
    cudaArray *tiled_tex_array;
};



#endif // _TEST_RES_H_
