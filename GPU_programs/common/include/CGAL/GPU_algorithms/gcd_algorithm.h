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

#ifndef GCD_ALGORITHM_H
#define GCD_ALGORITHM_H

#include <gmp.h>
#include <vector>
#include <CGAL/GPU_algorithms/common_defs.h>

namespace CGAL {

namespace internal {

struct GPU_gcd {

protected:
    GPU_gcd();

public:
    static GPU_gcd& instance() {
        static GPU_gcd obj;
        return obj;
    }

    bool internal_compute(const MPZ_vector_1& fv, const MPZ_vector_1& gv,
             MPZ_vector_1& r, unsigned bits_);

    bool debug_run(int argc, char **argv);

    ~GPU_gcd();

protected:

    bool setup(unsigned deg_f_, unsigned deg_g_, unsigned bits_);

    bool run_gpu_part(const MPZ_vector_1& fv, const MPZ_vector_1& gv,
         MPZ_vector_1& r);

    bool reference_solution(const MPZ_vector_1& fv, const MPZ_vector_1& gv,
            const unsigned *const_mem, const unsigned *GCDlc,
            const unsigned *Mods, const unsigned *InvMods,
            unsigned *reference, const unsigned *GPU_res);

    bool RNS_conversion(const MPZ_vector_1& fv, const MPZ_vector_1& gv,
         unsigned *U, void *pconst, unsigned *Mods, unsigned *InvKs,
         unsigned *Mus);

    bool RNS_recover(const unsigned *R, const unsigned *Mods, 
            MPZ_vector_1& out);

    void launch_kernel(unsigned *const_mem, const unsigned *Mods,
            const unsigned *U, unsigned *R, bool& coprime);

    void mod_reduce_kernel_dispatch(const unsigned *devU,
        unsigned *devR, const unsigned *Mods);

    void QR_gcd_kernel_dispatch(unsigned *& mem_out,
        unsigned *devIn, unsigned *devOut, unsigned *devU,
            unsigned *page_locked_buf);
    void QR_gcd_kernel_data_padding();
    
    void mod_inverse_kernel1_dispatch(unsigned *devOut, 
        const unsigned *devIn, const unsigned *devMods,
        const unsigned *devGCDLcf);

    void MRC_kernel_dispatch(unsigned *devOut, const unsigned *devIn, 
        const unsigned *devMods, const unsigned *devInvDets);

    bool alloc_device_mem();
    void free_device_mem();

    void device_static_setup();
    void host_static_setup();

protected:
    unsigned nu, nv; // poly degrees, estimated bitsize of result
    unsigned max_nu, max_nv, n_moduli;  // padded poly degrees, # of moduli
    unsigned nu_ofs4, nv_ofs4;  // offsets for correct data padding
    unsigned limbs_f, limbs_g; // max # of words per coefficient 

    unsigned chunk_sz, n_blocks; // chunk size, # of blocks

    unsigned c_stride, word_size; // constant memory stride & word size
    unsigned data_padding, mods_padding; // in-data padding, moduli padding

    unsigned data_sz, prod_sz, gcd_sz; // total amount of input and output data
    unsigned block_io_sz, gcd_real_sz; // mem for interblock communication

    unsigned mem_size_in, mem_size_prod;
    unsigned mem_size_cpu_aux, mem_size_dev_aux;
    unsigned aux_cpu_mem, aux_dev_mem;
    // amount of memory allocated on host & device
    unsigned CPU_bytes_alloced, DEV_bytes_allocated;
    unsigned n_streams; // # of streams to use for concurrent kernel

    void *CPU_mem_ptr, *DEV_mem_ptr;

    char benchmark_file[64]; // benchmark file name
    bool no_cpu_run;    // whether to compute reference solution
    bool write_to_file; // whether to write benchmark results to file

};

} // namespace internal

} // namespace CGAL

#endif // GCD_ALGORITHM_H
