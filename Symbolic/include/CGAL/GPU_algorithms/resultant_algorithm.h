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

#ifndef _RESULTANT_ALGORITHM_H
#define _RESULTANT_ALGORITHM_H

#include <vector>
#include <gmp.h>
#include <CGAL/GPU_algorithms/common_defs.h>

namespace CGAL {

namespace internal {

struct GPU_resultant {

protected:

    GPU_resultant();

public:
    static GPU_resultant& instance() {
        static GPU_resultant obj;
        return obj;
    }

    bool internal_compute(const MPZ_vector_2& fv,
        const MPZ_vector_2& gv, MPZ_vector_1& r, unsigned deg_y1_,
        unsigned deg_y2_, unsigned deg_x1_, unsigned deg_x2_,
        unsigned low_deg_, unsigned high_deg_, unsigned bits_);

    bool debug_run(int argc, char **argv);

    ~GPU_resultant();
    
protected:

    bool setup(unsigned deg_y1_, unsigned deg_y2_,unsigned deg_x1_,
        unsigned deg_x2_, unsigned low_deg_, unsigned high_deg_,
        unsigned bits_);


    bool quick_run();
    bool run_gpu_part(const MPZ_vector_2& fv, const MPZ_vector_2& gv,
            MPZ_vector_1& r);

    bool alloc_device_mem();
    void free_device_mem();

    void device_static_setup();
    void host_static_setup();

    void launch_kernel(const unsigned *const_mem, const unsigned *Mods,
            const unsigned *U, unsigned *R);

    //! dispatching kernel calls
    void interpolate_data_padding();
    void interpolate_kernel_dispatch(unsigned *devR, const unsigned *devXs,
            unsigned *devInvDets);

    void resultant_data_padding();
    void resultant_kernel_dispatch(unsigned *devR, unsigned *devR2,
        unsigned *devU);

//     void mod_inverse_kernel1_dispatch(unsigned *devR, unsigned *devR2);
    void mod_inverse_kernel1_dispatch(unsigned *devR,
        unsigned *devR2, unsigned *devR3,
        unsigned *devR4);
    void mod_inverse_kernel2_dispatch(unsigned *devR, const unsigned *devU,
            const unsigned *devMods);
    void CRA_kernel_dispatch(unsigned *devOut, const unsigned *devIn,
        const unsigned *devMods, const unsigned *devInvDets);

    bool RNS_conversion(const MPZ_vector_2& fv,
        const MPZ_vector_2& gv, unsigned *U, void *pconst,
            unsigned *Mods, unsigned *InvKs, unsigned *Mus);

    bool RNS_recover(const unsigned *R, const unsigned *Mods,
             MPZ_vector_1& out);

    bool reference_solution(const unsigned *const_mem,
        const unsigned *U, const unsigned *Mods, const unsigned *InvMods,
        unsigned *reference, const unsigned *GPU_res);

protected:
    unsigned CUMP_N_MODULI, CUMP_N_BATCHES_PER_MOD;
    bool no_cpu_run;

    // nv = nr - nv
    unsigned nr, nv, nu, max_uv, // maximal of nu & nv
            deg_x1, deg_x2; // polynomial degrees
    unsigned bits;
    unsigned max_nr, nr_sz; // maximal nr value and
                    // max # of poly coefficients for all moduli
    unsigned prod_sz; // size of the result (product)
    unsigned low_deg, high_deg, n_real_pts; // estimated lower resultant degree
                    // and # of required interpolation points
    unsigned pts_overrun, // minimal # of excess points
        n_pts_failed;  // maximal # of points per modulus for which
                            // the algorithm failed

    unsigned batch_size, word_size; 
    unsigned out_data_padding, mods_padding; // interpolate data & moduli pad

    unsigned mem_size_nr, mem_size_prod;
    unsigned mem_size_cpu_aux, mem_size_dev_aux;
    unsigned aux_cpu_mem, aux_dev_mem;
    // amount of memory allocated on host & device
    unsigned CPU_bytes_alloced, DEV_bytes_allocated;

    void *CPU_mem_ptr, *DEV_mem_ptr;

    char benchmark_file[64]; // benchmark file name
    unsigned fg_idx;         // to save benchmarks with an index

    bool write_to_file; // write benchmark results to file
};

} // namespace internal

} // namespace CGAL

#endif // _RESULTANT_ALGORITHM_H
