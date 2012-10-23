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

#ifndef GPU_DEVICE_MANAGER_H
#define GPU_DEVICE_MANAGER_H

#include <vector>
#include <include/MOD_entry.h>

namespace CGAL {

namespace internal {

typedef void (*GPU_DEVICE_CALLBACK)(int);

struct GPU_device_manager {

protected:
    GPU_device_manager();

    void device_static_setup();

public:
    static GPU_device_manager& instance() {
        static GPU_device_manager obj;
        return obj;
    }

    void register_callback(GPU_DEVICE_CALLBACK func) {
        callbacks.push_back(func);
    }

    void push_CUDA_context();
    void pop_CUDA_context();

    void notify_callbacks();
    ~GPU_device_manager();

    std::vector< MOD_entry > mod_table;
    std::vector< double > mod_bits;

protected:

    std::vector< GPU_DEVICE_CALLBACK > callbacks;
};

} // namespace internal

} // namespace CGAL

#endif // GPU_DEVICE_MANAGER_H
