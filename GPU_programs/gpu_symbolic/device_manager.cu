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

// #include <include/cudart_version.h>

#if CUDART_VERSION < 2000
#undef __EXCEPTIONS // disable exceptions for CUDA 1.0 (bug workaround)
#endif

#include <iostream>
#include <signal.h>
#include <gmp.h>

#include <include/macros.h>
#include <include/device_manager.h>
#include <include/modular_arithm.h>

#include <cuda.h>

namespace CGAL {

namespace internal {

#if (CUDART_VERSION >= 2000)
cudaEvent_t e_start, e_end;
#endif

void sigint_handler(int sig) {
    printf("SIGINT received, bailing out..\n");

    signal(SIGINT, SIG_DFL);
    signal(SIGABRT, SIG_DFL);
//     signal(SIGKILL, SIG_DFL);
    signal(SIGSEGV, SIG_DFL);

    GPU_device_manager& obj = GPU_device_manager::instance();
    obj.notify_callbacks();

    raise(SIGINT);
}

GPU_device_manager::GPU_device_manager() {

    device_static_setup();

    signal(SIGABRT, sigint_handler);
//     signal(SIGKILL, sigint_handler);
    signal(SIGINT, sigint_handler);
    signal(SIGSEGV, sigint_handler);
}

void GPU_device_manager::notify_callbacks() {

    for(unsigned i = 0; i < callbacks.size(); i++) {
        (*callbacks[i])(0);
    }
}

GPU_device_manager::~GPU_device_manager() {
#if (CUDART_VERSION >= 2000)
    cudaEventDestroy(e_start);
    cudaEventDestroy(e_end);
#endif
}

//! protected by mutex ??
CUcontext global_ctx;

void GPU_device_manager::push_CUDA_context() {
//     printf("pushing cuda context: %d\n", global_ctx);
    CUMP_SAFE_CALL_drv(cuCtxPushCurrent(global_ctx));
}

void GPU_device_manager::pop_CUDA_context() {

//     unsigned ver;
//     cuCtxGetApiVersion(global_ctx, &ver);
    CUMP_SAFE_CALL_drv(cuCtxPopCurrent(&global_ctx));
//     printf("popping current context: %d\n", global_ctx);
}

void GPU_device_manager::device_static_setup() {

    int cnt;
    CUMP_SAFE_CALL(cudaGetDeviceCount(&cnt));
    if(cnt == 0) {                                                  
        fprintf(stderr, "\nSorry, your graphics card does not support CUDA, "
            "do you want to play Solitaire instead ?\n");
        exit(1);
    }                                                           

    int dev = 0;                                                
    cudaDeviceProp props;                                  
    CUMP_SAFE_CALL(cudaGetDeviceProperties(&props, dev));
    if(props.major < 1) {
        fprintf(stderr, "Device does not support CUDA.\n");
        exit(1);                                               
    }
    printf("\nFound CUDA compatible device: %s; SM version: %d.%d"
           "\nTotal device mem: %d Mb\n", props.name, props.major,
            props.minor, props.totalGlobalMem);

    CUMP_SAFE_CALL(cudaSetDevice(dev));

#if CUDART_VERSION >= 2000
    cudaEventCreate(&e_start);
    cudaEventCreate(&e_end);
#endif

// #if CUDART_VERSION >= 2020
//    CUMP_SAFE_CALL(cudaSetDeviceFlags(cudaDeviceMapHost));
// #endif

#if CUMP_USE_32BIT_MODULI_SET
    printf("Using 31-bit moduli set\n");    
    const char *fname = "moduli_set31";
#else
    printf("Using 24-bit moduli set\n");
    const char *fname = "moduli_set24";
#endif

    unsigned sz = load_moduli_set(mod_table, mod_bits, fname);
    if(sz == -1u) {
        printf("Unable to load moduli set\n");
        abort();
    }

    pop_CUDA_context(); // detach context from the calling thread
                        // to enable other threads use it
}

} // namespace internal

} // namespace CGAL
