// SPDX-FileCopyrightText: Copyright (c) 2026
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "cuBQL/builder/cuda/builder_common.h"
#include "cubql_allocators.h"
#include <cstdio>

namespace cuBQL {
namespace utils {
namespace impl {

using cuBQL::gpuBuilder_impl::AtomicBox;
using cuBQL::gpuBuilder_impl::atomic_grow;

// 1. Initializer kernel using the native .set_empty() API
template<typename T, int D>
__global__ void initGlobalBoxKernel(AtomicBox<box_t<T,D>> *globalBox)
{
    globalBox->set_empty();
}

// 2. Elemental kernel: Each thread processes one box and grows the global box atomically
template<typename T, int D>
__global__ void computeGlobalBoxKernel(AtomicBox<box_t<T,D>> *globalBox,
                                       const box_t<T,D>      *primBoxes,
                                       uint32_t                numPrims)
{
    const int primID = threadIdx.x + blockIdx.x * blockDim.x;
    if (primID >= numPrims) return;

    const box_t<T,D> box = primBoxes[primID];
    
    if (box.get_lower(0) <= box.get_upper(0)) {
        atomic_grow(*globalBox, box);
    }
}

// 3. Resolution kernel: Decodes the internal atomic/integer representation back to true floats
template<typename T, int D>
__global__ void resolveGlobalBoxKernel(box_t<T,D> *resolvedBox, const AtomicBox<box_t<T,D>> *globalBox)
{
    for (int d = 0; d < D; ++d) {
        resolvedBox->lower[d] = globalBox->get_lower(d);
        resolvedBox->upper[d] = globalBox->get_upper(d);
    }
}

} // namespace impl

/**
 * @brief High-performance GPU-parallel global geometric object bounding box calculation.
 * Safely decodes internal atomic representations before host synchronization.
 */
template <typename T, int D>
inline void computeGlobalBoxParallel(cuBQL::box_t<T, D>& out_globalBox,
                                     const cuBQL::box_t<T, D>* d_boxesA, 
                                     int numPrims, 
                                     cudaStream_t s, 
                                     cuBQL::DeviceMemoryResource& memResource) {
    
    printf("[TRACKING-HOST] -> Entered computeGlobalBoxParallel (Object Strategy). Prims: %d\n", numPrims);
    fflush(stdout);

    if (numPrims <= 0) {
        printf("[TRACKING-HOST] -> Early exit: numPrims <= 0\n");
        fflush(stdout);
        return;
    }

    // Allocate workspace for the atomic collector wrapper
    gpuBuilder_impl::AtomicBox<box_t<T, D>>* d_globalBox = nullptr;
    _ALLOC(d_globalBox, 1, s, memResource);

    // Allocate workspace for a standard layout box to handle safe host resolution
    box_t<T, D>* d_resolvedBox = nullptr;
    _ALLOC(d_resolvedBox, 1, s, memResource);

    if (!d_globalBox || !d_resolvedBox) {
        printf("[CRITICAL-HOST] -> Workspace allocation failed!\n");
        fflush(stdout);
        return;
    }

    // Initialize to empty on the device
    impl::initGlobalBoxKernel<T, D><<<1, 1, 0, s>>>(d_globalBox);

    // Compute bounds globally across all primitives
    int blockSize = 1024;
    int numBlocks = (numPrims + blockSize - 1) / blockSize;
    impl::computeGlobalBoxKernel<T, D><<<numBlocks, blockSize, 0, s>>>(d_globalBox, d_boxesA, numPrims);

    // Decode the atomic box into a standard box structure on the device
    printf("[TRACKING-HOST] -> Resolving internal atomic representation to standard float primitives...\n");
    fflush(stdout);
    impl::resolveGlobalBoxKernel<T, D><<<1, 1, 0, s>>>(d_resolvedBox, d_globalBox);

    // Enqueue standard box memory back to the host destination variable
    printf("[TRACKING-HOST] -> Enqueuing DeviceToHost MemcpyAsync...\n");
    fflush(stdout);
    CUBQL_CUDA_CALL(MemcpyAsync(&out_globalBox, d_resolvedBox, sizeof(box_t<T, D>), cudaMemcpyDeviceToHost, s));

    // Wait cleanly for stream execution to finalize
    printf("[TRACKING-HOST] -> Awaiting Stream Synchronize...\n");
    fflush(stdout);
    cudaError_t errSync = cudaStreamSynchronize(s);
    
    if (errSync != cudaSuccess) {
        printf("[CRITICAL-GPU RUNTIME ERROR] -> Stream synchronization failed: %s\n", cudaGetErrorString(errSync));
    }
    fflush(stdout);

    // Free device allocations
    _FREE(d_globalBox, s, memResource);
    _FREE(d_resolvedBox, s, memResource);
}

} // namespace utils
} // namespace cuBQL