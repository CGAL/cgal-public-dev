// SPDX-FileCopyrightText: Copyright (c) 2026
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "cuBQL/builder/cuda/builder_common.h"
#include <cooperative_groups.h>

namespace cuBQL {
namespace utils {

namespace impl {

// Float/Double safe atomic Min loop using CAS
template <typename T>
__device__ inline void atomicMinGeneric(T* address, T val) {
    if constexpr (sizeof(T) == 4) {
        int* address_as_int = (int*)address;
        int old = *address_as_int, assumed;
        while (val < __int_as_float(old)) {
            assumed = old;
            old = atomicCAS(address_as_int, assumed, __float_as_int(val));
            if (assumed == old) break;
        }
    } else if constexpr (sizeof(T) == 8) {
        unsigned long long extern* address_as_ull = (unsigned long long extern*)address;
        unsigned long long old = *address_as_ull, assumed;
        while (val < __longlong_as_double(old)) {
            assumed = old;
            old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val));
            if (assumed == old) break;
        }
    }
}

// Float/Double safe atomic Max loop using CAS
template <typename T>
__device__ inline void atomicMaxGeneric(T* address, T val) {
    if constexpr (sizeof(T) == 4) {
        int* address_as_int = (int*)address;
        int old = *address_as_int, assumed;
        while (val > __int_as_float(old)) {
            assumed = old;
            old = atomicCAS(address_as_int, assumed, __float_as_int(val));
            if (assumed == old) break;
        }
    } else if constexpr (sizeof(T) == 8) {
        unsigned long long extern* address_as_ull = (unsigned long long extern*)address;
        unsigned long long old = *address_as_ull, assumed;
        while (val > __longlong_as_double(old)) {
            assumed = old;
            old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val));
            if (assumed == old) break;
        }
    }
}

// Parallel reduction kernel using dynamic shared memory allocation
template <typename T, int D>
__global__ void globalBoxReductionKernel(const cuBQL::box_t<T, D>* d_boxes, 
                                         uint32_t numPrims, 
                                         cuBQL::box_t<T, D>* d_result) {
    // 1. Initialize thread-local bounding box to inversion boundaries
    cuBQL::box_t<T, D> localBox;
    for (int d = 0; d < D; ++d) {
        localBox.lower[d] = cuBQL::gpuBuilder_impl::empty_box_lower<T>();
        localBox.upper[d] = cuBQL::gpuBuilder_impl::empty_box_upper<T>();
    }

    uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    uint32_t stride = blockDim.x * gridDim.x;

    // 2. Grid-stride loop processing for maximum coalesced reads
    for (uint32_t i = tid; i < numPrims; i += stride) {
        cuBQL::box_t<T, D> b = d_boxes[i];
        if (b.get_lower(0) <= b.get_upper(0)) { // Skip invalid/empty source boxes
            for (int d = 0; d < D; ++d) {
                localBox.lower[d] = min(localBox.lower[d], b.lower[d]);
                localBox.upper[d] = max(localBox.upper[d], b.upper[d]);
            }
        }
    }

    // 3. Layout dynamic shared memory safely for dimensions
    extern __shared__ char sharedMemLayout[];
    T* s_min = (T*)sharedMemLayout;
    T* s_max = (T*)(sharedMemLayout + (D * blockDim.x * sizeof(T)));

    for (int d = 0; d < D; ++d) {
        s_min[d * blockDim.x + threadIdx.x] = localBox.lower[d];
        s_max[d * blockDim.x + threadIdx.x] = localBox.upper[d];
    }
    __syncthreads();

    // 4. Parallel block-level reduction loop
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) {
            for (int d = 0; d < D; ++d) {
                uint32_t baseIdx = d * blockDim.x + threadIdx.x;
                s_min[baseIdx] = min(s_min[baseIdx], s_min[baseIdx + s]);
                s_max[baseIdx] = max(s_max[baseIdx], s_max[baseIdx + s]);
            }
        }
        __syncthreads();
    }

    // 5. Block leader flushes the block results to global memory atomically
    if (threadIdx.x == 0) {
        for (int d = 0; d < D; ++d) {
            atomicMinGeneric<T>(&d_result->lower[d], s_min[d * blockDim.x]);
            atomicMaxGeneric<T>(&d_result->upper[d], s_max[d * blockDim.x]);
        }
    }
}

} // namespace impl

/**
 * @brief High-performance GPU-parallel global bounding box calculation.
 * @return cuBQL::box_t<T, D> The final tightly computed bounding box.
 */
template <typename T, int D>
inline cuBQL::box_t<T, D> computeGlobalBoxParallel(const cuBQL::box_t<T, D>* d_boxesA, 
                                                   int numPrims, 
                                                   cudaStream_t s, 
                                                   cuBQL::DeviceMemoryResource& memResource) {
    // Standard initialization structure
    cuBQL::box_t<T, D> h_globalBox;
    for (int d = 0; d < D; ++d) {
        h_globalBox.lower[d] = cuBQL::gpuBuilder_impl::empty_box_lower<T>();
        h_globalBox.upper[d] = cuBQL::gpuBuilder_impl::empty_box_upper<T>();
    }

    if (numPrims <= 0) return h_globalBox;

    // Allocate temporary device structure workspace
    cuBQL::box_t<T, D>* d_globalBox = nullptr;
    _ALLOC(d_globalBox, 1, s, memResource);
    CUBQL_CUDA_CALL(MemcpyAsync(d_globalBox, &h_globalBox, sizeof(cuBQL::box_t<T, D>), cudaMemcpyHostToDevice, s));

    // Calculate grid performance metrics
    int blockSize = 256; 
    int numBlocks = std::min(256, (int)((numPrims + blockSize - 1) / blockSize));
    size_t sharedMemSize = 2 * D * blockSize * sizeof(T);

    // Launch reduction pass
    impl::globalBoxReductionKernel<T, D><<<numBlocks, blockSize, sharedMemSize, s>>>(d_boxesA, numPrims, d_globalBox);

    // Synchronize stream and pull layout structures cleanly back to host
    cuBQL::box_t<T, D> globalBoxA;
    CUBQL_CUDA_CALL(MemcpyAsync(&globalBoxA, d_globalBox, sizeof(cuBQL::box_t<T, D>), cudaMemcpyDeviceToHost, s));
    CUBQL_CUDA_CALL(StreamSynchronize(s));

    // Free device allocations
    _FREE(d_globalBox, s, memResource);

    return globalBoxA;
}

} // namespace utils
} // namespace cuBQL