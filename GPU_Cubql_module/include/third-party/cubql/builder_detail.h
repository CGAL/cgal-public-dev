// SPDX-FileCopyrightText: Copyright (c) 2026
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "cuBQL/builder/cuda/builder_common.h"
#include "cuBQL/builder/cuda/refit.h"
//#include <cub/cub.cuh>

namespace cuBQL {
namespace builder_detail {

using namespace cuBQL::gpuBuilder_impl;

// Consolidated atomic float comparison helpers using CAS loops
__device__ inline float atomicMinFloat(float* address, float val) {
    int* address_as_int = (int*)address;
    int old = *address_as_int, assumed;
    while (val < __int_as_float(old)) {
        assumed = old;
        old = atomicCAS(address_as_int, assumed, __float_as_int(val));
        if (assumed == old) break;
    }
    return __int_as_float(old);
}

__device__ inline float atomicMaxFloat(float* address, float val) {
    int* address_as_int = (int*)address;
    int old = *address_as_int, assumed;
    while (val > __int_as_float(old)) {
        assumed = old;
        old = atomicCAS(address_as_int, assumed, __float_as_int(val));
        if (assumed == old) break;
    }
    return __int_as_float(old);
}

// Unified PrimState tracking layout
struct PrimState {
    union {
        struct {
            uint64_t primID : 31;
            uint64_t done   : 1;
            uint64_t nodeID : 32;
        };
        uint64_t bits;
    };
};

// Unified TempNode representation across tree level and forest builders
template <typename T, int D>
struct CUBQL_ALIGN(16) TempNode {
    using box_t = cuBQL::box_t<T, D>;
    union {
        struct {
            AtomicBox<box_t> centBounds;
            uint32_t         count;
            uint32_t         level; // Absolute depth level of this node
        } openBranch;
        struct {
            uint32_t offset;
            int      dim;
            uint32_t tieBreaker;
            float    pos;
        } openNode;
        struct {
            uint32_t offset;
            uint32_t count;
            uint32_t unused[2];
        } doneNode;
    };
};

__global__ inline void clearPrimStates(PrimState* primStates, int totalAllocationPrims) {
    const int primID = threadIdx.x + blockIdx.x * blockDim.x;
    if (primID >= totalAllocationPrims) return;
    primStates[primID].nodeID = (uint32_t)-1;
    primStates[primID].done   = true;
    primStates[primID].primID = primID;
}

} // namespace builder_detail
} // namespace cuBQL