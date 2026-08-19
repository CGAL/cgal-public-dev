// SPDX-FileCopyrightText: Copyright (c) 2026
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/execution_policy.h>
#include "cuBQL/builder/cuda/builder_common.h" // Keeping native box types

// 1. Sleek box intersection helper (unchanged, still a classic)
template<typename T, int D>
__device__ inline bool boxesIntersect(
    const cuBQL::box_t<T,D>& a, 
    const cuBQL::box_t<T,D>& b) 
{
    for (int i = 0; i < D; ++i) {
        if (a.get_lower(i) > b.get_upper(i) || a.get_upper(i) < b.get_lower(i)) {
            return false;
        }
    }
    return true;
}

// 2. Pass 1: Count Intersections (Now with 100% less struct digging!)
template<typename T, int D>
__global__ void countBoxIntersections(
    const cuBQL::box_t<T,D>* d_boxesA, uint32_t numA,
    const cuBQL::box_t<T,D>* d_boxesB, uint32_t numB,
    uint32_t* d_perThreadCounts) 
{
    uint32_t idxA = threadIdx.x + blockIdx.x * blockDim.x;
    if (idxA >= numA) return;

    // Direct memory loading, cache hits go BRRRRR
    const auto boxA = d_boxesA[idxA]; 
    uint32_t localIntersections = 0;

    for (uint32_t idxB = 0; idxB < numB; ++idxB) {
        if (boxesIntersect<T,D>(boxA, d_boxesB[idxB])) { 
            localIntersections++; 
        }
    }
    d_perThreadCounts[idxA] = localIntersections;
}

// 3. Pass 2: Fill Intersections (The ultra-clean scatter)
template<typename T, int D>
__global__ void fillBoxIntersections(
    const cuBQL::box_t<T,D>* d_boxesA, uint32_t numA,
    const cuBQL::box_t<T,D>* d_boxesB, uint32_t numB,
    const uint32_t* d_offsets,
    uint32_t* d_outPairsA, uint32_t* d_outPairsB) 
{
    uint32_t idxA = threadIdx.x + blockIdx.x * blockDim.x;
    if (idxA >= numA) return;

    const auto boxA = d_boxesA[idxA];
    uint32_t writePos = d_offsets[idxA];

    for (uint32_t idxB = 0; idxB < numB; ++idxB) {
        if (boxesIntersect<T,D>(boxA, d_boxesB[idxB])) { 
            d_outPairsA[writePos] = idxA;
            d_outPairsB[writePos] = idxB;
            writePos++;
        }
    }
}

// 4. The Grand Orchestrator (Simblified Edition)
template<typename T, int D>
uint32_t executeBoxCrossCheck(
    const cuBQL::box_t<T,D>* d_boxesA, uint32_t numA,
    const cuBQL::box_t<T,D>* d_boxesB, uint32_t numB,
    thrust::device_vector<uint32_t>& d_outPairsA,
    thrust::device_vector<uint32_t>& d_outPairsB,
    cudaStream_t s)
{
    if (numA == 0 || numB == 0) return 0;

    uint32_t threadsCross = 256;
    uint32_t blocksCross = (numA + threadsCross - 1) / threadsCross;

    auto exec_policy = thrust::cuda::par.on(s);

    // Pass 1: Count those overlaps!
    thrust::device_vector<uint32_t> d_perThreadCounts(numA, 0);
    countBoxIntersections<T,D><<<blocksCross, threadsCross, 0, s>>>(
        d_boxesA, numA,
        d_boxesB, numB,
        thrust::raw_pointer_cast(d_perThreadCounts.data())
    );

    // Prefix Scan magic
    thrust::device_vector<uint32_t> d_offsets(numA, 0);
    thrust::exclusive_scan(exec_policy, d_perThreadCounts.begin(), d_perThreadCounts.end(), d_offsets.begin());
    
    uint32_t lastCount = d_perThreadCounts.back();
    uint32_t lastOffset = d_offsets.back();
    uint32_t totalIntersections = lastOffset + lastCount;

    if (totalIntersections > 0) {
        d_outPairsA.resize(totalIntersections);
        d_outPairsB.resize(totalIntersections);

        // Pass 2: Fill those pair vectors!
        fillBoxIntersections<T,D><<<blocksCross, threadsCross, 0, s>>>(
            d_boxesA, numA,
            d_boxesB, numB,
            thrust::raw_pointer_cast(d_offsets.data()),
            thrust::raw_pointer_cast(d_outPairsA.data()),
            thrust::raw_pointer_cast(d_outPairsB.data())
        );
    }

    return totalIntersections;
}