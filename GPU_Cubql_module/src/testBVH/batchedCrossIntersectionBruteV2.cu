// SPDX-FileCopyrightText: Copyright (c) 2026
// SPDX-License-Identifier: Apache-2.0

#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/execution_policy.h>
#include <iostream>
#include "include/third-party/cubql/fixedBoxQueryv2.h"

#ifndef CUBQL_CUDA_CALL
#define CUBQL_CUDA_CALL( call )                                                              \
{                                                                                            \
    cudaError_t err = cuda##call;                                                             \
    if( cudaSuccess != err ) {                                                               \
        fprintf(stderr, "CUDA error in file '%s' in line %i : %s.\n",                        \
                __FILE__, __LINE__, cudaGetErrorString( err ) );                             \
                exit(EXIT_FAILURE);                                                          \
    }                                                                                        \
}
#endif

// --------------------------------------------------------------------
// NATIVE INLINE BOX INTERSECT
// --------------------------------------------------------------------
template<typename T, int D>
__device__ inline bool primitiveBoxOverlaps(
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

// Functor to calculate comparative combinations per localized pair sequence slot
struct CellPairChunkWorkFunctor {
    const uint32_t* nodeOffsetsA;
    const uint32_t* nodeOffsetsB;
    const uint32_t* intersectPairsA;
    const uint32_t* intersectPairsB;
    size_t globalOffset;

    __host__ __device__
    CellPairChunkWorkFunctor(const uint32_t* _offsetsA, const uint32_t* _offsetsB,
                             const uint32_t* _pairsA, const uint32_t* _pairsB, size_t _globalOffset)
        : nodeOffsetsA(_offsetsA), nodeOffsetsB(_offsetsB), 
          intersectPairsA(_pairsA), intersectPairsB(_pairsB), globalOffset(_globalOffset) {}

    __host__ __device__
    uint32_t operator()(const size_t i) const {
        size_t actualIdx = globalOffset + i;
        uint32_t cellA = intersectPairsA[actualIdx];
        uint32_t cellB = intersectPairsB[actualIdx];
        
        uint32_t sizeA = nodeOffsetsA[cellA + 1] - nodeOffsetsA[cellA];
        uint32_t sizeB = nodeOffsetsB[cellB + 1] - nodeOffsetsB[cellB];
        return sizeA * sizeB; 
    }
};

// --------------------------------------------------------------------
// KERNELS
// --------------------------------------------------------------------

__global__ void populatePrimitiveWorkMapKernel(
    uint2* d_threadToPrimMap,
    const uint32_t* d_intersectPairsA, const uint32_t* d_intersectPairsB,
    uint32_t chunkStart, uint32_t activePairsInChunk,
    const uint32_t* d_nodeOffsetsA, const uint32_t* d_nodeOffsetsB,
    const uint32_t* d_chunkWorkOffsets)
{
    uint32_t localPairIdx = threadIdx.x + blockIdx.x * blockDim.x;
    if (localPairIdx >= activePairsInChunk) return;

    uint32_t globalPairIdx = chunkStart + localPairIdx;
    uint32_t cellA = d_intersectPairsA[globalPairIdx];
    uint32_t cellB = d_intersectPairsB[globalPairIdx];

    uint32_t startA = d_nodeOffsetsA[cellA];
    uint32_t endA   = d_nodeOffsetsA[cellA + 1];
    uint32_t sizeA  = endA - startA;

    uint32_t startB = d_nodeOffsetsB[cellB];
    uint32_t endB   = d_nodeOffsetsB[cellB + 1];
    uint32_t sizeB  = endB - startB;

    if (sizeA == 0 || sizeB == 0) return;

    uint32_t baseWritePos = d_chunkWorkOffsets[localPairIdx];

    for (uint32_t i = 0; i < sizeA; ++i) {
        for (uint32_t j = 0; j < sizeB; ++j) {
            uint32_t threadWriteIdx = baseWritePos + (i * sizeB) + j;
            d_threadToPrimMap[threadWriteIdx] = make_uint2(startA + i, startB + j);
        }
    }
}

template <typename T, int D>
__global__ void evaluateFlatPrimitiveOverlapsKernel(
    uint32_t* d_matchFlags,
    const uint2* d_threadToPrimMap, uint32_t totalPrimitiveThreads,
    const uint32_t* d_outSortedPrimIDsA, const cuBQL::box_t<T,D>* d_boxesA,
    const uint32_t* d_outSortedPrimIDsB, const cuBQL::box_t<T,D>* d_boxesB)
{
    uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= totalPrimitiveThreads) return;

    uint2 lookup = d_threadToPrimMap[tid];
    
    uint32_t primIDA = d_outSortedPrimIDsA[lookup.x];
    uint32_t primIDB = d_outSortedPrimIDsB[lookup.y];

    if (primitiveBoxOverlaps<T, D>(d_boxesA[primIDA], d_boxesB[primIDB])) {
        d_matchFlags[tid] = 1;
    } else {
        d_matchFlags[tid] = 0;
    }
}

// --------------------------------------------------------------------
// HOST EXECUTABLE REGION
// --------------------------------------------------------------------
template <typename T, int D>
uint64_t executePrimitiveBoxCrossIntersectionCheck(
    int batchMultiplier,
    const thrust::device_vector<uint32_t>& d_intersectPairsA,
    const thrust::device_vector<uint32_t>& d_intersectPairsB,
    const uint32_t* d_outSortedPrimIDsA, const uint32_t* d_outNodeOffsetsA, const cuBQL::box_t<T,D>* d_boxesA,
    const uint32_t* d_outSortedPrimIDsB, const uint32_t* d_outNodeOffsetsB, const cuBQL::box_t<T,D>* d_boxesB,
    cudaStream_t s)
{
    uint32_t totalPairs = d_intersectPairsA.size();
    if (totalPairs == 0) return 0;

    if (batchMultiplier <= 0 || batchMultiplier > (int)totalPairs) {
        batchMultiplier = totalPairs;
    }

    uint64_t globalIntersectionCount = 0;

    thrust::device_vector<uint32_t> d_chunkSizes(batchMultiplier);
    thrust::device_vector<uint32_t> d_chunkOffsets(batchMultiplier);

    for (uint32_t i = 0; i < totalPairs; i += batchMultiplier) {
        uint32_t activePairsInChunk = std::min((uint32_t)batchMultiplier, totalPairs - i);

        CellPairChunkWorkFunctor chunkWorkFunctor(
            d_outNodeOffsetsA, d_outNodeOffsetsB, 
            thrust::raw_pointer_cast(d_intersectPairsA.data()), 
            thrust::raw_pointer_cast(d_intersectPairsB.data()), 
            i
        );
        
        thrust::transform(
            thrust::cuda::par.on(s),
            thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(activePairsInChunk),
            d_chunkSizes.begin(),
            chunkWorkFunctor
        );

        thrust::exclusive_scan(
            thrust::cuda::par.on(s),
            d_chunkSizes.begin(),
            d_chunkSizes.begin() + activePairsInChunk,
            d_chunkOffsets.begin()
        );

        uint32_t lastSize = 0;
        uint32_t lastOffset = 0;
        CUBQL_CUDA_CALL(MemcpyAsync(&lastSize, thrust::raw_pointer_cast(d_chunkSizes.data()) + activePairsInChunk - 1, sizeof(uint32_t), cudaMemcpyDeviceToHost, s));
        CUBQL_CUDA_CALL(MemcpyAsync(&lastOffset, thrust::raw_pointer_cast(d_chunkOffsets.data()) + activePairsInChunk - 1, sizeof(uint32_t), cudaMemcpyDeviceToHost, s));
        CUBQL_CUDA_CALL(StreamSynchronize(s));

        uint32_t totalPrimitiveThreads = lastSize + lastOffset;
        if (totalPrimitiveThreads == 0) continue;

        thrust::device_vector<uint2> d_threadToPrimMap(totalPrimitiveThreads);
        thrust::device_vector<uint32_t> d_matchFlags(totalPrimitiveThreads);

        uint32_t populateBlockSize = 128;
        uint32_t populateGridSize  = (activePairsInChunk + populateBlockSize - 1) / populateBlockSize;

        populatePrimitiveWorkMapKernel<<<populateGridSize, populateBlockSize, 0, s>>>(
            thrust::raw_pointer_cast(d_threadToPrimMap.data()),
            thrust::raw_pointer_cast(d_intersectPairsA.data()),
            thrust::raw_pointer_cast(d_intersectPairsB.data()),
            i, activePairsInChunk,
            d_outNodeOffsetsA, d_outNodeOffsetsB,
            thrust::raw_pointer_cast(d_chunkOffsets.data())
        );

        uint32_t evalBlockSize = 256;
        uint32_t evalGridSize  = (totalPrimitiveThreads + evalBlockSize - 1) / evalBlockSize;

        evaluateFlatPrimitiveOverlapsKernel<T, D><<<evalGridSize, evalBlockSize, 0, s>>>(
            thrust::raw_pointer_cast(d_matchFlags.data()),
            thrust::raw_pointer_cast(d_threadToPrimMap.data()),
            totalPrimitiveThreads,
            d_outSortedPrimIDsA, d_boxesA,
            d_outSortedPrimIDsB, d_boxesB
        );

        uint64_t chunkIntersections = thrust::reduce(
            thrust::cuda::par.on(s),
            d_matchFlags.begin(),
            d_matchFlags.end(),
            (uint64_t)0,
            thrust::plus<uint64_t>()
        );

        globalIntersectionCount += chunkIntersections;
    }

    CUBQL_CUDA_CALL(StreamSynchronize(s));
    return globalIntersectionCount;
}

// --------------------------------------------------------------------
// EXPLICIT INSTANTIATION FOR TEST PIPELINE
// --------------------------------------------------------------------
template uint64_t executePrimitiveBoxCrossIntersectionCheck<float, 3>(
    int batchMultiplier,
    const thrust::device_vector<uint32_t>& d_intersectPairsA,
    const thrust::device_vector<uint32_t>& d_intersectPairsB,
    const uint32_t* d_outSortedPrimIDsA, const uint32_t* d_outNodeOffsetsA, const cuBQL::box_t<float, 3>* d_boxesA,
    const uint32_t* d_outSortedPrimIDsB, const uint32_t* d_outNodeOffsetsB, const cuBQL::box_t<float, 3>* d_boxesB,
    cudaStream_t s
);