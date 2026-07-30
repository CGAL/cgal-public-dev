#ifndef BATCHED_CROSS_INTERSECTION_COMMON_H
#define BATCHED_CROSS_INTERSECTION_COMMON_H

#include <cuda_runtime.h>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include "cuBQL/bvh.h"
#include "samples/common/loadOBJ.h"
#include "include/third-party/cubql/fixedBoxQueryv2.h"

// Shared helper macro for CUDA error checking
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
// SHARED FUNCTORS
// --------------------------------------------------------------------

struct BatchSizeInChunkFunctor {
    const uint32_t* outPairsB;
    const uint32_t* outOffsetsB;
    const uint32_t* reverseMapB;
    uint32_t outMarkedCountB;
    uint32_t totalPrimsB;
    int chunkStartBatchIdx;
    int activeBatchesInChunk;

    __host__ __device__
    BatchSizeInChunkFunctor(const uint32_t* _pairsB, const uint32_t* _offsetsB, const uint32_t* _reverseMapB,
                            uint32_t _markedCount, uint32_t _totalPrims,
                            int _chunkStartBatchIdx, int _activeBatchesInChunk)
        : outPairsB(_pairsB), outOffsetsB(_offsetsB), reverseMapB(_reverseMapB), outMarkedCountB(_markedCount),
          totalPrimsB(_totalPrims), chunkStartBatchIdx(_chunkStartBatchIdx),
          activeBatchesInChunk(_activeBatchesInChunk) {}

    __host__ __device__
    int operator()(const int b) const {
        if (b >= activeBatchesInChunk) return 0;
        uint32_t bIdx = outPairsB[chunkStartBatchIdx + b];
        uint32_t markedIdxB = reverseMapB[bIdx];
        uint32_t startOffset = outOffsetsB[markedIdxB];
        uint32_t endOffset = (markedIdxB + 1 < outMarkedCountB) ? outOffsetsB[markedIdxB + 1] : totalPrimsB;
        return (int)(endOffset - startOffset);
    }
};

// --------------------------------------------------------------------
// SHARED KERNEL DECLARATIONS
// --------------------------------------------------------------------

__global__ void computeMaxChunkSizeKernel(
    const uint32_t* outPairsB,
    const uint32_t* outOffsetsB,
    const uint32_t* reverseMapB,
    uint32_t outMarkedCountB,
    uint32_t totalPrimsB,
    int totalBatches,
    int batchMultiplier,
    int numChunks,
    int* d_globalMax
);

__global__ void assembleChunkBuffersByBatchKernel(
    uint32_t* d_BIter,
    uint64_t* d_AIter,
    const uint32_t* outPairsA,
    const uint32_t* outPairsB,
    const uint32_t* reverseMapB,
    const uint32_t* outOffsetsB,
    const uint32_t* outPrimsFlatB,
    uint32_t outMarkedCountB,
    uint32_t totalPrimsB,
    int chunkStartBatchIdx,
    int activeBatchesInChunk,
    const int* d_chunkBatchOffsets 
);

__global__ void countAABBOverlapsKernel_Indirected(
    int *pairCounts, 
    cuBQL::bvh3f bvhA, 
    const cuBQL::Triangle *triA, 
    const cuBQL::Triangle *triB, 
    const uint32_t* d_BIter,       
    uint32_t startOffsetB,         
    int numPrimsB, 
    const uint64_t* d_AIter
);

__global__ void fillAABBOverlapsKernel_Indirected(
    int2 *candidatePairs, 
    const int *offsets, 
    cuBQL::bvh3f bvhA, 
    const cuBQL::Triangle *triA, 
    const cuBQL::Triangle *triB, 
    const uint32_t* d_BIter,       
    uint32_t startOffsetB,         
    int numPrimsB, 
    const uint64_t* d_AIter
);

#endif // BATCHED_CROSS_INTERSECTION_COMMON_H