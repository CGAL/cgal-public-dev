#include "batchedCrossIntersection.h"
#include "include/third-party/cubql/fixedBoxQueryv2.h"
#include <thrust/device_ptr.h>
#include <thrust/scan.h>
#include <thrust/execution_policy.h>
#include <thrust/fill.h>
#include <thrust/iterator/counting_iterator.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include "../custom_pipeline/GPUPredicatesCheck.h"
#include "TargetStatus.h"

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

// Functors for Chunk Size Extraction
struct GetSideChunkSizeFunctor {
    const uint32_t* outPairs;
    const uint32_t* outOffsets;
    uint32_t outMarkedCount;
    uint32_t totalPrims;
    int totalBatches;
    int batchMultiplier;

    __host__ __device__
    GetSideChunkSizeFunctor(const uint32_t* _pairs, const uint32_t* _offsets, 
                            uint32_t _markedCount, uint32_t _totalPrims, 
                            int _totalBatches, int _batchMultiplier) 
        : outPairs(_pairs), outOffsets(_offsets), outMarkedCount(_markedCount), 
          totalPrims(_totalPrims), totalBatches(_totalBatches), batchMultiplier(_batchMultiplier) {}

    __host__ __device__
    int operator()(const int chunkIdx) const {
        int i = chunkIdx * batchMultiplier; 
        if (i >= totalBatches) return 0;

        int activeBatchesInChunk = (i + batchMultiplier <= totalBatches) ? batchMultiplier : (totalBatches - i);
        int chunkTotal = 0;

        for (int b = 0; b < activeBatchesInChunk; ++b) {
            uint32_t bIdx = outPairs[i + b];
            uint32_t startOffset = outOffsets[bIdx];
            uint32_t endOffset = (bIdx + 1 < outMarkedCount) ? outOffsets[bIdx + 1] : totalPrims;
            chunkTotal += (int)(endOffset - startOffset);
        }
        return chunkTotal;
    }
};

struct SideBatchSizeInChunkFunctor {
    const uint32_t* outPairs;
    const uint32_t* outOffsets;
    uint32_t outMarkedCount;
    uint32_t totalPrims;
    int chunkStartBatchIdx;
    int activeBatchesInChunk;

    __host__ __device__
    SideBatchSizeInChunkFunctor(const uint32_t* _pairs, const uint32_t* _offsets,
                                uint32_t _markedCount, uint32_t _totalPrims,
                                int _chunkStartBatchIdx, int _activeBatchesInChunk)
        : outPairs(_pairs), outOffsets(_offsets), outMarkedCount(_markedCount),
          totalPrims(_totalPrims), chunkStartBatchIdx(_chunkStartBatchIdx),
          activeBatchesInChunk(_activeBatchesInChunk) {}

    __host__ __device__
    int operator()(const int b) const {
        if (b >= activeBatchesInChunk) return 0;
        uint32_t bIdx = outPairs[chunkStartBatchIdx + b];
        uint32_t startOffset = outOffsets[bIdx];
        uint32_t endOffset = (bIdx + 1 < outMarkedCount) ? outOffsets[bIdx + 1] : totalPrims;
        return (int)(endOffset - startOffset);
    }
};

__global__ void computeBatchPairCountsKernel(int* d_batchPairCounts, const int* d_sizesA, const int* d_sizesB, int count) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < count) {
        d_batchPairCounts[idx] = d_sizesA[idx] * d_sizesB[idx];
    }
}

// Computes how many 16x16 macro blocks are required to fill the comparison matrix for each batch
__global__ void computeBatchBlockCountsKernel(int* d_batchBlockCounts, const int* d_sizesA, const int* d_sizesB, int count) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < count) {
        int blocksA = (d_sizesA[idx] + 16 - 1) / 16;
        int blocksB = (d_sizesB[idx] + 16 - 1) / 16;
        d_batchBlockCounts[idx] = blocksA * blocksB;
    }
}

// --------------------------------------------------------------------
// OPTIMIZED KERNEL DEFINITIONS
// --------------------------------------------------------------------

__global__ void assembleSideBufferKernel_BlockPerBatch(
    uint32_t* d_Iter,
    const uint32_t* outPairs,
    const uint32_t* outOffsets,
    const uint32_t* outPrimsFlat,
    int chunkStartBatchIdx,
    int activeBatchesInChunk,
    const int* d_chunkBatchOffsets,
    const int* d_chunkBatchSizes
) {
    int b = blockIdx.x;
    if (b >= activeBatchesInChunk) return;

    int numPrims = d_chunkBatchSizes[b];
    if (numPrims == 0) return;

    uint32_t batchPrimsIdx = outPairs[chunkStartBatchIdx + b];
    uint32_t startOffset = outOffsets[batchPrimsIdx];
    int writeOffset = d_chunkBatchOffsets[b];

    for (int p = threadIdx.x; p < numPrims; p += blockDim.x) {
        d_Iter[writeOffset + p] = outPrimsFlat[startOffset + p];
    }
}

__global__ void generatePairToBatchMapKernel(
    int* d_pairToBatchMap, 
    const int* d_batchPairOffsets, 
    const int* d_batchPairCounts, 
    int activeBatchesInChunk
) {
    int b = blockIdx.x;
    if (b >= activeBatchesInChunk) return;

    int startOffset = d_batchPairOffsets[b];
    int totalPairs  = d_batchPairCounts[b];

    for (int i = threadIdx.x; i < totalPairs; i += blockDim.x) {
        d_pairToBatchMap[startOffset + i] = b;
    }
}

// OPTIMIZED EVALUATION KERNEL: Shared Memory Tiled Execution Layout
__global__ void evaluatePairsKernel_Tiled(
    int* d_pairValid,
    const cuBQL::Triangle* triA,
    const cuBQL::Triangle* triB,
    const uint32_t* d_IterA,
    const uint32_t* d_IterB,
    const int* d_chunkBatchOffsetsA,
    const int* d_chunkBatchOffsetsB,
    const int* d_chunkBatchSizesA,
    const int* d_chunkBatchSizesB, 
    const int* d_batchPairOffsets,
    const int* d_batchBlockOffsets,
    int activeBatchesInChunk,
    int totalBlocksToCheck
) {
    int blockIdx1D = blockIdx.x;
    if (blockIdx1D >= totalBlocksToCheck) return;

    // 1. Coarse-grained Block Binary Search to safely identify the batch ID uniform to this block
    __shared__ int shared_b;
    if (threadIdx.x == 0 && threadIdx.y == 0) {
        int low = 0;
        int high = activeBatchesInChunk - 1;
        int ans = 0;
        while (low <= high) {
            int mid = low + (high - low) / 2;
            if (d_batchBlockOffsets[mid] <= blockIdx1D) {
                ans = mid;
                low = mid + 1;
            } else {
                high = mid - 1;
            }
        }
        shared_b = ans;
    }
    __syncthreads();
    int b = shared_b;

    int numPrimsA = d_chunkBatchSizesA[b];
    int numPrimsB = d_chunkBatchSizesB[b];
    int blocksB = (numPrimsB + 16 - 1) / 16;

    int relBlockIdx = blockIdx1D - d_batchBlockOffsets[b];
    int blockRow = relBlockIdx / blocksB;
    int blockCol = relBlockIdx % blocksB;

    int ty = threadIdx.y;
    int tx = threadIdx.x;

    int localIdxA = blockRow * 16 + ty;
    int localIdxB = blockCol * 16 + tx;

    // Shared local caches typed to exact geometry layouts
    using BoxType = decltype(triA[0].bounds());
    __shared__ BoxType sharedBoundsA[16];
    __shared__ BoxType sharedBoundsB[16];

    // 2. Coalesced Shared Memory Cooperative Loading Phase
    int linearId = ty * 16 + tx; 
    if (linearId < 16) {
        int loadIdxA = blockRow * 16 + linearId;
        if (loadIdxA < numPrimsA) {
            uint32_t primIdA = d_IterA[d_chunkBatchOffsetsA[b] + loadIdxA];
            sharedBoundsA[linearId] = triA[primIdA].bounds();
        }
    }
    if (linearId >= 16 && linearId < 32) {
        int loadIdxB = blockCol * 16 + (linearId - 16);
        if (loadIdxB < numPrimsB) {
            uint32_t primIdB = d_IterB[d_chunkBatchOffsetsB[b] + loadIdxB];
            sharedBoundsB[linearId - 16] = triB[primIdB].bounds();
        }
    }
    __syncthreads();

    // 3. Complete Intersection Pass entirely inside local shared boundaries
    if (localIdxA < numPrimsA && localIdxB < numPrimsB) {
        bool overlaps = sharedBoundsA[ty].overlaps(sharedBoundsB[tx]);
        int gIdx = d_batchPairOffsets[b] + (localIdxA * numPrimsB) + localIdxB;
        d_pairValid[gIdx] = overlaps ? 1 : 0;
    }
}

// Collection Kernel: Retains linear mapping to perform compaction
__global__ void collectPairsKernel_BSFree(
    int2* candidatePairs,
    const int* d_pairValid,
    const int* d_pairOutputPositions,
    const uint32_t* d_IterA,
    const uint32_t* d_IterB,
    const int* d_chunkBatchOffsetsA,
    const int* d_chunkBatchOffsetsB,
    const int* d_chunkBatchSizesB,
    const int* d_batchPairOffsets,
    const int* d_pairToBatchMap,
    int totalPairsToCheck
) {
    int gIdx = threadIdx.x + blockIdx.x * blockDim.x;
    if (gIdx >= totalPairsToCheck || d_pairValid[gIdx] == 0) return;

    int wPos = d_pairOutputPositions[gIdx];

    int b = d_pairToBatchMap[gIdx];
    int relPairIdx = gIdx - d_batchPairOffsets[b];

    int numPrimsB = d_chunkBatchSizesB[b];
    int idxA = relPairIdx / numPrimsB;
    int idxB = relPairIdx % numPrimsB;

    uint32_t primIdA = d_IterA[d_chunkBatchOffsetsA[b] + idxA];
    uint32_t primIdB = d_IterB[d_chunkBatchOffsetsB[b] + idxB];

    candidatePairs[wPos] = make_int2((int)primIdA, (int)primIdB);
}

// --------------------------------------------------------------------
// HOST EXECUTABLE REGION
// --------------------------------------------------------------------
uint64_t executeBatchedCrossIntersectionLoop(
    int batchMultiplier,
    int totalBatches,
    const thrust::device_vector<uint32_t>& d_outPairsA,
    const thrust::device_vector<uint32_t>& d_outPairsB,
    const thrust::device_vector<uint32_t>& d_markedNodeIndicesA,
    const thrust::device_vector<uint32_t>& d_markedNodeIndicesB,
    const thrust::device_vector<uint32_t>& d_outOffsetsB,
    const thrust::device_vector<uint32_t>& d_outPrimsFlatB,
    const thrust::device_vector<uint32_t>& d_nodeDescendantCountsB, 
    const thrust::device_vector<uint32_t>& d_outOffsetsA,
    const thrust::device_vector<uint32_t>& d_outPrimsFlatA,
    const thrust::device_vector<uint32_t>& d_nodeDescendantCountsA, 
    uint32_t h_outMarkedCountA,
    uint32_t h_outMarkedCountB,
    const cuBQL::Triangle* dMeshA,
    const cuBQL::Triangle* dMeshB,
    std::vector<int2>& hGreenPairs,
    std::vector<int2>& hYellowPairs,
    IntersectionTimeTracker& tracker 
) {
    double tTotalStart = cuBQL::getCurrentTime();
    double tAllocStart = tTotalStart;

    uint64_t finalCandidatePairs = 0;

    const uint32_t* ptr_outPairsA     = thrust::raw_pointer_cast(d_outPairsA.data());
    const uint32_t* ptr_outPairsB     = thrust::raw_pointer_cast(d_outPairsB.data());
    const uint32_t* ptr_outOffsetsA   = thrust::raw_pointer_cast(d_outOffsetsA.data());
    const uint32_t* ptr_outOffsetsB   = thrust::raw_pointer_cast(d_outOffsetsB.data());
    const uint32_t* ptr_outPrimsFlatA = thrust::raw_pointer_cast(d_outPrimsFlatA.data());
    const uint32_t* ptr_outPrimsFlatB = thrust::raw_pointer_cast(d_outPrimsFlatB.data());

    if (totalBatches > 0) {
        batchMultiplier = std::min(batchMultiplier, totalBatches);
    } else {
        batchMultiplier = 1; 
    }

    int maxChunkPrimsA = 0;
    int maxChunkPrimsB = 0;
    if (totalBatches > 0) {
        int numChunks = (totalBatches + batchMultiplier - 1) / batchMultiplier;
        auto chunkIdxIterator = thrust::make_counting_iterator(0);
        
        maxChunkPrimsA = thrust::reduce(
            thrust::device,
            thrust::make_transform_iterator(chunkIdxIterator, GetSideChunkSizeFunctor(ptr_outPairsA, ptr_outOffsetsA, h_outMarkedCountA, d_outPrimsFlatA.size(), totalBatches, batchMultiplier)),
            thrust::make_transform_iterator(chunkIdxIterator, GetSideChunkSizeFunctor(ptr_outPairsA, ptr_outOffsetsA, h_outMarkedCountA, d_outPrimsFlatA.size(), totalBatches, batchMultiplier)) + numChunks,
            0, thrust::maximum<int>()
        );

        maxChunkPrimsB = thrust::reduce(
            thrust::device,
            thrust::make_transform_iterator(chunkIdxIterator, GetSideChunkSizeFunctor(ptr_outPairsB, ptr_outOffsetsB, h_outMarkedCountB, d_outPrimsFlatB.size(), totalBatches, batchMultiplier)),
            thrust::make_transform_iterator(chunkIdxIterator, GetSideChunkSizeFunctor(ptr_outPairsB, ptr_outOffsetsB, h_outMarkedCountB, d_outPrimsFlatB.size(), totalBatches, batchMultiplier)) + numChunks,
            0, thrust::maximum<int>()
        );
    }

    uint32_t* d_IterA = nullptr;
    uint32_t* d_IterB = nullptr;
    
    int* d_chunkBatchSizesA   = nullptr;
    int* d_chunkBatchOffsetsA = nullptr;
    int* d_chunkBatchSizesB   = nullptr;
    int* d_chunkBatchOffsetsB = nullptr;
    
    int* d_batchPairCounts    = nullptr;
    int* d_batchPairOffsets   = nullptr;

    // Grid tracking maps for tiled block scanning boundaries
    int* d_batchBlockCounts   = nullptr;
    int* d_batchBlockOffsets  = nullptr;

    if (maxChunkPrimsA > 0) CUBQL_CUDA_CALL(Malloc(&d_IterA, maxChunkPrimsA * sizeof(uint32_t)));
    if (maxChunkPrimsB > 0) CUBQL_CUDA_CALL(Malloc(&d_IterB, maxChunkPrimsB * sizeof(uint32_t)));

    CUBQL_CUDA_CALL(Malloc(&d_chunkBatchSizesA, batchMultiplier * sizeof(int)));
    CUBQL_CUDA_CALL(Malloc(&d_chunkBatchOffsetsA, batchMultiplier * sizeof(int)));
    CUBQL_CUDA_CALL(Malloc(&d_chunkBatchSizesB, batchMultiplier * sizeof(int)));
    CUBQL_CUDA_CALL(Malloc(&d_chunkBatchOffsetsB, batchMultiplier * sizeof(int)));
    
    CUBQL_CUDA_CALL(Malloc(&d_batchPairCounts, batchMultiplier * sizeof(int)));
    CUBQL_CUDA_CALL(Malloc(&d_batchPairOffsets, batchMultiplier * sizeof(int)));

    CUBQL_CUDA_CALL(Malloc(&d_batchBlockCounts, batchMultiplier * sizeof(int)));
    CUBQL_CUDA_CALL(Malloc(&d_batchBlockOffsets, batchMultiplier * sizeof(int)));

    tracker.preallocateTimeMs = (cuBQL::getCurrentTime() - tAllocStart) * 1000.0;

    cudaEvent_t evComputeStart, evComputeEnd;
    CUBQL_CUDA_CALL(EventCreate(&evComputeStart));
    CUBQL_CUDA_CALL(EventCreate(&evComputeEnd));

    for (int i = 0; i < totalBatches; i += batchMultiplier) {
        
        double tAssemblyStart = cuBQL::getCurrentTime();
        int activeBatchesInChunk = std::min(batchMultiplier, totalBatches - i);
        auto localBatchIdxIterator = thrust::make_counting_iterator(0);

        // Scan Side A Layout
        auto transformA = thrust::make_transform_iterator(localBatchIdxIterator, SideBatchSizeInChunkFunctor(ptr_outPairsA, ptr_outOffsetsA, h_outMarkedCountA, d_outPrimsFlatA.size(), i, activeBatchesInChunk));
        thrust::copy_n(thrust::device, transformA, activeBatchesInChunk, thrust::device_ptr<int>(d_chunkBatchSizesA));
        thrust::exclusive_scan(thrust::device, thrust::device_ptr<int>(d_chunkBatchSizesA), thrust::device_ptr<int>(d_chunkBatchSizesA) + activeBatchesInChunk, thrust::device_ptr<int>(d_chunkBatchOffsetsA));

        // Scan Side B Layout
        auto transformB = thrust::make_transform_iterator(localBatchIdxIterator, SideBatchSizeInChunkFunctor(ptr_outPairsB, ptr_outOffsetsB, h_outMarkedCountB, d_outPrimsFlatB.size(), i, activeBatchesInChunk));
        thrust::copy_n(thrust::device, transformB, activeBatchesInChunk, thrust::device_ptr<int>(d_chunkBatchSizesB));
        thrust::exclusive_scan(thrust::device, thrust::device_ptr<int>(d_chunkBatchSizesB), thrust::device_ptr<int>(d_chunkBatchSizesB) + activeBatchesInChunk, thrust::device_ptr<int>(d_chunkBatchOffsetsB));

        int lastSizeA = 0, lastOffsetA = 0;
        CUBQL_CUDA_CALL(MemcpyAsync(&lastSizeA, d_chunkBatchSizesA + activeBatchesInChunk - 1, sizeof(int), cudaMemcpyDeviceToHost));
        CUBQL_CUDA_CALL(Memcpy(&lastOffsetA, d_chunkBatchOffsetsA + activeBatchesInChunk - 1, sizeof(int), cudaMemcpyDeviceToHost));
        int totalChunkPrimsA = lastSizeA + lastOffsetA;

        int lastSizeB = 0, lastOffsetB = 0;
        CUBQL_CUDA_CALL(MemcpyAsync(&lastSizeB, d_chunkBatchSizesB + activeBatchesInChunk - 1, sizeof(int), cudaMemcpyDeviceToHost));
        CUBQL_CUDA_CALL(Memcpy(&lastOffsetB, d_chunkBatchOffsetsB + activeBatchesInChunk - 1, sizeof(int), cudaMemcpyDeviceToHost));
        int totalChunkPrimsB = lastSizeB + lastOffsetB;

        if (totalChunkPrimsA == 0 || totalChunkPrimsB == 0) continue;

        int assemblyBlockSz = 128; 
        assembleSideBufferKernel_BlockPerBatch<<<activeBatchesInChunk, assemblyBlockSz>>>(
            d_IterA, ptr_outPairsA, ptr_outOffsetsA, ptr_outPrimsFlatA, 
            i, activeBatchesInChunk, d_chunkBatchOffsetsA, d_chunkBatchSizesA
        );
        assembleSideBufferKernel_BlockPerBatch<<<activeBatchesInChunk, assemblyBlockSz>>>(
            d_IterB, ptr_outPairsB, ptr_outOffsetsB, ptr_outPrimsFlatB, 
            i, activeBatchesInChunk, d_chunkBatchOffsetsB, d_chunkBatchSizesB
        );

        int blockSz = 256;
        int gridSetup = (activeBatchesInChunk + blockSz - 1) / blockSz;
        computeBatchPairCountsKernel<<<gridSetup, blockSz>>>(d_batchPairCounts, d_chunkBatchSizesA, d_chunkBatchSizesB, activeBatchesInChunk);
        thrust::exclusive_scan(thrust::device, thrust::device_ptr<int>(d_batchPairCounts), thrust::device_ptr<int>(d_batchPairCounts) + activeBatchesInChunk, thrust::device_ptr<int>(d_batchPairOffsets));

        // Scan Block configurations required for Tiled mapping arrays
        computeBatchBlockCountsKernel<<<gridSetup, blockSz>>>(d_batchBlockCounts, d_chunkBatchSizesA, d_chunkBatchSizesB, activeBatchesInChunk);
        thrust::exclusive_scan(thrust::device, thrust::device_ptr<int>(d_batchBlockCounts), thrust::device_ptr<int>(d_batchBlockCounts) + activeBatchesInChunk, thrust::device_ptr<int>(d_batchBlockOffsets));

        int lastPairCount = 0, lastPairOffset = 0;
        CUBQL_CUDA_CALL(MemcpyAsync(&lastPairCount, d_batchPairCounts + activeBatchesInChunk - 1, sizeof(int), cudaMemcpyDeviceToHost));
        CUBQL_CUDA_CALL(Memcpy(&lastPairOffset, d_batchPairOffsets + activeBatchesInChunk - 1, sizeof(int), cudaMemcpyDeviceToHost));
        int totalPairsToCheck = lastPairCount + lastPairOffset;

        int lastBlockCount = 0, lastBlockOffset = 0;
        CUBQL_CUDA_CALL(MemcpyAsync(&lastBlockCount, d_batchBlockCounts + activeBatchesInChunk - 1, sizeof(int), cudaMemcpyDeviceToHost));
        CUBQL_CUDA_CALL(Memcpy(&lastBlockOffset, d_batchBlockOffsets + activeBatchesInChunk - 1, sizeof(int), cudaMemcpyDeviceToHost));
        int totalBlocksToCheck = lastBlockCount + lastBlockOffset;

        tracker.assemblyPhaseMs += (cuBQL::getCurrentTime() - tAssemblyStart) * 1000.0;

        if (totalPairsToCheck == 0) continue;

        int* d_pairValid = nullptr;
        int* d_pairOutputPositions = nullptr;
        int* d_pairToBatchMap = nullptr; 
        
        CUBQL_CUDA_CALL(Malloc(&d_pairValid, totalPairsToCheck * sizeof(int)));
        CUBQL_CUDA_CALL(Malloc(&d_pairOutputPositions, totalPairsToCheck * sizeof(int)));
        CUBQL_CUDA_CALL(Malloc(&d_pairToBatchMap, totalPairsToCheck * sizeof(int)));

        CUBQL_CUDA_CALL(EventRecord(evComputeStart, 0));

        // Generate the flat combination-to-batch indexing map (retained exclusively for linear selection kernel)
        generatePairToBatchMapKernel<<<activeBatchesInChunk, 128>>>(
            d_pairToBatchMap, d_batchPairOffsets, d_batchPairCounts, activeBatchesInChunk
        );

        // Execute Tiled evaluation matrix configuration (16x16 2D Blocks)
        dim3 tiledBlockSize(16, 16);
        int tiledGridSize = totalBlocksToCheck;

        if (tiledGridSize > 0) {
            evaluatePairsKernel_Tiled<<<tiledGridSize, tiledBlockSize>>>(
                d_pairValid, dMeshA, dMeshB, d_IterA, d_IterB, 
                d_chunkBatchOffsetsA, d_chunkBatchOffsetsB, 
                d_chunkBatchSizesA, d_chunkBatchSizesB, 
                d_batchPairOffsets, d_batchBlockOffsets, 
                activeBatchesInChunk, totalBlocksToCheck
            );
        }

        thrust::exclusive_scan(thrust::device, thrust::device_ptr<int>(d_pairValid), thrust::device_ptr<int>(d_pairValid) + totalPairsToCheck, thrust::device_ptr<int>(d_pairOutputPositions));

        int lastValid = 0, lastValidOffset = 0;
        CUBQL_CUDA_CALL(MemcpyAsync(&lastValid, d_pairValid + totalPairsToCheck - 1, sizeof(int), cudaMemcpyDeviceToHost));
        CUBQL_CUDA_CALL(Memcpy(&lastValidOffset, d_pairOutputPositions + totalPairsToCheck - 1, sizeof(int), cudaMemcpyDeviceToHost));
        int totalChunkPairs = lastValid + lastValidOffset;

        if (totalChunkPairs == 0) {
            CUBQL_CUDA_CALL(Free(d_pairValid));
            CUBQL_CUDA_CALL(Free(d_pairOutputPositions));
            CUBQL_CUDA_CALL(Free(d_pairToBatchMap));
            CUBQL_CUDA_CALL(EventRecord(evComputeEnd, 0));
            CUBQL_CUDA_CALL(EventSynchronize(evComputeEnd));
            float elapsedMs = 0.0f;
            CUBQL_CUDA_CALL(EventElapsedTime(&elapsedMs, evComputeStart, evComputeEnd));
            tracker.executionPhaseMs += elapsedMs;
            continue;
        }

        finalCandidatePairs += totalChunkPairs;

        int2* d_candidatePairs = nullptr;
        CUBQL_CUDA_CALL(Malloc(&d_candidatePairs, totalChunkPairs * sizeof(int2)));

        int evaluationBlockSize = 256;
        int evaluationGridSize = (totalPairsToCheck + evaluationBlockSize - 1) / evaluationBlockSize;

        collectPairsKernel_BSFree<<<evaluationGridSize, evaluationBlockSize>>>(
            d_candidatePairs, d_pairValid, d_pairOutputPositions, d_IterA, d_IterB, 
            d_chunkBatchOffsetsA, d_chunkBatchOffsetsB, d_chunkBatchSizesB, 
            d_batchPairOffsets, d_pairToBatchMap, totalPairsToCheck
        );

        CUBQL_CUDA_CALL(EventRecord(evComputeEnd, 0));
        CUBQL_CUDA_CALL(EventSynchronize(evComputeEnd));
        float chunkComputeMs = 0.0f;
        CUBQL_CUDA_CALL(EventElapsedTime(&chunkComputeMs, evComputeStart, evComputeEnd));
        tracker.executionPhaseMs += chunkComputeMs;

        CUBQL_CUDA_CALL(Free(d_pairValid));
        CUBQL_CUDA_CALL(Free(d_pairOutputPositions));
        CUBQL_CUDA_CALL(Free(d_pairToBatchMap));

        // --- COMPACTION PHASE ---
        double tEvalStart = cuBQL::getCurrentTime();
        int* d_pairStatuses = nullptr;
        CUBQL_CUDA_CALL(Malloc(&d_pairStatuses, totalChunkPairs * sizeof(int)));

        double internalPredicateTimeDummy = 0.0;
        evaluateAndCompactPairs(d_candidatePairs, d_pairStatuses, dMeshA, dMeshB, totalChunkPairs, internalPredicateTimeDummy);

        tracker.fineEvaluationPhaseMs += (cuBQL::getCurrentTime() - tEvalStart) * 1000.0;

        double tEvalStartTwo = cuBQL::getCurrentTime();

        thrust::device_ptr<int2> dev_evaluated(d_candidatePairs);
        thrust::device_ptr<int> dev_statuses(d_pairStatuses);
        thrust::device_vector<int2> dev_green(totalChunkPairs);
        thrust::device_vector<int2> dev_yellow(totalChunkPairs);

        auto green_end = thrust::copy_if(thrust::device, dev_evaluated, dev_evaluated + totalChunkPairs, dev_statuses, thrust::device_ptr<int2>(thrust::raw_pointer_cast(dev_green.data())), IsTargetPairStatus{(int)PAIR_GREEN});
        auto yellow_end = thrust::copy_if(thrust::device, dev_evaluated, dev_evaluated + totalChunkPairs, dev_statuses, thrust::device_ptr<int2>(thrust::raw_pointer_cast(dev_yellow.data())), IsTargetPairStatus{(int)PAIR_YELLOW});

        int totalGreen = green_end - thrust::device_ptr<int2>(thrust::raw_pointer_cast(dev_green.data()));
        int totalYellow = yellow_end - thrust::device_ptr<int2>(thrust::raw_pointer_cast(dev_yellow.data()));

        if (totalGreen > 0) {
            size_t oldSize = hGreenPairs.size();
            hGreenPairs.resize(oldSize + totalGreen);
            CUBQL_CUDA_CALL(Memcpy(hGreenPairs.data() + oldSize, thrust::raw_pointer_cast(dev_green.data()), totalGreen * sizeof(int2), cudaMemcpyDeviceToHost));
        }
        if (totalYellow > 0) {
            size_t oldSize = hYellowPairs.size();
            hYellowPairs.resize(oldSize + totalYellow);
            CUBQL_CUDA_CALL(Memcpy(hYellowPairs.data() + oldSize, thrust::raw_pointer_cast(dev_yellow.data()), totalYellow * sizeof(int2), cudaMemcpyDeviceToHost));
        }

        CUBQL_CUDA_CALL(Free(d_candidatePairs));
        CUBQL_CUDA_CALL(Free(d_pairStatuses));

        tracker.DownloadAndClean += (cuBQL::getCurrentTime() - tEvalStartTwo) * 1000.0;
    }

    double tCleanupStart = cuBQL::getCurrentTime();
    CUBQL_CUDA_CALL(EventDestroy(evComputeStart));
    CUBQL_CUDA_CALL(EventDestroy(evComputeEnd));

    if (d_IterA)              CUBQL_CUDA_CALL(Free(d_IterA));
    if (d_IterB)              CUBQL_CUDA_CALL(Free(d_IterB));
    CUBQL_CUDA_CALL(Free(d_chunkBatchSizesA));
    CUBQL_CUDA_CALL(Free(d_chunkBatchOffsetsA));
    CUBQL_CUDA_CALL(Free(d_chunkBatchSizesB));
    CUBQL_CUDA_CALL(Free(d_chunkBatchOffsetsB));
    CUBQL_CUDA_CALL(Free(d_batchPairCounts));
    CUBQL_CUDA_CALL(Free(d_batchPairOffsets));
    CUBQL_CUDA_CALL(Free(d_batchBlockCounts));
    CUBQL_CUDA_CALL(Free(d_batchBlockOffsets));

    tracker.cleanupTimeMs = (cuBQL::getCurrentTime() - tCleanupStart) * 1000.0;
    tracker.totalLoopTimeMs = (cuBQL::getCurrentTime() - tTotalStart) * 1000.0;

    std::cout << "--------------------------------------------------\n\n";
    return finalCandidatePairs;
}