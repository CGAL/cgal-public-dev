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
#include "batchedCrossIntersection.h"
#include "../custom_pipeline/GPUPredicatesCheck.h"
#include "TargetStatus.h"

// Helper macro matching your project's naming convention
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

// Functor matching your batch mapping strategy to quickly populate single batch sizes directly on GPU
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
// KERNEL DEFINITIONS
// --------------------------------------------------------------------

// Coalesced parallel maximum chunk space calculator 
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
) {
    // Each block manages exactly one chunk
    int chunkIdx = blockIdx.x;
    if (chunkIdx >= numChunks) return;

    int startBatch = chunkIdx * batchMultiplier;
    int activeBatchesInChunk = (startBatch + batchMultiplier <= totalBatches) 
                               ? batchMultiplier 
                               : (totalBatches - startBatch);

    int threadSum = 0;

    // Grid-stride loop within the block for continuous, coalesced memory reads
    for (int b = threadIdx.x; b < activeBatchesInChunk; b += blockDim.x) {
        uint32_t bIdx = outPairsB[startBatch + b];
        uint32_t markedIdxB = reverseMapB[bIdx];
        uint32_t startOffset = outOffsetsB[markedIdxB];
        uint32_t endOffset = (markedIdxB + 1 < outMarkedCountB) ? outOffsetsB[markedIdxB + 1] : totalPrimsB;
        threadSum += (int)(endOffset - startOffset);
    }

    // Shared memory reduction allocation dynamically passed from the host
    extern __shared__ int sdata[];
    int tid = threadIdx.x;
    sdata[tid] = threadSum;
    __syncthreads();

    // Perform intra-block tree reduction
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sdata[tid] += sdata[tid + s];
        }
        __syncthreads();
    }

    // Accumulate the block total into the final maximum bounds using atomics
    if (tid == 0) {
        atomicMax(d_globalMax, sdata[0]);
    }
}

// Upgraded Parallel Batch-centric Assembly Kernel (Zero Binary Search)
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
) {
    int b = threadIdx.x + blockIdx.x * blockDim.x;
    if (b >= activeBatchesInChunk) return;

    int currentBatchArrayIdx = chunkStartBatchIdx + b;
    
    uint32_t bIdx = outPairsB[currentBatchArrayIdx];
    uint32_t markedIdxB = reverseMapB[bIdx];
    
    uint32_t startOffsetB = outOffsetsB[markedIdxB];
    uint32_t endOffsetB = (markedIdxB + 1 < outMarkedCountB) ? outOffsetsB[markedIdxB + 1] : totalPrimsB;
    int numPrims = (int)(endOffsetB - startOffsetB);
    
    if (numPrims <= 0) return;

    int writeSandboxOffset = d_chunkBatchOffsets[b];

    // Mesh A bypasses secondary lookup tables and reads directly from outPairsA
    uint64_t startNodeIdxA = (uint64_t)(outPairsA[currentBatchArrayIdx]);

    for (int p = 0; p < numPrims; ++p) {
        int targetWriteIdx = writeSandboxOffset + p;
        d_BIter[targetWriteIdx] = outPrimsFlatB[startOffsetB + p];
        d_AIter[targetWriteIdx] = startNodeIdxA;
    }
}

__global__ void countAABBOverlapsKernel_Indirected(
    int *pairCounts, 
    cuBQL::bvh3f bvhA, 
    const cuBQL::Triangle *triA, 
    const cuBQL::Triangle *triB, 
    const uint32_t* d_BIter,       
    uint32_t startOffsetB,         
    int numPrimsB, 
    const uint64_t* d_AIter)       
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= numPrimsB) return;

    uint32_t actualPrimIdB = d_BIter[startOffsetB + tid];
    uint64_t startNodeIdxA = d_AIter[startOffsetB + tid]; 
    cuBQL::Triangle b = triB[actualPrimIdB];
    cuBQL::box3f query = b.bounds();
    
    int count = 0;
    cuBQL::fixedBoxQueryv2::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        for (uint32_t i = 0; i < num; i++) {
            if (triA[ids[i]].bounds().overlaps(query)) { count++; }
        }
        return CUBQL_CONTINUE_TRAVERSAL;
    }, bvhA, query, startNodeIdxA);
    
    pairCounts[tid] = count;
}

__global__ void fillAABBOverlapsKernel_Indirected(
    int2 *candidatePairs, 
    const int *offsets, 
    cuBQL::bvh3f bvhA, 
    const cuBQL::Triangle *triA, 
    const cuBQL::Triangle *triB, 
    const uint32_t* d_BIter,       
    uint32_t startOffsetB,         
    int numPrimsB, 
    const uint64_t* d_AIter)       
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= numPrimsB) return;

    int wPos = offsets[tid];
    uint32_t actualPrimIdB = d_BIter[startOffsetB + tid];
    uint64_t startNodeIdxA = d_AIter[startOffsetB + tid]; 
    cuBQL::Triangle b = triB[actualPrimIdB];
    cuBQL::box3f query = b.bounds();
    
    cuBQL::fixedBoxQueryv2::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        for (uint32_t i = 0; i < num; i++) {
            if (triA[ids[i]].bounds().overlaps(query)) { 
                candidatePairs[wPos++] = make_int2((int)ids[i], (int)actualPrimIdB); 
            }
        }
        return CUBQL_CONTINUE_TRAVERSAL;
    }, bvhA, query, startNodeIdxA);
}

// --------------------------------------------------------------------
// HOST EXECUTABLE REGION
// --------------------------------------------------------------------
uint64_t executeBatchedCrossIntersectionLoop(
    int batchMultiplier,
    int totalBatches,
    const thrust::device_vector<uint32_t>& d_outPairsA,
    const thrust::device_vector<uint32_t>& d_outPairsB,
    const thrust::device_vector<uint32_t>& d_reverseMapB,
    const thrust::device_vector<uint32_t>& d_markedNodeIndicesB, 
    const thrust::device_vector<uint32_t>& d_outOffsetsB,
    const thrust::device_vector<uint32_t>& d_outPrimsFlatB,
    const thrust::device_vector<uint32_t>& d_nodeDescendantCountsB, 
    uint32_t h_outMarkedCountB,
    const cuBQL::bvh3f& bvhA,
    const cuBQL::Triangle* dMeshA,
    const cuBQL::Triangle* dMeshB,
    std::vector<int2>& hGreenPairs,
    std::vector<int2>& hYellowPairs,
    IntersectionTimeTracker& tracker 
) {
    double tTotalStart = cuBQL::getCurrentTime();
    double tAllocStart = tTotalStart;

    uint64_t finalCandidatePairs = 0;

    // 1. Unpack Raw Resource Pointers
    const uint32_t* ptr_outPairsA           = thrust::raw_pointer_cast(d_outPairsA.data());
    const uint32_t* ptr_outPairsB           = thrust::raw_pointer_cast(d_outPairsB.data());
    const uint32_t* ptr_reverseMapB         = thrust::raw_pointer_cast(d_reverseMapB.data());
    const uint32_t* ptr_outOffsetsB         = thrust::raw_pointer_cast(d_outOffsetsB.data());
    const uint32_t* ptr_outPrimsFlatB       = thrust::raw_pointer_cast(d_outPrimsFlatB.data());

    if (totalBatches > 0) {
        batchMultiplier = std::min(batchMultiplier, totalBatches);
    } else {
        batchMultiplier = 1; 
    }

    // 2. Compute Exact Maximum Chunk Space Bounds (Upgraded Parallel Pass)
    int maxChunkPrims = 0;
    if (totalBatches > 0) {
        int numChunks = (totalBatches + batchMultiplier - 1) / batchMultiplier;
        
        int* d_globalMax = nullptr;
        CUBQL_CUDA_CALL(Malloc(&d_globalMax, sizeof(int)));
        CUBQL_CUDA_CALL(Memset(d_globalMax, 0, sizeof(int)));

        int blockSize = 256; 
        int gridSize = numChunks;
        size_t sharedMemSize = blockSize * sizeof(int);

        // Executes cooperatively with total GPU parallelism over wide chunks
        computeMaxChunkSizeKernel<<<gridSize, blockSize, sharedMemSize>>>(
            ptr_outPairsB, ptr_outOffsetsB, ptr_reverseMapB, 
            h_outMarkedCountB, d_outPrimsFlatB.size(), 
            totalBatches, batchMultiplier, numChunks, d_globalMax
        );

        CUBQL_CUDA_CALL(Memcpy(&maxChunkPrims, d_globalMax, sizeof(int), cudaMemcpyDeviceToHost));
        CUBQL_CUDA_CALL(Free(d_globalMax));
    }

    // 3. Preallocate Memory Sandboxes
    uint32_t* d_BIter = nullptr;
    uint64_t* d_AIter = nullptr;
    int* d_pairCounts = nullptr;
    int* d_offsets    = nullptr;
    int* d_chunkBatchSizes   = nullptr;
    int* d_chunkBatchOffsets = nullptr;

    if (maxChunkPrims > 0) {
        CUBQL_CUDA_CALL(Malloc(&d_BIter, maxChunkPrims * sizeof(uint32_t)));
        CUBQL_CUDA_CALL(Malloc(&d_AIter, maxChunkPrims * sizeof(uint64_t)));
        CUBQL_CUDA_CALL(Malloc(&d_pairCounts, maxChunkPrims * sizeof(int)));
        CUBQL_CUDA_CALL(Malloc(&d_offsets, maxChunkPrims * sizeof(int)));
        CUBQL_CUDA_CALL(Malloc(&d_chunkBatchSizes, batchMultiplier * sizeof(int)));
        CUBQL_CUDA_CALL(Malloc(&d_chunkBatchOffsets, batchMultiplier * sizeof(int)));
    }
    tracker.preallocateTimeMs = (cuBQL::getCurrentTime() - tAllocStart) * 1000.0;

    cudaEvent_t evComputeStart, evComputeEnd;
    CUBQL_CUDA_CALL(EventCreate(&evComputeStart));
    CUBQL_CUDA_CALL(EventCreate(&evComputeEnd));

    // 4. Coarse-Grained Execution Chunk Loop
    for (int i = 0; i < totalBatches; i += batchMultiplier) {
        
        // --- ASSEMBLY TIMING START ---
        double tAssemblyStart = cuBQL::getCurrentTime();
        int activeBatchesInChunk = std::min(batchMultiplier, totalBatches - i);

        auto localBatchIdxIterator = thrust::make_counting_iterator(0);
        auto batchSizeTransformIterator = thrust::make_transform_iterator(
            localBatchIdxIterator,
            BatchSizeInChunkFunctor(ptr_outPairsB, ptr_outOffsetsB, ptr_reverseMapB, h_outMarkedCountB, d_outPrimsFlatB.size(), i, activeBatchesInChunk)
        );

        thrust::device_ptr<int> dev_batchSizes(d_chunkBatchSizes);
        thrust::copy_n(thrust::device, batchSizeTransformIterator, activeBatchesInChunk, dev_batchSizes);

        thrust::device_ptr<int> dev_batchOffsets(d_chunkBatchOffsets);
        thrust::exclusive_scan(thrust::device, dev_batchSizes, dev_batchSizes + activeBatchesInChunk, dev_batchOffsets);

        int lastBatchSize = 0;
        int lastBatchOffset = 0;
        CUBQL_CUDA_CALL(MemcpyAsync(&lastBatchSize, d_chunkBatchSizes + activeBatchesInChunk - 1, sizeof(int), cudaMemcpyDeviceToHost));
        CUBQL_CUDA_CALL(Memcpy(&lastBatchOffset, d_chunkBatchOffsets + activeBatchesInChunk - 1, sizeof(int), cudaMemcpyDeviceToHost));
        int totalChunkPrims = lastBatchSize + lastBatchOffset;

        if (totalChunkPrims == 0) continue;

        int assembleBlockSize = 64;
        int assembleGridSize  = (activeBatchesInChunk + assembleBlockSize - 1) / assembleBlockSize;

        assembleChunkBuffersByBatchKernel<<<assembleGridSize, assembleBlockSize>>>(
            d_BIter, d_AIter, ptr_outPairsA, ptr_outPairsB, ptr_reverseMapB,
            ptr_outOffsetsB, ptr_outPrimsFlatB, h_outMarkedCountB, d_outPrimsFlatB.size(),
            i, activeBatchesInChunk, d_chunkBatchOffsets
        );

        tracker.assemblyPhaseMs += (cuBQL::getCurrentTime() - tAssemblyStart) * 1000.0;
        // --- ASSEMBLY TIMING END ---

        // --- CUDA EXECUTION TIMING START ---
        CUBQL_CUDA_CALL(EventRecord(evComputeStart, 0));

        int kernelBlockSize = 128;
        int kernelGridSize  = (totalChunkPrims + kernelBlockSize - 1) / kernelBlockSize;

        countAABBOverlapsKernel_Indirected<<<kernelGridSize, kernelBlockSize>>>(
            d_pairCounts, bvhA, const_cast<cuBQL::Triangle*>(dMeshA), const_cast<cuBQL::Triangle*>(dMeshB), 
            d_BIter, 0, totalChunkPrims, d_AIter
        );

        thrust::device_ptr<int> dev_counts(d_pairCounts);
        thrust::device_ptr<int> dev_offsets(d_offsets);
        thrust::exclusive_scan(thrust::device, dev_counts, dev_counts + totalChunkPrims, dev_offsets);

        int lastCount = 0;
        int lastOffset = 0;
        CUBQL_CUDA_CALL(Memcpy(&lastCount, d_pairCounts + totalChunkPrims - 1, sizeof(int), cudaMemcpyDeviceToHost));
        CUBQL_CUDA_CALL(Memcpy(&lastOffset, d_offsets + totalChunkPrims - 1, sizeof(int), cudaMemcpyDeviceToHost));
        int totalChunkPairs = lastCount + lastOffset;

        if (totalChunkPairs == 0) {
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

        fillAABBOverlapsKernel_Indirected<<<kernelGridSize, kernelBlockSize>>>(
            d_candidatePairs, d_offsets, bvhA, const_cast<cuBQL::Triangle*>(dMeshA), const_cast<cuBQL::Triangle*>(dMeshB), 
            d_BIter, 0, totalChunkPrims, d_AIter
        );

        CUBQL_CUDA_CALL(EventRecord(evComputeEnd, 0));
        CUBQL_CUDA_CALL(EventSynchronize(evComputeEnd));
        
        float chunkComputeMs = 0.0f;
        CUBQL_CUDA_CALL(EventElapsedTime(&chunkComputeMs, evComputeStart, evComputeEnd));
        tracker.executionPhaseMs += chunkComputeMs;
        // --- CUDA EXECUTION TIMING END ---

        // ====================================================================
        // NEW GEOMETRIC EVALUATION & TRUST COMPACTION PHASE
        // ====================================================================
        double tEvalStart = cuBQL::getCurrentTime();

        int* d_pairStatuses = nullptr;
        CUBQL_CUDA_CALL(Malloc(&d_pairStatuses, totalChunkPairs * sizeof(int)));

        double internalPredicateTimeDummy = 0.0;
        evaluateAndCompactPairs(
            d_candidatePairs, 
            d_pairStatuses, 
            dMeshA, 
            dMeshB, 
            totalChunkPairs, 
            internalPredicateTimeDummy
        );

        tracker.fineEvaluationPhaseMs += (cuBQL::getCurrentTime() - tEvalStart) * 1000.0;

        double tEvalStartTwo = cuBQL::getCurrentTime();

        thrust::device_ptr<int2> dev_evaluated(d_candidatePairs);
        thrust::device_ptr<int> dev_statuses(d_pairStatuses);

        thrust::device_vector<int2> dev_green(totalChunkPairs);
        thrust::device_vector<int2> dev_yellow(totalChunkPairs);

        thrust::device_ptr<int2> dev_green_out(thrust::raw_pointer_cast(dev_green.data()));
        thrust::device_ptr<int2> dev_yellow_out(thrust::raw_pointer_cast(dev_yellow.data()));

        auto green_end = thrust::copy_if(thrust::device, dev_evaluated, dev_evaluated + totalChunkPairs, dev_statuses, dev_green_out, IsTargetPairStatus{(int)PAIR_GREEN});
        auto yellow_end = thrust::copy_if(thrust::device, dev_evaluated, dev_evaluated + totalChunkPairs, dev_statuses, dev_yellow_out, IsTargetPairStatus{(int)PAIR_YELLOW});

        int totalGreen = green_end - dev_green_out;
        int totalYellow = yellow_end - dev_yellow_out;

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

        tracker.DownloadAndClean += (cuBQL::getCurrentTime() - tEvalStartTwo ) * 1000.0;
    }

    // --- CLEANUP TIMING START ---
    double tCleanupStart = cuBQL::getCurrentTime();

    CUBQL_CUDA_CALL(EventDestroy(evComputeStart));
    CUBQL_CUDA_CALL(EventDestroy(evComputeEnd));

    if (d_BIter)              CUBQL_CUDA_CALL(Free(d_BIter));
    if (d_AIter)              CUBQL_CUDA_CALL(Free(d_AIter));
    if (d_pairCounts)         CUBQL_CUDA_CALL(Free(d_pairCounts));
    if (d_offsets)            CUBQL_CUDA_CALL(Free(d_offsets));
    if (d_chunkBatchSizes)    CUBQL_CUDA_CALL(Free(d_chunkBatchSizes));
    if (d_chunkBatchOffsets)  CUBQL_CUDA_CALL(Free(d_chunkBatchOffsets));

    tracker.cleanupTimeMs = (cuBQL::getCurrentTime() - tCleanupStart) * 1000.0;
    tracker.totalLoopTimeMs = (cuBQL::getCurrentTime() - tTotalStart) * 1000.0;

    std::cout << "--------------------------------------------------\n\n";
    return finalCandidatePairs;
}