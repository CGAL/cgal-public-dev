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

// Functor to calculate the total primitive size of a grouped CHUNK starting at index chunkIdx
struct GetChunkSizeFunctor {
    const uint32_t* outPairsB;
    const uint32_t* outOffsetsB;
    uint32_t outMarkedCountB;
    uint32_t totalPrimsB;
    int totalBatches;
    int batchMultiplier;

    __host__ __device__
    GetChunkSizeFunctor(const uint32_t* _pairsB, const uint32_t* _offsetsB, 
                        uint32_t _markedCount, uint32_t _totalPrims, 
                        int _totalBatches, int _batchMultiplier) 
        : outPairsB(_pairsB), outOffsetsB(_offsetsB), outMarkedCountB(_markedCount), 
          totalPrimsB(_totalPrims), totalBatches(_totalBatches), batchMultiplier(_batchMultiplier) {}

    __host__ __device__
    int operator()(const int chunkIdx) const {
        int i = chunkIdx * batchMultiplier; // Map chunk sequence index to the true batch array index
        if (i >= totalBatches) return 0;

        int activeBatchesInChunk = (i + batchMultiplier <= totalBatches) ? batchMultiplier : (totalBatches - i);
        int chunkTotal = 0;

        for (int b = 0; b < activeBatchesInChunk; ++b) {
            uint32_t bIdx = outPairsB[i + b];
            uint32_t startOffset = outOffsetsB[bIdx];
            uint32_t endOffset = (bIdx + 1 < outMarkedCountB) ? outOffsetsB[bIdx + 1] : totalPrimsB;
            chunkTotal += (int)(endOffset - startOffset);
        }
        return chunkTotal;
    }
};

// --------------------------------------------------------------------
// KERNEL DEFINITIONS
// --------------------------------------------------------------------
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
    const thrust::device_vector<uint32_t>& d_markedNodeIndicesA,
    const thrust::device_vector<uint32_t>& d_markedNodeIndicesB,
    const thrust::device_vector<uint32_t>& d_outOffsetsB,
    const thrust::device_vector<uint32_t>& d_outPrimsFlatB,
    const thrust::device_vector<uint32_t>& d_nodeDescendantCountsB, 
    uint32_t h_outMarkedCountB,
    const cuBQL::bvh3f& bvhA,
    const cuBQL::Triangle* dMeshA,
    const cuBQL::Triangle* dMeshB,
    IntersectionTimeTracker& tracker 
) {
    double tTotalStart = cuBQL::getCurrentTime();
    double tAllocStart = tTotalStart;

    uint64_t finalCandidatePairs = 0;

    // 1. Unpack Raw Resource Pointers
    const uint32_t* ptr_outPairsA           = thrust::raw_pointer_cast(d_outPairsA.data());
    const uint32_t* ptr_outPairsB           = thrust::raw_pointer_cast(d_outPairsB.data());
    const uint32_t* ptr_markedNodeIndicesA  = thrust::raw_pointer_cast(d_markedNodeIndicesA.data());
    const uint32_t* ptr_outOffsetsB         = thrust::raw_pointer_cast(d_outOffsetsB.data());
    const uint32_t* ptr_outPrimsFlatB       = thrust::raw_pointer_cast(d_outPrimsFlatB.data());

    // Defensive Clamping Optimization: Protect against huge multiplier values causing integer shifts/overflows
    if (totalBatches > 0) {
        batchMultiplier = std::min(batchMultiplier, totalBatches);
    } else {
        batchMultiplier = 1; 
    }

    // 2. Compute Exact Maximum Chunk Space Bounds via Parallel Device Window Reduction
    int maxChunkPrims = 0;
    if (totalBatches > 0) {
        int numChunks = (totalBatches + batchMultiplier - 1) / batchMultiplier;
        
        auto chunkIdxIterator = thrust::make_counting_iterator(0);
        auto chunkSizeTransformIterator = thrust::make_transform_iterator(
            chunkIdxIterator,
            GetChunkSizeFunctor(ptr_outPairsB, ptr_outOffsetsB, h_outMarkedCountB, d_outPrimsFlatB.size(), totalBatches, batchMultiplier)
        );

        maxChunkPrims = thrust::reduce(
            thrust::device,
            chunkSizeTransformIterator,
            chunkSizeTransformIterator + numChunks,
            0,
            thrust::maximum<int>()
        );
    }

    // 3. Preallocate Memory Sandboxes tightly scaled to peak structural limits
    uint32_t* d_BIter = nullptr;
    uint64_t* d_AIter = nullptr;
    int* d_pairCounts = nullptr;
    int* d_offsets    = nullptr;

    if (maxChunkPrims > 0) {
        CUBQL_CUDA_CALL(Malloc(&d_BIter, maxChunkPrims * sizeof(uint32_t)));
        CUBQL_CUDA_CALL(Malloc(&d_AIter, maxChunkPrims * sizeof(uint64_t)));
        CUBQL_CUDA_CALL(Malloc(&d_pairCounts, maxChunkPrims * sizeof(int)));
        CUBQL_CUDA_CALL(Malloc(&d_offsets, maxChunkPrims * sizeof(int)));
    }
    tracker.preallocateTimeMs = (cuBQL::getCurrentTime() - tAllocStart) * 1000.0;

    // Initialize CUDA events for microsecond-accurate hardware runtime tracking
    cudaEvent_t evComputeStart, evComputeEnd;
    CUBQL_CUDA_CALL(EventCreate(&evComputeStart));
    CUBQL_CUDA_CALL(EventCreate(&evComputeEnd));

    // 4. Coarse-Grained Execution Chunk Loop
    for (int i = 0; i < totalBatches; i += batchMultiplier) {
        
        // --- ASSEMBLY TIMING START ---
        double tAssemblyStart = cuBQL::getCurrentTime();

        int activeBatchesInChunk = std::min(batchMultiplier, totalBatches - i);
        int totalChunkPrims = 0;

        // Assembly Step: Pack B primitives and build A search roots across the current batch chunk
        for (int b = 0; b < activeBatchesInChunk; ++b) {
            uint32_t bIdx;
            CUBQL_CUDA_CALL(Memcpy(&bIdx, ptr_outPairsB + i + b, sizeof(uint32_t), cudaMemcpyDeviceToHost));

            uint32_t startOffset;
            CUBQL_CUDA_CALL(Memcpy(&startOffset, ptr_outOffsetsB + bIdx, sizeof(uint32_t), cudaMemcpyDeviceToHost));

            uint32_t endOffset;
            if (bIdx + 1 < h_outMarkedCountB) {
                CUBQL_CUDA_CALL(Memcpy(&endOffset, ptr_outOffsetsB + bIdx + 1, sizeof(uint32_t), cudaMemcpyDeviceToHost));
            } else {
                endOffset = d_outPrimsFlatB.size();
            }

            int numPrims = (int)(endOffset - startOffset);
            if (numPrims <= 0) continue;

            // Gather B Primitives into the Sandbox via continuous Device-to-Device transfer
            CUBQL_CUDA_CALL(MemcpyAsync(
                d_BIter + totalChunkPrims,
                ptr_outPrimsFlatB + startOffset,
                numPrims * sizeof(uint32_t),
                cudaMemcpyDeviceToDevice
            ));

            // Fetch matched Node Root A values for the entire primitive range
            uint32_t pairAVal;
            CUBQL_CUDA_CALL(Memcpy(&pairAVal, ptr_outPairsA + i + b, sizeof(uint32_t), cudaMemcpyDeviceToHost));
            
            uint32_t startNodeIdxA_32;
            CUBQL_CUDA_CALL(Memcpy(&startNodeIdxA_32, ptr_markedNodeIndicesA + pairAVal, sizeof(uint32_t), cudaMemcpyDeviceToHost));
            uint64_t startNodeIdxA = (uint64_t)startNodeIdxA_32;

            // Broadcast Node Root A straight across the chunk-local window
            thrust::device_ptr<uint64_t> thrust_AIter(d_AIter);
            thrust::fill_n(
                thrust::device,
                thrust_AIter + totalChunkPrims, 
                numPrims,
                startNodeIdxA
            );

            totalChunkPrims += numPrims;
        }

        tracker.assemblyPhaseMs += (cuBQL::getCurrentTime() - tAssemblyStart) * 1000.0;
        // --- ASSEMBLY TIMING END ---

        // If the entire grouped chunk contains zero elements, pass immediately
        if (totalChunkPrims == 0) continue;

        // --- CUDA EXECUTION TIMING START ---
        CUBQL_CUDA_CALL(EventRecord(evComputeStart, 0));

        // --- PHASE 1: Execution Count (Clean 0 Offset Baseline) ---
        int kernelBlockSize = 128;
        int kernelGridSize  = (totalChunkPrims + kernelBlockSize - 1) / kernelBlockSize;

        countAABBOverlapsKernel_Indirected<<<kernelGridSize, kernelBlockSize>>>(
            d_pairCounts, bvhA, const_cast<cuBQL::Triangle*>(dMeshA), const_cast<cuBQL::Triangle*>(dMeshB), 
            d_BIter, 0, totalChunkPrims, d_AIter
        );

        // --- PHASE 2: Reusable Device Exclusive Prefix Scan ---
        thrust::device_ptr<int> dev_counts(d_pairCounts);
        thrust::device_ptr<int> dev_offsets(d_offsets);
        thrust::exclusive_scan(thrust::device, dev_counts, dev_counts + totalChunkPrims, dev_offsets);

        // --- PHASE 3: Read Combined Output Window Bounds ---
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

        // --- PHASE 4: Local Write out and Sync ---
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

        // Free allocation inside this step
        CUBQL_CUDA_CALL(Free(d_candidatePairs));
        // --- CUDA EXECUTION TIMING END ---
    }

    // --- CLEANUP TIMING START ---
    double tCleanupStart = cuBQL::getCurrentTime();

    CUBQL_CUDA_CALL(EventDestroy(evComputeStart));
    CUBQL_CUDA_CALL(EventDestroy(evComputeEnd));

    // Final Cleanup of reusable chunk sandbox allocations
    if (d_BIter)      CUBQL_CUDA_CALL(Free(d_BIter));
    if (d_AIter)      CUBQL_CUDA_CALL(Free(d_AIter));
    if (d_pairCounts) CUBQL_CUDA_CALL(Free(d_pairCounts));
    if (d_offsets)    CUBQL_CUDA_CALL(Free(d_offsets));

    tracker.cleanupTimeMs = (cuBQL::getCurrentTime() - tCleanupStart) * 1000.0;
    // --- CLEANUP TIMING END ---

    tracker.totalLoopTimeMs = (cuBQL::getCurrentTime() - tTotalStart) * 1000.0;
    std::cout << "--------------------------------------------------\n\n";

    return finalCandidatePairs;
}