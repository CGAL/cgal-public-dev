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
        int i = chunkIdx * batchMultiplier; 
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

// Upgraded Parallel GPU Assembly Kernel
__global__ void assembleChunkBuffersKernel(
    uint32_t* d_BIter,
    uint64_t* d_AIter,
    const uint32_t* outPairsA,
    const uint32_t* outPairsB,
    const uint32_t* markedNodeIndicesA,
    const uint32_t* outOffsetsB,
    const uint32_t* outPrimsFlatB,
    uint32_t outMarkedCountB,
    uint32_t totalPrimsB,
    int chunkStartBatchIdx,
    int activeBatchesInChunk,
    const int* d_chunkBatchOffsets,  // Maps where each batch's primitive block starts inside the chunk sandbox
    int totalChunkPrims
) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= totalChunkPrims) return;

    // Binary search to find which batch this specific thread/primitive slot belongs to
    int low = 0;
    int high = activeBatchesInChunk - 1;
    int b = 0;
    while (low <= high) {
        int mid = (low + high) / 2;
        if (d_chunkBatchOffsets[mid] <= tid) {
            b = mid;
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }

    // Determine target offsets and identity details for the mapped batch completely on GPU
    int batchLocalPrimIdx = tid - d_chunkBatchOffsets[b];
    int currentBatchArrayIdx = chunkStartBatchIdx + b;

    uint32_t bIdx = outPairsB[currentBatchArrayIdx];
    uint32_t startOffsetB = outOffsetsB[bIdx];

    // Read flat primitive tracking identity and assign to target iteration sandbox buffer
    d_BIter[tid] = outPrimsFlatB[startOffsetB + batchLocalPrimIdx];

    // Broadcast A-tree search root values out to matched thread indices
    uint32_t pairAVal = outPairsA[currentBatchArrayIdx];
    d_AIter[tid] = (uint64_t)(markedNodeIndicesA[pairAVal]);
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
    int* d_chunkBatchOffsets = nullptr; // GPU offset lookup layout map for assembly execution

    if (maxChunkPrims > 0) {
        CUBQL_CUDA_CALL(Malloc(&d_BIter, maxChunkPrims * sizeof(uint32_t)));
        CUBQL_CUDA_CALL(Malloc(&d_AIter, maxChunkPrims * sizeof(uint64_t)));
        CUBQL_CUDA_CALL(Malloc(&d_pairCounts, maxChunkPrims * sizeof(int)));
        CUBQL_CUDA_CALL(Malloc(&d_offsets, maxChunkPrims * sizeof(int)));
        CUBQL_CUDA_CALL(Malloc(&d_chunkBatchOffsets, batchMultiplier * sizeof(int)));
    }
    tracker.preallocateTimeMs = (cuBQL::getCurrentTime() - tAllocStart) * 1000.0;

    // Allocate small staging pinned host structures to parse local chunk batch metadata tracking safely
    std::vector<uint32_t> h_localPairsB(batchMultiplier);
    std::vector<int> h_localBatchOffsets(batchMultiplier);
    std::vector<uint32_t> h_localOffsetsB; 

    cudaEvent_t evComputeStart, evComputeEnd;
    CUBQL_CUDA_CALL(EventCreate(&evComputeStart));
    CUBQL_CUDA_CALL(EventCreate(&evComputeEnd));

    // 4. Coarse-Grained Execution Chunk Loop
    for (int i = 0; i < totalBatches; i += batchMultiplier) {
        
        // --- ASSEMBLY TIMING START ---
        double tAssemblyStart = cuBQL::getCurrentTime();

        int activeBatchesInChunk = std::min(batchMultiplier, totalBatches - i);
        
        // Block-transfer active batch identifiers straight from GPU to host tracking vectors in one swoop
        CUBQL_CUDA_CALL(Memcpy(
            h_localPairsB.data(), 
            ptr_outPairsB + i, 
            activeBatchesInChunk * sizeof(uint32_t), 
            cudaMemcpyDeviceToHost
        ));

        // Fetch corresponding bounding index values via bulk transfer
        uint32_t maxBIdxUsed = 0;
        for(int b = 0; b < activeBatchesInChunk; ++b) {
            maxBIdxUsed = std::max(maxBIdxUsed, h_localPairsB[b]);
        }
        
        h_localOffsetsB.resize(maxBIdxUsed + 2);
        CUBQL_CUDA_CALL(Memcpy(
            h_localOffsetsB.data(), 
            ptr_outOffsetsB, 
            (maxBIdxUsed + 2) * sizeof(uint32_t), 
            cudaMemcpyDeviceToHost
        ));

        // Rapid host sequence setup to find global limits and thread map configurations
        int totalChunkPrims = 0;
        for (int b = 0; b < activeBatchesInChunk; ++b) {
            h_localBatchOffsets[b] = totalChunkPrims;
            uint32_t bIdx = h_localPairsB[b];
            uint32_t startOffset = h_localOffsetsB[bIdx];
            uint32_t endOffset = (bIdx + 1 < h_outMarkedCountB) ? h_localOffsetsB[bIdx + 1] : d_outPrimsFlatB.size();
            
            totalChunkPrims += std::max(0, (int)(endOffset - startOffset));
        }

        if (totalChunkPrims == 0) continue;

        // Move the computed offset layouts back to the GPU once to configure the assembly execution map
        CUBQL_CUDA_CALL(MemcpyAsync(
            d_chunkBatchOffsets, 
            h_localBatchOffsets.data(), 
            activeBatchesInChunk * sizeof(int), 
            cudaMemcpyHostToDevice
        ));

        // Launch the parallel initialization kernel to layout data structures simultaneously on the hardware
        int assembleBlockSize = 256;
        int assembleGridSize  = (totalChunkPrims + assembleBlockSize - 1) / assembleBlockSize;

        assembleChunkBuffersKernel<<<assembleGridSize, assembleBlockSize>>>(
            d_BIter, d_AIter, ptr_outPairsA, ptr_outPairsB, ptr_markedNodeIndicesA,
            ptr_outOffsetsB, ptr_outPrimsFlatB, h_outMarkedCountB, d_outPrimsFlatB.size(),
            i, activeBatchesInChunk, d_chunkBatchOffsets, totalChunkPrims
        );

        tracker.assemblyPhaseMs += (cuBQL::getCurrentTime() - tAssemblyStart) * 1000.0;
        // --- ASSEMBLY TIMING END ---

        // --- CUDA EXECUTION TIMING START ---
        CUBQL_CUDA_CALL(EventRecord(evComputeStart, 0));

        // --- PHASE 1: Execution Count ---
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

        CUBQL_CUDA_CALL(Free(d_candidatePairs));
        // --- CUDA EXECUTION TIMING END ---
    }

    // --- CLEANUP TIMING START ---
    double tCleanupStart = cuBQL::getCurrentTime();

    CUBQL_CUDA_CALL(EventDestroy(evComputeStart));
    CUBQL_CUDA_CALL(EventDestroy(evComputeEnd));

    if (d_BIter)              CUBQL_CUDA_CALL(Free(d_BIter));
    if (d_AIter)              CUBQL_CUDA_CALL(Free(d_AIter));
    if (d_pairCounts)         CUBQL_CUDA_CALL(Free(d_pairCounts));
    if (d_offsets)            CUBQL_CUDA_CALL(Free(d_offsets));
    if (d_chunkBatchOffsets)  CUBQL_CUDA_CALL(Free(d_chunkBatchOffsets));

    tracker.cleanupTimeMs = (cuBQL::getCurrentTime() - tCleanupStart) * 1000.0;
    // --- CLEANUP TIMING END ---

    tracker.totalLoopTimeMs = (cuBQL::getCurrentTime() - tTotalStart) * 1000.0;
    std::cout << "--------------------------------------------------\n\n";

    return finalCandidatePairs;
}