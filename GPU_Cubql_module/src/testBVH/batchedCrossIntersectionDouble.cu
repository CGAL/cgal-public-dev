#include <thrust/device_ptr.h>
#include <thrust/scan.h>
#include <thrust/execution_policy.h>
#include <thrust/fill.h>
#include <thrust/copy.h>
#include <thrust/iterator/counting_iterator.h>
#include <iostream>
#include <vector>
#include <algorithm>

#include "../custom_pipeline/GPUPredicatesCheckDoubleAssisted.h"
#include "batchedCrossIntersectionCommon.h"
#include "TargetStatus.h"
#include "../src/CPU/YellowFilter.h"
#include "batchedCrossIntersectionDouble.h"

uint64_t executeBatchedCrossIntersectionLoopDouble(
    Mesh & meshAcpu, Mesh & meshBcpu,
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
    const cuBQL::box3f* d_boxesA,
    const cuBQL::box3f* d_boxesB,
    const double3* d_vertsA,
    const uint3* d_indicesA,
    const double3* d_vertsB,
    const uint3* d_indicesB,
    tbb::concurrent_vector<int2> & finalExactPairs,
    IntersectionTimeTracker& tracker,
    Point3 m_centerA, Point3 m_centerB, double3 m_rotA, double3 m_transA, double3 m_rotB, double3 m_transB,
    cudaStream_t stream
) {
    double tTotalStart = cuBQL::getCurrentTime();
    double tAllocStart = tTotalStart;

    uint64_t finalCandidatePairs = 0;

    // 1. Unpack Raw Resource Pointers
    const uint32_t* ptr_outPairsA     = thrust::raw_pointer_cast(d_outPairsA.data());
    const uint32_t* ptr_outPairsB     = thrust::raw_pointer_cast(d_outPairsB.data());
    const uint32_t* ptr_reverseMapB   = thrust::raw_pointer_cast(d_reverseMapB.data());
    const uint32_t* ptr_outOffsetsB   = thrust::raw_pointer_cast(d_outOffsetsB.data());
    const uint32_t* ptr_outPrimsFlatB = thrust::raw_pointer_cast(d_outPrimsFlatB.data());

    if (totalBatches > 0) {
        batchMultiplier = std::min(batchMultiplier, totalBatches);
    } else {
        batchMultiplier = 1; 
    }

    // 2. Compute Exact Maximum Chunk Space Bounds
    int maxChunkPrims = 0;
    if (totalBatches > 0) {
        int numChunks = (totalBatches + batchMultiplier - 1) / batchMultiplier;
        
        int* d_globalMax = nullptr;
        CUBQL_CUDA_CALL(MallocAsync(&d_globalMax, sizeof(int), stream));
        CUBQL_CUDA_CALL(MemsetAsync(d_globalMax, 0, sizeof(int), stream));

        int blockSize = 256; 
        int gridSize = numChunks;
        size_t sharedMemSize = blockSize * sizeof(int);

        computeMaxChunkSizeKernel<<<gridSize, blockSize, sharedMemSize, stream>>>(
            ptr_outPairsB, ptr_outOffsetsB, ptr_reverseMapB, 
            h_outMarkedCountB, d_outPrimsFlatB.size(), 
            totalBatches, batchMultiplier, numChunks, d_globalMax
        );

        CUBQL_CUDA_CALL(MemcpyAsync(&maxChunkPrims, d_globalMax, sizeof(int), cudaMemcpyDeviceToHost, stream));
        CUBQL_CUDA_CALL(StreamSynchronize(stream));
        CUBQL_CUDA_CALL(FreeAsync(d_globalMax, stream));
    }

    // 3. Preallocate Memory Sandboxes
    uint32_t* d_BIter = nullptr;
    uint64_t* d_AIter = nullptr;
    int* d_pairCounts = nullptr;
    int* d_offsets    = nullptr;
    int* d_chunkBatchSizes   = nullptr;
    int* d_chunkBatchOffsets = nullptr;

    if (maxChunkPrims > 0) {
        CUBQL_CUDA_CALL(MallocAsync(&d_BIter, maxChunkPrims * sizeof(uint32_t), stream));
        CUBQL_CUDA_CALL(MallocAsync(&d_AIter, maxChunkPrims * sizeof(uint64_t), stream));
        CUBQL_CUDA_CALL(MallocAsync(&d_pairCounts, maxChunkPrims * sizeof(int), stream));
        CUBQL_CUDA_CALL(MallocAsync(&d_offsets, maxChunkPrims * sizeof(int), stream));
        CUBQL_CUDA_CALL(MallocAsync(&d_chunkBatchSizes, batchMultiplier * sizeof(int), stream));
        CUBQL_CUDA_CALL(MallocAsync(&d_chunkBatchOffsets, batchMultiplier * sizeof(int), stream));
    }
    tracker.preallocateTimeMs = (cuBQL::getCurrentTime() - tAllocStart) * 1000.0;

    // Use multiple distinct events to track sub-phases purely asynchronously
    cudaEvent_t evAabbStart, evAabbEnd, evEvalStart, evEvalEnd;
    CUBQL_CUDA_CALL(EventCreate(&evAabbStart));
    CUBQL_CUDA_CALL(EventCreate(&evAabbEnd));
    CUBQL_CUDA_CALL(EventCreate(&evEvalStart));
    CUBQL_CUDA_CALL(EventCreate(&evEvalEnd));

    std::vector<int2> hGreenPairs;
    std::vector<int2> hYellowPairs;

    // 4. Coarse-Grained Execution Chunk Loop
    for (int i = 0; i < totalBatches; i += batchMultiplier) {
        tracker.numberOfBatchLoops++;
        
        double tAssemblyStart = cuBQL::getCurrentTime();
        int activeBatchesInChunk = std::min(batchMultiplier, totalBatches - i);

        auto localBatchIdxIterator = thrust::make_counting_iterator(0);
        auto batchSizeTransformIterator = thrust::make_transform_iterator(
            localBatchIdxIterator,
            BatchSizeInChunkFunctor(ptr_outPairsB, ptr_outOffsetsB, ptr_reverseMapB, h_outMarkedCountB, d_outPrimsFlatB.size(), i, activeBatchesInChunk)
        );

        thrust::device_ptr<int> dev_batchSizes(d_chunkBatchSizes);
        thrust::copy_n(thrust::cuda::par.on(stream), batchSizeTransformIterator, activeBatchesInChunk, dev_batchSizes);

        thrust::device_ptr<int> dev_batchOffsets(d_chunkBatchOffsets);
        thrust::exclusive_scan(thrust::cuda::par.on(stream), dev_batchSizes, dev_batchSizes + activeBatchesInChunk, dev_batchOffsets);

        int lastBatchSize = 0;
        int lastBatchOffset = 0;
        CUBQL_CUDA_CALL(MemcpyAsync(&lastBatchSize, d_chunkBatchSizes + activeBatchesInChunk - 1, sizeof(int), cudaMemcpyDeviceToHost, stream));
        CUBQL_CUDA_CALL(MemcpyAsync(&lastBatchOffset, d_chunkBatchOffsets + activeBatchesInChunk - 1, sizeof(int), cudaMemcpyDeviceToHost, stream));
        
        // Required sync to get sizes for the chunk
        CUBQL_CUDA_CALL(StreamSynchronize(stream));
        
        int totalChunkPrims = lastBatchSize + lastBatchOffset;
        if (totalChunkPrims == 0) continue;

        int assembleBlockSize = 64;
        int assembleGridSize  = (activeBatchesInChunk + assembleBlockSize - 1) / assembleBlockSize;

        assembleChunkBuffersByBatchKernel<<<assembleGridSize, assembleBlockSize, 0, stream>>>(
            d_BIter, d_AIter, ptr_outPairsA, ptr_outPairsB, ptr_reverseMapB,
            ptr_outOffsetsB, ptr_outPrimsFlatB, h_outMarkedCountB, d_outPrimsFlatB.size(),
            i, activeBatchesInChunk, d_chunkBatchOffsets
        );

        tracker.assemblyPhaseMs += (cuBQL::getCurrentTime() - tAssemblyStart) * 1000.0;

        int kernelBlockSize = 128;
        int kernelGridSize  = (totalChunkPrims + kernelBlockSize - 1) / kernelBlockSize;

        // Count Overlaps
        countAABBOverlapsKernel_Indirected_boxes<<<kernelGridSize, kernelBlockSize, 0, stream>>>(
            d_pairCounts, bvhA, d_boxesA, d_boxesB, 
            d_BIter, 0, totalChunkPrims, d_AIter
        );

        thrust::device_ptr<int> dev_counts(d_pairCounts);
        thrust::device_ptr<int> dev_offsets(d_offsets);
        thrust::exclusive_scan(thrust::cuda::par.on(stream), dev_counts, dev_counts + totalChunkPrims, dev_offsets);

        int lastCount = 0;
        int lastOffset = 0;
        CUBQL_CUDA_CALL(MemcpyAsync(&lastCount, d_pairCounts + totalChunkPrims - 1, sizeof(int), cudaMemcpyDeviceToHost, stream));
        CUBQL_CUDA_CALL(MemcpyAsync(&lastOffset, d_offsets + totalChunkPrims - 1, sizeof(int), cudaMemcpyDeviceToHost, stream));
        
        // Required sync to allocate the candidate array size
        CUBQL_CUDA_CALL(StreamSynchronize(stream));
        
        int totalChunkPairs = lastCount + lastOffset;
        if (totalChunkPairs == 0) continue;
        
        finalCandidatePairs += totalChunkPairs;

        // ----------------------------------------------------------------
        // AABB CANDIDATE GENERATION
        // ----------------------------------------------------------------
        int2* d_candidatePairs = nullptr;
        CUBQL_CUDA_CALL(MallocAsync(&d_candidatePairs, totalChunkPairs * sizeof(int2), stream));

        CUBQL_CUDA_CALL(EventRecord(evAabbStart, stream));
        fillAABBOverlapsKernel_Indirected_Boxes<<<kernelGridSize, kernelBlockSize, 0, stream>>>(
            d_candidatePairs, d_offsets, bvhA, d_boxesA, d_boxesB, 
            d_BIter, 0, totalChunkPrims, d_AIter
        );
        CUBQL_CUDA_CALL(EventRecord(evAabbEnd, stream)); // NO SYNC HERE

        // ----------------------------------------------------------------
        // DIRECT DOUBLE-PRECISION GPU PREDICATES EVALUATION
        // ----------------------------------------------------------------
        int* d_pairStatuses = nullptr;
        CUBQL_CUDA_CALL(MallocAsync(&d_pairStatuses, totalChunkPairs * sizeof(int), stream));

        CUBQL_CUDA_CALL(EventRecord(evEvalStart, stream));
        evaluateAndCompactPairsDoubleAssisted(
            d_candidatePairs,
            d_pairStatuses,
            d_vertsA,
            d_indicesA,
            d_vertsB,
            d_indicesB,
            totalChunkPairs,
            stream
        );
        CUBQL_CUDA_CALL(EventRecord(evEvalEnd, stream)); // NO SYNC HERE

        // ----------------------------------------------------------------
        // PAIR COMPACTION PASS (GREEN -> FINAL, YELLOW -> CPU CGAL)
        // ----------------------------------------------------------------
        double tEvalStartTwo = cuBQL::getCurrentTime();

        thrust::device_ptr<int2> dev_evaluated(d_candidatePairs);
        thrust::device_ptr<int> dev_statuses(d_pairStatuses);

        // MallocAsync avoids device-wide driver stalls
        int2* d_green_raw = nullptr;
        int2* d_yellow_raw = nullptr;
        CUBQL_CUDA_CALL(MallocAsync(&d_green_raw, totalChunkPairs * sizeof(int2), stream));
        CUBQL_CUDA_CALL(MallocAsync(&d_yellow_raw, totalChunkPairs * sizeof(int2), stream));

        thrust::device_ptr<int2> dev_green_out(d_green_raw);
        thrust::device_ptr<int2> dev_yellow_out(d_yellow_raw);

        auto green_end  = thrust::copy_if(thrust::cuda::par.on(stream), dev_evaluated, dev_evaluated + totalChunkPairs, dev_statuses, dev_green_out,  IsTargetPairStatus{(int)PAIR_GREEN});
        auto yellow_end = thrust::copy_if(thrust::cuda::par.on(stream), dev_evaluated, dev_evaluated + totalChunkPairs, dev_statuses, dev_yellow_out, IsTargetPairStatus{(int)PAIR_YELLOW});

        int totalGreen  = green_end - dev_green_out;
        int totalYellow = yellow_end - dev_yellow_out;

        if (totalGreen > 0) {
            size_t oldSize = hGreenPairs.size();
            hGreenPairs.resize(oldSize + totalGreen);
            CUBQL_CUDA_CALL(MemcpyAsync(hGreenPairs.data() + oldSize, d_green_raw, totalGreen * sizeof(int2), cudaMemcpyDeviceToHost, stream));
        }

        if (totalYellow > 0) {
            size_t oldSize = hYellowPairs.size();
            hYellowPairs.resize(oldSize + totalYellow);
            CUBQL_CUDA_CALL(MemcpyAsync(hYellowPairs.data() + oldSize, d_yellow_raw, totalYellow * sizeof(int2), cudaMemcpyDeviceToHost, stream));
        }

        // ----------------------------------------------------------------
        // FINALIZE CHUNK, SYNC, AND RECORD EXACT TIMERS
        // ----------------------------------------------------------------
        CUBQL_CUDA_CALL(StreamSynchronize(stream)); // Only ONE sync per chunk for the entire execution

        float aabbMs = 0.0f;
        float evalMs = 0.0f;
        CUBQL_CUDA_CALL(EventElapsedTime(&aabbMs, evAabbStart, evAabbEnd));
        CUBQL_CUDA_CALL(EventElapsedTime(&evalMs, evEvalStart, evEvalEnd));
        
        tracker.executionPhaseMs += aabbMs;
        tracker.fineEvaluationPhaseMs += evalMs;

        CUBQL_CUDA_CALL(FreeAsync(d_candidatePairs, stream));
        CUBQL_CUDA_CALL(FreeAsync(d_pairStatuses, stream));
        CUBQL_CUDA_CALL(FreeAsync(d_green_raw, stream));
        CUBQL_CUDA_CALL(FreeAsync(d_yellow_raw, stream));

        tracker.DownloadAndClean += (cuBQL::getCurrentTime() - tEvalStartTwo) * 1000.0;
    }

    // --------------------------------------------------------------------
    // CPU EXACT PREDICATES PASS (YELLOW FILTERING VIA TBB / CGAL)
    // --------------------------------------------------------------------
    double tCPUPredicates = cuBQL::getCurrentTime();

    tracker.confirmedGreenPairs  = hGreenPairs.size();
    tracker.confirmedYellowPairs = hYellowPairs.size();

    finalExactPairs = tbb::concurrent_vector<int2>(hGreenPairs.begin(), hGreenPairs.end());
    
    // CGAL TBB Filter executes directly on all Yellow pairs
    filterYellowPairsTBB(meshAcpu, meshBcpu, hYellowPairs.data(), hYellowPairs.size(), finalExactPairs, m_centerA, m_rotA, m_transA, m_centerB, m_rotB, m_transB);

    tracker.CPUPredicates = (cuBQL::getCurrentTime() - tCPUPredicates) * 1000.0;

    double tCleanupStart = cuBQL::getCurrentTime();

    CUBQL_CUDA_CALL(EventDestroy(evAabbStart));
    CUBQL_CUDA_CALL(EventDestroy(evAabbEnd));
    CUBQL_CUDA_CALL(EventDestroy(evEvalStart));
    CUBQL_CUDA_CALL(EventDestroy(evEvalEnd));

    if (d_BIter)              CUBQL_CUDA_CALL(FreeAsync(d_BIter, stream));
    if (d_AIter)              CUBQL_CUDA_CALL(FreeAsync(d_AIter, stream));
    if (d_pairCounts)         CUBQL_CUDA_CALL(FreeAsync(d_pairCounts, stream));
    if (d_offsets)            CUBQL_CUDA_CALL(FreeAsync(d_offsets, stream));
    if (d_chunkBatchSizes)    CUBQL_CUDA_CALL(FreeAsync(d_chunkBatchSizes, stream));
    if (d_chunkBatchOffsets)  CUBQL_CUDA_CALL(FreeAsync(d_chunkBatchOffsets, stream));

    CUBQL_CUDA_CALL(StreamSynchronize(stream));

    tracker.cleanupTimeMs = (cuBQL::getCurrentTime() - tCleanupStart) * 1000.0;
    tracker.totalLoopTimeMs = (cuBQL::getCurrentTime() - tTotalStart) * 1000.0;

    std::cout << "--------------------------------------------------\n\n";
    return finalCandidatePairs;
}