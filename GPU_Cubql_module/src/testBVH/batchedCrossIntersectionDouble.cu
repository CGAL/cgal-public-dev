#include <thrust/device_ptr.h>
#include <thrust/scan.h>
#include <thrust/execution_policy.h>
#include <thrust/fill.h>
#include <thrust/copy.h>
#include <thrust/iterator/counting_iterator.h>
#include <iostream>
#include <vector>
#include <algorithm>

#include "batchedCrossIntersectionDouble.h"
#include "batchedCrossIntersectionCommon.h"
#include "../custom_pipeline/GPUPredicatesCheckV2.h"
#include "TargetStatus.h"
#include "../src/CPU/YellowFilter.h"

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
    const cuBQL::Triangle* dMeshA,
    const cuBQL::Triangle* dMeshB,
    const float2 *triAMetrics,
    const float2 *triBMetrics,
    tbb::concurrent_vector<int2> & finalExactPairs,
    IntersectionTimeTracker& tracker,
    cudaStream_t stream,
    const TriangleDouble* dMeshADouble,
    const TriangleDouble* dMeshBDouble
) {
    double tTotalStart = cuBQL::getCurrentTime();
    double tAllocStart = tTotalStart;

    uint64_t finalCandidatePairs = 0;
    bool useDoublePass = (dMeshADouble != nullptr && dMeshBDouble != nullptr);

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

    // 2. Compute Exact Maximum Chunk Space Bounds (Parallel Pass)
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

    cudaEvent_t evComputeStart, evComputeEnd;
    CUBQL_CUDA_CALL(EventCreate(&evComputeStart));
    CUBQL_CUDA_CALL(EventCreate(&evComputeEnd));

    std::vector<int2> hGreenPairs;
    std::vector<int2> hOrangePairs; // Pairs that remained unresolved after GPU passes
    int totalFloatYellowCount = 0;

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

        CUBQL_CUDA_CALL(EventRecord(evComputeStart, stream));

        int kernelBlockSize = 128;
        int kernelGridSize  = (totalChunkPrims + kernelBlockSize - 1) / kernelBlockSize;

        countAABBOverlapsKernel_Indirected<<<kernelGridSize, kernelBlockSize, 0, stream>>>(
            d_pairCounts, bvhA, const_cast<cuBQL::Triangle*>(dMeshA), const_cast<cuBQL::Triangle*>(dMeshB), 
            d_BIter, 0, totalChunkPrims, d_AIter
        );

        thrust::device_ptr<int> dev_counts(d_pairCounts);
        thrust::device_ptr<int> dev_offsets(d_offsets);
        thrust::exclusive_scan(thrust::cuda::par.on(stream), dev_counts, dev_counts + totalChunkPrims, dev_offsets);

        int lastCount = 0;
        int lastOffset = 0;
        CUBQL_CUDA_CALL(MemcpyAsync(&lastCount, d_pairCounts + totalChunkPrims - 1, sizeof(int), cudaMemcpyDeviceToHost, stream));
        CUBQL_CUDA_CALL(MemcpyAsync(&lastOffset, d_offsets + totalChunkPrims - 1, sizeof(int), cudaMemcpyDeviceToHost, stream));
        CUBQL_CUDA_CALL(StreamSynchronize(stream));
        
        int totalChunkPairs = lastCount + lastOffset;

        if (totalChunkPairs == 0) {
            CUBQL_CUDA_CALL(EventRecord(evComputeEnd, stream));
            CUBQL_CUDA_CALL(EventSynchronize(evComputeEnd));
            float elapsedMs = 0.0f;
            CUBQL_CUDA_CALL(EventElapsedTime(&elapsedMs, evComputeStart, evComputeEnd));
            tracker.executionPhaseMs += elapsedMs;
            continue;
        }
        
        finalCandidatePairs += totalChunkPairs;

        int2* d_candidatePairs = nullptr;
        CUBQL_CUDA_CALL(MallocAsync(&d_candidatePairs, totalChunkPairs * sizeof(int2), stream));

        fillAABBOverlapsKernel_Indirected<<<kernelGridSize, kernelBlockSize, 0, stream>>>(
            d_candidatePairs, d_offsets, bvhA, const_cast<cuBQL::Triangle*>(dMeshA), const_cast<cuBQL::Triangle*>(dMeshB), 
            d_BIter, 0, totalChunkPrims, d_AIter
        );

        CUBQL_CUDA_CALL(EventRecord(evComputeEnd, stream));
        CUBQL_CUDA_CALL(EventSynchronize(evComputeEnd));
        
        float chunkComputeMs = 0.0f;
        CUBQL_CUDA_CALL(EventElapsedTime(&chunkComputeMs, evComputeStart, evComputeEnd));
        tracker.executionPhaseMs += chunkComputeMs;

        // ----------------------------------------------------------------
        // STAGE A: GPU FLOAT PREDICATES EVALUATION
        // ----------------------------------------------------------------
        double tEvalStart = cuBQL::getCurrentTime();

        int* d_pairStatuses = nullptr;
        CUBQL_CUDA_CALL(MallocAsync(&d_pairStatuses, totalChunkPairs * sizeof(int), stream));

        evaluateAndCompactPairsV2(
            d_candidatePairs, 
            d_pairStatuses, 
            dMeshA, 
            dMeshB, 
            triAMetrics,
            triBMetrics,
            totalChunkPairs, 
            stream
        );

        tracker.fineEvaluationPhaseMs += (cuBQL::getCurrentTime() - tEvalStart) * 1000.0;

        double tEvalStartTwo = cuBQL::getCurrentTime();

        thrust::device_ptr<int2> dev_evaluated(d_candidatePairs);
        thrust::device_ptr<int> dev_statuses(d_pairStatuses);

        thrust::device_vector<int2> dev_green(totalChunkPairs);
        thrust::device_vector<int2> dev_yellow(totalChunkPairs);

        thrust::device_ptr<int2> dev_green_out(thrust::raw_pointer_cast(dev_green.data()));
        thrust::device_ptr<int2> dev_yellow_out(thrust::raw_pointer_cast(dev_yellow.data()));

        auto green_end = thrust::copy_if(thrust::cuda::par.on(stream), dev_evaluated, dev_evaluated + totalChunkPairs, dev_statuses, dev_green_out, IsTargetPairStatus{(int)PAIR_GREEN});
        auto yellow_end = thrust::copy_if(thrust::cuda::par.on(stream), dev_evaluated, dev_evaluated + totalChunkPairs, dev_statuses, dev_yellow_out, IsTargetPairStatus{(int)PAIR_YELLOW});

        int totalGreen = green_end - dev_green_out;
        int totalYellow = yellow_end - dev_yellow_out;
        totalFloatYellowCount += totalYellow;

        if (totalGreen > 0) {
            size_t oldSize = hGreenPairs.size();
            hGreenPairs.resize(oldSize + totalGreen);
            CUBQL_CUDA_CALL(MemcpyAsync(hGreenPairs.data() + oldSize, thrust::raw_pointer_cast(dev_green.data()), totalGreen * sizeof(int2), cudaMemcpyDeviceToHost, stream));
        }

        // ----------------------------------------------------------------
        // STAGE B: GPU DOUBLE PREDICATES EVALUATION (OPTIONAL PASS)
        // ----------------------------------------------------------------
        if (useDoublePass && totalYellow > 0) {
            double tDoubleStart = cuBQL::getCurrentTime();

            int* d_yellowStatuses = nullptr;
            CUBQL_CUDA_CALL(MallocAsync(&d_yellowStatuses, totalYellow * sizeof(int), stream));

            // Launch 64-bit interval evaluation on yellow candidates
            evaluateAndCompactPairsDouble(
                thrust::raw_pointer_cast(dev_yellow.data()),
                d_yellowStatuses,
                dMeshADouble,
                dMeshBDouble,
                totalYellow,
                stream
            );

            thrust::device_ptr<int2> dev_yellow_in(thrust::raw_pointer_cast(dev_yellow.data()));
            thrust::device_ptr<int> dev_dbl_statuses(d_yellowStatuses);

            thrust::device_vector<int2> dev_double_green(totalYellow);
            thrust::device_vector<int2> dev_orange(totalYellow);

            thrust::device_ptr<int2> dev_dbl_green_out(thrust::raw_pointer_cast(dev_double_green.data()));
            thrust::device_ptr<int2> dev_orange_out(thrust::raw_pointer_cast(dev_orange.data()));

            auto dbl_green_end = thrust::copy_if(thrust::cuda::par.on(stream), dev_yellow_in, dev_yellow_in + totalYellow, dev_dbl_statuses, dev_dbl_green_out, IsTargetPairStatus{(int)PAIR_GREEN});
            auto orange_end    = thrust::copy_if(thrust::cuda::par.on(stream), dev_yellow_in, dev_yellow_in + totalYellow, dev_dbl_statuses, dev_orange_out, IsTargetPairStatus{(int)PAIR_YELLOW});

            int totalDblGreen = dbl_green_end - dev_dbl_green_out;
            int totalOrange   = orange_end - dev_orange_out;

            if (totalDblGreen > 0) {
                size_t oldSize = hGreenPairs.size();
                hGreenPairs.resize(oldSize + totalDblGreen);
                CUBQL_CUDA_CALL(MemcpyAsync(hGreenPairs.data() + oldSize, thrust::raw_pointer_cast(dev_double_green.data()), totalDblGreen * sizeof(int2), cudaMemcpyDeviceToHost, stream));
            }
            if (totalOrange > 0) {
                size_t oldSize = hOrangePairs.size();
                hOrangePairs.resize(oldSize + totalOrange);
                CUBQL_CUDA_CALL(MemcpyAsync(hOrangePairs.data() + oldSize, thrust::raw_pointer_cast(dev_orange.data()), totalOrange * sizeof(int2), cudaMemcpyDeviceToHost, stream));
            }

            CUBQL_CUDA_CALL(FreeAsync(d_yellowStatuses, stream));
            tracker.gpuDoublePredicatesMs += (cuBQL::getCurrentTime() - tDoubleStart) * 1000.0;
        } 
        else if (!useDoublePass && totalYellow > 0) {
            // No double array passed: Float yellow pairs directly become orange pairs for CPU fallback
            size_t oldSize = hOrangePairs.size();
            hOrangePairs.resize(oldSize + totalYellow);
            CUBQL_CUDA_CALL(MemcpyAsync(hOrangePairs.data() + oldSize, thrust::raw_pointer_cast(dev_yellow.data()), totalYellow * sizeof(int2), cudaMemcpyDeviceToHost, stream));
        }

        CUBQL_CUDA_CALL(StreamSynchronize(stream));

        CUBQL_CUDA_CALL(FreeAsync(d_candidatePairs, stream));
        CUBQL_CUDA_CALL(FreeAsync(d_pairStatuses, stream));

        tracker.DownloadAndClean += (cuBQL::getCurrentTime() - tEvalStartTwo) * 1000.0;
    }

    // --------------------------------------------------------------------
    // STAGE C: CPU EXACT PREDICATES PASS (ORANGE LIST FILTERING)
    // --------------------------------------------------------------------
    double tCPUPredicates = cuBQL::getCurrentTime();

    tracker.confirmedGreenPairs = hGreenPairs.size();
    tracker.confirmedYellowPairs = totalFloatYellowCount;
    tracker.confirmedOrangePairs = hOrangePairs.size();

    finalExactPairs = tbb::concurrent_vector<int2>(hGreenPairs.begin(), hGreenPairs.end());
    
    // CGAL TBB Filter executes strictly on remaining Orange pairs
    filterYellowPairsTBB(meshAcpu, meshBcpu, hOrangePairs.data(), hOrangePairs.size(), finalExactPairs);

    tracker.CPUPredicates = (cuBQL::getCurrentTime() - tCPUPredicates) * 1000.0;

    double tCleanupStart = cuBQL::getCurrentTime();

    CUBQL_CUDA_CALL(EventDestroy(evComputeStart));
    CUBQL_CUDA_CALL(EventDestroy(evComputeEnd));

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