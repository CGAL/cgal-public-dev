#include <thrust/device_ptr.h>
#include <thrust/scan.h>
#include <thrust/execution_policy.h>
#include <thrust/fill.h>
#include <thrust/copy.h>
#include <thrust/iterator/counting_iterator.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>

#include "../custom_pipeline/GPUPredicatesCheckDoubleAssisted.h"
#include "../custom_pipeline/GPUPredicatesCheckDouble.h" 
#include "../custom_pipeline/GPUPredicatesCheckShewchukFloat.h" 
#include "../custom_pipeline/GPUPredicatesCheckBigInteger.h" 
#include "../custom_pipeline/GPUPredicatesTwoLap.h"
#include "batchedCrossIntersectionCommon.h"
#include "TargetStatus.h"
#include "../src/CPU/YellowFilter.h"
#include "batchedCrossIntersectionDouble.h"


uint64_t executeBatchedCrossIntersectionLoopDouble(Mesh& meshAcpu,
                                                   Mesh& meshBcpu,
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
                                                   int2*& outFinalExactPairs, 
                                                   size_t& outFinalCount,     
                                                   IntersectionTimeTracker& tracker,
                                                   Point3 m_centerA,
                                                   Point3 m_centerB,
                                                   double3 m_rotA,
                                                   double3 m_transA,
                                                   double3 m_rotB,
                                                   double3 m_transB,
                                                   int exactPredicateComputeMode,
                                                   cudaStream_t stream) {
  double tTotalStart = cuBQL::getCurrentTime();
  double tAllocStart = tTotalStart;

  uint64_t finalCandidatePairs = 0;

  // 1. Unpack Raw Resource Pointers
  const uint32_t* ptr_outPairsA = thrust::raw_pointer_cast(d_outPairsA.data());
  const uint32_t* ptr_outPairsB = thrust::raw_pointer_cast(d_outPairsB.data());
  const uint32_t* ptr_reverseMapB = thrust::raw_pointer_cast(d_reverseMapB.data());
  const uint32_t* ptr_outOffsetsB = thrust::raw_pointer_cast(d_outOffsetsB.data());
  const uint32_t* ptr_outPrimsFlatB = thrust::raw_pointer_cast(d_outPrimsFlatB.data());

  if(totalBatches > 0) {
    batchMultiplier = std::min(batchMultiplier, totalBatches);
  } else {
    batchMultiplier = 1;
  }

  // 2. Compute Exact Maximum Chunk Space Bounds
  int maxChunkPrims = 0;
  if(totalBatches > 0) {
    int numChunks = (totalBatches + batchMultiplier - 1) / batchMultiplier;

    int* d_globalMax = nullptr;
    CUBQL_CUDA_CALL(MallocAsync(&d_globalMax, sizeof(int), stream));
    CUBQL_CUDA_CALL(MemsetAsync(d_globalMax, 0, sizeof(int), stream));

    int blockSize = 256;
    int gridSize = numChunks;
    size_t sharedMemSize = blockSize * sizeof(int);

    computeMaxChunkSizeKernel<<<gridSize, blockSize, sharedMemSize, stream>>>(
        ptr_outPairsB, ptr_outOffsetsB, ptr_reverseMapB, h_outMarkedCountB, d_outPrimsFlatB.size(), totalBatches,
        batchMultiplier, numChunks, d_globalMax);

    CUBQL_CUDA_CALL(MemcpyAsync(&maxChunkPrims, d_globalMax, sizeof(int), cudaMemcpyDeviceToHost, stream));
    CUBQL_CUDA_CALL(StreamSynchronize(stream));
    CUBQL_CUDA_CALL(FreeAsync(d_globalMax, stream));
  }

  // 3. Preallocate Memory Sandboxes
  uint32_t* d_BIter = nullptr;
  uint64_t* d_AIter = nullptr;
  int* d_pairCounts = nullptr;
  int* d_offsets = nullptr;
  int* d_chunkBatchSizes = nullptr;
  int* d_chunkBatchOffsets = nullptr;

  if(maxChunkPrims > 0) {
    CUBQL_CUDA_CALL(MallocAsync(&d_BIter, maxChunkPrims * sizeof(uint32_t), stream));
    CUBQL_CUDA_CALL(MallocAsync(&d_AIter, maxChunkPrims * sizeof(uint64_t), stream));
    CUBQL_CUDA_CALL(MallocAsync(&d_pairCounts, maxChunkPrims * sizeof(int), stream));
    CUBQL_CUDA_CALL(MallocAsync(&d_offsets, maxChunkPrims * sizeof(int), stream));
    CUBQL_CUDA_CALL(MallocAsync(&d_chunkBatchSizes, batchMultiplier * sizeof(int), stream));
    CUBQL_CUDA_CALL(MallocAsync(&d_chunkBatchOffsets, batchMultiplier * sizeof(int), stream));
  }
  tracker.preallocateTimeMs = (cuBQL::getCurrentTime() - tAllocStart) * 1000.0;

cudaEvent_t evAabbStart, evAabbEnd, evEvalStart, evEvalEnd, evDoubleStart, evDoubleEnd;
cudaEvent_t evCompact1Start, evCompact1End, evCompact2Start, evCompact2End;

CUBQL_CUDA_CALL(EventCreate(&evAabbStart));
CUBQL_CUDA_CALL(EventCreate(&evAabbEnd));
CUBQL_CUDA_CALL(EventCreate(&evEvalStart));
CUBQL_CUDA_CALL(EventCreate(&evEvalEnd));
CUBQL_CUDA_CALL(EventCreate(&evDoubleStart)); 
CUBQL_CUDA_CALL(EventCreate(&evDoubleEnd));   
CUBQL_CUDA_CALL(EventCreate(&evCompact1Start));
CUBQL_CUDA_CALL(EventCreate(&evCompact1End));
CUBQL_CUDA_CALL(EventCreate(&evCompact2Start));
CUBQL_CUDA_CALL(EventCreate(&evCompact2End));
  // Dynamic Raw Host Buffer Management (0.00 ms zero-init cost)
  int2* h_finalBuffer = nullptr;
  size_t hostCapacity = 0;
  size_t hostCount = 0;

  auto ensureHostCapacity = [&](size_t neededCapacity) {
    if(neededCapacity > hostCapacity) {
      size_t newCap = std::max(neededCapacity, hostCapacity == 0 ? (size_t)1048576 : hostCapacity * 2);
      int2* newBuf = static_cast<int2*>(std::realloc(h_finalBuffer, newCap * sizeof(int2)));
      if(!newBuf) {
        std::cerr << "Host Memory Allocation Failed!\n";
        std::exit(EXIT_FAILURE);
      }
      h_finalBuffer = newBuf;
      hostCapacity = newCap;
    }
  };

  // 4. Coarse-Grained Execution Chunk Loop
  for(int i = 0; i < totalBatches; i += batchMultiplier) {
    tracker.numberOfBatchLoops++;

    double tAssemblyStart = cuBQL::getCurrentTime();
    int activeBatchesInChunk = std::min(batchMultiplier, totalBatches - i);

    auto localBatchIdxIterator = thrust::make_counting_iterator(0);
    auto batchSizeTransformIterator = thrust::make_transform_iterator(
        localBatchIdxIterator,
        BatchSizeInChunkFunctor(ptr_outPairsB, ptr_outOffsetsB, ptr_reverseMapB, h_outMarkedCountB,
                                d_outPrimsFlatB.size(), i, activeBatchesInChunk));

    thrust::device_ptr<int> dev_batchSizes(d_chunkBatchSizes);
    thrust::copy_n(thrust::cuda::par.on(stream), batchSizeTransformIterator, activeBatchesInChunk, dev_batchSizes);

    thrust::device_ptr<int> dev_batchOffsets(d_chunkBatchOffsets);
    thrust::exclusive_scan(thrust::cuda::par.on(stream), dev_batchSizes, dev_batchSizes + activeBatchesInChunk,
                           dev_batchOffsets);

    int lastBatchSize = 0;
    int lastBatchOffset = 0;
    CUBQL_CUDA_CALL(MemcpyAsync(&lastBatchSize, d_chunkBatchSizes + activeBatchesInChunk - 1, sizeof(int),
                                cudaMemcpyDeviceToHost, stream));
    CUBQL_CUDA_CALL(MemcpyAsync(&lastBatchOffset, d_chunkBatchOffsets + activeBatchesInChunk - 1, sizeof(int),
                                cudaMemcpyDeviceToHost, stream));

    CUBQL_CUDA_CALL(StreamSynchronize(stream));

    int totalChunkPrims = lastBatchSize + lastBatchOffset;
    if(totalChunkPrims == 0)
      continue;

    int assembleBlockSize = 64;
    int assembleGridSize = (activeBatchesInChunk + assembleBlockSize - 1) / assembleBlockSize;

    assembleChunkBuffersByBatchKernel<<<assembleGridSize, assembleBlockSize, 0, stream>>>(
        d_BIter, d_AIter, ptr_outPairsA, ptr_outPairsB, ptr_reverseMapB, ptr_outOffsetsB, ptr_outPrimsFlatB,
        h_outMarkedCountB, d_outPrimsFlatB.size(), i, activeBatchesInChunk, d_chunkBatchOffsets);

    tracker.assemblyPhaseMs += (cuBQL::getCurrentTime() - tAssemblyStart) * 1000.0;

    int kernelBlockSize = 128;
    int kernelGridSize = (totalChunkPrims + kernelBlockSize - 1) / kernelBlockSize;

    CUBQL_CUDA_CALL(EventRecord(evAabbStart, stream));

    countAABBOverlapsKernel_Indirected_boxes<<<kernelGridSize, kernelBlockSize, 0, stream>>>(
        d_pairCounts, bvhA, d_boxesA, d_boxesB, d_BIter, 0, totalChunkPrims, d_AIter);

    thrust::device_ptr<int> dev_counts(d_pairCounts);
    thrust::device_ptr<int> dev_offsets(d_offsets);
    thrust::exclusive_scan(thrust::cuda::par.on(stream), dev_counts, dev_counts + totalChunkPrims, dev_offsets);

    int lastCount = 0;
    int lastOffset = 0;
    CUBQL_CUDA_CALL(
        MemcpyAsync(&lastCount, d_pairCounts + totalChunkPrims - 1, sizeof(int), cudaMemcpyDeviceToHost, stream));
    CUBQL_CUDA_CALL(
        MemcpyAsync(&lastOffset, d_offsets + totalChunkPrims - 1, sizeof(int), cudaMemcpyDeviceToHost, stream));

    CUBQL_CUDA_CALL(StreamSynchronize(stream));

    int totalChunkPairs = lastCount + lastOffset;
    if(totalChunkPairs == 0)
      continue;

    finalCandidatePairs += totalChunkPairs;

    // ----------------------------------------------------------------
    // AABB CANDIDATE GENERATION
    // ----------------------------------------------------------------
    int2* d_candidatePairs = nullptr;
    CUBQL_CUDA_CALL(MallocAsync(&d_candidatePairs, totalChunkPairs * sizeof(int2), stream));

    
    fillAABBOverlapsKernel_Indirected_Boxes<<<kernelGridSize, kernelBlockSize, 0, stream>>>(
        d_candidatePairs, d_offsets, bvhA, d_boxesA, d_boxesB, d_BIter, 0, totalChunkPrims, d_AIter);
    CUBQL_CUDA_CALL(EventRecord(evAabbEnd, stream));

    // ----------------------------------------------------------------
    // DIRECT DOUBLE-PRECISION GPU PREDICATES EVALUATION
    // ----------------------------------------------------------------
    int* d_pairStatuses = nullptr;
    CUBQL_CUDA_CALL(MallocAsync(&d_pairStatuses, totalChunkPairs * sizeof(int), stream));

    CUBQL_CUDA_CALL(EventRecord(evEvalStart, stream));
    evaluateAndCompactPairsDoubleAssisted(d_candidatePairs, d_pairStatuses, d_vertsA, d_indicesA, d_vertsB, d_indicesB,
                                          totalChunkPairs, stream);
    CUBQL_CUDA_CALL(EventRecord(evEvalEnd, stream));

   // ----------------------------------------------------------------
    // PAIR COMPACTION PASS 1: FLOAT PREDICATES -> GREEN & YELLOW
    // ----------------------------------------------------------------
   CUBQL_CUDA_CALL(EventRecord(evCompact1Start, stream)); // START PASS 1

    thrust::device_ptr<int2> dev_evaluated(d_candidatePairs);
    thrust::device_ptr<int> dev_statuses(d_pairStatuses);

    int2* d_green_raw = nullptr;
    int2* d_yellow_raw = nullptr;
    CUBQL_CUDA_CALL(MallocAsync(&d_green_raw, totalChunkPairs * sizeof(int2), stream));
    CUBQL_CUDA_CALL(MallocAsync(&d_yellow_raw, totalChunkPairs * sizeof(int2), stream));

    thrust::device_ptr<int2> dev_green_out(d_green_raw);
    thrust::device_ptr<int2> dev_yellow_out(d_yellow_raw);

    auto green_end = thrust::copy_if(thrust::cuda::par.on(stream), dev_evaluated, dev_evaluated + totalChunkPairs,
                                     dev_statuses, dev_green_out, IsTargetPairStatus{(int)PAIR_GREEN});
    auto yellow_end = thrust::copy_if(thrust::cuda::par.on(stream), dev_evaluated, dev_evaluated + totalChunkPairs,
                                      dev_statuses, dev_yellow_out, IsTargetPairStatus{(int)PAIR_YELLOW});

    int totalGreen = green_end - dev_green_out;
    int totalYellow = yellow_end - dev_yellow_out;

    tracker.confirmedGreenPairs += totalGreen;
    tracker.confirmedYellowPairs += totalYellow;

    ensureHostCapacity(hostCount + totalGreen + totalYellow);

    if(totalGreen > 0) {
      CUBQL_CUDA_CALL(MemcpyAsync(h_finalBuffer + hostCount, d_green_raw, totalGreen * sizeof(int2),
                                  cudaMemcpyDeviceToHost, stream));
      hostCount += totalGreen;
    }

    CUBQL_CUDA_CALL(EventRecord(evCompact1End, stream)); // END PASS 1

    // ----------------------------------------------------------------
    // DIRECT DOUBLE-PRECISION PREDICATES EVALUATION (ON YELLOW LIST)
    // ----------------------------------------------------------------
    int2* d_double_green_raw = nullptr;
    int2* d_orange_raw = nullptr;
    int* d_double_statuses = nullptr;
    std::vector<int2> hChunkOrange;

    if(totalYellow > 0) {
      if(exactPredicateComputeMode > 0) {
        CUBQL_CUDA_CALL(MallocAsync(&d_double_statuses, totalYellow * sizeof(int), stream));
        CUBQL_CUDA_CALL(MallocAsync(&d_double_green_raw, totalYellow * sizeof(int2), stream));
        CUBQL_CUDA_CALL(MallocAsync(&d_orange_raw, totalYellow * sizeof(int2), stream));

        CUBQL_CUDA_CALL(EventRecord(evDoubleStart, stream));
        if (exactPredicateComputeMode == 1)
        {
          evaluateAndCompactPairsDouble(d_yellow_raw, d_double_statuses, d_vertsA, d_indicesA, d_vertsB, d_indicesB,
                                      totalYellow, stream);
        }
        else if (exactPredicateComputeMode == 2)
        {
          evaluateGeometricPairsKernelBigInt(d_yellow_raw, d_double_statuses, d_vertsA, d_indicesA, d_vertsB, d_indicesB,
                                      totalYellow, stream);

        }
        else if (exactPredicateComputeMode == 3)
        {

          evaluateAndCompactPairsShewchukFloat(d_yellow_raw, d_double_statuses, d_vertsA, d_indicesA, d_vertsB, d_indicesB,
                                      totalYellow, stream);

        }
        else
        {
          evaluateTwoLapPairs(d_yellow_raw, d_double_statuses, d_vertsA, d_indicesA, d_vertsB, d_indicesB,
                                      totalYellow, stream);

        }
        CUBQL_CUDA_CALL(EventRecord(evDoubleEnd, stream));

        // PAIR COMPACTION PASS 2: DOUBLE PREDICATES -> GREEN & ORANGE
        CUBQL_CUDA_CALL(EventRecord(evCompact2Start, stream)); // START PASS 2

        thrust::device_ptr<int2> dev_yellow_in(d_yellow_raw);
        thrust::device_ptr<int> dev_dbl_statuses(d_double_statuses);
        thrust::device_ptr<int2> dev_dbl_green_out(d_double_green_raw);
        thrust::device_ptr<int2> dev_orange_out(d_orange_raw);

        auto dbl_green_end = thrust::copy_if(thrust::cuda::par.on(stream), dev_yellow_in, dev_yellow_in + totalYellow,
                                             dev_dbl_statuses, dev_dbl_green_out, IsTargetPairStatus{(int)PAIR_GREEN});
        auto orange_end = thrust::copy_if(thrust::cuda::par.on(stream), dev_yellow_in, dev_yellow_in + totalYellow,
                                          dev_dbl_statuses, dev_orange_out, IsTargetPairStatus{(int)PAIR_YELLOW});

        int totalDoubleGreen = dbl_green_end - dev_dbl_green_out;
        int totalOrange = orange_end - dev_orange_out;

        tracker.confirmedGreenPairs += totalDoubleGreen;
        tracker.confirmedOrangePairs += totalOrange;

        ensureHostCapacity(hostCount + totalDoubleGreen + totalOrange);

        if(totalDoubleGreen > 0) {
          CUBQL_CUDA_CALL(MemcpyAsync(h_finalBuffer + hostCount, d_double_green_raw, totalDoubleGreen * sizeof(int2),
                                      cudaMemcpyDeviceToHost, stream));
          hostCount += totalDoubleGreen;
        }

        if(totalOrange > 0) {
          hChunkOrange.resize(totalOrange);
          CUBQL_CUDA_CALL(
              MemcpyAsync(hChunkOrange.data(), d_orange_raw, totalOrange * sizeof(int2), cudaMemcpyDeviceToHost, stream));
        }

        CUBQL_CUDA_CALL(EventRecord(evCompact2End, stream)); // END PASS 2
      } else {
        tracker.confirmedOrangePairs += totalYellow;
        hChunkOrange.resize(totalYellow);
        CUBQL_CUDA_CALL(
            MemcpyAsync(hChunkOrange.data(), d_yellow_raw, totalYellow * sizeof(int2), cudaMemcpyDeviceToHost, stream));
      }
    }

  CUBQL_CUDA_CALL(StreamSynchronize(stream)); // Single stream sync for chunk completion

    // ----------------------------------------------------------------
    // 1. CPU EXACT FALLBACK (ORANGE LIST)
    // ----------------------------------------------------------------
    if(!hChunkOrange.empty()) {
      double tCpuStart = cuBQL::getCurrentTime();
      size_t confirmedOranges =
          filterYellowPairsTBB(meshAcpu, meshBcpu, hChunkOrange.data(), hChunkOrange.size(), h_finalBuffer + hostCount,
                               m_centerA, m_rotA, m_transA, m_centerB, m_rotB, m_transB);
      tracker.CPUPredicates += (cuBQL::getCurrentTime() - tCpuStart) * 1000.0;
      hostCount += confirmedOranges;
    }

    // ----------------------------------------------------------------
    // 2. METRICS RECORDING
    // ----------------------------------------------------------------
    float aabbMs = 0.0f, evalMs = 0.0f, doubleMs = 0.0f;
    float compact1Ms = 0.0f, compact2Ms = 0.0f;

    CUBQL_CUDA_CALL(EventElapsedTime(&aabbMs, evAabbStart, evAabbEnd));
    CUBQL_CUDA_CALL(EventElapsedTime(&evalMs, evEvalStart, evEvalEnd));
    CUBQL_CUDA_CALL(EventElapsedTime(&compact1Ms, evCompact1Start, evCompact1End));

    if(totalYellow > 0 && exactPredicateComputeMode > 0) {
      CUBQL_CUDA_CALL(EventElapsedTime(&doubleMs, evDoubleStart, evDoubleEnd));
      CUBQL_CUDA_CALL(EventElapsedTime(&compact2Ms, evCompact2Start, evCompact2End));
    }

    tracker.executionPhaseMs += aabbMs;
    tracker.fineEvaluationPhaseMs += evalMs;
    tracker.gpuDoublePredicatesMs += doubleMs;
    tracker.DownloadAndClean += (compact1Ms + compact2Ms);

    // ----------------------------------------------------------------
    // 3. RESTORED PER-CHUNK GPU MEMORY DEALLOCATION (FIXES OOM)
    // ----------------------------------------------------------------
    double tCleanStart = cuBQL::getCurrentTime();

    CUBQL_CUDA_CALL(FreeAsync(d_candidatePairs, stream));
    CUBQL_CUDA_CALL(FreeAsync(d_pairStatuses, stream));
    CUBQL_CUDA_CALL(FreeAsync(d_green_raw, stream));
    CUBQL_CUDA_CALL(FreeAsync(d_yellow_raw, stream));

    if(totalYellow > 0 && exactPredicateComputeMode > 0) {
      CUBQL_CUDA_CALL(FreeAsync(d_double_statuses, stream));
      CUBQL_CUDA_CALL(FreeAsync(d_double_green_raw, stream));
      CUBQL_CUDA_CALL(FreeAsync(d_orange_raw, stream));
    }

    tracker.DownloadAndClean += (cuBQL::getCurrentTime() - tCleanStart) * 1000.0;
  }

  // Pass ownership directly out to caller (0.00 ms overhead)
  outFinalExactPairs = h_finalBuffer;
  outFinalCount = hostCount;

  double tCleanupStart = cuBQL::getCurrentTime();

  CUBQL_CUDA_CALL(EventDestroy(evAabbStart));
  CUBQL_CUDA_CALL(EventDestroy(evAabbEnd));
  CUBQL_CUDA_CALL(EventDestroy(evEvalStart));
  CUBQL_CUDA_CALL(EventDestroy(evEvalEnd));
  CUBQL_CUDA_CALL(EventDestroy(evDoubleStart)); 
  CUBQL_CUDA_CALL(EventDestroy(evDoubleEnd));   
  CUBQL_CUDA_CALL(EventDestroy(evCompact1Start));
  CUBQL_CUDA_CALL(EventDestroy(evCompact1End));
  CUBQL_CUDA_CALL(EventDestroy(evCompact2Start));
  CUBQL_CUDA_CALL(EventDestroy(evCompact2End));
  
  if(d_BIter)
    CUBQL_CUDA_CALL(FreeAsync(d_BIter, stream));
  if(d_AIter)
    CUBQL_CUDA_CALL(FreeAsync(d_AIter, stream));
  if(d_pairCounts)
    CUBQL_CUDA_CALL(FreeAsync(d_pairCounts, stream));
  if(d_offsets)
    CUBQL_CUDA_CALL(FreeAsync(d_offsets, stream));
  if(d_chunkBatchSizes)
    CUBQL_CUDA_CALL(FreeAsync(d_chunkBatchSizes, stream));
  if(d_chunkBatchOffsets)
    CUBQL_CUDA_CALL(FreeAsync(d_chunkBatchOffsets, stream));

  CUBQL_CUDA_CALL(StreamSynchronize(stream));

  tracker.cleanupTimeMs = (cuBQL::getCurrentTime() - tCleanupStart) * 1000.0;
  tracker.totalLoopTimeMs = (cuBQL::getCurrentTime() - tTotalStart) * 1000.0;

  std::cout << "--------------------------------------------------\n\n";
  return finalCandidatePairs;
}