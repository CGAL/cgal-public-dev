#define CUBQL_GPU_BUILDER_IMPLEMENTATION 1

// Core cuBQL headers
#include "cuBQL/bvh.h"
#include "cuBQL/builder/cuda.h"

// Your custom experimental builder layout
#include "include/third-party/cubql/sm_builder_v3.h"
#include "include/third-party/cubql/sm_builder_v4.h"
#include "include/third-party/cubql/refit_forest.h"

// Custom traversal
#include "include/third-party/cubql/fixedBoxQueryv2.h"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/scan.h>
#include <thrust/execution_policy.h>
#include <thrust/system/cuda/execution_policy.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/unique.h>

#include <vector>
#include <algorithm>
#include <iostream> // For deep debug printings
#include "samples/common/loadOBJ.h"

// Include modular execution targets
#include "DualTreeStep.h"
#include "rapidDescendKernel.h"
// #include "batchedCrossIntersection.h"

#include "crossCheckNew.h"
#include "kernelsTestBVHAlternative.h"
#include "prune_pipeline.h"
#include "batchedCrossIntersectionDouble.h"

#include "global_box.h"
#include "utils.h"



#ifndef _ALLOC
#define _ALLOC(ptr, count, stream, memResource)                                                                        \
  CUBQL_CUDA_CALL(MallocAsync((void**)&(ptr), (size_t)(count) * sizeof(*(ptr)), stream))
#endif

#ifndef _FREE
#define _FREE(ptr, stream, memResource) CUBQL_CUDA_CALL(FreeAsync((ptr), stream))
#endif





extern "C" void kernelsTestBVHV3(Mesh& meshAcpu,
                                 Mesh& meshBcpu,
                                 const double3* hVertsA,
                                 int numVertsA,
                                 const uint3* hIndicesA,
                                 int numTrianglesA,
                                 int maxCellSizeA,
                                 const double3* hVertsB,
                                 int numVertsB,
                                 const uint3* hIndicesB,
                                 int numTrianglesB,
                                 int maxCellSizeB,
                                 int batchMultiplier,
                                 int mode,
                                 int leafThreshold,
                                 ExecutionStats& stats,
                                 int2*& outFinalExactPairs,      
    size_t& outFinalCount) {
  if(numTrianglesA <= 0 || numTrianglesB <= 0) {
    return;
  }

  // std::cout << "\n==================================================" << std::endl;
  // std::cout << " [DEEP DEBUG] => Entering kernelsTestBVHV2 Pipeline..." << std::endl;
  // std::cout << "==================================================" << std::endl;

  double tPipelineStart = cuBQL::getCurrentTime();
  cudaStream_t stream = 0;
  cuBQL::DeviceMemoryResource memResource;
  cuBQL::BuildConfig buildConfig;
  buildConfig.makeLeafThreshold = leafThreshold;

  // --------------------------------------------------------------------
  // INITIAL RAW ALLOCATION & COPY TRACKING
  // --------------------------------------------------------------------
  double tAllocStart = cuBQL::getCurrentTime();

  // Allocate and stream indexed buffers for Mesh A
  double3* dVertsA = nullptr;
  uint3* dIndicesA = nullptr;
  CUBQL_CUDA_CALL(Malloc(&dVertsA, numVertsA * sizeof(float3)));
  CUBQL_CUDA_CALL(Memcpy(dVertsA, hVertsA, numVertsA * sizeof(float3), cudaMemcpyHostToDevice));
  CUBQL_CUDA_CALL(Malloc(&dIndicesA, numTrianglesA * sizeof(uint3)));
  CUBQL_CUDA_CALL(Memcpy(dIndicesA, hIndicesA, numTrianglesA * sizeof(uint3), cudaMemcpyHostToDevice));

  // Allocate and stream indexed buffers for Mesh B
  double3* dVertsB = nullptr;
  uint3* dIndicesB = nullptr;
  CUBQL_CUDA_CALL(Malloc(&dVertsB, numVertsB * sizeof(float3)));
  CUBQL_CUDA_CALL(Memcpy(dVertsB, hVertsB, numVertsB * sizeof(float3), cudaMemcpyHostToDevice));
  CUBQL_CUDA_CALL(Malloc(&dIndicesB, numTrianglesB * sizeof(uint3)));
  CUBQL_CUDA_CALL(Memcpy(dIndicesB, hIndicesB, numTrianglesB * sizeof(uint3), cudaMemcpyHostToDevice));



cuBQL::box3f* dBoxesA;
cuBQL::box3f* dBoxesB;
  // BVH Bounding Box Buffers
  CUBQL_CUDA_CALL(MallocAsync((void**)&dBoxesA, numTrianglesA * sizeof(cuBQL::box3f), stream));
  CUBQL_CUDA_CALL(MallocAsync((void**)&dBoxesB, numTrianglesB * sizeof(cuBQL::box3f), stream));

  // Generate initial bounding boxes directly from vertices and indices
  generateBoxesTrisKernel<<<cuBQL::divRoundUp(numTrianglesA, 256), 256, 0, stream>>>(dBoxesA, dVertsA,
                                                                                         dIndicesA, numTrianglesA);

  generateBoxesTrisKernel<<<cuBQL::divRoundUp(numTrianglesB, 256), 256, 0, stream>>>(dBoxesB, dVertsB,
                                                                                         dIndicesB, numTrianglesB);



  cudaDeviceSynchronize();
  double tAllocEnd = cuBQL::getCurrentTime();
  stats.initialAllocAndCopyMs = (tAllocEnd - tAllocStart) * 1000.0;

  cuBQL::box3f globalBoxA;
  cuBQL::box3f globalBoxB;

  // --------------------------------------------------------------------
  // THRUST DEVICE ALLOCATION / INITIALIZATION OVERHEAD TRACKING
  // --------------------------------------------------------------------
  double tThrustInitStart = cuBQL::getCurrentTime();

  thrust::device_vector<uint32_t> d_outPairsA;
  thrust::device_vector<uint32_t> d_outPairsB;

  thrust::device_vector<uint32_t> d_outOffsetsB;
  thrust::device_vector<uint32_t> d_outPrimsFlatB;


  cudaDeviceSynchronize();
  double tThrustInitEnd = cuBQL::getCurrentTime();
  stats.thrustInitOverheadMs = (tThrustInitEnd - tThrustInitStart) * 1000.0;

  // --------------------------------------------------------------------
  // INJECTED V3 WORKFLOW: MESH A INITIALIZATION
  // --------------------------------------------------------------------
  CUBQL_CUDA_CALL(StreamSynchronize(stream));
  double tInitAStart = cuBQL::getCurrentTime();

  // Compute global bounding box for Mesh A
  cuBQL::utils::computeGlobalBoxParallel<float, 3>(globalBoxA, dBoxesA, numTrianglesA, stream, memResource);

  cuBQL::box3f* outNodeBoxesA = nullptr;
  uint32_t* outSortedPrimIDsA = nullptr;
  uint32_t* outNodeOffsetsA = nullptr;
  uint32_t outTotalActiveCellsA = 0;

  cuBQL::gpuBuilder_v3::test_speedrun_initialization_linear(dBoxesA, numTrianglesA, (uint32_t)maxCellSizeA, globalBoxA,
                                                            outNodeBoxesA, outSortedPrimIDsA, outNodeOffsetsA,
                                                            outTotalActiveCellsA, stream, memResource);
  CUBQL_CUDA_CALL(StreamSynchronize(stream));
  stats.buildRefitMeshAMs = (cuBQL::getCurrentTime() - tInitAStart) * 1000.0;

  // --------------------------------------------------------------------
  // INJECTED V3 WORKFLOW: MESH B INITIALIZATION
  // --------------------------------------------------------------------
  CUBQL_CUDA_CALL(StreamSynchronize(stream));
  double tInitBStart = cuBQL::getCurrentTime();

  // Compute global bounding box for Mesh B
  cuBQL::utils::computeGlobalBoxParallel<float, 3>(globalBoxB, dBoxesB, numTrianglesB, stream, memResource);

  cuBQL::box3f* outNodeBoxesB = nullptr;
  uint32_t* outSortedPrimIDsB = nullptr;
  uint32_t* outNodeOffsetsB = nullptr;
  uint32_t outTotalActiveCellsB = 0;

  cuBQL::gpuBuilder_v3::test_speedrun_initialization_linear(dBoxesB, numTrianglesB, (uint32_t)maxCellSizeB, globalBoxB,
                                                            outNodeBoxesB, outSortedPrimIDsB, outNodeOffsetsB,
                                                            outTotalActiveCellsB, stream, memResource);
  CUBQL_CUDA_CALL(StreamSynchronize(stream));
  stats.buildRefitMeshBMs = (cuBQL::getCurrentTime() - tInitBStart) * 1000.0;

  // std::cout << "\n--------------------------------------------------" << std::endl;
  // std::cout << " [LINEAR INITIALIZATION METRICS OVERVIEW]" << std::endl;
  // std::cout << " -> Tree A: Total Compounded Active Cells = " << outTotalActiveCellsA << std::endl;
  // std::cout << " -> Tree B: Total Compounded Active Cells = " << outTotalActiveCellsB << std::endl;
  // std::cout << "--------------------------------------------------" << std::endl;

  // Map to historical baseline counts expected by downstream processing before modification
  uint32_t h_outMarkedCountA = outTotalActiveCellsA;
  uint32_t h_outMarkedCountB = outTotalActiveCellsB;

  // --------------------------------------------------------------------
  // CROSS CHECK
  // --------------------------------------------------------------------
  double tCrossStart = cuBQL::getCurrentTime();

  uint32_t totalIntersections = executeBoxCrossCheck<float, 3>(outNodeBoxesA, h_outMarkedCountA, outNodeBoxesB,
                                                               h_outMarkedCountB, d_outPairsA, d_outPairsB, stream);
  CUBQL_CUDA_CALL(StreamSynchronize(stream));

  double tCrossEnd = cuBQL::getCurrentTime();
  stats.gpuCrossCheckEngineMs = (tCrossEnd - tCrossStart) * 1000.0;

  uint64_t totalPossiblePairs = (uint64_t)h_outMarkedCountA * h_outMarkedCountB;
  double intersectionPercentage =
      totalPossiblePairs > 0 ? ((double)totalIntersections / totalPossiblePairs) * 100.0 : 0.0;

  // std::cout << "\n==================================================" << std::endl;
  // std::cout << "          LAUNCHING CROSS-CHECK KERNELS           " << std::endl;
  // std::cout << "==================================================" << std::endl;
  // std::cout << " -> Intersecting Pairs Detected     : " << totalIntersections << " pairs" << std::endl;
  // std::cout << " -> True Active Matrix Evaluated    : " << totalPossiblePairs << " pairings" << std::endl;
  // std::cout << " -> Percentage of Total Overlaps    : " << intersectionPercentage << "%" << std::endl;

  // Save layout baselines prior to pruning
  uint32_t initialCellsA = outTotalActiveCellsA;
  uint32_t initialCellsB = outTotalActiveCellsB;

  // --------------------------------------------------------------------
  // CALCULATE CURRENT PRIMITIVE SUMS FOR PRUNING
  // --------------------------------------------------------------------
  uint32_t firstOffsetA = 0, lastOffsetA = 0;
  CUBQL_CUDA_CALL(MemcpyAsync(&firstOffsetA, outNodeOffsetsA, sizeof(uint32_t), cudaMemcpyDeviceToHost, stream));
  CUBQL_CUDA_CALL(
      MemcpyAsync(&lastOffsetA, outNodeOffsetsA + h_outMarkedCountA, sizeof(uint32_t), cudaMemcpyDeviceToHost, stream));

  uint32_t firstOffsetB = 0, lastOffsetB = 0;
  CUBQL_CUDA_CALL(MemcpyAsync(&firstOffsetB, outNodeOffsetsB, sizeof(uint32_t), cudaMemcpyDeviceToHost, stream));
  CUBQL_CUDA_CALL(
      MemcpyAsync(&lastOffsetB, outNodeOffsetsB + h_outMarkedCountB, sizeof(uint32_t), cudaMemcpyDeviceToHost, stream));
  CUBQL_CUDA_CALL(StreamSynchronize(stream));

  int currentPrimsNumA = (int)(lastOffsetA - firstOffsetA);
  int currentPrimsNumB = (int)(lastOffsetB - firstOffsetB);

  // std::cout << " [POPULATED CELLS CHECK SUMMARY]" << std::endl;
  // std::cout << " -> Tree A Active Nodes: " << initialCellsA << " | Total Verified Prims = " << currentPrimsNumA
  //           << std::endl;
  // std::cout << " -> Tree B Active Nodes: " << initialCellsB << " | Total Verified Prims = " << currentPrimsNumB
  //           << std::endl;
  // std::cout << "--------------------------------------------------\n" << std::endl;

  // --------------------------------------------------------------------
  // PARALLEL STREAM COMPACTION & PRUNING ALGORITHM (POST-BUILD)
  // --------------------------------------------------------------------
  // std::cout << " \n[PIPELINE PRUNING] => Dispatched grid streaming parallel re-index compaction (POST-BUILD)..."
  //           << std::endl;

  double pruneStartA = cuBQL::getCurrentTime();
  parallelPruneAndReindexAll(thrust::raw_pointer_cast(d_outPairsA.data()), totalIntersections, outSortedPrimIDsA,
                             outNodeOffsetsA, outTotalActiveCellsA, currentPrimsNumA, stream, memResource);
  CUBQL_CUDA_CALL(StreamSynchronize(stream));
  double pruneMsA = (cuBQL::getCurrentTime() - pruneStartA) * 1000.0;
  // std::cout << " -> Mesh A Structural Pruning (Post)  : " << pruneMsA << " ms" << std::endl;
  // std::cout << "    * Active Cells Surviving: " << initialCellsA << " -> " << outTotalActiveCellsA << std::endl;

  double pruneStartB = cuBQL::getCurrentTime();
  parallelPruneAndReindexAll(thrust::raw_pointer_cast(d_outPairsB.data()), totalIntersections, outSortedPrimIDsB,
                             outNodeOffsetsB, outTotalActiveCellsB, currentPrimsNumB, stream, memResource);
  CUBQL_CUDA_CALL(StreamSynchronize(stream));
  double pruneMsB = (cuBQL::getCurrentTime() - pruneStartB) * 1000.0;
  // std::cout << " -> Mesh B Structural Pruning (Post)  : " << pruneMsB << " ms" << std::endl;
  // std::cout << "    * Active Cells Surviving: " << initialCellsB << " -> " << outTotalActiveCellsB << std::endl;


  // --------------------------------------------------------------------
  // GPU BUILDER V4 FOREST EXPANSION
  // --------------------------------------------------------------------
  // std::cout << " [FOREST EXPANSION] => Launching Level-by-Level Parallel Sub-Tree Compilation..." << std::endl;

  const uint32_t maxNodesA = 2u * (uint32_t)numTrianglesA + (outTotalActiveCellsA + (outTotalActiveCellsA & 1u)) + 2u;
  const uint32_t maxNodesB = 2u * (uint32_t)numTrianglesB + (outTotalActiveCellsB + (outTotalActiveCellsB & 1u)) + 2u;

  thrust::device_vector<uint32_t> d_nodeDescendantCountsA(maxNodesA, 0);
  thrust::device_vector<uint32_t> d_nodeDescendantCountsB(maxNodesB, 0);

  // --- BUILD FOREST MESH A ---
  double tForestAStart = cuBQL::getCurrentTime();
  cuBQL::BinaryBVH<float, 3> bvhA;
  cuBQL::gpuBuilder_v4::build_forest<float, 3>(
      bvhA, dBoxesA, currentPrimsNumA, numTrianglesA, outTotalActiveCellsA, outSortedPrimIDsA, outNodeOffsetsA,
      buildConfig, thrust::raw_pointer_cast(d_nodeDescendantCountsA.data()), stream, memResource);

  cuBQL::cuda_forest::refit_forest<float, 3>(bvhA, dBoxesA, outTotalActiveCellsA, stream, memResource);
  CUBQL_CUDA_CALL(StreamSynchronize(stream));
  double forestAMs = (cuBQL::getCurrentTime() - tForestAStart) * 1000.0;
  stats.buildRefitMeshAMs += forestAMs;
  // std::cout << " -> Mesh A Forest BVH Expansion     : " << forestAMs << " ms | Nodes: " << bvhA.numNodes << std::endl;

  // --- BUILD FOREST MESH B ---
  double tForestBStart = cuBQL::getCurrentTime();
  cuBQL::BinaryBVH<float, 3> bvhB;
  cuBQL::gpuBuilder_v4::build_forest<float, 3>(
      bvhB, dBoxesB, currentPrimsNumB, numTrianglesB, outTotalActiveCellsB, outSortedPrimIDsB, outNodeOffsetsB,
      buildConfig, thrust::raw_pointer_cast(d_nodeDescendantCountsB.data()), stream, memResource);
  cuBQL::cuda_forest::refit_forest<float, 3>(bvhB, dBoxesB, outTotalActiveCellsB, stream, memResource);
  CUBQL_CUDA_CALL(StreamSynchronize(stream));
  double forestBMs = (cuBQL::getCurrentTime() - tForestBStart) * 1000.0;
  stats.buildRefitMeshBMs += forestBMs;
  // std::cout << " -> Mesh B Forest BVH Expansion     : " << forestBMs << " ms | Nodes: " << bvhB.numNodes << std::endl;


  // --------------------------------------------------------------------
  // PREPARE DOWNSTREAM MAPPINGS USING COMPACTED SURVIVAL COUNTS
  // --------------------------------------------------------------------
  uint32_t finalActiveCellsA = outTotalActiveCellsA;
  uint32_t finalActiveCellsB = outTotalActiveCellsB;

  thrust::device_vector<uint32_t> d_markedNodeIndicesA(finalActiveCellsA);
  thrust::sequence(thrust::device.on(stream), d_markedNodeIndicesA.begin(), d_markedNodeIndicesA.end());

  thrust::device_vector<uint32_t> d_markedNodeIndicesB(finalActiveCellsB);
  thrust::sequence(thrust::device.on(stream), d_markedNodeIndicesB.begin(), d_markedNodeIndicesB.end());

  uint32_t maxPossibleNodesB = bvhB.numNodes;
  thrust::device_vector<uint32_t> d_reverseMapB(maxPossibleNodesB, 0);

  // std::cout << " [OK] Sequences initialization complete. Proceeding to Dual Tree Step..." << std::endl;
  // std::cout << "==================================================\n" << std::endl;

  // --------------------------------------------------------------------
  // DUAL-TREE TRAVERSAL STEP OVERHEAD PHASE
  // --------------------------------------------------------------------
  double tDualStepStart = cuBQL::getCurrentTime();

  if(mode > 0) {
    executeDualTreeStep(mode, maxCellSizeA, maxCellSizeB, d_outPairsA, d_outPairsB, d_markedNodeIndicesA,
                        d_markedNodeIndicesB, d_nodeDescendantCountsA, d_nodeDescendantCountsB, finalActiveCellsA,
                        finalActiveCellsB, bvhA, bvhB);
    cudaDeviceSynchronize();
  }

  double tDualStepEnd = cuBQL::getCurrentTime();
  stats.dualTreeStepMs = (mode > 0) ? (tDualStepEnd - tDualStepStart) * 1000.0 : 0.0;

  int totalBatches = d_outPairsA.size();
  IntersectionTimeTracker tracker;
  uint64_t finalCandidatePairs = 0;

  // --------------------------------------------------------------------
  // CONDITIONAL PIPELINE EXECUTION BASED ON MODE
  // --------------------------------------------------------------------
  double tGpuBfsStart = cuBQL::getCurrentTime();

  executeRapidDescentBFS(bvhB, currentPrimsNumB, d_markedNodeIndicesB, d_nodeDescendantCountsB, finalActiveCellsB,
                         d_outOffsetsB, d_outPrimsFlatB);

  if(finalActiveCellsB > 0) {
    int blockSize = 256;
    int gridSize = (finalActiveCellsB + blockSize - 1) / blockSize;
    populateReverseMapBKernel<<<gridSize, blockSize, 0, stream>>>(thrust::raw_pointer_cast(d_reverseMapB.data()),
                                                                  thrust::raw_pointer_cast(d_markedNodeIndicesB.data()),
                                                                  finalActiveCellsB);
  }

  cudaDeviceSynchronize();
  double tGpuBfsEnd = cuBQL::getCurrentTime();

    Point3 m_centerA = Point3{0.0, 0.0, 0.0};
    Point3 m_centerB = Point3{0.0, 0.0, 0.0};
    double3 m_rotA    = double3{0.0, 0.0, 0.0};
    double3 m_transA  = double3{0.0, 0.0, 0.0};
    double3 m_rotB    = double3{0.0, 0.0, 0.0};
    double3 m_transB  = double3{0.0, 0.0, 0.0};

  finalCandidatePairs = executeBatchedCrossIntersectionLoopDouble(
        meshAcpu, meshBcpu, batchMultiplier, totalBatches, d_outPairsA, d_outPairsB, d_reverseMapB,
        d_markedNodeIndicesB, d_outOffsetsB, d_outPrimsFlatB, d_nodeDescendantCountsB, finalActiveCellsB, bvhA,
        dBoxesA, // <-- ADDED: Precomputed boxes for Mesh A
        dBoxesB, // <-- ADDED: Precomputed boxes for Mesh B
        dVertsA, dIndicesA, dVertsB, dIndicesB, outFinalExactPairs, outFinalCount, tracker, m_centerA,
        m_centerB, m_rotA, m_transA, m_rotB, m_transB, stream);



//   if(activateAsyncDownload == 0) {

//     finalCandidatePairs = executeBatchedCrossIntersectionLoopV2(
//         meshAcpu, meshBcpu, batchMultiplier, totalBatches, d_outPairsA, d_outPairsB, d_reverseMapB,
//         d_markedNodeIndicesB, d_outOffsetsB, d_outPrimsFlatB, d_nodeDescendantCountsB, finalActiveCellsB, bvhA, dMeshA,
//         dMeshB, dMeshMetricsA, dMeshMetricsB, finalExactPairs, tracker, stream);

//   } else {
//     finalCandidatePairs = executeBatchedCrossIntersectionLoopV3(
//         meshAcpu, meshBcpu, batchMultiplier, totalBatches, d_outPairsA, d_outPairsB, d_reverseMapB,
//         d_markedNodeIndicesB, d_outOffsetsB, d_outPrimsFlatB, d_nodeDescendantCountsB, finalActiveCellsB, bvhA, dMeshA,
//         dMeshB, dMeshMetricsA, dMeshMetricsB, finalExactPairs, tracker, stream, activateAsyncDownload);
//   }
  // --------------------------------------------------------------------
  // EXPLICIT CLEANUP & RECOVERY METRIC TRACKING
  // --------------------------------------------------------------------

  //std::cout << "Starting to clean up" << std::endl;

  double tCleanupStart = cuBQL::getCurrentTime();

  //CUBQL_CUDA_CALL(Free(dMeshA));
  CUBQL_CUDA_CALL(Free(dBoxesA));
  //CUBQL_CUDA_CALL(Free(dMeshB));
  CUBQL_CUDA_CALL(Free(dBoxesB));

  // Release the intermediate indexed topology allocations
  CUBQL_CUDA_CALL(Free(dVertsA));
  CUBQL_CUDA_CALL(Free(dIndicesA));
  CUBQL_CUDA_CALL(Free(dVertsB));
  CUBQL_CUDA_CALL(Free(dIndicesB));

  if(outNodeBoxesA)
    _FREE(outNodeBoxesA, stream, memResource);
  if(outSortedPrimIDsA)
    _FREE(outSortedPrimIDsA, stream, memResource);
  if(outNodeOffsetsA)
    _FREE(outNodeOffsetsA, stream, memResource);
  if(outNodeBoxesB)
    _FREE(outNodeBoxesB, stream, memResource);
  if(outSortedPrimIDsB)
    _FREE(outSortedPrimIDsB, stream, memResource);
  if(outNodeOffsetsB)
    _FREE(outNodeOffsetsB, stream, memResource);

  cuBQL::cuda::free(bvhA, stream, memResource);
  cuBQL::cuda::free(bvhB, stream, memResource);

  d_markedNodeIndicesA.shrink_to_fit();
  d_nodeDescendantCountsA.shrink_to_fit();
  d_markedNodeIndicesB.shrink_to_fit();
  d_nodeDescendantCountsB.shrink_to_fit();
  d_reverseMapB.shrink_to_fit();
  d_outPairsA.shrink_to_fit();
  d_outPairsB.shrink_to_fit();
  d_outOffsetsB.shrink_to_fit();
  d_outPrimsFlatB.shrink_to_fit();

  cudaDeviceSynchronize();
  double tCleanupEnd = cuBQL::getCurrentTime();
  stats.finalCleanupSyncMs = (tCleanupEnd - tCleanupStart) * 1000.0;

  double tPipelineEnd = cuBQL::getCurrentTime();

  stats.meshATotalNodes = bvhA.numNodes;
  stats.meshAExtractedTargets = h_outMarkedCountA;
  stats.meshBTotalNodes = bvhB.numNodes;
  stats.meshBExtractedTargets = h_outMarkedCountB;
  stats.totalIntersections = totalIntersections;
  stats.totalPossiblePairs = totalPossiblePairs;
  stats.intersectionPercentage = intersectionPercentage;
  stats.gpuCrossCheckEngineMs = (tCrossEnd - tCrossStart) * 1000.0;
  stats.parallelDfsDescentBMs = (tGpuBfsEnd - tGpuBfsStart) * 1000.0;
  stats.GPUTotalTime = (tPipelineEnd - tPipelineStart) * 1000.0;
  stats.totalCrissCrossBatches = totalBatches;
  stats.finalAabbCandidatePairs = finalCandidatePairs;
  stats.confirmedGreenPairs = stats.loopTracker.confirmedGreenPairs;
  stats.confirmedYellowPairs = stats.loopTracker.confirmedYellowPairs;
  stats.loopTracker = tracker;
}