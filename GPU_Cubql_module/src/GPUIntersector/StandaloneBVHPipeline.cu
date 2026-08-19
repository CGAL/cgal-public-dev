#define CUBQL_GPU_BUILDER_IMPLEMENTATION 1

// Core cuBQL headers
#include "cuBQL/bvh.h"
#include "cuBQL/builder/cuda.h"

// Grid forest and partitioner modules
#include "third-party/cubql/grid_partitioner.h"
#include "third-party/cubql/forest_builder.h"
#include "third-party/cubql/refit_forest.h"

// Traversal modules
#include "third-party/cubql/fixedBoxQueryv2.h"

// CUDA & Thrust Headers
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/scan.h>
#include <thrust/execution_policy.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>

// Standard Library Headers
#include <vector>
#include <algorithm>
#include <iostream>

// Mesh Loader & Traversal Engines
#include "loadOBJ.h"
#include "traversal/DualTreeStep.h"
#include "traversal/rapidDescendKernel.h"
#include "traversal/crossCheckNoBVH.h"
#include "traversal/prune_pipeline.h"
#include "single_tree/SingleTreeBatchIntersector.h"

// Utilities & Geometry
#include "common/ExecutionStats.h"
#include "common/global_box.h"
#include "common/utils.h"

#ifndef _ALLOC
#define _ALLOC(ptr, count, stream, memResource)                                                                        \
  CUBQL_CUDA_CALL(MallocAsync((void**)&(ptr), (size_t)(count) * sizeof(*(ptr)), stream))
#endif

#ifndef _FREE
#define _FREE(ptr, stream, memResource) CUBQL_CUDA_CALL(FreeAsync((ptr), stream))
#endif

/**
 * @brief Main execution entry point for the hybrid Grid-Forest BVH Mesh Intersection Pipeline.
 *
 * Performs end-to-end GPU spatial partitioning, coarse AABB overlap cross-checking,
 * post-build structural pruning, forest sub-tree expansion, and batched single-tree
 * traversal with exact GPU predicates.
 */
extern "C" void runGridForestIntersectionPipeline(Mesh& meshA,
                                                  Mesh& meshB,
                                                  const double3* h_vertsA,
                                                  int numVertsA,
                                                  const uint3* h_indicesA,
                                                  int numTrianglesA,
                                                  int maxCellSizeA,
                                                  const double3* h_vertsB,
                                                  int numVertsB,
                                                  const uint3* h_indicesB,
                                                  int numTrianglesB,
                                                  int maxCellSizeB,
                                                  int batchMultiplier,
                                                  int numberOfDualTreeSteps,
                                                  int gpuPredicateMode,
                                                  int leafThreshold,
                                                  ExecutionStats& stats,
                                                  int2*& outFinalExactPairs,
                                                  size_t& outFinalCount,
                                                  Point3 centerA,
                                                  Point3 centerB,
                                                  double3 rotA,
                                                  double3 transA,
                                                  double3 rotB,
                                                  double3 transB) {
  if(numTrianglesA <= 0 || numTrianglesB <= 0) {
    return;
  }

  const double tPipelineStart = cuBQL::getCurrentTime();
  cudaStream_t stream = 0;
  cuBQL::DeviceMemoryResource memResource;
  cuBQL::BuildConfig buildConfig;
  buildConfig.makeLeafThreshold = leafThreshold;

  // ============================================================================
  // 1. DEVICE ALLOCATION & MESH DATA TRANSFER
  // ============================================================================
  const double tAllocStart = cuBQL::getCurrentTime();

  double3* d_vertsA = nullptr;
  uint3* d_indicesA = nullptr;
  _ALLOC(d_vertsA, numVertsA, stream, memResource);
  _ALLOC(d_indicesA, numTrianglesA, stream, memResource);
  CUBQL_CUDA_CALL(MemcpyAsync(d_vertsA, h_vertsA, numVertsA * sizeof(double3), cudaMemcpyHostToDevice, stream));
  CUBQL_CUDA_CALL(MemcpyAsync(d_indicesA, h_indicesA, numTrianglesA * sizeof(uint3), cudaMemcpyHostToDevice, stream));

  double3* d_vertsB = nullptr;
  uint3* d_indicesB = nullptr;
  _ALLOC(d_vertsB, numVertsB, stream, memResource);
  _ALLOC(d_indicesB, numTrianglesB, stream, memResource);
  CUBQL_CUDA_CALL(MemcpyAsync(d_vertsB, h_vertsB, numVertsB * sizeof(double3), cudaMemcpyHostToDevice, stream));
  CUBQL_CUDA_CALL(MemcpyAsync(d_indicesB, h_indicesB, numTrianglesB * sizeof(uint3), cudaMemcpyHostToDevice, stream));

  // ============================================================================
  // 2. WORLD-SPACE VERTEX TRANSFORMATIONS & AABB GENERATION
  // ============================================================================
  const Mat3x3 rotMatrixA = makeRotationMatrixDeg(rotA.x, rotA.y, rotA.z);
  const double3 centerPosA = make_double3(centerA.x(), centerA.y(), centerA.z());
  launchTransformVertices(d_vertsA, d_vertsA, numVertsA, rotMatrixA, centerPosA, transA, stream);

  const Mat3x3 rotMatrixB = makeRotationMatrixDeg(rotB.x, rotB.y, rotB.z);
  const double3 centerPosB = make_double3(centerB.x(), centerB.y(), centerB.z());
  launchTransformVertices(d_vertsB, d_vertsB, numVertsB, rotMatrixB, centerPosB, transB, stream);

  cuBQL::box3f* d_boxesA = nullptr;
  cuBQL::box3f* d_boxesB = nullptr;
  _ALLOC(d_boxesA, numTrianglesA, stream, memResource);
  _ALLOC(d_boxesB, numTrianglesB, stream, memResource);

  generateBoxesTrisKernel<<<cuBQL::divRoundUp(numTrianglesA, 256), 256, 0, stream>>>(d_boxesA, d_vertsA, d_indicesA,
                                                                                     numTrianglesA);
  generateBoxesTrisKernel<<<cuBQL::divRoundUp(numTrianglesB, 256), 256, 0, stream>>>(d_boxesB, d_vertsB, d_indicesB,
                                                                                     numTrianglesB);

  CUBQL_CUDA_CALL(StreamSynchronize(stream));
  stats.initialAllocAndCopyMs = (cuBQL::getCurrentTime() - tAllocStart) * 1000.0;

  // ============================================================================
  // 3. LEVEL-CUT SPATIAL PARTITIONING (GRID FOREST PREPARATION)
  // ============================================================================
  const double tThrustStart = cuBQL::getCurrentTime();
  thrust::device_vector<uint32_t> d_outPairsA;
  thrust::device_vector<uint32_t> d_outPairsB;
  thrust::device_vector<uint32_t> d_outOffsetsB;
  thrust::device_vector<uint32_t> d_outPrimsFlatB;
  CUBQL_CUDA_CALL(StreamSynchronize(stream));
  stats.thrustInitOverheadMs = (cuBQL::getCurrentTime() - tThrustStart) * 1000.0;

  // Partition Mesh A
  const double tInitAStart = cuBQL::getCurrentTime();
  cuBQL::box3f globalBoxA;
  cuBQL::utils::computeGlobalBoxParallel<float, 3>(globalBoxA, d_boxesA, numTrianglesA, stream, memResource);

  cuBQL::box3f* d_nodeBoxesA = nullptr;
  uint32_t* d_sortedPrimIDsA = nullptr;
  uint32_t* d_nodeOffsetsA = nullptr;
  uint32_t totalActiveCellsA = 0;

  cuBQL::ext::grid_forest::partitioner::partitionPrimitivesToLevelCut(
      d_boxesA, numTrianglesA, (uint32_t)maxCellSizeA, globalBoxA, d_nodeBoxesA, d_sortedPrimIDsA, d_nodeOffsetsA,
      totalActiveCellsA, stream, memResource);

  CUBQL_CUDA_CALL(StreamSynchronize(stream));
  stats.buildRefitMeshAMs = (cuBQL::getCurrentTime() - tInitAStart) * 1000.0;

  // Partition Mesh B
  const double tInitBStart = cuBQL::getCurrentTime();
  cuBQL::box3f globalBoxB;
  cuBQL::utils::computeGlobalBoxParallel<float, 3>(globalBoxB, d_boxesB, numTrianglesB, stream, memResource);

  cuBQL::box3f* d_nodeBoxesB = nullptr;
  uint32_t* d_sortedPrimIDsB = nullptr;
  uint32_t* d_nodeOffsetsB = nullptr;
  uint32_t totalActiveCellsB = 0;

  cuBQL::ext::grid_forest::partitioner::partitionPrimitivesToLevelCut(
      d_boxesB, numTrianglesB, (uint32_t)maxCellSizeB, globalBoxB, d_nodeBoxesB, d_sortedPrimIDsB, d_nodeOffsetsB,
      totalActiveCellsB, stream, memResource);

  CUBQL_CUDA_CALL(StreamSynchronize(stream));
  stats.buildRefitMeshBMs = (cuBQL::getCurrentTime() - tInitBStart) * 1000.0;

  // ============================================================================
  // 4. COARSE BOUNDING BOX CROSS-CHECK
  // ============================================================================
  const double tCrossStart = cuBQL::getCurrentTime();
  const uint32_t initialMarkedCountA = totalActiveCellsA;
  const uint32_t initialMarkedCountB = totalActiveCellsB;

  const uint32_t totalIntersections = executeBoxCrossCheck<float, 3>(
      d_nodeBoxesA, initialMarkedCountA, d_nodeBoxesB, initialMarkedCountB, d_outPairsA, d_outPairsB, stream);

  CUBQL_CUDA_CALL(StreamSynchronize(stream));
  const double tCrossEnd = cuBQL::getCurrentTime();
  stats.gpuCrossCheckEngineMs = (tCrossEnd - tCrossStart) * 1000.0;

  const uint64_t totalPossiblePairs = (uint64_t)initialMarkedCountA * initialMarkedCountB;
  const double intersectionPercentage =
      totalPossiblePairs > 0 ? ((double)totalIntersections / totalPossiblePairs) * 100.0 : 0.0;

  // ============================================================================
  // 5. POST-BUILD PARALLEL STREAM COMPACTION & PRUNING
  // ============================================================================
  uint32_t firstOffsetA = 0, lastOffsetA = 0;
  CUBQL_CUDA_CALL(MemcpyAsync(&firstOffsetA, d_nodeOffsetsA, sizeof(uint32_t), cudaMemcpyDeviceToHost, stream));
  CUBQL_CUDA_CALL(MemcpyAsync(&lastOffsetA, d_nodeOffsetsA + initialMarkedCountA, sizeof(uint32_t),
                              cudaMemcpyDeviceToHost, stream));

  uint32_t firstOffsetB = 0, lastOffsetB = 0;
  CUBQL_CUDA_CALL(MemcpyAsync(&firstOffsetB, d_nodeOffsetsB, sizeof(uint32_t), cudaMemcpyDeviceToHost, stream));
  CUBQL_CUDA_CALL(MemcpyAsync(&lastOffsetB, d_nodeOffsetsB + initialMarkedCountB, sizeof(uint32_t),
                              cudaMemcpyDeviceToHost, stream));
  CUBQL_CUDA_CALL(StreamSynchronize(stream));

  // Remove 'const' from currentPrimsNumA and currentPrimsNumB
  int currentPrimsNumA = (int)(lastOffsetA - firstOffsetA);
  int currentPrimsNumB = (int)(lastOffsetB - firstOffsetB);

  parallelPruneAndReindexAll(thrust::raw_pointer_cast(d_outPairsA.data()), totalIntersections, d_sortedPrimIDsA,
                             d_nodeOffsetsA, totalActiveCellsA, currentPrimsNumA, stream, memResource);

  parallelPruneAndReindexAll(thrust::raw_pointer_cast(d_outPairsB.data()), totalIntersections, d_sortedPrimIDsB,
                             d_nodeOffsetsB, totalActiveCellsB, currentPrimsNumB, stream, memResource);
  CUBQL_CUDA_CALL(StreamSynchronize(stream));

  // ============================================================================
  // 6. SUB-TREE FOREST BVH CONSTRUCTION
  // ============================================================================
  const uint32_t maxNodesA = 2u * (uint32_t)numTrianglesA + (totalActiveCellsA + (totalActiveCellsA & 1u)) + 2u;
  const uint32_t maxNodesB = 2u * (uint32_t)numTrianglesB + (totalActiveCellsB + (totalActiveCellsB & 1u)) + 2u;

  thrust::device_vector<uint32_t> d_nodeDescendantCountsA(maxNodesA, 0);
  thrust::device_vector<uint32_t> d_nodeDescendantCountsB(maxNodesB, 0);

  // Build Sub-BVH Forest for Mesh A
  const double tForestAStart = cuBQL::getCurrentTime();
  cuBQL::BinaryBVH<float, 3> bvhA;
  cuBQL::ext::grid_forest::forest_builder::build_forest<float, 3>(
      bvhA, d_boxesA, currentPrimsNumA, numTrianglesA, totalActiveCellsA, d_sortedPrimIDsA, d_nodeOffsetsA, buildConfig,
      thrust::raw_pointer_cast(d_nodeDescendantCountsA.data()), stream, memResource);
  cuBQL::cuda_forest::refit_forest<float, 3>(bvhA, d_boxesA, totalActiveCellsA, stream, memResource);
  CUBQL_CUDA_CALL(StreamSynchronize(stream));
  stats.buildRefitMeshAMs += (cuBQL::getCurrentTime() - tForestAStart) * 1000.0;

  // Build Sub-BVH Forest for Mesh B
  const double tForestBStart = cuBQL::getCurrentTime();
  cuBQL::BinaryBVH<float, 3> bvhB;
  cuBQL::ext::grid_forest::forest_builder::build_forest<float, 3>(
      bvhB, d_boxesB, currentPrimsNumB, numTrianglesB, totalActiveCellsB, d_sortedPrimIDsB, d_nodeOffsetsB, buildConfig,
      thrust::raw_pointer_cast(d_nodeDescendantCountsB.data()), stream, memResource);
  cuBQL::cuda_forest::refit_forest<float, 3>(bvhB, d_boxesB, totalActiveCellsB, stream, memResource);
  CUBQL_CUDA_CALL(StreamSynchronize(stream));
  stats.buildRefitMeshBMs += (cuBQL::getCurrentTime() - tForestBStart) * 1000.0;

  // ============================================================================
  // 7. DUAL-TREE REFINEMENT & RAPID DESCENT BFS
  // ============================================================================
  uint32_t finalActiveCellsA = totalActiveCellsA;
  uint32_t finalActiveCellsB = totalActiveCellsB;

  thrust::device_vector<uint32_t> d_markedNodeIndicesA(finalActiveCellsA);
  thrust::sequence(thrust::device.on(stream), d_markedNodeIndicesA.begin(), d_markedNodeIndicesA.end());

  thrust::device_vector<uint32_t> d_markedNodeIndicesB(finalActiveCellsB);
  thrust::sequence(thrust::device.on(stream), d_markedNodeIndicesB.begin(), d_markedNodeIndicesB.end());

  thrust::device_vector<uint32_t> d_reverseMapB(bvhB.numNodes, 0);

  const double tDualStepStart = cuBQL::getCurrentTime();
  if(numberOfDualTreeSteps > 0) {
    executeDualTreeStep(numberOfDualTreeSteps, maxCellSizeA, maxCellSizeB, d_outPairsA, d_outPairsB, d_markedNodeIndicesA,
                        d_markedNodeIndicesB, d_nodeDescendantCountsA, d_nodeDescendantCountsB, finalActiveCellsA,
                        finalActiveCellsB, bvhA, bvhB);
  }

  CUBQL_CUDA_CALL(StreamSynchronize(stream));
  stats.dualTreeStepMs = (numberOfDualTreeSteps > 0) ? (cuBQL::getCurrentTime() - tDualStepStart) * 1000.0 : 0.0;

  const int totalBatches = d_outPairsA.size();
  const double tGpuBfsStart = cuBQL::getCurrentTime();

  executeRapidDescentBFS(bvhB, currentPrimsNumB, d_markedNodeIndicesB, d_nodeDescendantCountsB, finalActiveCellsB,
                         d_outOffsetsB, d_outPrimsFlatB);

  if(finalActiveCellsB > 0) {
    const int blockSize = 256;
    const int gridSize = (finalActiveCellsB + blockSize - 1) / blockSize;
    populateReverseMapBKernel<<<gridSize, blockSize, 0, stream>>>(thrust::raw_pointer_cast(d_reverseMapB.data()),
                                                                  thrust::raw_pointer_cast(d_markedNodeIndicesB.data()),
                                                                  finalActiveCellsB);
  }
  CUBQL_CUDA_CALL(StreamSynchronize(stream));
  const double tGpuBfsEnd = cuBQL::getCurrentTime();

  // ============================================================================
  // 8. BATCHED SINGLE-TREE TRAVERSAL & EXACT GPU PREDICATES
  // ============================================================================
  IntersectionTimeTracker tracker;
  const uint64_t finalCandidatePairs = executeSingleTreeBatchedTraversalWithPredicates(
      meshA, meshB, batchMultiplier, totalBatches, d_outPairsA, d_outPairsB, d_reverseMapB, d_markedNodeIndicesB,
      d_outOffsetsB, d_outPrimsFlatB, d_nodeDescendantCountsB, finalActiveCellsB, bvhA, d_boxesA, d_boxesB, d_vertsA,
      d_indicesA, d_vertsB, d_indicesB, outFinalExactPairs, outFinalCount, tracker, centerA, centerB, rotA, transA,
      rotB, transB, gpuPredicateMode, stream);

  // ============================================================================
  // 9. RESOURCE CLEANUP & STATS AGGREGATION
  // ============================================================================
  const double tCleanupStart = cuBQL::getCurrentTime();

  _FREE(d_boxesA, stream, memResource);
  _FREE(d_boxesB, stream, memResource);
  _FREE(d_vertsA, stream, memResource);
  _FREE(d_indicesA, stream, memResource);
  _FREE(d_vertsB, stream, memResource);
  _FREE(d_indicesB, stream, memResource);

  if(d_nodeBoxesA)
    _FREE(d_nodeBoxesA, stream, memResource);
  if(d_sortedPrimIDsA)
    _FREE(d_sortedPrimIDsA, stream, memResource);
  if(d_nodeOffsetsA)
    _FREE(d_nodeOffsetsA, stream, memResource);
  if(d_nodeBoxesB)
    _FREE(d_nodeBoxesB, stream, memResource);
  if(d_sortedPrimIDsB)
    _FREE(d_sortedPrimIDsB, stream, memResource);
  if(d_nodeOffsetsB)
    _FREE(d_nodeOffsetsB, stream, memResource);

  cuBQL::cuda::free(bvhA, stream, memResource);
  cuBQL::cuda::free(bvhB, stream, memResource);

  CUBQL_CUDA_CALL(StreamSynchronize(stream));
  stats.finalCleanupSyncMs = (cuBQL::getCurrentTime() - tCleanupStart) * 1000.0;

  // Record Pipeline Execution Statistics
  stats.meshATotalNodes = bvhA.numNodes;
  stats.meshAExtractedTargets = initialMarkedCountA;
  stats.meshBTotalNodes = bvhB.numNodes;
  stats.meshBExtractedTargets = initialMarkedCountB;
  stats.totalIntersections = totalIntersections;
  stats.totalPossiblePairs = totalPossiblePairs;
  stats.intersectionPercentage = intersectionPercentage;
  stats.gpuCrossCheckEngineMs = stats.gpuCrossCheckEngineMs;
  stats.parallelDfsDescentBMs = (tGpuBfsEnd - tGpuBfsStart) * 1000.0;
  stats.GPUTotalTime = (cuBQL::getCurrentTime() - tPipelineStart) * 1000.0;
  stats.totalCrissCrossBatches = totalBatches;
  stats.finalAabbCandidatePairs = finalCandidatePairs;
  stats.confirmedGreenPairs = stats.loopTracker.confirmedGreenPairs;
  stats.confirmedYellowPairs = stats.loopTracker.confirmedYellowPairs;
  stats.loopTracker = tracker;
}