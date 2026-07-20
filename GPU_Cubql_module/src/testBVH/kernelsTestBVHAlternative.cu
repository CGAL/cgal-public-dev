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
#include "batchedCrossIntersectionV2.h"
#include "batchedCrossIntersectionV3.h"
#include "crossCheckNew.h"
#include "kernelsTestBVHAlternative.h"
#include "prune_pipeline.h"

#include "global_box.h"

#ifndef _ALLOC
#define _ALLOC(ptr, count, stream, memResource)                                                                        \
  CUBQL_CUDA_CALL(MallocAsync((void**)&(ptr), (size_t)(count) * sizeof(*(ptr)), stream))
#endif

#ifndef _FREE
#define _FREE(ptr, stream, memResource) CUBQL_CUDA_CALL(FreeAsync((ptr), stream))
#endif


// UPGRADED: Instantly builds triangles AND evaluates tight precision error bounds
__global__ void assembleTrianglesKernel(cuBQL::Triangle* dMesh,
                                        float2* dMetrics,
                                        const float3* dVerts,
                                        const float* dVertErrors, // Can now be nullptr
                                        const uint3* dIndices,
                                        int numTriangles) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if(idx >= numTriangles)
    return;

  uint3 triIdx = dIndices[idx];

  float3 p0 = dVerts[triIdx.x];
  float3 p1 = dVerts[triIdx.y];
  float3 p2 = dVerts[triIdx.z];

  // 1. Map components directly into cuBQL fields
  dMesh[idx].a.x = p0.x;
  dMesh[idx].a.y = p0.y;
  dMesh[idx].a.z = p0.z;
  dMesh[idx].b.x = p1.x;
  dMesh[idx].b.y = p1.y;
  dMesh[idx].b.z = p1.z;
  dMesh[idx].c.x = p2.x;
  dMesh[idx].c.y = p2.y;
  dMesh[idx].c.z = p2.z;

  // 2. Compute L (Absolute Spatial Scale)
  float l0 = fmaxf(fmaxf(fabsf(p0.x), fabsf(p0.y)), fabsf(p0.z));
  float l1 = fmaxf(fmaxf(fabsf(p1.x), fabsf(p1.y)), fabsf(p1.z));
  float l2 = fmaxf(fmaxf(fabsf(p2.x), fabsf(p2.y)), fabsf(p2.z));
  float L = fmaxf(fmaxf(l0, l1), l2);

  // 3. Compute E (Maximum Edge Axis Extent)
  float ex = fmaxf(fmaxf(fabsf(p0.x - p1.x), fabsf(p1.x - p2.x)), fabsf(p2.x - p0.x));
  float ey = fmaxf(fmaxf(fabsf(p0.y - p1.y), fabsf(p1.y - p2.y)),
                   fabsf(p2.y - p0.y)); // corrected typographical index '1.y' to 'p1.y'
  float ez = fmaxf(fmaxf(fabsf(p0.z - p1.z), fabsf(p1.z - p2.z)), fabsf(p2.z - p0.z));
  float E = fmaxf(fmaxf(ex, ey), ez);

  // 4. Absorb Host-Generated Rounding Error Bounds (or compute them dynamically)
  float e0, e1, e2;
  if(dVertErrors != nullptr) {
    e0 = dVertErrors[triIdx.x];
    e1 = dVertErrors[triIdx.y];
    e2 = dVertErrors[triIdx.z];
  } else {
    // 2^-24 upper bound constant (conservative rounding to 5.9604645e-8f)
    const float EPS_BOUND = 5.9604645e-8f;
    e0 = l0 * EPS_BOUND;
    e1 = l1 * EPS_BOUND;
    e2 = l2 * EPS_BOUND;
  }
  float maxVertexError = fmaxf(fmaxf(e0, e1), e2);

  // Pad spatial parameters conservatively to encompass downcast drift
  L += maxVertexError;
  E += (2.0f * maxVertexError);

  // 5. Write out packed parameters to companion buffer
  dMetrics[idx] = make_float2(L, E * E);
}

__global__ void generateBoxes(cuBQL::box3f* boxes, const cuBQL::Triangle* tris, int N) {
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if(i < N) {
    boxes[i] = tris[i].bounds();
  }
}

__global__ void
populateReverseMapBKernel(uint32_t* d_reverseMapB, const uint32_t* d_markedNodeIndicesB, uint32_t h_outMarkedCountB) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if(idx < h_outMarkedCountB) {
    uint32_t directBvhNodeId = d_markedNodeIndicesB[idx];
    d_reverseMapB[directBvhNodeId] = idx;
  }
}

extern "C" void kernelsTestBVHV2(Mesh& meshAcpu,
                                 Mesh& meshBcpu,
                                 const float3* hVertsA,
                                 int numVertsA,
                                 const uint3* hIndicesA,
                                 const float* hVertErrorsA,
                                 int numTrianglesA,
                                 int maxCellSizeA,
                                 const float3* hVertsB,
                                 int numVertsB,
                                 const uint3* hIndicesB,
                                 const float* hVertErrorsB,
                                 int numTrianglesB,
                                 int maxCellSizeB,
                                 int batchMultiplier,
                                 int mode,
                                 int leafThreshold,
                                 int activateAsyncDownload,
                                 ExecutionStats& stats,
                                 tbb::concurrent_vector<int2>& finalExactPairs) {
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
  float3* dVertsA = nullptr;
  uint3* dIndicesA = nullptr;
  CUBQL_CUDA_CALL(Malloc(&dVertsA, numVertsA * sizeof(float3)));
  CUBQL_CUDA_CALL(Memcpy(dVertsA, hVertsA, numVertsA * sizeof(float3), cudaMemcpyHostToDevice));
  CUBQL_CUDA_CALL(Malloc(&dIndicesA, numTrianglesA * sizeof(uint3)));
  CUBQL_CUDA_CALL(Memcpy(dIndicesA, hIndicesA, numTrianglesA * sizeof(uint3), cudaMemcpyHostToDevice));

  // Allocate and stream indexed buffers for Mesh B
  float3* dVertsB = nullptr;
  uint3* dIndicesB = nullptr;
  CUBQL_CUDA_CALL(Malloc(&dVertsB, numVertsB * sizeof(float3)));
  CUBQL_CUDA_CALL(Memcpy(dVertsB, hVertsB, numVertsB * sizeof(float3), cudaMemcpyHostToDevice));
  CUBQL_CUDA_CALL(Malloc(&dIndicesB, numTrianglesB * sizeof(uint3)));
  CUBQL_CUDA_CALL(Memcpy(dIndicesB, hIndicesB, numTrianglesB * sizeof(uint3), cudaMemcpyHostToDevice));


  // Instantiate standard unrolled triangle buffers on device
  cuBQL::Triangle* dMeshA;
  CUBQL_CUDA_CALL(Malloc(&dMeshA, numTrianglesA * sizeof(cuBQL::Triangle)));

  cuBQL::Triangle* dMeshB;
  CUBQL_CUDA_CALL(Malloc(&dMeshB, numTrianglesB * sizeof(cuBQL::Triangle)));

  float* dVertErrorsA = nullptr;
  if(hVertErrorsA != nullptr) {
    CUBQL_CUDA_CALL(Malloc(&dVertErrorsA, numVertsA * sizeof(float)));
    CUBQL_CUDA_CALL(Memcpy(dVertErrorsA, hVertErrorsA, numVertsA * sizeof(float), cudaMemcpyHostToDevice));
  }

  float* dVertErrorsB = nullptr;
  if(hVertErrorsB != nullptr) {
    CUBQL_CUDA_CALL(Malloc(&dVertErrorsB, numVertsB * sizeof(float)));
    CUBQL_CUDA_CALL(Memcpy(dVertErrorsB, hVertErrorsB, numVertsB * sizeof(float), cudaMemcpyHostToDevice));
  }

  // 2. Allocate the Triangle Metrics Companion Buffers
  float2* dMeshMetricsA = nullptr;
  CUBQL_CUDA_CALL(Malloc(&dMeshMetricsA, numTrianglesA * sizeof(float2)));

  float2* dMeshMetricsB = nullptr;
  CUBQL_CUDA_CALL(Malloc(&dMeshMetricsB, numTrianglesB * sizeof(float2)));

  assembleTrianglesKernel<<<cuBQL::divRoundUp(numTrianglesA, 256), 256, 0, stream>>>(
      dMeshA, dMeshMetricsA, dVertsA, dVertErrorsA, dIndicesA, numTrianglesA);
  assembleTrianglesKernel<<<cuBQL::divRoundUp(numTrianglesB, 256), 256, 0, stream>>>(
      dMeshB, dMeshMetricsB, dVertsB, dVertErrorsB, dIndicesB, numTrianglesB);

  cuBQL::box3f* dBoxesA;
  CUBQL_CUDA_CALL(Malloc(&dBoxesA, numTrianglesA * sizeof(cuBQL::box3f)));
  generateBoxes<<<cuBQL::divRoundUp(numTrianglesA, 256), 256, 0, stream>>>(dBoxesA, dMeshA, numTrianglesA);

  cuBQL::box3f* dBoxesB;
  CUBQL_CUDA_CALL(Malloc(&dBoxesB, numTrianglesB * sizeof(cuBQL::box3f)));
  generateBoxes<<<cuBQL::divRoundUp(numTrianglesB, 256), 256, 0, stream>>>(dBoxesB, dMeshB, numTrianglesB);

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
  thrust::device_vector<uint32_t> d_outOffsetsA;
  thrust::device_vector<uint32_t> d_outPrimsFlatA;

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

  std::cout << "\n--------------------------------------------------" << std::endl;
  std::cout << " [LINEAR INITIALIZATION METRICS OVERVIEW]" << std::endl;
  std::cout << " -> Tree A: Total Compounded Active Cells = " << outTotalActiveCellsA << std::endl;
  std::cout << " -> Tree B: Total Compounded Active Cells = " << outTotalActiveCellsB << std::endl;
  std::cout << "--------------------------------------------------" << std::endl;

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

  std::cout << "\n==================================================" << std::endl;
  std::cout << "          LAUNCHING CROSS-CHECK KERNELS           " << std::endl;
  std::cout << "==================================================" << std::endl;
  std::cout << " -> Intersecting Pairs Detected     : " << totalIntersections << " pairs" << std::endl;
  std::cout << " -> True Active Matrix Evaluated    : " << totalPossiblePairs << " pairings" << std::endl;
  std::cout << " -> Percentage of Total Overlaps    : " << intersectionPercentage << "%" << std::endl;

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

  std::cout << " [POPULATED CELLS CHECK SUMMARY]" << std::endl;
  std::cout << " -> Tree A Active Nodes: " << initialCellsA << " | Total Verified Prims = " << currentPrimsNumA
            << std::endl;
  std::cout << " -> Tree B Active Nodes: " << initialCellsB << " | Total Verified Prims = " << currentPrimsNumB
            << std::endl;
  std::cout << "--------------------------------------------------\n" << std::endl;

  // --------------------------------------------------------------------
  // PARALLEL STREAM COMPACTION & PRUNING ALGORITHM (POST-BUILD)
  // --------------------------------------------------------------------
  std::cout << " \n[PIPELINE PRUNING] => Dispatched grid streaming parallel re-index compaction (POST-BUILD)..."
            << std::endl;

  double pruneStartA = cuBQL::getCurrentTime();
  parallelPruneAndReindexAll(thrust::raw_pointer_cast(d_outPairsA.data()), totalIntersections, outSortedPrimIDsA,
                             outNodeOffsetsA, outTotalActiveCellsA, currentPrimsNumA, stream, memResource);
  CUBQL_CUDA_CALL(StreamSynchronize(stream));
  double pruneMsA = (cuBQL::getCurrentTime() - pruneStartA) * 1000.0;
  std::cout << " -> Mesh A Structural Pruning (Post)  : " << pruneMsA << " ms" << std::endl;
  std::cout << "    * Active Cells Surviving: " << initialCellsA << " -> " << outTotalActiveCellsA << std::endl;

  double pruneStartB = cuBQL::getCurrentTime();
  parallelPruneAndReindexAll(thrust::raw_pointer_cast(d_outPairsB.data()), totalIntersections, outSortedPrimIDsB,
                             outNodeOffsetsB, outTotalActiveCellsB, currentPrimsNumB, stream, memResource);
  CUBQL_CUDA_CALL(StreamSynchronize(stream));
  double pruneMsB = (cuBQL::getCurrentTime() - pruneStartB) * 1000.0;
  std::cout << " -> Mesh B Structural Pruning (Post)  : " << pruneMsB << " ms" << std::endl;
  std::cout << "    * Active Cells Surviving: " << initialCellsB << " -> " << outTotalActiveCellsB << std::endl;


  // --------------------------------------------------------------------
  // GPU BUILDER V4 FOREST EXPANSION
  // --------------------------------------------------------------------
  std::cout << " [FOREST EXPANSION] => Launching Level-by-Level Parallel Sub-Tree Compilation..." << std::endl;

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
  std::cout << " -> Mesh A Forest BVH Expansion     : " << forestAMs << " ms | Nodes: " << bvhA.numNodes << std::endl;

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
  std::cout << " -> Mesh B Forest BVH Expansion     : " << forestBMs << " ms | Nodes: " << bvhB.numNodes << std::endl;


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
                        finalActiveCellsB, bvhA, bvhB, dMeshA, dMeshB);
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

  if(activateAsyncDownload == 0) {

    finalCandidatePairs = executeBatchedCrossIntersectionLoopV2(
        meshAcpu, meshBcpu, batchMultiplier, totalBatches, d_outPairsA, d_outPairsB, d_reverseMapB,
        d_markedNodeIndicesB, d_outOffsetsB, d_outPrimsFlatB, d_nodeDescendantCountsB, finalActiveCellsB, bvhA, dMeshA,
        dMeshB, dMeshMetricsA, dMeshMetricsB, finalExactPairs, tracker, stream);

  } else {
    finalCandidatePairs = executeBatchedCrossIntersectionLoopV3(
        meshAcpu, meshBcpu, batchMultiplier, totalBatches, d_outPairsA, d_outPairsB, d_reverseMapB,
        d_markedNodeIndicesB, d_outOffsetsB, d_outPrimsFlatB, d_nodeDescendantCountsB, finalActiveCellsB, bvhA, dMeshA,
        dMeshB, dMeshMetricsA, dMeshMetricsB, finalExactPairs, tracker, stream);
  }
  // --------------------------------------------------------------------
  // EXPLICIT CLEANUP & RECOVERY METRIC TRACKING
  // --------------------------------------------------------------------

  std::cout << "Starting to clean up" << std::endl;

  double tCleanupStart = cuBQL::getCurrentTime();

  CUBQL_CUDA_CALL(Free(dMeshA));
  CUBQL_CUDA_CALL(Free(dBoxesA));
  CUBQL_CUDA_CALL(Free(dMeshB));
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
  d_outOffsetsA.shrink_to_fit();
  d_outPrimsFlatA.shrink_to_fit();

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