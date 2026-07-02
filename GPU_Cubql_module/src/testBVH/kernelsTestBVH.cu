#define CUBQL_GPU_BUILDER_IMPLEMENTATION 1



// Core cuBQL headers
#include "cuBQL/bvh.h"
#include "cuBQL/builder/cuda.h"

// Your custom experimental builder layout
#include "include/third-party/cubql/sm_builder_v2_2.h"
#include "include/third-party/cubql/sm_builder_v3.h"
#include "include/third-party/cubql/sm_builder_v4.h"

// Custom traversal 
#include "include/third-party/cubql/fixedBoxQueryv2.h"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/scan.h>
#include <thrust/execution_policy.h>
#include <thrust/system/cuda/execution_policy.h>

#include <vector>
#include <algorithm>
#include <iostream>
#include "samples/common/loadOBJ.h"

// Include modular execution targets
#include "crossCheck.h"
#include "DualTreeStep.h" 
#include "rapidDescendKernel.h"
#include "batchedCrossIntersection.h"
#include "crossCheckNew.h"

// Include the updated header matching this implementation
#include "kernelsTestBVH.h"

#include "batchedCrossIntersectionBrute.h"

#include "batchedCrossIntersectionBruteV2.h"




// --------------------------------------------------------------------
// LOCAL STREAM-SAFE REPLACEMENTS FOR INTERNAL CUBQL ALLOCATORS
// --------------------------------------------------------------------
#ifndef _ALLOC
#define _ALLOC(ptr, count, stream, memResource) \
    CUBQL_CUDA_CALL(MallocAsync((void**)&(ptr), (size_t)(count) * sizeof(*(ptr)), stream))
#endif

#ifndef _FREE
#define _FREE(ptr, stream, memResource) \
    CUBQL_CUDA_CALL(FreeAsync((ptr), stream))
#endif

// --------------------------------------------------------------------
// EXISTING KERNELS & HELPERS
// --------------------------------------------------------------------
__global__ void generateBoxes(cuBQL::box3f* boxes, const cuBQL::Triangle* tris, int N) {
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if(i < N) { boxes[i] = tris[i].bounds(); }
}

__global__ void populateReverseMapBKernel(uint32_t* d_reverseMapB, const uint32_t* d_markedNodeIndicesB, uint32_t h_outMarkedCountB) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx < h_outMarkedCountB) {
    uint32_t directBvhNodeId = d_markedNodeIndicesB[idx];
    d_reverseMapB[directBvhNodeId] = idx;
  }
}


template<typename T, int D>
void testFastBVH(const cuBQL::box_t<T,D>* d_boxesA, int numPrimsA, int maxCellSizeA,
                 const cuBQL::box_t<T,D>* d_boxesB, int numPrimsB, int maxCellSizeB,
                 cuBQL::box_t<T,D> globalBoxA, cuBQL::box_t<T,D> globalBoxB,
                 cudaStream_t s, cuBQL::DeviceMemoryResource& memResource, int batchMultiplier , int userleafThreshold)
{
  std::cout << "\n==================================================" << std::endl;
  std::cout << " [DEEP DEBUG] => Entering ZIMBLIFIED testFastBVH..." << std::endl;
  std::cout << "==================================================" << std::endl;

  // Check raw pointers and parameters right away
  std::cout << " -> Device Pointer A : " << d_boxesA << " | Prims: " << numPrimsA << " | MaxCell: " << maxCellSizeA << std::endl;
  std::cout << " -> Device Pointer B : " << d_boxesB << " | Prims: " << numPrimsB << " | MaxCell: " << maxCellSizeB << std::endl;
  
  if (d_boxesA == d_boxesB) {
      std::cout << " [CRITICAL WARNING] d_boxesA and d_boxesB point to the EXACT SAME memory address!" << std::endl;
  } else {
      std::cout << " [OK] Input pointers A and B are distinct memory blocks." << std::endl;
  }

  if (maxCellSizeA > 20 || maxCellSizeB > 20) {
      std::cerr << "[CRITICAL ERROR] numIterations is too high! Aborting testFastBVH." << std::endl;
      return;
  }

  double tStart = cuBQL::getCurrentTime();

  // --- MESH A RUN (Zimblified Edition) ---
  cuBQL::box_t<T, D>* outNodeBoxesA = nullptr;
  uint32_t* outSortedPrimIDsA = nullptr;
  uint32_t* outNodeOffsetsA = nullptr;
  uint32_t outTotalActiveCellsA = 0;

  // Profile Mesh A Initialization
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double tInitAStart = cuBQL::getCurrentTime();

  cuBQL::gpuBuilder_v3::test_speedrun_initialization_linear(
      d_boxesA, numPrimsA, (uint32_t)maxCellSizeA, globalBoxA,
      outNodeBoxesA, outSortedPrimIDsA, outNodeOffsetsA, outTotalActiveCellsA,
      s, memResource
  );
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double tInitAMs = (cuBQL::getCurrentTime() - tInitAStart) * 1000.0;

  // --- MESH B RUN (Zimblified Edition) ---
  cuBQL::box_t<T, D>* outNodeBoxesB = nullptr;
  uint32_t* outSortedPrimIDsB = nullptr;
  uint32_t* outNodeOffsetsB = nullptr;
  uint32_t outTotalActiveCellsB = 0;

  // Profile Mesh B Initialization
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double tInitBStart = cuBQL::getCurrentTime();

  cuBQL::gpuBuilder_v3::test_speedrun_initialization_linear(
      d_boxesB, numPrimsB, (uint32_t)maxCellSizeB, globalBoxB,
      outNodeBoxesB, outSortedPrimIDsB, outNodeOffsetsB, outTotalActiveCellsB,
      s, memResource
  );
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double tInitBMs = (cuBQL::getCurrentTime() - tInitBStart) * 1000.0;

  double tEnd = cuBQL::getCurrentTime();
  double elapsedMs = (tEnd - tStart) * 1000.0;

  // --------------------------------------------------------------------
  // VERIFY BUILDING METRICS READBACK (No structural junk, 100% Active)
  // --------------------------------------------------------------------
  std::cout << "\n--------------------------------------------------" << std::endl;
  std::cout << " [LINEAR INITIALIZATION METRICS OVERVIEW]" << std::endl;
  std::cout << " -> Tree A: Total Compounded Active Cells = " << outTotalActiveCellsA << std::endl;
  std::cout << " -> Tree B: Total Compounded Active Cells = " << outTotalActiveCellsB << std::endl;
  std::cout << "--------------------------------------------------" << std::endl;

  // Read back only the active, valid bounding boxes
  std::vector<cuBQL::box_t<T,D>> h_boxesA(outTotalActiveCellsA);
  std::vector<cuBQL::box_t<T,D>> h_boxesB(outTotalActiveCellsB);
  std::vector<uint32_t> h_offsetsA(outTotalActiveCellsA + 1);
  std::vector<uint32_t> h_offsetsB(outTotalActiveCellsB + 1);

  CUBQL_CUDA_CALL(Memcpy(h_boxesA.data(), outNodeBoxesA, outTotalActiveCellsA * sizeof(cuBQL::box_t<T,D>), cudaMemcpyDeviceToHost));
  CUBQL_CUDA_CALL(Memcpy(h_boxesB.data(), outNodeBoxesB, outTotalActiveCellsB * sizeof(cuBQL::box_t<T,D>), cudaMemcpyDeviceToHost));
  CUBQL_CUDA_CALL(Memcpy(h_offsetsA.data(), outNodeOffsetsA, (outTotalActiveCellsA + 1) * sizeof(uint32_t), cudaMemcpyDeviceToHost));
  CUBQL_CUDA_CALL(Memcpy(h_offsetsB.data(), outNodeOffsetsB, (outTotalActiveCellsB + 1) * sizeof(uint32_t), cudaMemcpyDeviceToHost));

  // Calculate true primitive allocations via offsets
  uint32_t sumPrimsA = h_offsetsA[outTotalActiveCellsA] - h_offsetsA[0];
  uint32_t sumPrimsB = h_offsetsB[outTotalActiveCellsB] - h_offsetsB[0];

  std::cout << "--------------------------------------------------" << std::endl;
  std::cout << " [POPULATED CELLS CHECK SUMMARY]" << std::endl;
  std::cout << " -> Tree A Active Nodes: " << outTotalActiveCellsA << " | Total Verified Prims = " << sumPrimsA << std::endl;
  std::cout << " -> Tree B Active Nodes: " << outTotalActiveCellsB << " | Total Verified Prims = " << sumPrimsB << std::endl;
  std::cout << "--------------------------------------------------\n" << std::endl;

  // --------------------------------------------------------------------
  // PARALLEL TEMP NODE CROSS CHECK
  // --------------------------------------------------------------------
  std::cout << "\n==================================================" << std::endl;
  std::cout << "          LAUNCHING CROSS-CHECK KERNELS           " << std::endl;
  std::cout << "==================================================" << std::endl;

  thrust::device_vector<uint32_t> d_intersectPairsA;
  thrust::device_vector<uint32_t> d_intersectPairsB;

  double crossCheckStart = cuBQL::getCurrentTime();

  uint32_t totalIntersections = executeBoxCrossCheck<T,D>(
      outNodeBoxesA, outTotalActiveCellsA,
      outNodeBoxesB, outTotalActiveCellsB,
      d_intersectPairsA, d_intersectPairsB, s
  );
  CUBQL_CUDA_CALL(StreamSynchronize(s));

  double crossCheckEnd = cuBQL::getCurrentTime();
  double crossCheckMs = (crossCheckEnd - crossCheckStart) * 1000.0;
  
  // --------------------------------------------------------------------
  // NEW: GPU BUILDER V4 FOREST EXPANSION & PERFORMANCE BENCHMARK
  // --------------------------------------------------------------------
  std::cout << " [FOREST EXPANSION] => Launching Level-by-Level Parallel Sub-Tree Compilation..." << std::endl;

  // 1. Calculate EXACT maxNodes using the same formula as build_forest to prevent out-of-bounds Memset
  const uint32_t numInitA = outTotalActiveCellsA + (outTotalActiveCellsA & 1u);
  const uint32_t maxNodesA = 2u * (uint32_t)numPrimsA + numInitA + 2u;

  const uint32_t numInitB = outTotalActiveCellsB + (outTotalActiveCellsB & 1u);
  const uint32_t maxNodesB = 2u * (uint32_t)numPrimsB + numInitB + 2u;

  // 2. Safely allocate target-level diagnostics with the correct uniform sizing
  thrust::device_vector<uint32_t> d_nodeDescendantCountsA(maxNodesA, 0);
  thrust::device_vector<uint32_t> d_nodeDescendantCountsB(maxNodesB, 0);

  // 3. Configure build options explicitly to avoid uninitialized defaults
  cuBQL::BuildConfig buildConfig; 
  buildConfig.makeLeafThreshold = userleafThreshold; // Set explicit leaf configuration threshold

  // --- Compile Forest A ---
  cuBQL::BinaryBVH<T, D> bvhA;
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double forestAStart = cuBQL::getCurrentTime();

  cuBQL::gpuBuilder_v4::build_forest<T, D>(
      bvhA, d_boxesA, numPrimsA, outTotalActiveCellsA,
      outSortedPrimIDsA, outNodeOffsetsA, buildConfig, thrust::raw_pointer_cast(d_nodeDescendantCountsA.data()), s, memResource
  );
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double forestAMs = (cuBQL::getCurrentTime() - forestAStart) * 1000.0;

  // --- Compile Forest B ---
  cuBQL::BinaryBVH<T, D> bvhB;
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double forestBStart = cuBQL::getCurrentTime();

  cuBQL::gpuBuilder_v4::build_forest<T, D>(
      bvhB, d_boxesB, numPrimsB, outTotalActiveCellsB,
      outSortedPrimIDsB, outNodeOffsetsB, buildConfig, thrust::raw_pointer_cast(d_nodeDescendantCountsB.data()), s, memResource
  );
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double forestBMs = (cuBQL::getCurrentTime() - forestBStart) * 1000.0;

  // --------------------------------------------------------------------
  // FINAL PERFORMANCE & RESULTS READOUT
  // --------------------------------------------------------------------
  uint64_t totalValidActivePairs = (uint64_t)outTotalActiveCellsA * outTotalActiveCellsB;
  double overlapPercentage = 0.0;
  if (totalValidActivePairs > 0) {
      overlapPercentage = (double)totalIntersections / (double)totalValidActivePairs * 100.0;
  }

  std::cout << "\n==================================================" << std::endl;
  std::cout << "               FINAL BENCHMARK SUMMARY            " << std::endl;
  std::cout << "==================================================" << std::endl;
  std::cout << " -> Mesh A Speedrun Init Speed      : " << tInitAMs << " ms" << std::endl;
  std::cout << " -> Mesh B Speedrun Init Speed      : " << tInitBMs << " ms" << std::endl;
  std::cout << " -> Cell Cross-Check Execution      : " << crossCheckMs << " ms" << std::endl;
  std::cout << " ------------------------------------------------" << std::endl;
  std::cout << " -> Mesh A Forest BVH Expansion     : " << forestAMs << " ms | Nodes: " << bvhA.numNodes << std::endl;
  std::cout << " -> Mesh B Forest BVH Expansion     : " << forestBMs << " ms | Nodes: " << bvhB.numNodes << std::endl;
  std::cout << " ------------------------------------------------" << std::endl;
  std::cout << " -> Intersecting Nodes Detected     : " << totalIntersections << " pairs" << std::endl;
  std::cout << " -> True Active Matrix Evaluated    : " << totalValidActivePairs << " pairings" << std::endl;
  std::cout << " -> Percentage of Node Overlaps     : " << overlapPercentage << "%" << std::endl;
  std::cout << "==================================================\n" << std::endl;

  if (outNodeBoxesA)     _FREE(outNodeBoxesA, s, memResource);
  if (outSortedPrimIDsA) _FREE(outSortedPrimIDsA, s, memResource);
  if (outNodeOffsetsA)   _FREE(outNodeOffsetsA, s, memResource);
  if (outNodeBoxesB)     _FREE(outNodeBoxesB, s, memResource);
  if (outSortedPrimIDsB) _FREE(outSortedPrimIDsB, s, memResource);
  if (outNodeOffsetsB)   _FREE(outNodeOffsetsB, s, memResource);
}
// --------------------------------------------------------------------
// MAIN ENTRY TESTING PIPELINE
// --------------------------------------------------------------------
extern "C" void kernelsTestBVH(const cuBQL::Triangle* hMeshA, int numTrianglesA, int maxCellSizeA,
                               const cuBQL::Triangle* hMeshB, int numTrianglesB, int maxCellSizeB,
                               int batchMultiplier,
                               int mode, 
                               int leafThreshold,
                               ExecutionStats& stats, 
                               std::vector<int2>& hGreenPairs,  
                               std::vector<int2>& hYellowPairs)
{
  if(numTrianglesA <= 0 || numTrianglesB <= 0) {
    return;
  }

  double tPipelineStart = cuBQL::getCurrentTime();
  cudaStream_t stream = 0;
  cuBQL::DeviceMemoryResource memResource;
  cuBQL::BuildConfig buildConfig;
  buildConfig.makeLeafThreshold = leafThreshold;

  // --------------------------------------------------------------------
  // INITIAL RAW ALLOCATION & COPY TRACKING
  // --------------------------------------------------------------------
  double tAllocStart = cuBQL::getCurrentTime();

  cuBQL::Triangle* dMeshA;
  CUBQL_CUDA_CALL(Malloc(&dMeshA, numTrianglesA * sizeof(cuBQL::Triangle)));
  CUBQL_CUDA_CALL(Memcpy(dMeshA, hMeshA, numTrianglesA * sizeof(cuBQL::Triangle), cudaMemcpyDefault));

  cuBQL::Triangle* dMeshB;
  CUBQL_CUDA_CALL(Malloc(&dMeshB, numTrianglesB * sizeof(cuBQL::Triangle)));
  CUBQL_CUDA_CALL(Memcpy(dMeshB, hMeshB, numTrianglesB * sizeof(cuBQL::Triangle), cudaMemcpyDefault));

  cudaDeviceSynchronize();
  double tAllocEnd = cuBQL::getCurrentTime();
  stats.initialAllocAndCopyMs = (tAllocEnd - tAllocStart) * 1000.0;

  // --------------------------------------------------------------------
  // THRUST DEVICE ALLOCATION / INITIALIZATION OVERHEAD TRACKING
  // --------------------------------------------------------------------
  double tThrustInitStart = cuBQL::getCurrentTime();

  uint32_t maxPossibleNodesA = 2 * numTrianglesA;
  uint32_t maxPossibleNodesB = 2 * numTrianglesB;

  thrust::device_vector<uint32_t> d_markedNodeIndicesA(maxPossibleNodesA, 0);
  thrust::device_vector<uint32_t> d_nodeDescendantCountsA(maxPossibleNodesA, 0); 
  thrust::device_vector<uint32_t> d_markedNodeIndicesB(maxPossibleNodesB, 0);
  thrust::device_vector<uint32_t> d_nodeDescendantCountsB(maxPossibleNodesB, 0); 

  thrust::device_vector<uint32_t> d_reverseMapB(maxPossibleNodesB, 0);
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
  // HIERARCHY A: BUILD & REFIT
  // --------------------------------------------------------------------
  double tBuildAStart = cuBQL::getCurrentTime();
  cuBQL::box3f* dBoxesA;
  CUBQL_CUDA_CALL(Malloc(&dBoxesA, numTrianglesA * sizeof(cuBQL::box3f)));
  generateBoxes<<<cuBQL::divRoundUp(numTrianglesA, 256), 256>>>(dBoxesA, dMeshA, numTrianglesA);

  uint32_t h_outMarkedCountA = 0;
  cuBQL::bvh3f bvhA;
  cuBQL::gpuBuilder_v2_2::build_custom(
      bvhA, dBoxesA, numTrianglesA, buildConfig, (uint32_t)maxCellSizeA, thrust::raw_pointer_cast(d_markedNodeIndicesA.data()),
      thrust::raw_pointer_cast(d_nodeDescendantCountsA.data()), &h_outMarkedCountA, stream, memResource);
  cuBQL::cuda::refit(bvhA, dBoxesA, stream);
  
  cudaDeviceSynchronize();
  double tBuildAEnd = cuBQL::getCurrentTime();

  // --------------------------------------------------------------------
  // HIERARCHY B: BUILD & REFIT
  // --------------------------------------------------------------------
  double tBuildBStart = cuBQL::getCurrentTime();
  cuBQL::box3f* dBoxesB;
  CUBQL_CUDA_CALL(Malloc(&dBoxesB, numTrianglesB * sizeof(cuBQL::box3f)));
  generateBoxes<<<cuBQL::divRoundUp(numTrianglesB, 256), 256>>>(dBoxesB, dMeshB, numTrianglesB);

  uint32_t h_outMarkedCountB = 0;
  cuBQL::bvh3f bvhB;
  cuBQL::gpuBuilder_v2_2::build_custom(
      bvhB, dBoxesB, numTrianglesB, buildConfig, (uint32_t)maxCellSizeB, thrust::raw_pointer_cast(d_markedNodeIndicesB.data()),
      thrust::raw_pointer_cast(d_nodeDescendantCountsB.data()), &h_outMarkedCountB, stream, memResource);
  cuBQL::cuda::refit(bvhB, dBoxesB, stream);
  
  cudaDeviceSynchronize();
  double tBuildBEnd = cuBQL::getCurrentTime();

  // --------------------------------------------------------------------
  // DISPATCH EXPERIMENTAL INJECTED V3 WORKFLOW TEST (FIXED)
  // --------------------------------------------------------------------
//   cuBQL::box3f hostBoxA, hostBoxB;
//   CUBQL_CUDA_CALL(Memcpy(&hostBoxA, &(bvhA.nodes[0].bounds), sizeof(cuBQL::box3f), cudaMemcpyDeviceToHost));
//   CUBQL_CUDA_CALL(Memcpy(&hostBoxB, &(bvhB.nodes[0].bounds), sizeof(cuBQL::box3f), cudaMemcpyDeviceToHost));
//   cudaDeviceSynchronize();

//   testFastBVH<float, 3>(
//       dBoxesA, numTrianglesA, maxCellSizeA,
//       dBoxesB, numTrianglesB, maxCellSizeB,
//       hostBoxA, hostBoxB,
//       stream, memResource, batchMultiplier, leafThreshold
//   );

  // --------------------------------------------------------------------
  // EXTRACTED GPU PARALLEL CRISS-CROSS INTERSECTION MODULE
  // --------------------------------------------------------------------
  double tCrossStart = cuBQL::getCurrentTime();

  uint32_t totalIntersections = executeCrissCrossIntersection(
      bvhA, d_markedNodeIndicesA, h_outMarkedCountA,
      bvhB, d_markedNodeIndicesB, h_outMarkedCountB,
      d_outPairsA, d_outPairsB
  );

  cudaDeviceSynchronize(); 
  double tCrossEnd = cuBQL::getCurrentTime();

  uint64_t totalPossiblePairs = (uint64_t)h_outMarkedCountA * h_outMarkedCountB;
  double intersectionPercentage = totalPossiblePairs > 0 ? ((double)totalIntersections / totalPossiblePairs) * 100.0 : 0.0;

  // --------------------------------------------------------------------
  // DUAL-TREE TRAVERSAL STEP OVERHEAD PHASE
  // --------------------------------------------------------------------
  double tDualStepStart = cuBQL::getCurrentTime();
  
  if (mode > 0) {
    executeDualTreeStep(
        mode, maxCellSizeA, maxCellSizeB, d_outPairsA, d_outPairsB,
        d_markedNodeIndicesA, d_markedNodeIndicesB, d_nodeDescendantCountsA, d_nodeDescendantCountsB,
        h_outMarkedCountA, h_outMarkedCountB, bvhA, bvhB, dMeshA, dMeshB
    );
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
  executeRapidDescentBFS(bvhB, numTrianglesB, d_markedNodeIndicesB, d_nodeDescendantCountsB, h_outMarkedCountB, d_outOffsetsB, d_outPrimsFlatB);

  if (h_outMarkedCountB > 0) {
    int blockSize = 256;
    int gridSize = (h_outMarkedCountB + blockSize - 1) / blockSize;
    populateReverseMapBKernel<<<gridSize, blockSize, 0, stream>>>(
        thrust::raw_pointer_cast(d_reverseMapB.data()), thrust::raw_pointer_cast(d_markedNodeIndicesB.data()), h_outMarkedCountB
    );
  }
  
  cudaDeviceSynchronize();
  double tGpuBfsEnd = cuBQL::getCurrentTime();

  finalCandidatePairs = executeBatchedCrossIntersectionLoop(
      batchMultiplier, totalBatches, d_outPairsA, d_outPairsB, d_reverseMapB,             
      d_markedNodeIndicesB, d_outOffsetsB, d_outPrimsFlatB, d_nodeDescendantCountsB,
      h_outMarkedCountB, bvhA, dMeshA, dMeshB, hGreenPairs, hYellowPairs, tracker        
  );

  // --------------------------------------------------------------------
  // EXPLICIT CLEANUP & RECOVERY METRIC TRACKING
  // --------------------------------------------------------------------
  double tCleanupStart = cuBQL::getCurrentTime();

  CUBQL_CUDA_CALL(Free(dMeshA));
  CUBQL_CUDA_CALL(Free(dBoxesA));
  CUBQL_CUDA_CALL(Free(dMeshB));
  CUBQL_CUDA_CALL(Free(dBoxesB));
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

  stats.meshATotalNodes         = bvhA.numNodes;
  stats.meshAExtractedTargets   = h_outMarkedCountA;
  stats.meshBTotalNodes         = bvhB.numNodes;
  stats.meshBExtractedTargets   = h_outMarkedCountB;
  stats.totalIntersections     = totalIntersections;
  stats.totalPossiblePairs     = totalPossiblePairs;
  stats.intersectionPercentage = intersectionPercentage;
  stats.buildRefitMeshAMs     = (tBuildAEnd - tBuildAStart) * 1000.0;
  stats.buildRefitMeshBMs     = (tBuildBEnd - tBuildBStart) * 1000.0;
  stats.gpuCrossCheckEngineMs = (tCrossEnd - tCrossStart) * 1000.0;
  stats.parallelDfsDescentBMs = (tGpuBfsEnd - tGpuBfsStart) * 1000.0; 
  stats.GPUTotalTime          = (tPipelineEnd - tPipelineStart) * 1000.0;
  stats.totalCrissCrossBatches  = totalBatches;
  stats.finalAabbCandidatePairs = finalCandidatePairs;
  stats.confirmedGreenPairs     = hGreenPairs.size();
  stats.confirmedYellowPairs    = hYellowPairs.size();
  stats.loopTracker             = tracker;
}