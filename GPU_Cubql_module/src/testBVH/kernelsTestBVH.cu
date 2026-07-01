#define CUBQL_GPU_BUILDER_IMPLEMENTATION 1



// Core cuBQL headers
#include "cuBQL/bvh.h"
#include "cuBQL/builder/cuda.h"

// Your custom experimental builder layout
#include "include/third-party/cubql/sm_builder_v2_2.h"
#include "include/third-party/cubql/sm_builder_v3.h"

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
                 cudaStream_t s, cuBQL::DeviceMemoryResource& memResource)
{
  std::cout << "\n==================================================" << std::endl;
  std::cout << " [DEEP DEBUG] => Entering testFastBVH..." << std::endl;
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

  // --- MESH A RUN ---
  size_t countA = std::max((size_t)(1 << maxCellSizeA) * 4, (size_t)2 * numPrimsA);
  uint32_t* d_dummyCountsA = nullptr;
  _ALLOC(d_dummyCountsA, countA, s, memResource);

  cuBQL::gpuBuilder_v3::FinalNormalNode<T,D>* outTempNodesA = nullptr;
  cuBQL::gpuBuilder_v3::PrimState* outPrimStatesA = nullptr;
  uint32_t outTotalAllocatedNodesA = 0;
  uint32_t outFirstActiveNodeIDA = 0;

  // Profile Mesh A Initialization
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double tInitAStart = cuBQL::getCurrentTime();

  cuBQL::gpuBuilder_v3::test_speedrun_initialization(
      d_boxesA, numPrimsA, (uint32_t)maxCellSizeA, globalBoxA, d_dummyCountsA,
      outTempNodesA, outPrimStatesA, outTotalAllocatedNodesA, outFirstActiveNodeIDA,
      s, memResource
  );
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double tInitAMs = (cuBQL::getCurrentTime() - tInitAStart) * 1000.0;

  // --- MESH B RUN ---
  size_t countB = std::max((size_t)(1 << maxCellSizeB) * 4, (size_t)2 * numPrimsB);
  uint32_t* d_dummyCountsB = nullptr;
  _ALLOC(d_dummyCountsB, countB, s, memResource);

  cuBQL::gpuBuilder_v3::FinalNormalNode<T,D>* outTempNodesB = nullptr;
  cuBQL::gpuBuilder_v3::PrimState* outPrimStatesB = nullptr;
  uint32_t outTotalAllocatedNodesB = 0;
  uint32_t outFirstActiveNodeIDB = 0;

  // Profile Mesh B Initialization
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double tInitBStart = cuBQL::getCurrentTime();

  cuBQL::gpuBuilder_v3::test_speedrun_initialization(
      d_boxesB, numPrimsB, (uint32_t)maxCellSizeB, globalBoxB, d_dummyCountsB,
      outTempNodesB, outPrimStatesB, outTotalAllocatedNodesB, outFirstActiveNodeIDB,
      s, memResource
  );
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double tInitBMs = (cuBQL::getCurrentTime() - tInitBStart) * 1000.0;

  double tEnd = cuBQL::getCurrentTime();
  double elapsedMs = (tEnd - tStart) * 1000.0;

  // --------------------------------------------------------------------
  // VERIFY BUILDING METRICS READBACK
  // --------------------------------------------------------------------
  std::cout << "\n--------------------------------------------------" << std::endl;
  std::cout << " [INITIALIZATION METRICS OVERVIEW]" << std::endl;
  std::cout << " -> Tree A: Total Allocated Nodes = " << outTotalAllocatedNodesA << " | First Active ID = " << outFirstActiveNodeIDA << std::endl;
  std::cout << " -> Tree B: Total Allocated Nodes = " << outTotalAllocatedNodesB << " | First Active ID = " << outFirstActiveNodeIDB << std::endl;
  
  uint32_t activeCountA = outTotalAllocatedNodesA - outFirstActiveNodeIDA;
  uint32_t activeCountB = outTotalAllocatedNodesB - outFirstActiveNodeIDB;
  std::cout << " -> Active Traversal Ranges: Tree A = " << activeCountA << " nodes | Tree B = " << activeCountB << " nodes" << std::endl;
  std::cout << "--------------------------------------------------" << std::endl;

  std::vector<cuBQL::gpuBuilder_v3::FinalNormalNode<T,D>> h_nodesA(outTotalAllocatedNodesA);
  std::vector<cuBQL::gpuBuilder_v3::FinalNormalNode<T,D>> h_nodesB(outTotalAllocatedNodesB);

  CUBQL_CUDA_CALL(Memcpy(h_nodesA.data(), outTempNodesA, outTotalAllocatedNodesA * sizeof(cuBQL::gpuBuilder_v3::FinalNormalNode<T,D>), cudaMemcpyDeviceToHost));
  CUBQL_CUDA_CALL(Memcpy(h_nodesB.data(), outTempNodesB, outTotalAllocatedNodesB * sizeof(cuBQL::gpuBuilder_v3::FinalNormalNode<T,D>), cudaMemcpyDeviceToHost));

  // --------------------------------------------------------------------
  // DETAILED NODE-BY-NODE METRICS PRINTING
  // --------------------------------------------------------------------
  std::cout << "\n==================================================" << std::endl;
  std::cout << "          TREE A SPATIAL NODE DETAIL MAP          " << std::endl;
  std::cout << "==================================================" << std::endl;
  uint32_t sumPrimsA = 0;
  uint32_t validNodesA = 0;
  for (uint32_t i = outFirstActiveNodeIDA; i < outTotalAllocatedNodesA; ++i) {
      auto& openBranch = h_nodesA[i].openBranch;
      if (openBranch.count > 0) {
          validNodesA++;
          sumPrimsA += openBranch.count;
       //   std::cout << " Node [#" << i << "] -> Prims: " << openBranch.count << "\n";
      //    std::cout << "   Bounds Min: [";
      //    for(int d=0; d<D; ++d) std::cout << openBranch.bounds.get_lower(d) << (d < D-1 ? ", " : "");
      //    std::cout << "]\n   Bounds Max: [";
       //   for(int d=0; d<D; ++d) std::cout << openBranch.bounds.get_upper(d) << (d < D-1 ? ", " : "");
       //   std::cout << "]\n" << std::endl;
      }
  }

  std::cout << "==================================================" << std::endl;
  std::cout << "          TREE B SPATIAL NODE DETAIL MAP          " << std::endl;
  std::cout << "==================================================" << std::endl;
  uint32_t sumPrimsB = 0;
  uint32_t validNodesB = 0;
  for (uint32_t i = outFirstActiveNodeIDB; i < outTotalAllocatedNodesB; ++i) {
      auto& openBranch = h_nodesB[i].openBranch;
      if (openBranch.count > 0) {
          validNodesB++;
          sumPrimsB += openBranch.count;
       //   std::cout << " Node [#" << i << "] -> Prims: " << openBranch.count << "\n";
       //   std::cout << "   Bounds Min: [";
        //  for(int d=0; d<D; ++d) std::cout << openBranch.bounds.get_lower(d) << (d < D-1 ? ", " : "");
        //  std::cout << "]\n   Bounds Max: [";
        //  for(int d=0; d<D; ++d) std::cout << openBranch.bounds.get_upper(d) << (d < D-1 ? ", " : "");
        //  std::cout << "]\n" << std::endl;
      }
  }

  std::cout << "--------------------------------------------------" << std::endl;
  std::cout << " [POPULATED CELLS CHECK SUMMARY]" << std::endl;
  std::cout << " -> Tree A Populated Nodes: " << validNodesA << " / " << activeCountA << " | Total Prims = " << sumPrimsA << std::endl;
  std::cout << " -> Tree B Populated Nodes: " << validNodesB << " / " << activeCountB << " | Total Prims = " << sumPrimsB << std::endl;
  std::cout << "--------------------------------------------------\n" << std::endl;

  if (outTotalAllocatedNodesA == outTotalAllocatedNodesB && sumPrimsA == sumPrimsB) {
      std::cout << " [DIAGNOSTIC WARNING] Tree structural signatures are identical! Checking specific node boundaries..." << std::endl;
      
      if (activeCountA > 0 && activeCountB > 0) {
          auto& boxA = h_nodesA[outFirstActiveNodeIDA].openBranch.bounds;
          auto& boxB = h_nodesB[outFirstActiveNodeIDB].openBranch.bounds;
          bool boundsIdentical = true;
          for(int d=0; d<D; ++d) {
              if (boxA.get_lower(d) != boxB.get_lower(d) || boxA.get_upper(d) != boxB.get_upper(d)) {
                  boundsIdentical = false;
              }
          }
          if(boundsIdentical) {
              std::cout << " -> [ALERT] The first active layout nodes share EXACT structural bounding boxes!" << std::endl;
          } else {
              std::cout << " -> [OK] Structural footprints are distinct despite matching allocation lengths." << std::endl;
          }
      }
  }

  // --------------------------------------------------------------------
  // PARALLEL TEMP NODE CROSS CHECK & PERFORMANCE BENCHMARK
  // --------------------------------------------------------------------
  std::cout << "\n==================================================" << std::endl;
  std::cout << "          LAUNCHING CROSS-CHECK KERNELS           " << std::endl;
  std::cout << "==================================================" << std::endl;

  thrust::device_vector<uint32_t> d_intersectPairsA;
  thrust::device_vector<uint32_t> d_intersectPairsB;

  double crossCheckStart = cuBQL::getCurrentTime();

  uint32_t totalIntersections = executeNormalNodeCrossCheck<T,D>(
      outTempNodesA, outFirstActiveNodeIDA, outTotalAllocatedNodesA,
      outTempNodesB, outFirstActiveNodeIDB, outTotalAllocatedNodesB,
      d_intersectPairsA, d_intersectPairsB, s
  );
  CUBQL_CUDA_CALL(StreamSynchronize(s));

  double crossCheckEnd = cuBQL::getCurrentTime();
  double crossCheckMs = (crossCheckEnd - crossCheckStart) * 1000.0;

  uint64_t totalValidActivePairs = (uint64_t)validNodesA * validNodesB;

  double overlapPercentage = 0.0;
  if (totalValidActivePairs > 0) {
      overlapPercentage = (double)totalIntersections / (double)totalValidActivePairs * 100.0;
  }

  // Added initialization metrics outputs cleanly mapped alongside cross-check tracking
  std::cout << " -> Mesh A Speedrun Init Speed : " << tInitAMs << " ms" << std::endl;
  std::cout << " -> Mesh B Speedrun Init Speed : " << tInitBMs << " ms" << std::endl;
  std::cout << " -> Cross-Check Execution Speed : " << crossCheckMs << " ms" << std::endl;
  std::cout << " -> Intersecting Nodes Detected : " << totalIntersections << " pairs" << std::endl;
  std::cout << " -> True Active Matrix Evaluated: " << totalValidActivePairs << " pairings" << std::endl;
  std::cout << " -> Percentage of Node Overlaps : " << overlapPercentage << "%" << std::endl;
  std::cout << "==================================================\n" << std::endl;

  // --- CLEANUP ---
  if (outTempNodesA)  _FREE(outTempNodesA, s, memResource);
  if (outPrimStatesA) _FREE(outPrimStatesA, s, memResource);
  if (outTempNodesB)  _FREE(outTempNodesB, s, memResource);
  if (outPrimStatesB) _FREE(outPrimStatesB, s, memResource);
  if (d_dummyCountsA) _FREE(d_dummyCountsA, s, memResource);
  if (d_dummyCountsB) _FREE(d_dummyCountsB, s, memResource);
}
// --------------------------------------------------------------------
// MAIN ENTRY TESTING PIPELINE
// --------------------------------------------------------------------
extern "C" void kernelsTestBVH(const cuBQL::Triangle* hMeshA, int numTrianglesA, int maxCellSizeA,
                               const cuBQL::Triangle* hMeshB, int numTrianglesB, int maxCellSizeB,
                               int batchMultiplier,
                               int mode, 
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
  buildConfig.makeLeafThreshold = 1;

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
  cuBQL::box3f hostBoxA, hostBoxB;
  CUBQL_CUDA_CALL(Memcpy(&hostBoxA, &(bvhA.nodes[0].bounds), sizeof(cuBQL::box3f), cudaMemcpyDeviceToHost));
  CUBQL_CUDA_CALL(Memcpy(&hostBoxB, &(bvhB.nodes[0].bounds), sizeof(cuBQL::box3f), cudaMemcpyDeviceToHost));
  cudaDeviceSynchronize();

  testFastBVH<float, 3>(
      dBoxesA, numTrianglesA, maxCellSizeA,
      dBoxesB, numTrianglesB, maxCellSizeB,
      hostBoxA, hostBoxB,
      stream, memResource
  );

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