#define CUBQL_GPU_BUILDER_IMPLEMENTATION 1

// Core cuBQL headers
#include "cuBQL/bvh.h"
#include "cuBQL/builder/cuda.h"

// Your custom experimental builder layout
#include "include/third-party/cubql/sm_builder_v2.h"

// Custom traversal 
#include "include/third-party/cubql/fixedBoxQueryv2.h"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/scan.h>
#include <thrust/execution_policy.h>
#include <thrust/system/cuda/execution_policy.h>

#include <vector>
#include <algorithm>
#include "samples/common/loadOBJ.h"

// Include modular execution targets
//#include "crossCheck.h"
//#include "rapidDescendKernel.h"
#include "../custom_pipeline/GPUPredicatesCheck.h"
#include "../testBVH/crossCheck.h"
#include "DualTreeMeshIntersection.h"

// Include the updated header matching this implementation
//#include "kernelsTestBVH.h"

#include "../testBVH/ExecutionStats.h"
#include "DualTreeLauncher.h"

// --------------------------------------------------------------------
// EXISTING KERNELS & HELPERS
// --------------------------------------------------------------------
__global__ void generateBoxes(cuBQL::box3f* boxes, const cuBQL::Triangle* tris, int N) {
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if(i < N) { boxes[i] = tris[i].bounds(); }
}

// --------------------------------------------------------------------
// MAIN ENTRY TESTING PIPELINE
// --------------------------------------------------------------------
extern "C" void dualTreeKernelLaunch(const cuBQL::Triangle* hMeshA, int numTrianglesA, int maxCellSizeA,
                               const cuBQL::Triangle* hMeshB, int numTrianglesB, int maxCellSizeB,
                               int batchMultiplier, float fillRatio,
                               ExecutionStats& stats, 
                               std::vector<int2>& hGreenPairs,  
                               std::vector<int2>& hYellowPairs)
{
  if(numTrianglesA <= 0 || numTrianglesB <= 0) {
    return;
  }

  // Start tracking the comprehensive GPU workflow time
  double tPipelineStart = cuBQL::getCurrentTime();

  cudaStream_t stream = 0;
  cuBQL::DeviceMemoryResource memResource;
  cuBQL::BuildConfig buildConfig;
  buildConfig.makeLeafThreshold = 1;

  // --------------------------------------------------------------------
  // HIERARCHY A: BUILD & REFIT
  // --------------------------------------------------------------------
  cuBQL::Triangle* dMeshA;
  CUBQL_CUDA_CALL(Malloc(&dMeshA, numTrianglesA * sizeof(cuBQL::Triangle)));
  CUBQL_CUDA_CALL(Memcpy(dMeshA, hMeshA, numTrianglesA * sizeof(cuBQL::Triangle), cudaMemcpyDefault));

  double tBuildAStart = cuBQL::getCurrentTime();
  cuBQL::box3f* dBoxesA;
  CUBQL_CUDA_CALL(Malloc(&dBoxesA, numTrianglesA * sizeof(cuBQL::box3f)));
  generateBoxes<<<cuBQL::divRoundUp(numTrianglesA, 256), 256>>>(dBoxesA, dMeshA, numTrianglesA);

  uint32_t maxPossibleNodesA = 2 * numTrianglesA;
  thrust::device_vector<uint32_t> d_markedNodeIndicesA(maxPossibleNodesA, 0);
  thrust::device_vector<uint32_t> d_nodeDescendantCountsA(maxPossibleNodesA, 0); 
  uint32_t h_outMarkedCountA = 0;

  cuBQL::bvh3f bvhA;
  cuBQL::gpuBuilder_v2::build_custom(
      bvhA, dBoxesA, numTrianglesA, buildConfig, (uint32_t)maxCellSizeA, thrust::raw_pointer_cast(d_markedNodeIndicesA.data()),
      thrust::raw_pointer_cast(d_nodeDescendantCountsA.data()), &h_outMarkedCountA, stream, memResource);
  cuBQL::cuda::refit(bvhA, dBoxesA, stream);
  cudaDeviceSynchronize();
  double tBuildAEnd = cuBQL::getCurrentTime();

  // --------------------------------------------------------------------
  // HIERARCHY B: BUILD & REFIT
  // --------------------------------------------------------------------
  cuBQL::Triangle* dMeshB;
  CUBQL_CUDA_CALL(Malloc(&dMeshB, numTrianglesB * sizeof(cuBQL::Triangle)));
  CUBQL_CUDA_CALL(Memcpy(dMeshB, hMeshB, numTrianglesB * sizeof(cuBQL::Triangle), cudaMemcpyDefault));

  double tBuildBStart = cuBQL::getCurrentTime();
  cuBQL::box3f* dBoxesB;
  CUBQL_CUDA_CALL(Malloc(&dBoxesB, numTrianglesB * sizeof(cuBQL::box3f)));
  generateBoxes<<<cuBQL::divRoundUp(numTrianglesB, 256), 256>>>(dBoxesB, dMeshB, numTrianglesB);

  uint32_t maxPossibleNodesB = 2 * numTrianglesB;
  thrust::device_vector<uint32_t> d_markedNodeIndicesB(maxPossibleNodesB, 0);
  thrust::device_vector<uint32_t> d_nodeDescendantCountsB(maxPossibleNodesB, 0); 
  uint32_t h_outMarkedCountB = 0;

  cuBQL::bvh3f bvhB;
  cuBQL::gpuBuilder_v2::build_custom(
      bvhB, dBoxesB, numTrianglesB, buildConfig, (uint32_t)maxCellSizeB, thrust::raw_pointer_cast(d_markedNodeIndicesB.data()),
      thrust::raw_pointer_cast(d_nodeDescendantCountsB.data()), &h_outMarkedCountB, stream, memResource);
  cuBQL::cuda::refit(bvhB, dBoxesB, stream);
  cudaDeviceSynchronize();
  double tBuildBEnd = cuBQL::getCurrentTime();

  // --------------------------------------------------------------------
  // EXTRACTED GPU PARALLEL CRISS-CROSS INTERSECTION MODULE
  // --------------------------------------------------------------------
  double tCrossStart = cuBQL::getCurrentTime();
  
  thrust::device_vector<uint32_t> d_outPairsA;
  thrust::device_vector<uint32_t> d_outPairsB;

  uint32_t totalIntersections = executeCrissCrossIntersection(
      bvhA, d_markedNodeIndicesA, h_outMarkedCountA,
      bvhB, d_markedNodeIndicesB, h_outMarkedCountB,
      d_outPairsA, d_outPairsB
  );

  double tCrossEnd = cuBQL::getCurrentTime();

  uint64_t totalPossiblePairs = (uint64_t)h_outMarkedCountA * h_outMarkedCountB;
  double intersectionPercentage = totalPossiblePairs > 0 ? ((double)totalIntersections / totalPossiblePairs) * 100.0 : 0.0;


  // ====================================================================
  // ASYNC BATCHED CROSS-INTERSECTION COUNTING LOOP
  // ====================================================================
  int totalBatches = d_outPairsA.size();
  
  IntersectionTimeTracker tracker;


   uint64_t finalCandidatePairs = executeDualTreeTraversal(
    batchMultiplier,
    totalBatches,
    maxCellSizeA,      // Max total primitives across subtrees in A
    maxCellSizeB,      // Max total primitives across subtrees in B
    fillRatio, // Value between (0.0f, 1.0f] matching expected overlap complexity
    d_outPairsA,
    d_outPairsB,
    d_markedNodeIndicesA,
    d_markedNodeIndicesB,
    d_nodeDescendantCountsA, 
    d_nodeDescendantCountsB, 
    h_outMarkedCountA,
    h_outMarkedCountB,
    bvhA,
    bvhB,
    dMeshA,
    dMeshB,
    hGreenPairs,  
    hYellowPairs, 
    tracker
);

  // --------------------------------------------------------------------
  // MEMORY CLEANUP
  // --------------------------------------------------------------------
  CUBQL_CUDA_CALL(Free(dMeshA));
  CUBQL_CUDA_CALL(Free(dBoxesA));
  CUBQL_CUDA_CALL(Free(dMeshB));
  CUBQL_CUDA_CALL(Free(dBoxesB));
  cuBQL::cuda::free(bvhA, stream, memResource);
  cuBQL::cuda::free(bvhB, stream, memResource);

  double tPipelineEnd = cuBQL::getCurrentTime();

  // ====================================================================
  // EXPORT METRICS TO EXECUTIONSTATS
  // ====================================================================
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
  //stats.parallelDfsDescentBMs = (tGpuBfsEnd - tGpuBfsStart) * 1000.0;
  
  stats.GPUTotalTime          = (tPipelineEnd - tPipelineStart) * 1000.0;

  stats.totalCrissCrossBatches  = totalBatches;
  stats.finalAabbCandidatePairs = finalCandidatePairs;
  stats.confirmedGreenPairs     = hGreenPairs.size();
  stats.confirmedYellowPairs    = hYellowPairs.size();

  stats.loopTracker             = tracker;
}