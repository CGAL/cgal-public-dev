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
#include "crossCheck.h"
#include "rapidDescendKernel.h"
#include "batchedCrossIntersection.h"

// Include the updated header matching this implementation
#include "kernelsTestBVH.h"

#include "batchedCrossIntersectionBrute.h"

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
extern "C" void kernelsTestBVH(const cuBQL::Triangle* hMeshA, int numTrianglesA, int maxCellSizeA,
                               const cuBQL::Triangle* hMeshB, int numTrianglesB, int maxCellSizeB,
                               int batchMultiplier,
                               int mode, // <-- Added mode parameter
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

  // Explicit scopes or instant instantiations invoke heavy driver loads
  thrust::device_vector<uint32_t> d_markedNodeIndicesA(maxPossibleNodesA, 0);
  thrust::device_vector<uint32_t> d_nodeDescendantCountsA(maxPossibleNodesA, 0); 
  thrust::device_vector<uint32_t> d_markedNodeIndicesB(maxPossibleNodesB, 0);
  thrust::device_vector<uint32_t> d_nodeDescendantCountsB(maxPossibleNodesB, 0); 

  thrust::device_vector<uint32_t> d_outPairsA;
  thrust::device_vector<uint32_t> d_outPairsB;
  
  // Mesh B structural structures
  thrust::device_vector<uint32_t> d_outOffsetsB;   
  thrust::device_vector<uint32_t> d_outPrimsFlatB;

  // Mesh A structural structures (Allocated but empty if mode == 0)
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
  cuBQL::gpuBuilder_v2::build_custom(
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

  uint32_t totalIntersections = executeCrissCrossIntersection(
      bvhA, d_markedNodeIndicesA, h_outMarkedCountA,
      bvhB, d_markedNodeIndicesB, h_outMarkedCountB,
      d_outPairsA, d_outPairsB
  );

  cudaDeviceSynchronize(); // Ensure engine operations finish before timestamp
  double tCrossEnd = cuBQL::getCurrentTime();

  uint64_t totalPossiblePairs = (uint64_t)h_outMarkedCountA * h_outMarkedCountB;
  double intersectionPercentage = totalPossiblePairs > 0 ? ((double)totalIntersections / totalPossiblePairs) * 100.0 : 0.0;

  int totalBatches = d_outPairsA.size();
  IntersectionTimeTracker tracker;
  uint64_t finalCandidatePairs = 0;

  // --------------------------------------------------------------------
  // CONDITIONAL PIPELINE EXECUTION BASED ON MODE
  // --------------------------------------------------------------------
  double tGpuBfsStart = cuBQL::getCurrentTime();
  double tGpuBfsEnd = tGpuBfsStart;

  if (mode == 0) {
    // --------------------------------------------------------------------
    // MODE 0: RAPID DESCENT FOR MESH B ONLY & STANDARD INTERSECTION LOOP
    // --------------------------------------------------------------------
    executeRapidDescentBFS(
        bvhB, numTrianglesB, d_markedNodeIndicesB, d_nodeDescendantCountsB, 
        h_outMarkedCountB, d_outOffsetsB, d_outPrimsFlatB
    );
    
    cudaDeviceSynchronize();
    tGpuBfsEnd = cuBQL::getCurrentTime();

    finalCandidatePairs = executeBatchedCrossIntersectionLoop(
        batchMultiplier,
        totalBatches,
        d_outPairsA,
        d_outPairsB,
        d_markedNodeIndicesA,
        d_markedNodeIndicesB,
        d_outOffsetsB,
        d_outPrimsFlatB,
        d_nodeDescendantCountsB,
        h_outMarkedCountB,
        bvhA,
        dMeshA,
        dMeshB,
        hGreenPairs,   
        hYellowPairs,  
        tracker        
    );
  } 
  else {
    // --------------------------------------------------------------------
    // MODE 1: RAPID DESCENT FOR MESH B AND MESH A & SYMMETRIC INTERSECTION LOOP
    // --------------------------------------------------------------------
    // Rapid Descent on Mesh B
    executeRapidDescentBFS(
        bvhB, numTrianglesB, d_markedNodeIndicesB, d_nodeDescendantCountsB, 
        h_outMarkedCountB, d_outOffsetsB, d_outPrimsFlatB
    );

    // Rapid Descent on Mesh A
    executeRapidDescentBFS(
        bvhA, numTrianglesA, d_markedNodeIndicesA, d_nodeDescendantCountsA, 
        h_outMarkedCountA, d_outOffsetsA, d_outPrimsFlatA
    );
    
    cudaDeviceSynchronize();
    tGpuBfsEnd = cuBQL::getCurrentTime();

    finalCandidatePairs = executeBatchedCrossIntersectionLoop(
        batchMultiplier,
        totalBatches,
        d_outPairsA,
        d_outPairsB,
        d_markedNodeIndicesA,
        d_markedNodeIndicesB,
        d_outOffsetsB,
        d_outPrimsFlatB,
        d_nodeDescendantCountsB, 
        d_outOffsetsA,
        d_outPrimsFlatA,
        d_nodeDescendantCountsA, 
        h_outMarkedCountA,
        h_outMarkedCountB,
        dMeshA,
        dMeshB,
        hGreenPairs,
        hYellowPairs,
        tracker 
    );
  }

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

  // Manually forcing Thrust vectors to deallocate before pipeline stop
  d_markedNodeIndicesA.shrink_to_fit();
  d_nodeDescendantCountsA.shrink_to_fit();
  d_markedNodeIndicesB.shrink_to_fit();
  d_nodeDescendantCountsB.shrink_to_fit();
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
  stats.parallelDfsDescentBMs = (tGpuBfsEnd - tGpuBfsStart) * 1000.0; // Tracks combined BFS execution time
  
  stats.GPUTotalTime          = (tPipelineEnd - tPipelineStart) * 1000.0;

  stats.totalCrissCrossBatches  = totalBatches;
  stats.finalAabbCandidatePairs = finalCandidatePairs;
  stats.confirmedGreenPairs     = hGreenPairs.size();
  stats.confirmedYellowPairs    = hYellowPairs.size();

  stats.loopTracker             = tracker;
}