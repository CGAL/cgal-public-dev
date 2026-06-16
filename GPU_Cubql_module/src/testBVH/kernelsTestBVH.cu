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
#include <iostream>
#include <algorithm>
#include <iomanip>
#include "samples/common/loadOBJ.h"

// Include modular execution targets
#include "crossCheck.h"
#include "rapidDescendKernel.h"
#include "volumeSanityCheck.h" 
#include "batchedCrossIntersection.h"


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
                               int batchMultiplier) // Defaulted to 4, tweak as needed
{
  std::cout << "\n==================================================\n";
  std::cout << " RUNNING DUAL-MESH CROSS-INTERSECTION PIPELINE V2\n";
  std::cout << "==================================================\n";
  std::cout << "Mesh A: " << numTrianglesA << " tris (MaxCell: " << maxCellSizeA << ")\n";
  std::cout << "Mesh B: " << numTrianglesB << " tris (MaxCell: " << maxCellSizeB << ")\n";

  if(numTrianglesA <= 0 || numTrianglesB <= 0) {
    std::cerr << "Error: One or both input meshes contain no triangles.\n";
    return;
  }

  double tUploadStart = cuBQL::getCurrentTime();
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

  // --------------------------------------------------------------------
  // DEBUG: INSPECT PAIRING FRONTIER
  // --------------------------------------------------------------------
  thrust::host_vector<uint32_t> h_debugPairsA = d_outPairsA;
  thrust::host_vector<uint32_t> h_debugPairsB = d_outPairsB;
  
  std::cout << "\n--- CRISS-CROSS PAIRING FRONTIER DEBUG ---" << std::endl;
  std::cout << "Total Pairs Generated: " << h_debugPairsA.size() << std::endl;
  
  // Print the first 20 pairs to check for logical redundancy
  int limit = std::min((int)h_debugPairsA.size(), 20);
  for(int i = 0; i < limit; ++i) {
      std::cout << "Pair [" << i << "]: MeshA NodeIdx[" << h_debugPairsA[i] 
                << "] <-> MeshB NodeIdx[" << h_debugPairsB[i] << "]" << std::endl;
  }
  if (h_debugPairsA.size() > 20) std::cout << "... (truncated)" << std::endl;
  std::cout << "------------------------------------------\n" << std::endl;

  // --------------------------------------------------------------------
  // RAPID DESCENT: EXECUTED FOR MESH B ONLY
  // --------------------------------------------------------------------
  double tGpuBfsStart = cuBQL::getCurrentTime();
  thrust::device_vector<uint32_t> d_outOffsetsB;   
  thrust::device_vector<uint32_t> d_outPrimsFlatB; 

  executeRapidDescentBFS(
      bvhB, numTrianglesB, d_markedNodeIndicesB, d_nodeDescendantCountsB, 
      h_outMarkedCountB, d_outOffsetsB, d_outPrimsFlatB
  );
  double tGpuBfsEnd = cuBQL::getCurrentTime();

 // ====================================================================
  // SANITY CHECK: Verify Offset Mapping against Descendant Counts
  // ====================================================================
  std::cout << "Verifying d_outOffsetsB mapping against d_outPairsB..." << std::endl;

  // 1. Ensure we have the pairs on the host
  thrust::host_vector<uint32_t> h_outPairsB = d_outPairsB;
  thrust::host_vector<uint32_t> h_outOffsetsB = d_outOffsetsB;
  thrust::host_vector<uint32_t> h_markedNodeIndicesB = d_markedNodeIndicesB;
  thrust::host_vector<uint32_t> h_nodeDescendantCountsB = d_nodeDescendantCountsB;

  uint32_t totalBatchez = (uint32_t)h_outPairsB.size();
  bool offsetError = false;

  // 2. Validate consistency
  for (uint32_t i = 0; i < totalBatchez; ++i) {
      uint32_t batchBId = h_outPairsB[i];
      
      if (batchBId >= h_outMarkedCountB) {
          std::cerr << "OUT OF BOUNDS: batchBId " << batchBId 
                    << " exceeds h_outMarkedCountB (" << h_outMarkedCountB << ")" << std::endl;
          offsetError = true;
          break;
      }

      uint32_t nodeIdx = h_markedNodeIndicesB[batchBId];
      uint32_t expectedCount = h_nodeDescendantCountsB[nodeIdx];
      uint32_t actualCount = h_outOffsetsB[batchBId + 1] - h_outOffsetsB[batchBId];

      if (actualCount != expectedCount) {
          std::cerr << "OFFSET ERROR at pair index " << i 
                    << ": Batch " << batchBId 
                    << " (Node " << nodeIdx << ")"
                    << " expected " << expectedCount << " primitives, "
                    << " but offset range reserved " << actualCount << "!" << std::endl;
          offsetError = true;
          break; 
      }
  }

  // 3. Final feedback based on validation result
  if (offsetError) {
      std::cerr << "CRITICAL: Offset structure validation FAILED." << std::endl;
      // Depending on your project requirements, you might want to throw an exception:
      // throw std::runtime_error("Offset validation failed!");
  } else {
      std::cout << "SUCCESS: d_outOffsetsB is perfectly aligned with descendant counts." << std::endl;
  }

  // ====================================================================
  // GEOMETRIC EVALUATION & VOLUMETRIC SANITY CHECK FOR MESH B
  // ====================================================================
 // runVolumeSanityCheck(
 //     bvhB, h_outMarkedCountB, d_markedNodeIndicesB, d_nodeDescendantCountsB,
 //     d_outOffsetsB, d_outPrimsFlatB, hMeshB
 // );
// ====================================================================
  // ASYNC BATCHED CROSS-INTERSECTION COUNTING LOOP (EXTRACTED)
  // ====================================================================
  int totalBatches = d_outPairsA.size();
  
  // UPDATED: Instantiate the specialized profiling tracker structure
  IntersectionTimeTracker tracker;
  
  uint64_t finalCandidatePairs = executeBatchedCrossIntersectionLoop(
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
      tracker // UPDATED: Pass the tracker struct by reference
  );

  // --------------------------------------------------------------------
  // PERFORMANCE OUTPUT
  // --------------------------------------------------------------------
  std::cout << "--------------------------------------------------\n";
  std::cout << " STRUCTURE SUMMARY & PROPORTIONS\n";
  std::cout << "--------------------------------------------------\n";
  std::cout << "Mesh A Total Generated Nodes: " << bvhA.numNodes << "\n";
  std::cout << "Mesh A Extracted Targets (<" << maxCellSizeA << "): " << h_outMarkedCountA << "\n";
  std::cout << "Mesh B Total Generated Nodes: " << bvhB.numNodes << "\n";
  std::cout << "Mesh B Extracted Targets (<" << maxCellSizeB << "): " << h_outMarkedCountB << "\n\n";

  std::cout << "--------------------------------------------------\n";
  std::cout << " CRISS-CROSS BOUNDING BOX CROSS-CHECK\n";
  std::cout << "--------------------------------------------------\n";
  std::cout << "Intersection found: " << totalIntersections << " / " << totalPossiblePairs << "\n";
  std::cout << "Intersection Ratio: " << std::fixed << std::setprecision(4) << intersectionPercentage << "%\n\n";

  std::cout << "--------------------------------------------------\n";
  std::cout << " TIMING METRICS OVERVIEW\n";
  std::cout << "--------------------------------------------------\n";
  std::cout << "  |- Build + Refit (Mesh A):  " << (tBuildAEnd - tBuildAStart) * 1000.0 << " ms\n";
  std::cout << "  |- Build + Refit (Mesh B):  " << (tBuildBEnd - tBuildBStart) * 1000.0 << " ms\n";
  std::cout << "  |- GPU Cross-Check Engine:  " << (tCrossEnd - tCrossStart) * 1000.0 << " ms\n";
  std::cout << "  |- Parallel DFS Descent (B): " << (tGpuBfsEnd - tGpuBfsStart) * 1000.0 << " ms\n";
  std::cout << "--------------------------------------------------\n";

  // --------------------------------------------------------------------
  // DUAL-TREE BATCHING METRICS
  // --------------------------------------------------------------------
  std::cout << "--------------------------------------------------\n";
  std::cout << " DUAL-TREE DESCENT BATCHING METRICS\n";
  std::cout << "--------------------------------------------------\n";
  std::cout << "Total Criss-Cross Batches Processed: " << totalBatches << "\n";
  std::cout << "Final AABB Candidate Pairs Found:    " << finalCandidatePairs << "\n";
  std::cout << "--------------------------------------------------\n";

  // UPDATED: Print out the fine-grained execution phases breakdown (Sandbox Vs. CUDA)
  tracker.print();

  CUBQL_CUDA_CALL(Free(dMeshA));
  CUBQL_CUDA_CALL(Free(dBoxesA));
  CUBQL_CUDA_CALL(Free(dMeshB));
  CUBQL_CUDA_CALL(Free(dBoxesB));
  cuBQL::cuda::free(bvhA, stream, memResource);
  cuBQL::cuda::free(bvhB, stream, memResource);
}