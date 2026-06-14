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
#include <vector>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include "samples/common/loadOBJ.h"

// --------------------------------------------------------------------
// DEPENDENT TYPES & STRUCTS
// --------------------------------------------------------------------
typedef typename cuBQL::BinaryBVH<float, 3>::Node BvhNodeT;

// --------------------------------------------------------------------
// EXISTING KERNELS & HELPERS (UNCHANGED)
// --------------------------------------------------------------------
__global__ void countAABBOverlapsKernel(int *pairCounts, cuBQL::bvh3f bvhA, const cuBQL::Triangle *triA, const cuBQL::Triangle *triB, int currentBatchSize, uint64_t startNodeIdx) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= currentBatchSize) return;
    cuBQL::Triangle b = triB[tid];
    cuBQL::box3f query = b.bounds();
    int count = 0;
    cuBQL::fixedBoxQueryv2::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        for (uint32_t i = 0; i < num; i++) {
            if (triA[ids[i]].bounds().overlaps(query)) { count++; }
        }
        return CUBQL_CONTINUE_TRAVERSAL;
    }, bvhA, query, startNodeIdx);
    pairCounts[tid] = count;
}

__global__ void fillAABBOverlapsKernel(int2 *candidatePairs, const int *offsets, cuBQL::bvh3f bvhA, const cuBQL::Triangle *triA, const cuBQL::Triangle *triB, int currentBatchSize, uint64_t startNodeIdx) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= currentBatchSize) return;
    int wPos = offsets[tid];
    cuBQL::Triangle b = triB[tid];
    cuBQL::box3f query = b.bounds();
    cuBQL::fixedBoxQueryv2::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        for (uint32_t i = 0; i < num; i++) {
            if (triA[ids[i]].bounds().overlaps(query)) { candidatePairs[wPos++] = make_int2((int)ids[i], tid); }
        }
        return CUBQL_CONTINUE_TRAVERSAL;
    }, bvhA, query, startNodeIdx);
}

__global__ void generateBoxes(cuBQL::box3f* boxes, const cuBQL::Triangle* tris, int N) {
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if(i < N) { boxes[i] = tris[i].bounds(); }
}

__device__ inline bool boxBoxesIntersect(const cuBQL::box3f& a, const cuBQL::box3f& b) {
  return (a.lower[0] <= b.upper[0] && a.upper[0] >= b.lower[0]) &&
         (a.lower[1] <= b.upper[1] && a.upper[1] >= b.lower[1]) &&
         (a.lower[2] <= b.upper[2] && a.upper[2] >= b.lower[2]);
}

__global__ void intersectTargetNodesDirectly(
    const BvhNodeT* d_nodesA, const uint32_t* d_indicesA, uint32_t countA,
    const BvhNodeT* d_nodesB, const uint32_t* d_indicesB, uint32_t countB,
    uint32_t* d_intersectionsCount) 
{
  uint32_t idxA = threadIdx.x + blockIdx.x * blockDim.x;
  if (idxA >= countA) return;
  uint32_t nodeIDA = d_indicesA[idxA];
  cuBQL::box3f boxA = d_nodesA[nodeIDA].bounds;
  uint32_t localIntersections = 0;
  for (uint32_t idxB = 0; idxB < countB; ++idxB) {
    uint32_t nodeIDB = d_indicesB[idxB];
    cuBQL::box3f boxB = d_nodesB[nodeIDB].bounds;
    if (boxBoxesIntersect(boxA, boxB)) { localIntersections++; }
  }
  if (localIntersections > 0) { atomicAdd(d_intersectionsCount, localIntersections); }
}

// --------------------------------------------------------------------
// PHASE 1: RECURSIVE STRUCTURAL TREE CHECK (UNCHANGED CORE LOGIC)
// --------------------------------------------------------------------
bool verifyDescendantTree(uint32_t nodeID, 
                          const std::vector<BvhNodeT>& h_nodes, 
                          const thrust::host_vector<uint32_t>& h_counts) 
{
  const auto& node = h_nodes[nodeID];
  uint32_t expectedCount = h_counts[nodeID];

  if (node.admin.count > 0) {
    if (node.admin.count != expectedCount) {
      std::cerr << "[SANITY ERROR] Leaf Node " << nodeID 
                << " mismatch! Actual leaf count: " << node.admin.count 
                << ", Tracked map says: " << expectedCount << "\n";
      return false;
    }
    return true;
  }

  uint32_t leftChildID = node.admin.offset + 0;
  uint32_t rightChildID = node.admin.offset + 1;

  if (leftChildID >= h_nodes.size() || rightChildID >= h_nodes.size()) {
    std::cerr << "[SANITY ERROR] Inner Node " << nodeID 
              << " points to out-of-bounds children (" << leftChildID << ", " << rightChildID << ")\n";
    return false;
  }

  uint32_t leftCount = h_counts[leftChildID];
  uint32_t rightCount = h_counts[rightChildID];

  if (expectedCount != (leftCount + rightCount)) {
    std::cerr << "[SANITY ERROR] Hierarchy Mismatch at Node " << nodeID << "!\n"
              << "  |- Stored Count:     " << expectedCount << "\n"
              << "  |- Left Child (" << leftChildID << "): " << leftCount << "\n"
              << "  |- Right Child (" << rightChildID << "): " << rightCount << "\n"
              << "  |- Expected Sum:     " << (leftCount + rightCount) << "\n";
    return false;
  }

  return verifyDescendantTree(leftChildID, h_nodes, h_counts) && 
         verifyDescendantTree(rightChildID, h_nodes, h_counts);
}

// --------------------------------------------------------------------
// PHASE 2: EXTENDED SEED COHERENCY VALIDATOR (UNCHANGED CORE LOGIC)
// --------------------------------------------------------------------
bool verifyMarkedSeedsCoherency(const std::string& meshLabel,
                                uint32_t markedCount,
                                const thrust::host_vector<uint32_t>& h_markedIndices,
                                const thrust::host_vector<uint32_t>& h_nodeDescendantCounts)
{
  bool allSeedsValid = true;
  uint32_t zeroCountErrors = 0;

  for (uint32_t i = 0; i < markedCount; ++i) {
    uint32_t nodeID = h_markedIndices[i];
    uint32_t registeredDescendants = h_nodeDescendantCounts[nodeID];

    if (registeredDescendants == 0) {
      if (zeroCountErrors < 5) { // Log the first few issues
        std::cerr << "  |--> [SEED COHERENCY CRITICAL FAILURE] " << meshLabel << " | Seed Index Element Slot [" << i << "]\n"
                  << "        |- Points to Tree nodeID    : " << nodeID << "\n"
                  << "        |- h_nodeDescendantCounts[" << nodeID << "] value : " << registeredDescendants << " (!!! ERROR !!!)\n";
      }
      zeroCountErrors++;
      allSeedsValid = false;
    }
  }

  if (zeroCountErrors > 0) {
    std::cerr << "  |--> [SUMMARY] " << meshLabel << " has " << zeroCountErrors 
              << " broken marked target seeds out of " << markedCount << " total.\n";
  }
  return allSeedsValid;
}

// --------------------------------------------------------------------
// EXTERNAL LINKAGE IMPORT
// --------------------------------------------------------------------
extern "C" void executeRapidDescentBFS(
    const cuBQL::bvh3f& bvh,
    int numTriangles,
    const thrust::device_vector<uint32_t>& d_markedNodeIndices,
    const thrust::device_vector<uint32_t>& d_nodeDescendantCounts,
    uint32_t h_outMarkedCount,
    thrust::device_vector<uint32_t>& d_outOffsets,
    thrust::device_vector<uint32_t>& d_outPrimsFlat);

// --------------------------------------------------------------------
// MAIN ENTRY TESTING PIPELINE
// --------------------------------------------------------------------
extern "C" void kernelsTestBVH(const cuBQL::Triangle* hMeshA, int numTrianglesA, int maxCellSizeA,
                               const cuBQL::Triangle* hMeshB, int numTrianglesB, int maxCellSizeB) {
  std::cout << "\n==================================================\n";
  std::cout << " RUNNING DUAL-MESH CROSS-INTERSECTION PIPELINE\n";
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
  // RUNNING THE RUNTIME SANITY CHECKS (WITH HOST TIMERS)
  // --------------------------------------------------------------------
  std::cout << "\n--------------------------------------------------\n";
  std::cout << " RUNNING GLOBAL DESCENDANT SANITY CHECKS\n";
  std::cout << "--------------------------------------------------\n";
  
  double tSanityStart = cuBQL::getCurrentTime();

  // Verify Mesh A
  std::vector<BvhNodeT> h_nodesVectorA(bvhA.numNodes);
  CUBQL_CUDA_CALL(Memcpy(h_nodesVectorA.data(), bvhA.nodes, bvhA.numNodes * sizeof(BvhNodeT), cudaMemcpyDefault));
  thrust::host_vector<uint32_t> h_countsHostA = d_nodeDescendantCountsA;
  thrust::host_vector<uint32_t> h_markedIndicesHostA = d_markedNodeIndicesA;
  
  bool meshASane = verifyDescendantTree(0, h_nodesVectorA, h_countsHostA);
  bool seedsASane = verifyMarkedSeedsCoherency("Mesh A", h_outMarkedCountA, h_markedIndicesHostA, h_countsHostA);
  std::cout << "  |- Mesh A Tree Consistency   : " << (meshASane ? "PASSED" : "FAILED") << "\n";
  std::cout << "  |- Mesh A Target Seed Check  : " << (seedsASane ? "PASSED SUCCESSFULLY" : "CRITICAL SEED COUNT MISMATCH ERROR") << "\n";

  // Verify Mesh B
  std::vector<BvhNodeT> h_nodesVectorB(bvhB.numNodes);
  CUBQL_CUDA_CALL(Memcpy(h_nodesVectorB.data(), bvhB.nodes, bvhB.numNodes * sizeof(BvhNodeT), cudaMemcpyDefault));
  thrust::host_vector<uint32_t> h_countsHostB = d_nodeDescendantCountsB;
  thrust::host_vector<uint32_t> h_markedIndicesHostB = d_markedNodeIndicesB;

  bool meshBSane = verifyDescendantTree(0, h_nodesVectorB, h_countsHostB);
  bool seedsBSane = verifyMarkedSeedsCoherency("Mesh B", h_outMarkedCountB, h_markedIndicesHostB, h_countsHostB);
  std::cout << "  |- Mesh B Tree Consistency   : " << (meshBSane ? "PASSED" : "FAILED") << "\n";
  std::cout << "  |- Mesh B Target Seed Check  : " << (seedsBSane ? "PASSED SUCCESSFULLY" : "CRITICAL SEED COUNT MISMATCH ERROR") << "\n";

  double tSanityEnd = cuBQL::getCurrentTime();
  std::cout << "  |- Host Verification Loop Time: " << (tSanityEnd - tSanityStart) * 1000.0 << " ms\n\n";

  if (!seedsASane || !seedsBSane) {
     std::cerr << "[CRITICAL TERMINATION] Halting pipeline. The builder failed to assign counts to targets.\n";
     exit(EXIT_FAILURE);
  }

  // --------------------------------------------------------------------
  // GPU PARALLEL CRISS-CROSS INTERSECTION TEST
  // --------------------------------------------------------------------
  double tCrossStart = cuBQL::getCurrentTime();
  thrust::device_vector<uint32_t> d_intersectionCounter(1, 0);

  if(h_outMarkedCountA > 0 && h_outMarkedCountB > 0) {
    uint32_t threadsCross = 256;
    uint32_t blocksCross = cuBQL::divRoundUp(h_outMarkedCountA, threadsCross);

    intersectTargetNodesDirectly<<<blocksCross, threadsCross>>>(
        bvhA.nodes, thrust::raw_pointer_cast(d_markedNodeIndicesA.data()), h_outMarkedCountA,
        bvhB.nodes, thrust::raw_pointer_cast(d_markedNodeIndicesB.data()), h_outMarkedCountB,
        thrust::raw_pointer_cast(d_intersectionCounter.data()));
    cudaDeviceSynchronize();
  }
  double tCrossEnd = cuBQL::getCurrentTime();

  uint32_t totalIntersections = d_intersectionCounter[0];
  uint64_t totalPossiblePairs = (uint64_t)h_outMarkedCountA * h_outMarkedCountB;
  double intersectionPercentage = totalPossiblePairs > 0 ? ((double)totalIntersections / totalPossiblePairs) * 100.0 : 0.0;

  // --------------------------------------------------------------------
  // RAPID DESCENT: NOW EXECUTED FOR MESH B ONLY
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
  // MODIFIED: GEOMETRIC EVALUATION & VOLUMETRIC SANITY CHECK FOR MESH B
  // ====================================================================
  if(h_outMarkedCountB > 0) {
    std::cout << "\n--------------------------------------------------------------------\n";
    std::cout << "--- EXTRACTED TARGET THRESHOLD NODES (GEOMETRIC EVALUATION - MESH B) ---\n";
    std::cout << "--------------------------------------------------------------------\n";
    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout.precision(5);

    // Compute scene global volume base from root node (Node 0) of Mesh B
    auto& rootNode = h_nodesVectorB[0];
    float rdx = rootNode.bounds.upper[0] - rootNode.bounds.lower[0];
    float rdy = rootNode.bounds.upper[1] - rootNode.bounds.lower[1];
    float rdz = rootNode.bounds.upper[2] - rootNode.bounds.lower[2];
    float globalSceneVolume = (rdx > 0.0f && rdy > 0.0f && rdz > 0.0f) ? (rdx * rdy * rdz) : 1.0f;

    // Gather Descent output tracking details onto host from Mesh B vectors
    thrust::host_vector<uint32_t> h_outOffsetsB = d_outOffsetsB;
    thrust::host_vector<uint32_t> h_outPrimsFlatB = d_outPrimsFlatB;

    for(size_t i = 0; i < h_outMarkedCountB; ++i) {
      uint32_t nodeID = h_markedIndicesHostB[i];
      uint32_t expectedPrimsCount = h_countsHostB[nodeID];

      // 1. Compute original node metadata bounding box volume ratio
      auto& node = h_nodesVectorB[nodeID];
      float dx = node.bounds.upper[0] - node.bounds.lower[0];
      float dy = node.bounds.upper[1] - node.bounds.lower[1];
      float dz = node.bounds.upper[2] - node.bounds.lower[2];
      float nodeVolume = (dx > 0.0f && dy > 0.0f && dz > 0.0f) ? (dx * dy * dz) : 0.0f;
      float metadataVolumePercentage = (nodeVolume / globalSceneVolume) * 100.0f;

      // 2. Extract slice bounds from flattened parallel harvest outputs
      uint32_t startOffset = h_outOffsetsB[i];
      uint32_t endOffset = (i + 1 < h_outMarkedCountB) ? h_outOffsetsB[i + 1] : (uint32_t)h_outPrimsFlatB.size();
      uint32_t harvestedPrimsCount = endOffset - startOffset;

      cuBQL::box3f primitiveHarvestBBox; // Empty initialization boundary limits
      primitiveHarvestBBox.lower[0] = primitiveHarvestBBox.lower[1] = primitiveHarvestBBox.lower[2] = INFINITY;
      primitiveHarvestBBox.upper[0] = primitiveHarvestBBox.upper[1] = primitiveHarvestBBox.upper[2] = -INFINITY;

      for (uint32_t pIdx = startOffset; pIdx < endOffset; ++pIdx) {
        uint32_t actualPrimID = h_outPrimsFlatB[pIdx];
        cuBQL::box3f primBox = hMeshB[actualPrimID].bounds();
        
        // Manual element level box union extension
        for(int axis=0; axis<3; ++axis) {
          primitiveHarvestBBox.lower[axis] = std::min(primitiveHarvestBBox.lower[axis], primBox.lower[axis]);
          primitiveHarvestBBox.upper[axis] = std::max(primitiveHarvestBBox.upper[axis], primBox.upper[axis]);
        }
      }

      // 3. Compute volume calculated directly from the combined primitive collection
      float bdx = primitiveHarvestBBox.upper[0] - primitiveHarvestBBox.lower[0];
      float bdy = primitiveHarvestBBox.upper[1] - primitiveHarvestBBox.lower[1];
      float bdz = primitiveHarvestBBox.upper[2] - primitiveHarvestBBox.lower[2];
      float bfsVolume = (bdx > 0.0f && bdy > 0.0f && bdz > 0.0f) ? (bdx * bdy * bdz) : 0.0f;
      float gpuBfsVolumePercentage = (bfsVolume / globalSceneVolume) * 100.0f;

      std::cout << "  Item [" << i << "] -> Node ID: " << nodeID << "\n";
      std::cout << "       |- Primitives Validation -> Expected (Tracked Map): " << expectedPrimsCount
                << " | GPU BFS Harvest: " << harvestedPrimsCount << "\n";
      std::cout << "       |- Node Metadata Volume Ratio:       " << metadataVolumePercentage << "%\n";
      std::cout << "       |- GPU Parallel BFS Descent Ratio:   " << gpuBfsVolumePercentage << "%\n";
      
      if (harvestedPrimsCount != expectedPrimsCount) {
         std::cout << "       |----> [WARNING] Harvest balance mismatch flagged!\n";
      }
      std::cout << "\n";
    }
    std::cout << "--------------------------------------------------------------------\n";
  }

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
  std::cout << "  |- Host Tree Sanity Check:  " << (tSanityEnd - tSanityStart) * 1000.0 << " ms\n";
  std::cout << "  |- GPU Cross-Check Engine:  " << (tCrossEnd - tCrossStart) * 1000.0 << " ms\n";
  std::cout << "  |- Parallel DFS Descent (B): " << (tGpuBfsEnd - tGpuBfsStart) * 1000.0 << " ms\n";
  std::cout << "--------------------------------------------------\n";

  CUBQL_CUDA_CALL(Free(dMeshA));
  CUBQL_CUDA_CALL(Free(dBoxesA));
  CUBQL_CUDA_CALL(Free(dMeshB));
  CUBQL_CUDA_CALL(Free(dBoxesB));
  cuBQL::cuda::free(bvhA, stream, memResource);
  cuBQL::cuda::free(bvhB, stream, memResource);
}