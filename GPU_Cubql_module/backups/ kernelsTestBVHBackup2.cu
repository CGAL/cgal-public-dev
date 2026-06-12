#define CUBQL_GPU_BUILDER_IMPLEMENTATION 1

// Core cuBQL headers
#include "cuBQL/bvh.h"
#include "cuBQL/builder/cuda.h"

// Your custom experimental builder layout
#include "include/third-party/cubql/sm_builder_v2.h"

// Custom traversal 
#include "include/third-party/cubql/fixedBoxQueryv2.h"

#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include "samples/common/loadOBJ.h"

// --------------------------------------------------------------------
// UPDATED EXPERIMENTAL TRAVERSAL KERNELS (INDIRECT PRIMITIVE FETCHING)
// --------------------------------------------------------------------
__global__ void countAABBOverlapsIndirectKernel(int *pairCounts, cuBQL::bvh3f bvhA, const cuBQL::Triangle* triA, const cuBQL::Triangle* triB, const uint32_t* batchPrimIDsB, int currentBatchSize, uint64_t startNodeIdx) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= currentBatchSize) return;

    // Fetch triangle indirectly using gathered sub-tree indices from Mesh B
    uint32_t primIDB = batchPrimIDsB[tid];
    cuBQL::Triangle b = triB[primIDB];
    cuBQL::box3f query = b.bounds();
    int count = 0;

    cuBQL::fixedBoxQueryv2::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        for (uint32_t i = 0; i < num; i++) {
            if (triA[ids[i]].bounds().overlaps(query)) {
                count++;
            }
        }
        return CUBQL_CONTINUE_TRAVERSAL;
    }, bvhA, query, startNodeIdx);

    pairCounts[tid] = count;
}

__global__ void fillAABBOverlapsIndirectKernel(int2 *candidatePairs, const int *offsets, cuBQL::bvh3f bvhA, const cuBQL::Triangle* triA, const cuBQL::Triangle* triB, const uint32_t* batchPrimIDsB, int currentBatchSize, uint64_t startNodeIdx) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= currentBatchSize) return;

    int wPos = offsets[tid];
    uint32_t primIDB = batchPrimIDsB[tid];
    cuBQL::Triangle b = triB[primIDB];
    cuBQL::box3f query = b.bounds();

    cuBQL::fixedBoxQueryv2::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        for (uint32_t i = 0; i < num; i++) {
            if (triA[ids[i]].bounds().overlaps(query)) {
                // Store actual global triangle indices relative to each mesh array
                candidatePairs[wPos++] = make_int2((int)ids[i], (int)primIDB);
            }
        }
        return CUBQL_CONTINUE_TRAVERSAL;
    }, bvhA, query, startNodeIdx);
}

// --------------------------------------------------------------------
// CORE HIERARCHY MAINTENANCE KERNELS
// --------------------------------------------------------------------
__global__ void generateBoxes(cuBQL::box3f* boxes, const cuBQL::Triangle* tris, int N) {
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if(i < N) {
    boxes[i] = tris[i].bounds();
  }
}

typedef typename cuBQL::BinaryBVH<float, 3>::Node BvhNodeT;

struct FrontierElement
{
  uint32_t nodeID;     
  uint32_t targetSlot; 
};

__device__ inline bool boxBoxesIntersect(const cuBQL::box3f& a, const cuBQL::box3f& b) {
  return (a.lower[0] <= b.upper[0] && a.upper[0] >= b.lower[0]) &&
         (a.lower[1] <= b.upper[1] && a.upper[1] >= b.lower[1]) &&
         (a.lower[2] <= b.upper[2] && a.upper[2] >= b.lower[2]);
}

__device__ inline void atomicMinFloat(float* addr, float val) {
  int* addr_as_int = (int*)addr;
  int old = *addr_as_int;
  int assumed;
  do {
    assumed = old;
    old = atomicCAS(addr_as_int, assumed, __float_as_int(fminf(__int_as_float(assumed), val)));
  } while(assumed != old);
}

__device__ inline void atomicMaxFloat(float* addr, float val) {
  int* addr_as_int = (int*)addr;
  int old = *addr_as_int;
  int assumed;
  do {
    assumed = old;
    old = atomicCAS(addr_as_int, assumed, __float_as_int(fmaxf(__int_as_float(assumed), val)));
  } while(assumed != old);
}

__global__ void processFrontierLayer(const FrontierElement* d_inFrontier,
                                     uint32_t inCount,
                                     FrontierElement* d_outFrontier,
                                     uint32_t* d_outCount,
                                     const BvhNodeT* d_nodes,
                                     uint32_t numNodes,
                                     const uint32_t* d_primIDs,
                                     uint32_t numPrims,
                                     cuBQL::box3f* d_targetBoxes,
                                     uint32_t* d_targetCounts,
                                     const cuBQL::box3f* d_primitiveBoxes) {
  uint32_t idx = threadIdx.x + blockIdx.x * blockDim.x;
  if(idx >= inCount)
    return;

  FrontierElement curr = d_inFrontier[idx];
  BvhNodeT node = d_nodes[curr.nodeID];

  if(node.admin.count > 0) {
    uint32_t baseOffset = node.admin.offset;
    atomicAdd(&d_targetCounts[curr.targetSlot], (uint32_t)node.admin.count);

    for(uint32_t p = 0; p < node.admin.count; ++p) {
      uint32_t actualPrimID = d_primIDs[baseOffset + p];
      cuBQL::box3f primBox = d_primitiveBoxes[actualPrimID];

      float* targetLower = (float*)&d_targetBoxes[curr.targetSlot].lower;
      float* targetUpper = (float*)&d_targetBoxes[curr.targetSlot].upper;

      atomicMinFloat(&targetLower[0], primBox.lower[0]);
      atomicMinFloat(&targetLower[1], primBox.lower[1]);
      atomicMinFloat(&targetLower[2], primBox.lower[2]);

      atomicMaxFloat(&targetUpper[0], primBox.upper[0]);
      atomicMaxFloat(&targetUpper[1], primBox.upper[1]);
      atomicMaxFloat(&targetUpper[2], primBox.upper[2]);
    }
  } else {
    uint32_t leftChild = node.admin.offset;
    if(leftChild < numNodes) {
      uint32_t slot1 = atomicAdd(d_outCount, 1);
      d_outFrontier[slot1] = {leftChild, curr.targetSlot};

      uint32_t slot2 = atomicAdd(d_outCount, 1);
      d_outFrontier[slot2] = {leftChild + 1, curr.targetSlot};
    }
  }
}

__global__ void intersectFrontiers(const cuBQL::box3f* d_boxesA, uint32_t countA,
                                   const BvhNodeT* d_nodesB, const uint32_t* d_indicesB, uint32_t countB,
                                   uint32_t* d_intersectionsCount) {
  uint32_t idxA = threadIdx.x + blockIdx.x * blockDim.x;
  if(idxA >= countA) return;

  cuBQL::box3f boxA = d_boxesA[idxA];
  uint32_t localIntersections = 0;

  for(uint32_t idxB = 0; idxB < countB; ++idxB) {
    uint32_t nodeIDB = d_indicesB[idxB];
    cuBQL::box3f boxB = d_nodesB[nodeIDB].bounds;

    if(boxBoxesIntersect(boxA, boxB)) {
      localIntersections++;
    }
  }

  if(localIntersections > 0) {
    atomicAdd(d_intersectionsCount, localIntersections);
  }
}

__global__ void collectIntersectingPairsKernel(const cuBQL::box3f* d_boxesA, uint32_t countA,
                                               const BvhNodeT* d_nodesB, const uint32_t* d_indicesB, uint32_t countB,
                                               int2* d_pairs, uint32_t* d_globalCounter) {
  uint32_t idxA = threadIdx.x + blockIdx.x * blockDim.x;
  if(idxA >= countA) return;

  cuBQL::box3f boxA = d_boxesA[idxA];

  for(uint32_t idxB = 0; idxB < countB; ++idxB) {
    uint32_t nodeIDB = d_indicesB[idxB];
    cuBQL::box3f boxB = d_nodesB[nodeIDB].bounds;

    if(boxBoxesIntersect(boxA, boxB)) {
      uint32_t slot = atomicAdd(d_globalCounter, 1);
      d_pairs[slot] = make_int2((int)idxA, (int)idxB);
    }
  }
}

__global__ void initTargetBoxes(cuBQL::box3f* boxes, uint32_t* counts, uint32_t count) {
  uint32_t idx = threadIdx.x + blockIdx.x * blockDim.x;
  if(idx < count) {
    counts[idx] = 0;
    boxes[idx].lower[0] = 1e30f;
    boxes[idx].lower[1] = 1e30f;
    boxes[idx].lower[2] = 1e30f;
    boxes[idx].upper[0] = -1e30f;
    boxes[idx].upper[1] = -1e30f;
    boxes[idx].upper[2] = -1e30f;
  }
}

// Host helper function to unroll and harvest all primitive IDs underneath any sub-hierarchy node of B
void gatherPrimitivesFromSubtree(uint32_t nodeID, const BvhNodeT* h_nodes, const uint32_t* h_primIDs, std::vector<uint32_t>& outPrims) {
    const BvhNodeT& node = h_nodes[nodeID];
    if (node.admin.count > 0) {
        for (uint32_t i = 0; i < node.admin.count; ++i) {
            outPrims.push_back(h_primIDs[node.admin.offset + i]);
        }
    } else {
        // Internal node traversal rule for BinaryBVH structures
        uint32_t leftChild = node.admin.offset;
        gatherPrimitivesFromSubtree(leftChild, h_nodes, h_primIDs, outPrims);
        gatherPrimitivesFromSubtree(leftChild + 1, h_nodes, h_primIDs, outPrims);
    }
}

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

  cuBQL::Triangle* dMeshA;
  CUBQL_CUDA_CALL(Malloc(&dMeshA, numTrianglesA * sizeof(cuBQL::Triangle)));
  CUBQL_CUDA_CALL(Memcpy(dMeshA, hMeshA, numTrianglesA * sizeof(cuBQL::Triangle), cudaMemcpyDefault));

  cuBQL::Triangle* dMeshB;
  CUBQL_CUDA_CALL(Malloc(&dMeshB, numTrianglesB * sizeof(cuBQL::Triangle)));
  CUBQL_CUDA_CALL(Memcpy(dMeshB, hMeshB, numTrianglesB * sizeof(cuBQL::Triangle), cudaMemcpyDefault));
  double tUploadEnd = cuBQL::getCurrentTime();

  // --------------------------------------------------------------------
  // HIERARCHY A: BUILD & REFIT
  // --------------------------------------------------------------------
  double tBuildAStart = cuBQL::getCurrentTime();
  cuBQL::box3f* dBoxesA;
  CUBQL_CUDA_CALL(Malloc(&dBoxesA, numTrianglesA * sizeof(cuBQL::box3f)));
  generateBoxes<<<cuBQL::divRoundUp(numTrianglesA, 256), 256>>>(dBoxesA, dMeshA, numTrianglesA);

  uint32_t maxPossibleNodesA = 2 * numTrianglesA;
  thrust::device_vector<uint32_t> d_markedNodeIndicesA(maxPossibleNodesA, 0);
  thrust::device_vector<uint32_t> d_markedNodeDescendantCountsA(maxPossibleNodesA, 0);
  uint32_t h_outMarkedCountA = 0;

  cuBQL::bvh3f bvhA;
  cuBQL::gpuBuilder_v2::build_custom(
      bvhA, dBoxesA, numTrianglesA, buildConfig, (uint32_t)maxCellSizeA, thrust::raw_pointer_cast(d_markedNodeIndicesA.data()),
      thrust::raw_pointer_cast(d_markedNodeDescendantCountsA.data()), &h_outMarkedCountA, stream, memResource);
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

  uint32_t maxPossibleNodesB = 2 * numTrianglesB;
  thrust::device_vector<uint32_t> d_markedNodeIndicesB(maxPossibleNodesB, 0);
  thrust::device_vector<uint32_t> d_markedNodeDescendantCountsB(maxPossibleNodesB, 0);
  uint32_t h_outMarkedCountB = 0;

  cuBQL::bvh3f bvhB;
  cuBQL::gpuBuilder_v2::build_custom(
      bvhB, dBoxesB, numTrianglesB, buildConfig, (uint32_t)maxCellSizeB, thrust::raw_pointer_cast(d_markedNodeIndicesB.data()),
      thrust::raw_pointer_cast(d_markedNodeDescendantCountsB.data()), &h_outMarkedCountB, stream, memResource);
  cuBQL::cuda::refit(bvhB, dBoxesB, stream);
  cudaDeviceSynchronize();
  double tBuildBEnd = cuBQL::getCurrentTime();

  // --------------------------------------------------------------------
  // RAPID DESCENT: EXECUTED FOR MESH A ONLY
  // --------------------------------------------------------------------
  double tGpuBfsStart = cuBQL::getCurrentTime();
  thrust::device_vector<cuBQL::box3f> d_bfsCollectedBoxesA(h_outMarkedCountA);
  thrust::device_vector<uint32_t> d_bfsCollectedCountsA(h_outMarkedCountA);

  if(h_outMarkedCountA > 0) {
    initTargetBoxes<<<cuBQL::divRoundUp(h_outMarkedCountA, (uint32_t)64), 64>>>(
        thrust::raw_pointer_cast(d_bfsCollectedBoxesA.data()), thrust::raw_pointer_cast(d_bfsCollectedCountsA.data()),
        h_outMarkedCountA);

    uint32_t frontierAllocationSize = std::max((uint32_t)numTrianglesA, (uint32_t)32768);
    thrust::device_vector<FrontierElement> d_frontierA(frontierAllocationSize);
    thrust::device_vector<FrontierElement> d_frontierB(frontierAllocationSize);
    thrust::device_vector<uint32_t> d_atomicCounterPool(1, 0);

    std::vector<FrontierElement> h_seedFrontier(h_outMarkedCountA);
    std::vector<uint32_t> h_tempIndicesA(h_outMarkedCountA);
    thrust::copy(d_markedNodeIndicesA.begin(), d_markedNodeIndicesA.begin() + h_outMarkedCountA, h_tempIndicesA.begin());

    for(uint32_t i = 0; i < h_outMarkedCountA; ++i) {
      h_seedFrontier[i] = {h_tempIndicesA[i], (uint32_t)i};
    }
    thrust::copy(h_seedFrontier.begin(), h_seedFrontier.end(), d_frontierA.begin());

    uint32_t currentFrontierSize = h_outMarkedCountA;
    FrontierElement* p_inBuffer = thrust::raw_pointer_cast(d_frontierA.data());
    FrontierElement* p_outBuffer = thrust::raw_pointer_cast(d_frontierB.data());

    while(currentFrontierSize > 0) {
      d_atomicCounterPool[0] = 0;
      uint32_t threadsPerBlock = 256;
      uint32_t blockCount = cuBQL::divRoundUp(currentFrontierSize, (uint32_t)threadsPerBlock);

      processFrontierLayer<<<blockCount, threadsPerBlock>>>(
          p_inBuffer, currentFrontierSize, p_outBuffer, thrust::raw_pointer_cast(d_atomicCounterPool.data()), bvhA.nodes,
          bvhA.numNodes, bvhA.primIDs, bvhA.numPrims, thrust::raw_pointer_cast(d_bfsCollectedBoxesA.data()),
          thrust::raw_pointer_cast(d_bfsCollectedCountsA.data()), dBoxesA);
      cudaDeviceSynchronize();

      currentFrontierSize = d_atomicCounterPool[0];
      std::swap(p_inBuffer, p_outBuffer);
    }
  }
  double tGpuBfsEnd = cuBQL::getCurrentTime();

  // --------------------------------------------------------------------
  // GPU PARALLEL CRISS-CROSS-APPLE INTERSECTION TEST
  // --------------------------------------------------------------------
  double tCrossStart = cuBQL::getCurrentTime();
  thrust::device_vector<uint32_t> d_intersectionCounter(1, 0);

  if(h_outMarkedCountA > 0 && h_outMarkedCountB > 0) {
    uint32_t threadsCross = 256;
    uint32_t blocksCross = cuBQL::divRoundUp(h_outMarkedCountA, threadsCross);

    intersectFrontiers<<<blocksCross, threadsCross>>>(
        thrust::raw_pointer_cast(d_bfsCollectedBoxesA.data()), h_outMarkedCountA,
        bvhB.nodes, thrust::raw_pointer_cast(d_markedNodeIndicesB.data()), h_outMarkedCountB,
        thrust::raw_pointer_cast(d_intersectionCounter.data()));
    cudaDeviceSynchronize();
  }
  double tCrossEnd = cuBQL::getCurrentTime();

  uint32_t totalIntersections = d_intersectionCounter[0];
  uint64_t totalPossiblePairs = (uint64_t)h_outMarkedCountA * h_outMarkedCountB;
  double intersectionPercentage = totalPossiblePairs > 0 ? ((double)totalIntersections / totalPossiblePairs) * 100.0 : 0.0;

  // --------------------------------------------------------------------
  // UPDATED SUBTREE BATCH EXPERIMENTATION PIPELINE
  // --------------------------------------------------------------------
  std::cout << "\n--------------------------------------------------\n";
  std::cout << " LAUNCHING TRACKED BATCH SUB-TREE TESTS\n";
  std::cout << "--------------------------------------------------\n";

  double tBatchExperimentTotal = 0.0;
  uint64_t totalFillsAccumulated = 0;

  if (totalIntersections > 0) {
    thrust::device_vector<int2> d_intersectingPairs(totalIntersections);
    thrust::device_vector<uint32_t> d_pairGatherCounter(1, 0);

    uint32_t threadsGather = 256;
    uint32_t blocksGather = cuBQL::divRoundUp(h_outMarkedCountA, threadsGather);

    collectIntersectingPairsKernel<<<blocksGather, threadsGather>>>(
        thrust::raw_pointer_cast(d_bfsCollectedBoxesA.data()), h_outMarkedCountA,
        bvhB.nodes, thrust::raw_pointer_cast(d_markedNodeIndicesB.data()), h_outMarkedCountB,
        thrust::raw_pointer_cast(d_intersectingPairs.data()), thrust::raw_pointer_cast(d_pairGatherCounter.data()));
    cudaDeviceSynchronize();

    std::vector<int2> h_pairs(totalIntersections);
    thrust::copy(d_intersectingPairs.begin(), d_intersectingPairs.end(), h_pairs.begin());

    std::vector<uint32_t> h_indicesA(h_outMarkedCountA);
    std::vector<uint32_t> h_indicesB(h_outMarkedCountB);
    thrust::copy(d_markedNodeIndicesA.begin(), d_markedNodeIndicesA.begin() + h_outMarkedCountA, h_indicesA.begin());
    thrust::copy(d_markedNodeIndicesB.begin(), d_markedNodeIndicesB.begin() + h_outMarkedCountB, h_indicesB.begin());

    // Pull down tree structural buffers from B to allow full traversal unrolling on host
    std::vector<BvhNodeT> h_nodesB(bvhB.numNodes);
    std::vector<uint32_t> h_primIDsB(bvhB.numPrims);
    CUBQL_CUDA_CALL(Memcpy(h_nodesB.data(), bvhB.nodes, bvhB.numNodes * sizeof(BvhNodeT), cudaMemcpyDeviceToHost));
    CUBQL_CUDA_CALL(Memcpy(h_primIDsB.data(), bvhB.primIDs, bvhB.numPrims * sizeof(uint32_t), cudaMemcpyDeviceToHost));

    // Profile each structural pair interaction as an isolated batch
    for (uint32_t batchID = 0; batchID < totalIntersections; ++batchID) {
        int2 structuralPair = h_pairs[batchID];
        uint32_t localFrontierIdxA = structuralPair.x; 
        uint32_t localFrontierIdxB = structuralPair.y; 

        uint32_t nodeIDA = h_indicesA[localFrontierIdxA];
        uint32_t nodeIDB = h_indicesB[localFrontierIdxB];

        // Gather all primitive indices underneath macro-cell nodeIDB
        std::vector<uint32_t> h_batchPrimIDsB;
        gatherPrimitivesFromSubtree(nodeIDB, h_nodesB.data(), h_primIDsB.data(), h_batchPrimIDsB);
        
        uint32_t activeBatchSizeB = h_batchPrimIDsB.size();
        if (activeBatchSizeB == 0) continue;

        double tBatchRunStart = cuBQL::getCurrentTime();

        // Upload collected primitive lists for this specific batch execution
        thrust::device_vector<uint32_t> d_batchPrimIDsB = h_batchPrimIDsB;
        thrust::device_vector<int> d_pairCounts(activeBatchSizeB, 0);
        
        countAABBOverlapsIndirectKernel<<<cuBQL::divRoundUp(activeBatchSizeB, 256u), 256>>>(
            thrust::raw_pointer_cast(d_pairCounts.data()), bvhA, dMeshA, dMeshB, 
            thrust::raw_pointer_cast(d_batchPrimIDsB.data()), activeBatchSizeB, (uint64_t)nodeIDA
        );
        cudaDeviceSynchronize();

        thrust::device_vector<int> d_offsets(activeBatchSizeB);
        thrust::exclusive_scan(d_pairCounts.begin(), d_pairCounts.end(), d_offsets.begin());
        
        int totalPairsInBatch = d_offsets[activeBatchSizeB - 1] + d_pairCounts[activeBatchSizeB - 1];
        totalFillsAccumulated += totalPairsInBatch;

        if (totalPairsInBatch > 0) {
            thrust::device_vector<int2> d_candidatePairs(totalPairsInBatch);

            fillAABBOverlapsIndirectKernel<<<cuBQL::divRoundUp(activeBatchSizeB, 256u), 256>>>(
                thrust::raw_pointer_cast(d_candidatePairs.data()), thrust::raw_pointer_cast(d_offsets.data()),
                bvhA, dMeshA, dMeshB, thrust::raw_pointer_cast(d_batchPrimIDsB.data()), activeBatchSizeB, (uint64_t)nodeIDA
            );
            cudaDeviceSynchronize();
        }
        
        double tBatchRunEnd = cuBQL::getCurrentTime();
        tBatchExperimentTotal += (tBatchRunEnd - tBatchRunStart) * 1000.0;
    }
  }

  // --------------------------------------------------------------------
  // PROPORTIONS & PERFORMANCE OUTPUT
  // --------------------------------------------------------------------
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
  std::cout << " LEAF BATCH EXPERIMENT DIAGNOSTICS\n";
  std::cout << "--------------------------------------------------\n";
  std::cout << "Total Active Batches Evaluated: " << totalIntersections << "\n";
  std::cout << "Total Intersecting Leaf Pairs:  " << totalFillsAccumulated << "\n";
  std::cout << "Total Leaf Traversal Compute:   " << tBatchExperimentTotal << " ms\n\n";

  std::cout << "--------------------------------------------------\n";
  std::cout << " TIMING METRICS OVERVIEW\n";
  std::cout << "--------------------------------------------------\n";
  std::cout << "  |- Upload time:              " << (tUploadEnd - tUploadStart) * 1000.0 << " ms\n";
  std::cout << "  |- Build + Refit (Mesh A):  " << (tBuildAEnd - tBuildAStart) * 1000.0 << " ms\n";
  std::cout << "  |- Build + Refit (Mesh B):  " << (tBuildBEnd - tBuildBStart) * 1000.0 << " ms\n";
  std::cout << "  |- Parallel BFS Descent (A): " << (tGpuBfsEnd - tGpuBfsStart) * 1000.0 << " ms\n";
  std::cout << "  |- GPU Cross-Check Engine:  " << (tCrossEnd - tCrossStart) * 1000.0 << " ms\n";
  std::cout << "--------------------------------------------------\n";

  // Free allocations
  CUBQL_CUDA_CALL(Free(dMeshA));
  CUBQL_CUDA_CALL(Free(dBoxesA));
  CUBQL_CUDA_CALL(Free(dMeshB));
  CUBQL_CUDA_CALL(Free(dBoxesB));
  cuBQL::cuda::free(bvhA, stream, memResource);
  cuBQL::cuda::free(bvhB, stream, memResource);
}