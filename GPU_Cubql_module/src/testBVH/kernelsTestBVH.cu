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

#include <thrust/sort.h>
#include <thrust/unique.h>

#include "prune_pipeline.h"

#include "global_box.h"


// --------------------------------------------------------------------
// LOCAL STREAM-SAFE REPLACEMENTS FOR INTERNAL CUBQL ALLOCATORS
// --------------------------------------------------------------------
#ifndef _ALLOC
#define _ALLOC(ptr, count, stream, memResource)                                                                        \
  CUBQL_CUDA_CALL(MallocAsync((void**)&(ptr), (size_t)(count) * sizeof(*(ptr)), stream))
#endif

#ifndef _FREE
#define _FREE(ptr, stream, memResource) CUBQL_CUDA_CALL(FreeAsync((ptr), stream))
#endif

// --------------------------------------------------------------------
// EXISTING KERNELS & HELPERS
// --------------------------------------------------------------------
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


template <typename T, int D>
void testFastBVH(const cuBQL::box_t<T, D>* d_boxesA,
                 int numPrimsA,
                 int maxCellSizeA,
                 const cuBQL::box_t<T, D>* d_boxesB,
                 int numPrimsB,
                 int maxCellSizeB,
                 cuBQL::box_t<T, D> globalBoxA,
                 cuBQL::box_t<T, D> globalBoxB,
                 cudaStream_t s,
                 cuBQL::DeviceMemoryResource& memResource,
                 int batchMultiplier,
                 int userleafThreshold) {
  std::cout << "\n==================================================" << std::endl;
  std::cout << " [DEEP DEBUG] => Entering ZIMBLIFIED testFastBVH..." << std::endl;
  std::cout << "==================================================" << std::endl;

  // Check raw pointers and parameters right away
  std::cout << " -> Device Pointer A : " << d_boxesA << " | Prims: " << numPrimsA << " | MaxCell: " << maxCellSizeA
            << std::endl;
  std::cout << " -> Device Pointer B : " << d_boxesB << " | Prims: " << numPrimsB << " | MaxCell: " << maxCellSizeB
            << std::endl;

  if(d_boxesA == d_boxesB) {
    std::cout << " [CRITICAL WARNING] d_boxesA and d_boxesB point to the EXACT SAME memory address!" << std::endl;
  } else {
    std::cout << " [OK] Input pointers A and B are distinct memory blocks." << std::endl;
  }

  if(maxCellSizeA > 20 || maxCellSizeB > 20) {
    std::cerr << "[CRITICAL ERROR] numIterations is too high! Aborting testFastBVH." << std::endl;
    return;
  }

  // --------------------------------------------------------------------
  // PARALLEL CENTROID BOUNDING BOX CALCULATION & TIMING (TEMPORAL VARIABLES)
  // --------------------------------------------------------------------
  std::cout << "\n -> [COMPUTE] Calculating parallel centroid bounding boxes into temporal variables..." << std::endl;

  // === MESH A PROCESSING ===
  cuBQL::box_t<T, D> globalBoxANew;
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double tGBoxAStart = cuBQL::getCurrentTime();
  cuBQL::utils::computeGlobalBoxParallel<T, D>(globalBoxANew, d_boxesA, numPrimsA, s, memResource);
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double tGBoxAMs = (cuBQL::getCurrentTime() - tGBoxAStart) * 1000.0;
  std::cout << "    * Mesh A Centroid Box Time: " << tGBoxAMs << " ms" << std::endl;

  // Added: Print Dimensions for Mesh A (Old vs New)
  std::cout << "    [DIMENSIONS MESH A - OLD]:" << std::endl;
  for(int d = 0; d < D; ++d) {
    T size = globalBoxA.upper[d] - globalBoxA.lower[d];
    std::cout << "        Dim " << d << " -> Bounds: [" << globalBoxA.lower[d] << ", " << globalBoxA.upper[d]
              << "] | Size: " << size << std::endl;
  }
  std::cout << "    [DIMENSIONS MESH A - NEW]:" << std::endl;
  for(int d = 0; d < D; ++d) {
    T size = globalBoxANew.upper[d] - globalBoxANew.lower[d];
    std::cout << "        Dim " << d << " -> Bounds: [" << globalBoxANew.lower[d] << ", " << globalBoxANew.upper[d]
              << "] | Size: " << size << std::endl;
  }
  std::cout << "--------------------------------------------------" << std::endl;

  // === MESH B PROCESSING ===
  cuBQL::box_t<T, D> globalBoxBNew;
  double tGBoxBStart = cuBQL::getCurrentTime();
  cuBQL::utils::computeGlobalBoxParallel<T, D>(globalBoxBNew, d_boxesB, numPrimsB, s, memResource);
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double tGBoxBMs = (cuBQL::getCurrentTime() - tGBoxBStart) * 1000.0;
  std::cout << "    * Mesh B Centroid Box Time: " << tGBoxBMs << " ms" << std::endl;

  // Added: Print Dimensions for Mesh B (Old vs New)
  std::cout << "    [DIMENSIONS MESH B - OLD]:" << std::endl;
  for(int d = 0; d < D; ++d) {
    T size = globalBoxB.upper[d] - globalBoxB.lower[d];
    std::cout << "        Dim " << d << " -> Bounds: [" << globalBoxB.lower[d] << ", " << globalBoxB.upper[d]
              << "] | Size: " << size << std::endl;
  }
  std::cout << "    [DIMENSIONS MESH B - NEW]:" << std::endl;
  for(int d = 0; d < D; ++d) {
    T size = globalBoxBNew.upper[d] - globalBoxBNew.lower[d];
    std::cout << "        Dim " << d << " -> Bounds: [" << globalBoxBNew.lower[d] << ", " << globalBoxBNew.upper[d]
              << "] | Size: " << size << std::endl;
  }
  std::cout << "==================================================" << std::endl;

  double tStart = cuBQL::getCurrentTime();

  // --- MESH A RUN (Zimblified Edition) ---
  cuBQL::box_t<T, D>* outNodeBoxesA = nullptr;
  uint32_t* outSortedPrimIDsA = nullptr;
  uint32_t* outNodeOffsetsA = nullptr;
  uint32_t outTotalActiveCellsA = 0;

  // Profile Mesh A Initialization
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double tInitAStart = cuBQL::getCurrentTime();

  // globalBoxANew
  cuBQL::gpuBuilder_v3::test_speedrun_initialization_linear(d_boxesA, numPrimsA, (uint32_t)maxCellSizeA, globalBoxANew,
                                                            outNodeBoxesA, outSortedPrimIDsA, outNodeOffsetsA,
                                                            outTotalActiveCellsA, s, memResource);
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

  // globalBoxBNew
  cuBQL::gpuBuilder_v3::test_speedrun_initialization_linear(d_boxesB, numPrimsB, (uint32_t)maxCellSizeB, globalBoxBNew,
                                                            outNodeBoxesB, outSortedPrimIDsB, outNodeOffsetsB,
                                                            outTotalActiveCellsB, s, memResource);
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
  std::vector<cuBQL::box_t<T, D>> h_boxesA(outTotalActiveCellsA);
  std::vector<cuBQL::box_t<T, D>> h_boxesB(outTotalActiveCellsB);
  std::vector<uint32_t> h_offsetsA(outTotalActiveCellsA + 1);
  std::vector<uint32_t> h_offsetsB(outTotalActiveCellsB + 1);

  CUBQL_CUDA_CALL(Memcpy(h_boxesA.data(), outNodeBoxesA, outTotalActiveCellsA * sizeof(cuBQL::box_t<T, D>),
                         cudaMemcpyDeviceToHost));
  CUBQL_CUDA_CALL(Memcpy(h_boxesB.data(), outNodeBoxesB, outTotalActiveCellsB * sizeof(cuBQL::box_t<T, D>),
                         cudaMemcpyDeviceToHost));
  CUBQL_CUDA_CALL(Memcpy(h_offsetsA.data(), outNodeOffsetsA, (outTotalActiveCellsA + 1) * sizeof(uint32_t),
                         cudaMemcpyDeviceToHost));
  CUBQL_CUDA_CALL(Memcpy(h_offsetsB.data(), outNodeOffsetsB, (outTotalActiveCellsB + 1) * sizeof(uint32_t),
                         cudaMemcpyDeviceToHost));

  // Calculate true primitive allocations via offsets
  uint32_t sumPrimsA = h_offsetsA[outTotalActiveCellsA] - h_offsetsA[0];
  uint32_t sumPrimsB = h_offsetsB[outTotalActiveCellsB] - h_offsetsB[0];

  std::cout << "--------------------------------------------------" << std::endl;
  std::cout << " [POPULATED CELLS CHECK SUMMARY]" << std::endl;
  std::cout << " -> Tree A Active Nodes: " << outTotalActiveCellsA << " | Total Verified Prims = " << sumPrimsA
            << std::endl;
  std::cout << " -> Tree B Active Nodes: " << outTotalActiveCellsB << " | Total Verified Prims = " << sumPrimsB
            << std::endl;
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

  uint32_t totalIntersections =
      executeBoxCrossCheck<T, D>(outNodeBoxesA, outTotalActiveCellsA, outNodeBoxesB, outTotalActiveCellsB,
                                 d_intersectPairsA, d_intersectPairsB, s);
  CUBQL_CUDA_CALL(StreamSynchronize(s));

  double crossCheckEnd = cuBQL::getCurrentTime();
  double crossCheckMs = (crossCheckEnd - crossCheckStart) * 1000.0;

  // Save historical metrics prior to post-build pruning actions
  uint32_t initialCellsA = outTotalActiveCellsA;
  uint32_t initialPrimsA = sumPrimsA;
  uint32_t initialCellsB = outTotalActiveCellsB;
  uint32_t initialPrimsB = sumPrimsB;

  int currentPrimsNumA = (int)sumPrimsA;
  int currentPrimsNumB = (int)sumPrimsB;

// --------------------------------------------------------------------
  // PARALLEL STREAM COMPACTION & PRUNING ALGORITHM (EXECUTION DEFERRED TO POST-BUILD)
  // --------------------------------------------------------------------
  std::cout << " [PIPELINE PRUNING] => Dispatched grid streaming parallel re-index compaction (POST-BUILD)..."
            << std::endl;

  // Track vector sizing state BEFORE reindexing
  size_t sizeBeforeA = d_intersectPairsA.size();
  size_t sizeBeforeB = d_intersectPairsB.size();
  std::cout << " -> [SIZE BEFORE] d_intersectPairsA size: " << sizeBeforeA 
            << " | d_intersectPairsB size: " << sizeBeforeB << std::endl;

  // Track and execute Compaction Pipeline for Mesh A
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double pruneStartA = cuBQL::getCurrentTime();
  parallelPruneAndReindexAll(thrust::raw_pointer_cast(d_intersectPairsA.data()), totalIntersections, outSortedPrimIDsA,
                             outNodeOffsetsA, outTotalActiveCellsA, currentPrimsNumA, s, memResource);
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double pruneMsA = (cuBQL::getCurrentTime() - pruneStartA) * 1000.0;

  // Track and execute Compaction Pipeline for Mesh B
  double pruneStartB = cuBQL::getCurrentTime();
  parallelPruneAndReindexAll(thrust::raw_pointer_cast(d_intersectPairsB.data()), totalIntersections, outSortedPrimIDsB,
                             outNodeOffsetsB, outTotalActiveCellsB, currentPrimsNumB, s, memResource);
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double pruneMsB = (cuBQL::getCurrentTime() - pruneStartB) * 1000.0;

  // Track vector sizing state AFTER reindexing
  size_t sizeAfterA = d_intersectPairsA.size();
  size_t sizeAfterB = d_intersectPairsB.size();
  std::cout << " -> [SIZE AFTER]  d_intersectPairsA size: " << sizeAfterA 
            << " | d_intersectPairsB size: " << sizeAfterB << std::endl;

  // --------------------------------------------------------------------
  // EXTRA DEBUG: PRINT INTERSECTING CELL INDEX RANGES
  // --------------------------------------------------------------------
  std::cout << "\n -> [PAIR RANGE CHECK] Analyzing intersection pair arrays (Total: " << totalIntersections << ")..." << std::endl;
  
  if (totalIntersections > 0) {
    // Extract Mesh A overlap bounds
    auto mmA = thrust::minmax_element(thrust::device, d_intersectPairsA.begin(), d_intersectPairsA.end());
    uint32_t minA = *mmA.first;
    uint32_t maxA = *mmA.second;
    std::cout << "    * d_intersectPairsA => Min Cell ID: " << minA << " | Max Cell ID: " << maxA << std::endl;

    // Extract Mesh B overlap bounds
    auto mmB = thrust::minmax_element(thrust::device, d_intersectPairsB.begin(), d_intersectPairsB.end());
    uint32_t minB = *mmB.first;
    uint32_t maxB = *mmB.second;
    std::cout << "    * d_intersectPairsB => Min Cell ID: " << minB << " | Max Cell ID: " << maxB << std::endl;
  } else {
    std::cout << "    * [NOTICE] No intersecting pairs found. Skipping bounds check." << std::endl;
  }
  std::cout << "--------------------------------------------------" << std::endl;

  // --------------------------------------------------------------------
  // --------------------------------------------------------------------
  // GPU BUILDER V4 FOREST EXPANSION
  // --------------------------------------------------------------------
  std::cout << " [FOREST EXPANSION] => Launching Level-by-Level Parallel Sub-Tree Compilation..." << std::endl;

  // Calculate EXACT maxNodes using the PRUNED current counts to match internal allocations exactly
  // Change host sizing to use the full stride maximum to prevent early device clipping
  const uint32_t numInitA = outTotalActiveCellsA + (outTotalActiveCellsA & 1u);
  const uint32_t maxNodesA = 2u * (uint32_t)numPrimsA + numInitA + 2u; // Use numPrimsA here!

  const uint32_t numInitB = outTotalActiveCellsB + (outTotalActiveCellsB & 1u);
  const uint32_t maxNodesB = 2u * (uint32_t)numPrimsB + numInitB + 2u; // Use numPrimsB here!

  // Allocate target-level diagnostics with the corrected uniform sizing
  thrust::device_vector<uint32_t> d_nodeDescendantCountsA(maxNodesA, 0);
  thrust::device_vector<uint32_t> d_nodeDescendantCountsB(maxNodesB, 0);

  // Configure build options
  cuBQL::BuildConfig buildConfig;
  buildConfig.makeLeafThreshold = userleafThreshold;

  // --- Compile Forest A ---
  // FIX: Parameter sequence changed to match (bvh, boxes, currentPrimsNumA, numPrimsA, outTotalActiveCellsA...)
  cuBQL::BinaryBVH<T, D> bvhA;
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double forestAStart = cuBQL::getCurrentTime();

  cuBQL::gpuBuilder_v4::build_forest<T, D>(bvhA, d_boxesA, currentPrimsNumA, numPrimsA, outTotalActiveCellsA,
                                           outSortedPrimIDsA, outNodeOffsetsA, buildConfig,
                                           thrust::raw_pointer_cast(d_nodeDescendantCountsA.data()), s, memResource);
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double forestAMs = (cuBQL::getCurrentTime() - forestAStart) * 1000.0;

  // --- Compile Forest B ---
  // FIX: Parameter sequence changed to match (bvh, boxes, currentPrimsNumB, numPrimsB, outTotalActiveCellsB...)
  cuBQL::BinaryBVH<T, D> bvhB;
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double forestBStart = cuBQL::getCurrentTime();

  cuBQL::gpuBuilder_v4::build_forest<T, D>(bvhB, d_boxesB, currentPrimsNumB, numPrimsB, outTotalActiveCellsB,
                                           outSortedPrimIDsB, outNodeOffsetsB, buildConfig,
                                           thrust::raw_pointer_cast(d_nodeDescendantCountsB.data()), s, memResource);
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double forestBMs = (cuBQL::getCurrentTime() - forestBStart) * 1000.0;

  // --------------------------------------------------------------------
  // STRUCTURAL INTEGRITY & DESCENDANT COUNT CROSS-CHECKS
  // --------------------------------------------------------------------
  std::cout << "\n==================================================" << std::endl;
  std::cout << " [RUNNING METRIC CHECKS: FOREST STRUCTURE & PASSES]" << std::endl;
  std::cout << "==================================================" << std::endl;

  // --- CHECK FOR TREE A ---
  std::vector<typename cuBQL::BinaryBVH<T, D>::Node> h_bvhNodesA(bvhA.numNodes);
  std::vector<uint32_t> h_descendantsA(bvhA.numNodes);

  CUBQL_CUDA_CALL(Memcpy(h_bvhNodesA.data(), bvhA.nodes, bvhA.numNodes * sizeof(typename cuBQL::BinaryBVH<T, D>::Node),
                         cudaMemcpyDeviceToHost));
  CUBQL_CUDA_CALL(Memcpy(h_descendantsA.data(), thrust::raw_pointer_cast(d_nodeDescendantCountsA.data()),
                         bvhA.numNodes * sizeof(uint32_t), cudaMemcpyDeviceToHost));

  uint32_t totalRootDescendantsA = 0;
  for(uint32_t i = 0; i < outTotalActiveCellsA; ++i) {
    totalRootDescendantsA += h_descendantsA[i];
  }

  std::cout << " [TREE A DATA] Total Nodes in Array: " << bvhA.numNodes << " | Root Cells: " << outTotalActiveCellsA
            << std::endl;
  std::cout << "   -> TEST 1 (Sum of Roots): Total Prims accounted for at Root Level = " << totalRootDescendantsA
            << " (Expected Surviving: " << currentPrimsNumA << ")" << std::endl;

  uint32_t verifiedTreePrimsA = 0;
  uint32_t brokenLinksA = 0;
  uint32_t leafCountA = 0;

  for(uint32_t rootID = 0; rootID < outTotalActiveCellsA; ++rootID) {
    std::vector<uint32_t> stack;
    stack.push_back(rootID);

    while(!stack.empty()) {
      uint32_t currID = stack.back();
      stack.pop_back();

      if(currID >= bvhA.numNodes) {
        brokenLinksA++;
        continue;
      }

      auto& node = h_bvhNodesA[currID];
      uint32_t childOffset = node.admin.offset;
      uint32_t primCount = node.admin.count;

      if(childOffset == (uint32_t)-1 || primCount > 0) {
        verifiedTreePrimsA += primCount;
        leafCountA++;
      } else {
        uint32_t leftChild = childOffset;
        uint32_t rightChild = childOffset + 1;

        if(leftChild < bvhA.numNodes && rightChild < bvhA.numNodes) {
          uint32_t sumChildren = h_descendantsA[leftChild] + h_descendantsA[rightChild];
          if(h_descendantsA[currID] != sumChildren && h_descendantsA[currID] != 0) {
            static int alertsA = 0;
            if(alertsA++ < 5) {
              std::cout << "      [MISMATCH ALERT A] Node " << currID << " says " << h_descendantsA[currID]
                        << " descendants, but children sum to " << sumChildren << std::endl;
            }
          }
        }
        stack.push_back(rightChild);
        stack.push_back(leftChild);
      }
    }
  }

  std::cout << "   -> TEST 2 (Full Hierarchy Walk): Traversed " << leafCountA << " leaves. Found " << verifiedTreePrimsA
            << " primitives inside the tree structures. (Broken links: " << brokenLinksA << ")" << std::endl;
  if(verifiedTreePrimsA != (uint32_t)currentPrimsNumA) {
    std::cout << "   !!! [CRITICAL] Tree A layout is structurally dropping primitives during top-down splits !!!"
              << std::endl;
  }

  // --- CHECK FOR TREE B ---
  std::vector<typename cuBQL::BinaryBVH<T, D>::Node> h_bvhNodesB(bvhB.numNodes);
  std::vector<uint32_t> h_descendantsB(bvhB.numNodes);

  CUBQL_CUDA_CALL(Memcpy(h_bvhNodesB.data(), bvhB.nodes, bvhB.numNodes * sizeof(typename cuBQL::BinaryBVH<T, D>::Node),
                         cudaMemcpyDeviceToHost));
  CUBQL_CUDA_CALL(Memcpy(h_descendantsB.data(), thrust::raw_pointer_cast(d_nodeDescendantCountsB.data()),
                         bvhB.numNodes * sizeof(uint32_t), cudaMemcpyDeviceToHost));

  uint32_t totalRootDescendantsB = 0;
  for(uint32_t i = 0; i < outTotalActiveCellsB; ++i) {
    totalRootDescendantsB += h_descendantsB[i];
  }

  std::cout << " [TREE B DATA] Total Nodes in Array: " << bvhB.numNodes << " | Root Cells: " << outTotalActiveCellsB
            << std::endl;
  std::cout << "   -> TEST 1 (Sum of Roots): Total Prims accounted for at Root Level = " << totalRootDescendantsB
            << " (Expected Surviving: " << currentPrimsNumB << ")" << std::endl;

  uint32_t verifiedTreePrimsB = 0;
  uint32_t brokenLinksB = 0;
  uint32_t leafCountB = 0;

  for(uint32_t rootID = 0; rootID < outTotalActiveCellsB; ++rootID) {
    std::vector<uint32_t> stack;
    stack.push_back(rootID);

    while(!stack.empty()) {
      uint32_t currID = stack.back();
      stack.pop_back();

      if(currID >= bvhB.numNodes) {
        brokenLinksB++;
        continue;
      }

      auto& node = h_bvhNodesB[currID];
      uint32_t childOffset = node.admin.offset;
      uint32_t primCount = node.admin.count;

      if(childOffset == (uint32_t)-1 || primCount > 0) {
        verifiedTreePrimsB += primCount;
        leafCountB++;
      } else {
        uint32_t leftChild = childOffset;
        uint32_t rightChild = childOffset + 1;

        if(leftChild < bvhB.numNodes && rightChild < bvhB.numNodes) {
          uint32_t sumChildren = h_descendantsB[leftChild] + h_descendantsB[rightChild];
          if(h_descendantsB[currID] != sumChildren && h_descendantsB[currID] != 0) {
            static int alertsB = 0;
            if(alertsB++ < 5) {
              std::cout << "      [MISMATCH ALERT B] Node " << currID << " says " << h_descendantsB[currID]
                        << " descendants, but children sum to " << sumChildren << std::endl;
            }
          }
        }
        stack.push_back(rightChild);
        stack.push_back(leftChild);
      }
    }
  }

  std::cout << "   -> TEST 2 (Full Hierarchy Walk): Traversed " << leafCountB << " leaves. Found " << verifiedTreePrimsB
            << " primitives inside the tree structures. (Broken links: " << brokenLinksB << ")" << std::endl;
  if(verifiedTreePrimsB != (uint32_t)currentPrimsNumB) {
    std::cout << "   !!! [CRITICAL] Tree B layout is structurally dropping primitives during top-down splits !!!"
              << std::endl;
  }
  std::cout << "==================================================\n" << std::endl;

  // --------------------------------------------------------------------
  // FINAL PERFORMANCE & RESULTS READOUT
  // --------------------------------------------------------------------
  uint64_t totalValidActivePairs = (uint64_t)initialCellsA * initialCellsB;
  double overlapPercentage = 0.0;
  if(totalValidActivePairs > 0) {
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
  std::cout << " -> Mesh A Structural Pruning (Post)  : " << pruneMsA << " ms" << std::endl;
  std::cout << "    * Active Cells Surviving: " << initialCellsA << " -> " << outTotalActiveCellsA
            << " (Dropped: " << (initialCellsA - outTotalActiveCellsA) << ")" << std::endl;
  std::cout << "    * Pack Prim Elements Size: " << initialPrimsA << " -> " << currentPrimsNumA << std::endl;
  std::cout << " ------------------------------------------------" << std::endl;
  std::cout << " -> Mesh B Structural Pruning (Post)  : " << pruneMsB << " ms" << std::endl;
  std::cout << "    * Active Cells Surviving: " << initialCellsB << " -> " << outTotalActiveCellsB
            << " (Dropped: " << (initialCellsB - outTotalActiveCellsB) << ")" << std::endl;
  std::cout << "    * Pack Prim Elements Size: " << initialPrimsB << " -> " << currentPrimsNumB << std::endl;
  std::cout << " ------------------------------------------------" << std::endl;
  std::cout << " -> Intersecting Pairs Detected     : " << totalIntersections << " pairs" << std::endl;
  std::cout << " -> True Active Matrix Evaluated    : " << totalValidActivePairs << " pairings" << std::endl;
  std::cout << " -> Percentage of Total Overlaps    : " << overlapPercentage << "%" << std::endl;
  std::cout << "==================================================\n" << std::endl;

  // --------------------------------------------------------------------
  // ADDED: RAPID DESCENT BFS INJECTION & PRIMITIVE COUNT VERIFICATION
  // --------------------------------------------------------------------
  std::cout << " [INJECTED DEBUG] Running executeRapidDescentBFS on v4 Forest B..." << std::endl;

  thrust::device_vector<uint32_t> d_testMarkedNodeIndicesB(outTotalActiveCellsB);
  thrust::sequence(thrust::cuda::par.on(s), d_testMarkedNodeIndicesB.begin(), d_testMarkedNodeIndicesB.end(), 0u);

  thrust::device_vector<uint32_t> d_testOutOffsetsB;
  thrust::device_vector<uint32_t> d_testOutPrimsFlatB;

  double tBfsStart = cuBQL::getCurrentTime();
  executeRapidDescentBFS(bvhB, currentPrimsNumB, d_testMarkedNodeIndicesB, d_nodeDescendantCountsB,
                         outTotalActiveCellsB, d_testOutOffsetsB, d_testOutPrimsFlatB);
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double tBfsMs = (cuBQL::getCurrentTime() - tBfsStart) * 1000.0;

  size_t extractedPrimsB = d_testOutPrimsFlatB.size();
  std::cout << " -> [BFS RESULT] Extracted " << extractedPrimsB << " primitives over " << outTotalActiveCellsB
            << " roots in " << tBfsMs << " ms." << std::endl;

  if(extractedPrimsB == static_cast<size_t>(currentPrimsNumB)) {
    std::cout << " [SUCCESS] BFS Flat array matches expected surviving input size: " << extractedPrimsB
              << " == " << currentPrimsNumB << std::endl;
  } else {
    std::cerr << " [CRITICAL MISMATCH] BFS extraction size error! Expected Surviving: " << currentPrimsNumB
              << " | Got: " << extractedPrimsB << std::endl;
  }

  // --------------------------------------------------------------------
  // EXTRA DEBUG: EXTRACT MINIMUM AND MAXIMUM PRIMITIVE IDS
  // --------------------------------------------------------------------
  if (extractedPrimsB > 0) {
    auto result_pair = thrust::minmax_element(thrust::device, d_testOutPrimsFlatB.begin(), d_testOutPrimsFlatB.end());
    
    // Copy the values back to host from the device iterators safely
    uint32_t hostMinVal = *result_pair.first;
    uint32_t hostMaxVal = *result_pair.second;

    std::cout << " -> [VALUE RANGE] d_testOutPrimsFlatB bounds:" << std::endl;
    std::cout << "    * Smallest Primitive ID = " << hostMinVal << std::endl;
    std::cout << "    * Largest Primitive ID  = " << hostMaxVal << std::endl;
    
    if (hostMaxVal >= static_cast<uint32_t>(numPrimsB)) {
      std::cerr << "    [WARNING] Out-of-bounds primitive ID found! Max ID " << hostMaxVal 
                << " meets or exceeds raw asset total " << numPrimsB << std::endl;
    }
  } else {
    std::cout << " -> [VALUE RANGE] d_testOutPrimsFlatB is empty. No range to look up." << std::endl;
  }

  std::cout << "==================================================\n" << std::endl;

  // Memory Cleanup Phase
  if(outNodeBoxesA)
    _FREE(outNodeBoxesA, s, memResource);
  if(outSortedPrimIDsA)
    _FREE(outSortedPrimIDsA, s, memResource);
  if(outNodeOffsetsA)
    _FREE(outNodeOffsetsA, s, memResource);
  if(outNodeBoxesB)
    _FREE(outNodeBoxesB, s, memResource);
  if(outSortedPrimIDsB)
    _FREE(outSortedPrimIDsB, s, memResource);
  if(outNodeOffsetsB)
    _FREE(outNodeOffsetsB, s, memResource);
}
// --------------------------------------------------------------------
// MAIN ENTRY TESTING PIPELINE
// --------------------------------------------------------------------
extern "C" void kernelsTestBVH(const cuBQL::Triangle* hMeshA,
                               int numTrianglesA,
                               int maxCellSizeA,
                               const cuBQL::Triangle* hMeshB,
                               int numTrianglesB,
                               int maxCellSizeB,
                               int batchMultiplier,
                               int mode,
                               int leafThreshold,
                               ExecutionStats& stats,
                               std::vector<int2>& hGreenPairs,
                               std::vector<int2>& hYellowPairs) {
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
  cuBQL::gpuBuilder_v2_2::build_custom(bvhA, dBoxesA, numTrianglesA, buildConfig, (uint32_t)maxCellSizeA,
                                       thrust::raw_pointer_cast(d_markedNodeIndicesA.data()),
                                       thrust::raw_pointer_cast(d_nodeDescendantCountsA.data()), &h_outMarkedCountA,
                                       stream, memResource);
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
  cuBQL::gpuBuilder_v2_2::build_custom(bvhB, dBoxesB, numTrianglesB, buildConfig, (uint32_t)maxCellSizeB,
                                       thrust::raw_pointer_cast(d_markedNodeIndicesB.data()),
                                       thrust::raw_pointer_cast(d_nodeDescendantCountsB.data()), &h_outMarkedCountB,
                                       stream, memResource);
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

  testFastBVH<float, 3>(dBoxesA, numTrianglesA, maxCellSizeA, dBoxesB, numTrianglesB, maxCellSizeB, hostBoxA, hostBoxB,
                        stream, memResource, batchMultiplier, leafThreshold);

  // --------------------------------------------------------------------
  // EXTRACTED GPU PARALLEL CRISS-CROSS INTERSECTION MODULE
  // --------------------------------------------------------------------
  double tCrossStart = cuBQL::getCurrentTime();

  uint32_t totalIntersections =
      executeCrissCrossIntersection(bvhA, d_markedNodeIndicesA, h_outMarkedCountA, bvhB, d_markedNodeIndicesB,
                                    h_outMarkedCountB, d_outPairsA, d_outPairsB);

  cudaDeviceSynchronize();
  double tCrossEnd = cuBQL::getCurrentTime();

  uint64_t totalPossiblePairs = (uint64_t)h_outMarkedCountA * h_outMarkedCountB;
  double intersectionPercentage =
      totalPossiblePairs > 0 ? ((double)totalIntersections / totalPossiblePairs) * 100.0 : 0.0;

  // --------------------------------------------------------------------
  // DUAL-TREE TRAVERSAL STEP OVERHEAD PHASE
  // --------------------------------------------------------------------
  double tDualStepStart = cuBQL::getCurrentTime();

  if(mode > 0) {
    executeDualTreeStep(mode, maxCellSizeA, maxCellSizeB, d_outPairsA, d_outPairsB, d_markedNodeIndicesA,
                        d_markedNodeIndicesB, d_nodeDescendantCountsA, d_nodeDescendantCountsB, h_outMarkedCountA,
                        h_outMarkedCountB, bvhA, bvhB, dMeshA, dMeshB);
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
  executeRapidDescentBFS(bvhB, numTrianglesB, d_markedNodeIndicesB, d_nodeDescendantCountsB, h_outMarkedCountB,
                         d_outOffsetsB, d_outPrimsFlatB);

  if(h_outMarkedCountB > 0) {
    int blockSize = 256;
    int gridSize = (h_outMarkedCountB + blockSize - 1) / blockSize;
    populateReverseMapBKernel<<<gridSize, blockSize, 0, stream>>>(thrust::raw_pointer_cast(d_reverseMapB.data()),
                                                                  thrust::raw_pointer_cast(d_markedNodeIndicesB.data()),
                                                                  h_outMarkedCountB);
  }

  cudaDeviceSynchronize();
  double tGpuBfsEnd = cuBQL::getCurrentTime();

  finalCandidatePairs =
      executeBatchedCrossIntersectionLoop(batchMultiplier, totalBatches, d_outPairsA, d_outPairsB, d_reverseMapB,
                                          d_markedNodeIndicesB, d_outOffsetsB, d_outPrimsFlatB, d_nodeDescendantCountsB,
                                          h_outMarkedCountB, bvhA, dMeshA, dMeshB, hGreenPairs, hYellowPairs, tracker);

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

  stats.meshATotalNodes = bvhA.numNodes;
  stats.meshAExtractedTargets = h_outMarkedCountA;
  stats.meshBTotalNodes = bvhB.numNodes;
  stats.meshBExtractedTargets = h_outMarkedCountB;
  stats.totalIntersections = totalIntersections;
  stats.totalPossiblePairs = totalPossiblePairs;
  stats.intersectionPercentage = intersectionPercentage;
  stats.buildRefitMeshAMs = (tBuildAEnd - tBuildAStart) * 1000.0;
  stats.buildRefitMeshBMs = (tBuildBEnd - tBuildBStart) * 1000.0;
  stats.gpuCrossCheckEngineMs = (tCrossEnd - tCrossStart) * 1000.0;
  stats.parallelDfsDescentBMs = (tGpuBfsEnd - tGpuBfsStart) * 1000.0;
  stats.GPUTotalTime = (tPipelineEnd - tPipelineStart) * 1000.0;
  stats.totalCrissCrossBatches = totalBatches;
  stats.finalAabbCandidatePairs = finalCandidatePairs;
  stats.confirmedGreenPairs = hGreenPairs.size();
  stats.confirmedYellowPairs = hYellowPairs.size();
  stats.loopTracker = tracker;
}