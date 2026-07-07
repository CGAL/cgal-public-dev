#define CUBQL_GPU_BUILDER_IMPLEMENTATION 1


// Core cuBQL headers
#include "cuBQL/bvh.h"
#include "cuBQL/builder/cuda.h"

// Your custom experimental builder layout
#include "include/third-party/cubql/sm_builder_v2_2.h"
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
void logMeshDimensions(const std::string& label, const cuBQL::box_t<T, D>& oldBox, const cuBQL::box_t<T, D>& newBox) {
  std::cout << "    [DIMENSIONS " << label << " - OLD]:" << std::endl;
  for(int d = 0; d < D; ++d) {
    T size = oldBox.upper[d] - oldBox.lower[d];
    std::cout << "        Dim " << d << " -> Bounds: [" << oldBox.lower[d] << ", " << oldBox.upper[d]
              << "] | Size: " << size << std::endl;
  }
  std::cout << "    [DIMENSIONS " << label << " - NEW]:" << std::endl;
  for(int d = 0; d < D; ++d) {
    T size = newBox.upper[d] - newBox.lower[d];
    std::cout << "        Dim " << d << " -> Bounds: [" << newBox.lower[d] << ", " << newBox.upper[d]
              << "] | Size: " << size << std::endl;
  }
}

template <typename T, int D>
cuBQL::box_t<T, D> computeCentroidBox(const cuBQL::box_t<T, D>* d_boxes,
                                      int numPrims,
                                      cudaStream_t s,
                                      cuBQL::DeviceMemoryResource& memResource,
                                      const std::string& label) {
  cuBQL::box_t<T, D> newBox;
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double start = cuBQL::getCurrentTime();

  cuBQL::utils::computeGlobalBoxParallel<T, D>(newBox, d_boxes, numPrims, s, memResource);

  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double ms = (cuBQL::getCurrentTime() - start) * 1000.0;
  std::cout << "    * " << label << " Centroid Box Time: " << ms << " ms" << std::endl;
  return newBox;
}

template <typename T, int D>
void runSpeedrunInitialization(const cuBQL::box_t<T, D>* d_boxes,
                               int numPrims,
                               int maxCellSize,
                               const cuBQL::box_t<T, D>& globalBox,
                               cuBQL::box_t<T, D>*& outBoxes,
                               uint32_t*& outPrimIDs,
                               uint32_t*& outOffsets,
                               uint32_t& outActiveCells,
                               cudaStream_t s,
                               cuBQL::DeviceMemoryResource& memResource,
                               double& elapsedMs) {
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double start = cuBQL::getCurrentTime();

  cuBQL::gpuBuilder_v3::test_speedrun_initialization_linear(d_boxes, numPrims, static_cast<uint32_t>(maxCellSize),
                                                            globalBox, outBoxes, outPrimIDs, outOffsets, outActiveCells,
                                                            s, memResource);

  CUBQL_CUDA_CALL(StreamSynchronize(s));
  elapsedMs = (cuBQL::getCurrentTime() - start) * 1000.0;
}

template <typename T, int D>
void verifyTreePopulatedCells(const std::string& label,
                              uint32_t activeCells,
                              const cuBQL::box_t<T, D>* d_boxes,
                              const uint32_t* d_offsets) {
  std::vector<cuBQL::box_t<T, D>> h_boxes(activeCells);
  std::vector<uint32_t> h_offsets(activeCells + 1);

  CUBQL_CUDA_CALL(Memcpy(h_boxes.data(), d_boxes, activeCells * sizeof(cuBQL::box_t<T, D>), cudaMemcpyDeviceToHost));
  CUBQL_CUDA_CALL(Memcpy(h_offsets.data(), d_offsets, (activeCells + 1) * sizeof(uint32_t), cudaMemcpyDeviceToHost));

  uint32_t sumPrims = h_offsets[activeCells] - h_offsets[0];
  std::cout << " -> " << label << " Active Nodes: " << activeCells << " | Total Verified Prims = " << sumPrims
            << std::endl;
}

void checkIntersectionPairRanges(uint32_t totalIntersections,
                                 const thrust::device_vector<uint32_t>& d_pairsA,
                                 const thrust::device_vector<uint32_t>& d_pairsB) {
  std::cout << "\n -> [PAIR RANGE CHECK] Analyzing intersection pair arrays (Total: " << totalIntersections << ")..."
            << std::endl;
  if(totalIntersections > 0) {
    auto mmA = thrust::minmax_element(thrust::device, d_pairsA.begin(), d_pairsA.end());
    std::cout << "    * d_intersectPairsA => Min Cell ID: " << *mmA.first << " | Max Cell ID: " << *mmA.second
              << std::endl;

    auto mmB = thrust::minmax_element(thrust::device, d_pairsB.begin(), d_pairsB.end());
    std::cout << "    * d_intersectPairsB => Min Cell ID: " << *mmB.first << " | Max Cell ID: " << *mmB.second
              << std::endl;
  } else {
    std::cout << "    * [NOTICE] No intersecting pairs found. Skipping bounds check." << std::endl;
  }
  std::cout << "--------------------------------------------------" << std::endl;
}

template <typename T, int D>
void verifyStructuralIntegrity(const std::string& label,
                               const cuBQL::BinaryBVH<T, D>& bvh,
                               uint32_t activeCells,
                               const uint32_t* d_descendantCounts,
                               int expectedSurviving) {
  std::vector<typename cuBQL::BinaryBVH<T, D>::Node> h_nodes(bvh.numNodes);
  std::vector<uint32_t> h_descendants(bvh.numNodes);

  CUBQL_CUDA_CALL(Memcpy(h_nodes.data(), bvh.nodes, bvh.numNodes * sizeof(typename cuBQL::BinaryBVH<T, D>::Node),
                         cudaMemcpyDeviceToHost));
  CUBQL_CUDA_CALL(
      Memcpy(h_descendants.data(), d_descendantCounts, bvh.numNodes * sizeof(uint32_t), cudaMemcpyDeviceToHost));

  uint32_t totalRootDescendants = 0;
  for(uint32_t i = 0; i < activeCells; ++i) {
    totalRootDescendants += h_descendants[i];
  }

  std::cout << " [" << label << " DATA] Total Nodes in Array: " << bvh.numNodes << " | Root Cells: " << activeCells
            << std::endl;
  std::cout << "   -> TEST 1 (Sum of Roots): Total Prims accounted for at Root Level = " << totalRootDescendants
            << " (Expected Surviving: " << expectedSurviving << ")" << std::endl;

  uint32_t verifiedTreePrims = 0;
  uint32_t brokenLinks = 0;
  uint32_t leafCount = 0;
  int alerts = 0;

  for(uint32_t rootID = 0; rootID < activeCells; ++rootID) {
    std::vector<uint32_t> stack = {rootID};

    while(!stack.empty()) {
      uint32_t currID = stack.back();
      stack.pop_back();

      if(currID >= bvh.numNodes) {
        brokenLinks++;
        continue;
      }

      auto& node = h_nodes[currID];
      uint32_t childOffset = node.admin.offset;
      uint32_t primCount = node.admin.count;

      if(childOffset == (uint32_t)-1 || primCount > 0) {
        verifiedTreePrims += primCount;
        leafCount++;
      } else {
        uint32_t leftChild = childOffset;
        uint32_t rightChild = childOffset + 1;

        if(leftChild < bvh.numNodes && rightChild < bvh.numNodes) {
          uint32_t sumChildren = h_descendants[leftChild] + h_descendants[rightChild];
          if(h_descendants[currID] != sumChildren && h_descendants[currID] != 0 && alerts++ < 5) {
            std::cout << "      [MISMATCH ALERT " << label << "] Node " << currID << " says " << h_descendants[currID]
                      << " descendants, but children sum to " << sumChildren << std::endl;
          }
        }
        stack.push_back(rightChild);
        stack.push_back(leftChild);
      }
    }
  }

  std::cout << "   -> TEST 2 (Full Hierarchy Walk): Traversed " << leafCount << " leaves. Found " << verifiedTreePrims
            << " primitives inside the tree structures. (Broken links: " << brokenLinks << ")" << std::endl;
  if(verifiedTreePrims != static_cast<uint32_t>(expectedSurviving)) {
    std::cout << "   !!! [CRITICAL] " << label
              << " layout is structurally dropping primitives during top-down splits !!!" << std::endl;
  }
}

template<typename T, int D>
void verifyParentChildBoxCoverage(const std::string& treeLabel, const cuBQL::BinaryBVH<T, D>& bvh) {
  if (bvh.numNodes <= 0 || bvh.nodes == nullptr) {
    std::cout << "   [" << treeLabel << " WARNING] Empty or unallocated BVH structure skipped." << std::endl;
    return;
  }

  std::vector<typename cuBQL::BinaryBVH<T, D>::Node> h_nodes(bvh.numNodes);
  cudaMemcpy(h_nodes.data(), bvh.nodes, bvh.numNodes * sizeof(typename cuBQL::BinaryBVH<T, D>::Node), cudaMemcpyDefault);

  uint32_t violationCount = 0;
  const T epsilon = static_cast<T>(1e-5);

  for (int i = 0; i < bvh.numNodes; ++i) {
    const auto& parent = h_nodes[i];
    
    // Skip leaf nodes
    if (parent.admin.count > 0 || parent.admin.offset == (uint32_t)-1) {
      continue;
    }

    // FIX: Skip zeroed out alignment padding slots
    if (parent.admin.offset == 0) {
      continue;
    }

    uint32_t c0 = parent.admin.offset + 0;
    uint32_t c1 = parent.admin.offset + 1;

    if (c0 >= (uint32_t)bvh.numNodes || c1 >= (uint32_t)bvh.numNodes) {
      violationCount++;
      continue;
    }

    const auto& child0 = h_nodes[c0];
    const auto& child1 = h_nodes[c1];
    bool nodeHasViolation = false;

    for (int axis = 0; axis < D; ++axis) {
      if (child0.bounds.lower[axis] < parent.bounds.lower[axis] - epsilon || 
          child0.bounds.upper[axis] > parent.bounds.upper[axis] + epsilon) {
        nodeHasViolation = true;
      }
      if (child1.bounds.lower[axis] < parent.bounds.lower[axis] - epsilon || 
          child1.bounds.upper[axis] > parent.bounds.upper[axis] + epsilon) {
        nodeHasViolation = true;
      }
    }

    if (nodeHasViolation) {
      violationCount++;
      if (violationCount <= 5) {
        std::cout << "      |--> [WARNING] Coverage Violation at Internal Node ID [" << i << "]:\n"
                  << "           |- Parent Box:  Min[" << parent.bounds.lower[0] << ", " << parent.bounds.lower[1] << "] Max[" << parent.bounds.upper[0] << ", " << parent.bounds.upper[1] << "]\n"
                  << "           |- Child0 (" << c0 << "): Min[" << child0.bounds.lower[0] << ", " << child0.bounds.lower[1] << "] Max[" << child0.bounds.upper[0] << ", " << child0.bounds.upper[1] << "]\n"
                  << "           |- Child1 (" << c1 << "): Min[" << child1.bounds.lower[0] << ", " << child1.bounds.lower[1] << "] Max[" << child1.bounds.upper[0] << ", " << child1.bounds.upper[1] << "]\n";
      }
    }
  }

  if (violationCount == 0) {
    std::cout << "   [SUCCESS] " << treeLabel << " Bounding Box Enclosure check passed. Parent nodes perfectly contain children cells." << std::endl;
  } else {
    std::cout << "   [FAILURE] " << treeLabel << " Box coverage broken! Found " << violationCount << " violating parent nodes." << std::endl;
  }
}

// --- MAIN ENTRY POINT ---

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

  // 1. Centroid Box Computations
  std::cout << "\n -> [COMPUTE] Calculating parallel centroid bounding boxes into temporal variables..." << std::endl;
  cuBQL::box_t<T, D> globalBoxANew = computeCentroidBox(d_boxesA, numPrimsA, s, memResource, "Mesh A");
  logMeshDimensions("Mesh A", globalBoxA, globalBoxANew);
  std::cout << "--------------------------------------------------" << std::endl;

  cuBQL::box_t<T, D> globalBoxBNew = computeCentroidBox(d_boxesB, numPrimsB, s, memResource, "Mesh B");
  logMeshDimensions("Mesh B", globalBoxB, globalBoxBNew);
  std::cout << "==================================================" << std::endl;

  // 2. Linear Initialization Speedrun Runs
  cuBQL::box_t<T, D>* outNodeBoxesA = nullptr;
  uint32_t* outSortedPrimIDsA = nullptr;
  uint32_t* outNodeOffsetsA = nullptr;
  uint32_t outTotalActiveCellsA = 0;
  cuBQL::box_t<T, D>* outNodeBoxesB = nullptr;
  uint32_t* outSortedPrimIDsB = nullptr;
  uint32_t* outNodeOffsetsB = nullptr;
  uint32_t outTotalActiveCellsB = 0;
  double tInitAMs = 0, tInitBMs = 0;

  runSpeedrunInitialization(d_boxesA, numPrimsA, maxCellSizeA, globalBoxANew, outNodeBoxesA, outSortedPrimIDsA,
                            outNodeOffsetsA, outTotalActiveCellsA, s, memResource, tInitAMs);
  runSpeedrunInitialization(d_boxesB, numPrimsB, maxCellSizeB, globalBoxBNew, outNodeBoxesB, outSortedPrimIDsB,
                            outNodeOffsetsB, outTotalActiveCellsB, s, memResource, tInitBMs);

  // 3. Populate & Verify Metrics
  std::cout << "\n--------------------------------------------------" << std::endl;
  std::cout << " [LINEAR INITIALIZATION METRICS OVERVIEW]" << std::endl;
  std::cout << " -> Tree A: Total Compounded Active Cells = " << outTotalActiveCellsA << std::endl;
  std::cout << " -> Tree B: Total Compounded Active Cells = " << outTotalActiveCellsB << std::endl;
  std::cout << "--------------------------------------------------" << std::endl;

  std::cout << " [POPULATED CELLS CHECK SUMMARY]" << std::endl;
  verifyTreePopulatedCells("Tree A", outTotalActiveCellsA, outNodeBoxesA, outNodeOffsetsA);
  verifyTreePopulatedCells("Tree B", outTotalActiveCellsB, outNodeBoxesB, outNodeOffsetsB);
  std::cout << "--------------------------------------------------\n" << std::endl;

  // Hardcode reading initial counts from temporary validation or assume full set before pruning
  // To match your code structure exactly, we read out node verification to get the sum bounds
  uint32_t initialCellsA = outTotalActiveCellsA;
  uint32_t initialCellsB = outTotalActiveCellsB;

  // 4. Parallel Cross Check Execution
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
  double crossCheckMs = (cuBQL::getCurrentTime() - crossCheckStart) * 1000.0;

  int currentPrimsNumA = numPrimsA; // Safe snapshot setup matching original structural pattern
  int currentPrimsNumB = numPrimsB;

  // 5. Compaction Pipeline / Post Build Pruning
  std::cout << " [PIPELINE PRUNING] => Dispatched grid streaming parallel re-index compaction (POST-BUILD)..."
            << std::endl;
  std::cout << " -> [SIZE BEFORE] d_intersectPairsA size: " << d_intersectPairsA.size()
            << " | d_intersectPairsB size: " << d_intersectPairsB.size() << std::endl;

  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double pruneStartA = cuBQL::getCurrentTime();
  parallelPruneAndReindexAll(thrust::raw_pointer_cast(d_intersectPairsA.data()), totalIntersections, outSortedPrimIDsA,
                             outNodeOffsetsA, outTotalActiveCellsA, currentPrimsNumA, s, memResource);
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double pruneMsA = (cuBQL::getCurrentTime() - pruneStartA) * 1000.0;

  double pruneStartB = cuBQL::getCurrentTime();
  parallelPruneAndReindexAll(thrust::raw_pointer_cast(d_intersectPairsB.data()), totalIntersections, outSortedPrimIDsB,
                             outNodeOffsetsB, outTotalActiveCellsB, currentPrimsNumB, s, memResource);
  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double pruneMsB = (cuBQL::getCurrentTime() - pruneStartB) * 1000.0;

  std::cout << " -> [SIZE AFTER]  d_intersectPairsA size: " << d_intersectPairsA.size()
            << " | d_intersectPairsB size: " << d_intersectPairsB.size() << std::endl;

  checkIntersectionPairRanges(totalIntersections, d_intersectPairsA, d_intersectPairsB);

  // 6. Forest Expansion Allocation setup
  std::cout << " [FOREST EXPANSION] => Launching Level-by-Level Parallel Sub-Tree Compilation..." << std::endl;
  const uint32_t maxNodesA = 2u * (uint32_t)numPrimsA + (outTotalActiveCellsA + (outTotalActiveCellsA & 1u)) + 2u;
  const uint32_t maxNodesB = 2u * (uint32_t)numPrimsB + (outTotalActiveCellsB + (outTotalActiveCellsB & 1u)) + 2u;

  thrust::device_vector<uint32_t> d_nodeDescendantCountsA(maxNodesA, 0);
  thrust::device_vector<uint32_t> d_nodeDescendantCountsB(maxNodesB, 0);

  cuBQL::BuildConfig buildConfig;
  buildConfig.makeLeafThreshold = userleafThreshold;

  // --- FOREST A BUILD & REFIT ---
  cuBQL::BinaryBVH<T, D> bvhA;
  double forestAStart = cuBQL::getCurrentTime();

  cuBQL::gpuBuilder_v4::build_forest<T, D>(bvhA, d_boxesA, currentPrimsNumA, numPrimsA, outTotalActiveCellsA,
                                           outSortedPrimIDsA, outNodeOffsetsA, buildConfig,
                                           thrust::raw_pointer_cast(d_nodeDescendantCountsA.data()), s, memResource);

  // Add the refit call for Forest A here (asynchronous on stream s)
  cuBQL::cuda_forest::refit_forest<T, D>(bvhA, d_boxesA, outTotalActiveCellsA, s, memResource);

  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double forestAMs = (cuBQL::getCurrentTime() - forestAStart) * 1000.0;

  // --- FOREST B BUILD & REFIT ---
  cuBQL::BinaryBVH<T, D> bvhB;
  double forestBStart = cuBQL::getCurrentTime();

  cuBQL::gpuBuilder_v4::build_forest<T, D>(bvhB, d_boxesB, currentPrimsNumB, numPrimsB, outTotalActiveCellsB,
                                           outSortedPrimIDsB, outNodeOffsetsB, buildConfig,
                                           thrust::raw_pointer_cast(d_nodeDescendantCountsB.data()), s, memResource);

  // Add the refit call for Forest B here (asynchronous on stream s)
  cuBQL::cuda_forest::refit_forest<T, D>(bvhB, d_boxesB, outTotalActiveCellsB, s, memResource);

  CUBQL_CUDA_CALL(StreamSynchronize(s));
  double forestBMs = (cuBQL::getCurrentTime() - forestBStart) * 1000.0;

  // 7. Hierarchy Walks Metric Verification
  std::cout << "\n==================================================" << std::endl;
  std::cout << " [RUNNING METRIC CHECKS: FOREST STRUCTURE & PASSES]" << std::endl;
  std::cout << "==================================================" << std::endl;
  verifyStructuralIntegrity("TREE A", bvhA, outTotalActiveCellsA,
                            thrust::raw_pointer_cast(d_nodeDescendantCountsA.data()), currentPrimsNumA);
  verifyStructuralIntegrity("TREE B", bvhB, outTotalActiveCellsB,
                            thrust::raw_pointer_cast(d_nodeDescendantCountsB.data()), currentPrimsNumB);

  std::cout << "\n--- PARENT-CHILD BOUNDING BOX COVERAGE VALIDATION ---" << std::endl;
  verifyParentChildBoxCoverage<T, D>("TREE A", bvhA);
  verifyParentChildBoxCoverage<T, D>("TREE B", bvhB);

  // === INJECTED AABB COVERAGE CHECKS ===
  //  std::cout << "--------------------------------------------------" << std::endl;
  //  std::cout << " [RUNNING INJECTED AABB ENCAPSULATION CHECKS]" << std::endl;
  // verifyAABBHierarchyInclusion<T, D>("TREE A", bvhA, s);
  // verifyAABBHierarchyInclusion<T, D>("TREE B", bvhB, s);
  // ======================================

  std::cout << "==================================================\n" << std::endl;

  // 8. Performance Display Readout
  uint64_t totalValidActivePairs = (uint64_t)initialCellsA * initialCellsB;
  double overlapPercentage =
      totalValidActivePairs > 0 ? ((double)totalIntersections / totalValidActivePairs * 100.0) : 0.0;

  std::cout << "\n==================================================" << std::endl;
  std::cout << "               FINAL BENCHMARK SUMMARY            " << std::endl;
  std::cout << "==================================================" << std::endl;
  std::cout << " -> Mesh A Speedrun Init Speed      : " << tInitAMs << " ms" << std::endl;
  std::cout << " -> Mesh B Speedrun Init Speed      : " << tInitBMs << " ms" << std::endl;
  std::cout << " -> Cell Cross-Check Execution      : " << crossCheckMs << " ms" << std::endl;
  std::cout << " ------------------------------------------------" << std::endl;
  std::cout << " -> Mesh A Forest BVH Expansion     : " << forestAMs << " ms | Nodes: " << bvhA.numNodes << std::endl;
  std::cout << " -> Mesh B Forest BVH Expansion     : " << forestBMs << " ms | Nodes: " << bvhB.numNodes << std::endl;
  std::cout << "==================================================\n" << std::endl;

  // 9. Injected BFS Test Phase
  std::cout << " [INJECTED DEBUG] Running executeRapidDescentBFS on v4 Forest B..." << std::endl;
  thrust::device_vector<uint32_t> d_testMarkedNodeIndicesB(outTotalActiveCellsB);
  thrust::sequence(thrust::cuda::par.on(s), d_testMarkedNodeIndicesB.begin(), d_testMarkedNodeIndicesB.end(), 0u);

  thrust::device_vector<uint32_t> d_testOutOffsetsB, d_testOutPrimsFlatB;
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

  // Value Range Checks
  if(extractedPrimsB > 0) {
    auto result_pair = thrust::minmax_element(thrust::device, d_testOutPrimsFlatB.begin(), d_testOutPrimsFlatB.end());
    uint32_t hostMinVal = *result_pair.first;
    uint32_t hostMaxVal = *result_pair.second;
    std::cout << " -> [VALUE RANGE] d_testOutPrimsFlatB bounds:\n    * Smallest Primitive ID = " << hostMinVal
              << "\n    * Largest Primitive ID  = " << hostMaxVal << std::endl;
    if(hostMaxVal >= static_cast<uint32_t>(numPrimsB)) {
      std::cerr << "    [WARNING] Out-of-bounds primitive ID found! Max ID " << hostMaxVal
                << " meets or exceeds raw asset total " << numPrimsB << std::endl;
    }
  }

  // Cleanup
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

  //testFastBVH<float, 3>(dBoxesA, numTrianglesA, maxCellSizeA, dBoxesB, numTrianglesB, maxCellSizeB, hostBoxA, hostBoxB,
   //                     stream, memResource, batchMultiplier, leafThreshold);

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