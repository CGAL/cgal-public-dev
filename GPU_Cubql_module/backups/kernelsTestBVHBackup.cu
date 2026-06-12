#define CUBQL_GPU_BUILDER_IMPLEMENTATION 1

// Core cuBQL headers
#include "cuBQL/bvh.h"
#include "cuBQL/builder/cuda.h"

// Your custom experimental builder layout
#include "include/third-party/cubql/sm_builder_v2.h"

#include <thrust/device_vector.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include "samples/common/loadOBJ.h"

// --------------------------------------------------------------------
// KERNELS & HELPERS
// --------------------------------------------------------------------
__global__ void generateBoxes(cuBQL::box3f* boxes, const cuBQL::Triangle* tris, int N) {
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if(i < N) {
    boxes[i] = tris[i].bounds();
  }
}

typedef typename cuBQL::BinaryBVH<float, 3>::Node BvhNodeT;

// Pair mapping tracker for the dynamic GPU BFS frontier
struct FrontierElement
{
  uint32_t nodeID;     // Current node we are traversing in the global hierarchy
  uint32_t targetSlot; // The tracked root index [0 ... h_outMarkedCount-1] this branch belongs to
};

// --------------------------------------------------------------------
// GPU KERNELS FOR BREADTH-FIRST LAYER EXPANSION & PRIMITIVE HARVESTING
// --------------------------------------------------------------------
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

// --------------------------------------------------------------------
// MAIN ENTRY TESTING PIPELINE
// --------------------------------------------------------------------
extern "C" void kernelsTestBVH(const cuBQL::Triangle* hMesh, int numTriangles, int leafSize) {
  std::cout << "\n==================================================\n";
  std::cout << " RUNNING CLEANED GEOMETRIC VERIFICATION PIPELINE\n";
  std::cout << "==================================================\n";
  std::cout << "Input Mesh Size: " << numTriangles << " triangles.\n";

  if(numTriangles <= 0) {
    std::cerr << "Error: Input mesh contains no triangles.\n";
    return;
  }

  double tSetupStart = cuBQL::getCurrentTime();

  // 1. Allocate and upload mesh data
  cuBQL::Triangle* dMesh;
  CUBQL_CUDA_CALL(Malloc(&dMesh, numTriangles * sizeof(cuBQL::Triangle)));
  CUBQL_CUDA_CALL(Memcpy(dMesh, hMesh, numTriangles * sizeof(cuBQL::Triangle), cudaMemcpyDefault));

  double tSetupEnd = cuBQL::getCurrentTime();
  double totalSetupTimeMs = (tSetupEnd - tSetupStart) * 1000.0;

  double tBuildStart = cuBQL::getCurrentTime();

  cuBQL::box3f* dBoxes;
  CUBQL_CUDA_CALL(Malloc(&dBoxes, numTriangles * sizeof(cuBQL::box3f)));

  generateBoxes<<<cuBQL::divRoundUp(numTriangles, 256), 256>>>(dBoxes, dMesh, numTriangles);
  cudaDeviceSynchronize();

  // 4. Configure build parameters
  cuBQL::BuildConfig buildConfig;
  buildConfig.makeLeafThreshold = 1;

  // 5. Allocate tracking buffers
  uint32_t maxPossibleNodes = 2 * numTriangles;
  thrust::device_vector<uint32_t> d_markedNodeIndices(maxPossibleNodes, 0);
  thrust::device_vector<uint32_t> d_markedNodeDescendantCounts(maxPossibleNodes, 0);
  uint32_t h_outMarkedCount = 0;

  cuBQL::bvh3f bvh;
  cudaStream_t stream = 0;
  cuBQL::DeviceMemoryResource memResource;

  // 6. Invoke custom GPU builder
  uint32_t thresholdX = (uint32_t)leafSize;

  cuBQL::gpuBuilder_v2::build_custom(
      bvh, dBoxes, numTriangles, buildConfig, thresholdX, thrust::raw_pointer_cast(d_markedNodeIndices.data()),
      thrust::raw_pointer_cast(d_markedNodeDescendantCounts.data()), &h_outMarkedCount, stream, memResource);
  cudaDeviceSynchronize();
  double tBuildEnd = cuBQL::getCurrentTime();

  // 7. Perform standard post-hierarchy refit calculation loop
  double tRefitStart = cuBQL::getCurrentTime();
  cuBQL::cuda::refit(bvh, dBoxes, stream);
  cudaDeviceSynchronize();
  double tRefitEnd = cuBQL::getCurrentTime();

  // --------------------------------------------------------------------
  // COMPUTE SCENE VOLUME FROM REFITTED ROOT NODE (INDEX 0)
  // --------------------------------------------------------------------
  BvhNodeT h_rootNode;
  CUBQL_CUDA_CALL(Memcpy(&h_rootNode, &bvh.nodes[0], sizeof(BvhNodeT), cudaMemcpyDeviceToHost));

  float sceneDx = h_rootNode.bounds.upper[0] - h_rootNode.bounds.lower[0];
  float sceneDy = h_rootNode.bounds.upper[1] - h_rootNode.bounds.lower[1];
  float sceneDz = h_rootNode.bounds.upper[2] - h_rootNode.bounds.lower[2];
  float globalSceneVolume = sceneDx * sceneDy * sceneDz;
  if(globalSceneVolume <= 0.0f)
    globalSceneVolume = 1e-6f;

  double totalPureBvhTimeMs = (tBuildEnd - tBuildStart) * 1000.0;
  double totalRefitTimeMs = (tRefitEnd - tRefitStart) * 1000.0;

  // --------------------------------------------------------------------
  // GPU BREADTH-FIRST LAYER DESCEND PERFORMANCE RUN
  // --------------------------------------------------------------------
  double tGpuBfsStart = cuBQL::getCurrentTime();

  thrust::device_vector<cuBQL::box3f> d_bfsCollectedBoxes(h_outMarkedCount);
  thrust::device_vector<uint32_t> d_bfsCollectedCounts(h_outMarkedCount);

  if(h_outMarkedCount > 0) {
    initTargetBoxes<<<cuBQL::divRoundUp(h_outMarkedCount, (uint32_t)64), 64>>>(
        thrust::raw_pointer_cast(d_bfsCollectedBoxes.data()), thrust::raw_pointer_cast(d_bfsCollectedCounts.data()),
        h_outMarkedCount);

    uint32_t frontierAllocationSize = std::max((uint32_t)numTriangles, (uint32_t)32768);
    thrust::device_vector<FrontierElement> d_frontierA(frontierAllocationSize);
    thrust::device_vector<FrontierElement> d_frontierB(frontierAllocationSize);
    thrust::device_vector<uint32_t> d_atomicCounterPool(1, 0);

    std::vector<FrontierElement> h_seedFrontier(h_outMarkedCount);
    std::vector<uint32_t> h_tempIndices(h_outMarkedCount);
    thrust::copy(d_markedNodeIndices.begin(), d_markedNodeIndices.begin() + h_outMarkedCount, h_tempIndices.begin());

    for(uint32_t i = 0; i < h_outMarkedCount; ++i) {
      h_seedFrontier[i] = {h_tempIndices[i], (uint32_t)i};
    }
    thrust::copy(h_seedFrontier.begin(), h_seedFrontier.end(), d_frontierA.begin());

    uint32_t currentFrontierSize = h_outMarkedCount;
    FrontierElement* p_inBuffer = thrust::raw_pointer_cast(d_frontierA.data());
    FrontierElement* p_outBuffer = thrust::raw_pointer_cast(d_frontierB.data());

    while(currentFrontierSize > 0) {
      d_atomicCounterPool[0] = 0;
      uint32_t threadsPerBlock = 256;
      uint32_t blockCount = cuBQL::divRoundUp(currentFrontierSize, (uint32_t)threadsPerBlock);

      processFrontierLayer<<<blockCount, threadsPerBlock>>>(
          p_inBuffer, currentFrontierSize, p_outBuffer, thrust::raw_pointer_cast(d_atomicCounterPool.data()), bvh.nodes,
          bvh.numNodes, bvh.primIDs, bvh.numPrims, thrust::raw_pointer_cast(d_bfsCollectedBoxes.data()),
          thrust::raw_pointer_cast(d_bfsCollectedCounts.data()), dBoxes);
      cudaDeviceSynchronize();

      currentFrontierSize = d_atomicCounterPool[0];
      std::swap(p_inBuffer, p_outBuffer);
    }
  }
  double tGpuBfsEnd = cuBQL::getCurrentTime();
  double totalGpuBfsTimeMs = (tGpuBfsEnd - tGpuBfsStart) * 1000.0;

  // --- DEVICE-TO-HOST COPY TIMING ---
  double tTransferStart = cuBQL::getCurrentTime();
  std::vector<uint32_t> h_indices(h_outMarkedCount);
  std::vector<uint32_t> h_counts(h_outMarkedCount);
  std::vector<BvhNodeT> h_finalNodes(bvh.numNodes);

  std::vector<cuBQL::box3f> h_gpuBfsBoxes(h_outMarkedCount);
  std::vector<uint32_t> h_gpuBfsCounts(h_outMarkedCount);

  if(h_outMarkedCount > 0) {
    thrust::copy(d_markedNodeIndices.begin(), d_markedNodeIndices.begin() + h_outMarkedCount, h_indices.begin());
    thrust::copy(d_markedNodeDescendantCounts.begin(), d_markedNodeDescendantCounts.begin() + h_outMarkedCount,
                 h_counts.begin());
    thrust::copy(d_bfsCollectedBoxes.begin(), d_bfsCollectedBoxes.end(), h_gpuBfsBoxes.begin());
    thrust::copy(d_bfsCollectedCounts.begin(), d_bfsCollectedCounts.end(), h_gpuBfsCounts.begin());
    CUBQL_CUDA_CALL(Memcpy(h_finalNodes.data(), bvh.nodes, bvh.numNodes * sizeof(BvhNodeT), cudaMemcpyDeviceToHost));
  }

  cudaDeviceSynchronize();
  double tTransferEnd = cuBQL::getCurrentTime();
  double totalTransferTimeMs = (tTransferEnd - tTransferStart) * 1000.0;

  // --------------------------------------------------------------------
  // PERFORMANCE METRICS PRINT-OUT
  // --------------------------------------------------------------------
  std::cout << "--------------------------------------------------\n";
  std::cout << " FULL BREAKDOWN PERFORMANCE METRICS\n";
  std::cout << "--------------------------------------------------\n";
  std::cout << "  |- 1. Input Data & AABB Setup Cost: " << totalSetupTimeMs << " ms\n";
  std::cout << "  |- 2. Pure BVH Hierarchy Generation: " << totalPureBvhTimeMs << " ms\n";
  std::cout << "  |- 3. Bounding Box Bottom-Up Refit: " << totalRefitTimeMs << " ms\n";
  std::cout << "  |- 4. GPU Parallel BFS Level Descend: " << totalGpuBfsTimeMs << " ms\n";
  std::cout << "  |- 5. Device-To-Host Copy Overhead : " << totalTransferTimeMs << " ms\n";
  std::cout << "  |= TOTAL PIPELINE EXECUTION TIME   : "
            << (totalSetupTimeMs + totalPureBvhTimeMs + totalRefitTimeMs + totalGpuBfsTimeMs + totalTransferTimeMs)
            << " ms\n";
  std::cout << "--------------------------------------------------\n\n";

  std::cout << "Total nodes created in BVH array: " << bvh.numNodes << "\n";
  std::cout << "Total distinct nodes meeting tracking threshold (< " << thresholdX << "): " << h_outMarkedCount
            << "\n\n";

  // 8. Print Cleaned Results with GPU Parallel BFS Verification
  if(h_outMarkedCount > 0) {
    std::cout << "--- EXTRACTED TARGET THRESHOLD NODES (GEOMETRIC EVALUATION) ---\n";
    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout.precision(5);

    for(size_t i = 0; i < h_outMarkedCount; ++i) {
      uint32_t nodeID = h_indices[i];
      uint32_t expectedPrimsCount = h_counts[i];

      // Get node bounding box metadata volume percentage
      auto& node = h_finalNodes[nodeID];
      float dx = node.bounds.upper[0] - node.bounds.lower[0];
      float dy = node.bounds.upper[1] - node.bounds.lower[1];
      float dz = node.bounds.upper[2] - node.bounds.lower[2];
      float nodeVolume = (dx > 0.0f && dy > 0.0f && dz > 0.0f) ? (dx * dy * dz) : 0.0f;
      float metadataVolumePercentage = (nodeVolume / globalSceneVolume) * 100.0f;

      // Process the boundaries built by our Level-by-Layer Parallel GPU BFS Pass
      auto& bfsBox = h_gpuBfsBoxes[i];
      float bdx = bfsBox.upper[0] - bfsBox.lower[0];
      float bdy = bfsBox.upper[1] - bfsBox.lower[1];
      float bdz = bfsBox.upper[2] - bfsBox.lower[2];
      float bfsVolume = (bdx > 0.0f && bdy > 0.0f && bdz > 0.0f) ? (bdx * bdy * bdz) : 0.0f;
      float gpuBfsVolumePercentage = (bfsVolume / globalSceneVolume) * 100.0f;

      std::cout << "  Item [" << i << "] -> Node ID: " << nodeID << "\n";
      std::cout << "       |- Primitives Validation -> Expected: " << expectedPrimsCount
                << " | GPU BFS Harvest: " << h_gpuBfsCounts[i] << "\n";
      std::cout << "       |- Node Metadata Volume Ratio:       " << metadataVolumePercentage << "%\n";
      std::cout << "       |- GPU Parallel BFS Descent Ratio:   " << gpuBfsVolumePercentage << "%\n\n";
    }
    std::cout << "--------------------------------------------------------------------\n";
  }

  // Clean resources
  CUBQL_CUDA_CALL(Free(dMesh));
  CUBQL_CUDA_CALL(Free(dBoxes));
  cuBQL::cuda::free(bvh, stream, memResource);
}