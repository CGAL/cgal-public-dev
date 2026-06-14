#define CUBQL_GPU_BUILDER_IMPLEMENTATION 1

#include "cuBQL/bvh.h"
#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/copy.h>
#include <thrust/execution_policy.h>
#include <vector>
#include <algorithm>

typedef typename cuBQL::BinaryBVH<float, 3>::Node BvhNodeT;

struct FrontierElement
{
  uint32_t nodeID;     
  uint32_t writeCursor; 
};

// --------------------------------------------------------------------
// KERNELS
// --------------------------------------------------------------------

__global__ void populateMarkedNodeCapacities(const uint32_t* d_markedNodeIndices,
                                             const uint32_t* d_nodeDescendantCounts,
                                             uint32_t* d_extractedCounts,
                                             uint32_t totalMarked) 
{
  uint32_t idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx < totalMarked) {
    uint32_t nodeID = d_markedNodeIndices[idx];
    d_extractedCounts[idx] = d_nodeDescendantCounts[nodeID];
  }
}

__global__ void processFrontierFlatLayer(const FrontierElement* d_inFrontier,
                                         uint32_t inCount,
                                         FrontierElement* d_outFrontier,
                                         uint32_t* d_outCount,
                                         uint32_t maxQueueCapacity,
                                         const BvhNodeT* d_nodes,
                                         uint32_t numNodes,
                                         const uint32_t* d_primIDs,
                                         uint32_t* d_outPrimsFlat,
                                         uint32_t maxPrimsFlatCapacity, 
                                         const uint32_t* d_nodeDescendantCounts)
{
  uint32_t idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= inCount) return;

  FrontierElement curr = d_inFrontier[idx];
  if (curr.nodeID >= numNodes) return;

  BvhNodeT node = d_nodes[curr.nodeID];

  // Leaf Node Case
  if (node.admin.count > 0) {
    uint32_t baseOffset = node.admin.offset;
    for (uint32_t p = 0; p < node.admin.count; ++p) {
      uint32_t writeIdx = curr.writeCursor + p;
      if (writeIdx < maxPrimsFlatCapacity) {
        d_outPrimsFlat[writeIdx] = d_primIDs[baseOffset + p];
      }
    }
  } 
  // Inner Node Case
  else {
    uint32_t leftChild = node.admin.offset;
    if (leftChild < numNodes) {
      uint32_t rightChild = leftChild + 1;
      uint32_t leftSize = d_nodeDescendantCounts[leftChild];

      uint32_t slot1 = atomicAdd(d_outCount, 1);
      if (slot1 < maxQueueCapacity) {
        d_outFrontier[slot1] = {leftChild, curr.writeCursor};
      }

      uint32_t slot2 = atomicAdd(d_outCount, 1);
      if (slot2 < maxQueueCapacity) {
        d_outFrontier[slot2] = {rightChild, curr.writeCursor + leftSize};
      }
    }
  }
}

// --------------------------------------------------------------------
// EXPOSED ORCHESTRATION PIPELINE
// --------------------------------------------------------------------
extern "C" void executeRapidDescentBFS(
    const cuBQL::bvh3f& bvhA,
    int numTrianglesA,
    const thrust::device_vector<uint32_t>& d_markedNodeIndicesA,
    const thrust::device_vector<uint32_t>& d_nodeDescendantCountsA, 
    uint32_t h_outMarkedCountA,
    thrust::device_vector<uint32_t>& d_outOffsets,                 
    thrust::device_vector<uint32_t>& d_outPrimsFlat)               
{
  if (h_outMarkedCountA == 0) return;

  // STEP 1: PRE-CALCULATE LAYOUT OFFSETS
  thrust::device_vector<uint32_t> d_tempCapacities(h_outMarkedCountA);
  
  uint32_t threadsInit = 256;
  uint32_t blocksInit = (h_outMarkedCountA + threadsInit - 1) / threadsInit;
  
  populateMarkedNodeCapacities<<<blocksInit, threadsInit>>>(
      thrust::raw_pointer_cast(d_markedNodeIndicesA.data()),
      thrust::raw_pointer_cast(d_nodeDescendantCountsA.data()),
      thrust::raw_pointer_cast(d_tempCapacities.data()),
      h_outMarkedCountA);
  
  d_outOffsets.resize(h_outMarkedCountA + 1);
  thrust::exclusive_scan(d_tempCapacities.begin(), d_tempCapacities.end(), d_outOffsets.begin());

  uint32_t totalRequiredPrims = thrust::reduce(d_tempCapacities.begin(), d_tempCapacities.end(), (uint32_t)0, thrust::plus<uint32_t>());
  d_outOffsets[h_outMarkedCountA] = totalRequiredPrims;

  d_outPrimsFlat.resize(totalRequiredPrims);

  // STEP 2: ALLOCATE COMPREHENSIVE FRONTIER COVERS
  uint32_t frontierAllocationSize = std::max((uint32_t)bvhA.numNodes, (uint32_t)32768);
  
  thrust::device_vector<FrontierElement> d_frontierA(frontierAllocationSize);
  thrust::device_vector<FrontierElement> d_frontierB(frontierAllocationSize);
  thrust::device_vector<uint32_t> d_atomicCounterPool(1, 0);
  
  std::vector<FrontierElement> h_seedFrontier(h_outMarkedCountA);
  std::vector<uint32_t> h_tempIndicesA(h_outMarkedCountA);
  std::vector<uint32_t> h_tempOffsetsA(h_outMarkedCountA);

  thrust::copy(d_markedNodeIndicesA.begin(), d_markedNodeIndicesA.begin() + h_outMarkedCountA, h_tempIndicesA.begin());
  thrust::copy(d_outOffsets.begin(), d_outOffsets.begin() + h_outMarkedCountA, h_tempOffsetsA.begin());

  for (uint32_t i = 0; i < h_outMarkedCountA; ++i) {
    h_seedFrontier[i] = {h_tempIndicesA[i], h_tempOffsetsA[i]};
  }
  
  thrust::copy(h_seedFrontier.begin(), h_seedFrontier.end(), d_frontierA.begin());

  uint32_t currentFrontierSize = h_outMarkedCountA;
  FrontierElement* p_inBuffer = thrust::raw_pointer_cast(d_frontierA.data());
  FrontierElement* p_outBuffer = thrust::raw_pointer_cast(d_frontierB.data());

  // STEP 3: LEVEL-BY-LEVEL BFS DESCENT LOOP
  while (currentFrontierSize > 0) {
    uint32_t zeroReset = 0;
    cudaMemcpy(thrust::raw_pointer_cast(d_atomicCounterPool.data()), &zeroReset, sizeof(uint32_t), cudaMemcpyHostToDevice);

    uint32_t threadsPerBlock = 256;
    uint32_t blockCount = (currentFrontierSize + threadsPerBlock - 1) / threadsPerBlock;

    processFrontierFlatLayer<<<blockCount, threadsPerBlock>>>(
        p_inBuffer, 
        currentFrontierSize, 
        p_outBuffer, 
        thrust::raw_pointer_cast(d_atomicCounterPool.data()), 
        frontierAllocationSize, 
        bvhA.nodes,
        bvhA.numNodes, 
        bvhA.primIDs, 
        thrust::raw_pointer_cast(d_outPrimsFlat.data()),
        totalRequiredPrims, 
        thrust::raw_pointer_cast(d_nodeDescendantCountsA.data()));
    
    cudaMemcpy(&currentFrontierSize, thrust::raw_pointer_cast(d_atomicCounterPool.data()), sizeof(uint32_t), cudaMemcpyDeviceToHost);

    if (currentFrontierSize > frontierAllocationSize) {
      currentFrontierSize = frontierAllocationSize;
    }

    std::swap(p_inBuffer, p_outBuffer);
  }
}