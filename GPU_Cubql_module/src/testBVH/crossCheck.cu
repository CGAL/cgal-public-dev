#include "crossCheck.h"
#include <thrust/scan.h>

// Helper device check for bounding box overlaps
__device__ inline bool boxBoxesIntersect(const cuBQL::box3f& a, const cuBQL::box3f& b) {
  return (a.lower[0] <= b.upper[0] && a.upper[0] >= b.lower[0]) &&
         (a.lower[1] <= b.upper[1] && a.upper[1] >= b.lower[1]) &&
         (a.lower[2] <= b.upper[2] && a.upper[2] >= b.lower[2]);
}

// Pass 1: Count intersections mapped per thread context
__global__ void countTargetNodesIntersections(
    const BvhNodeT* d_nodesA, const uint32_t* d_indicesA, uint32_t countA,
    const BvhNodeT* d_nodesB, const uint32_t* d_indicesB, uint32_t countB,
    uint32_t* d_perThreadCounts) 
{
  uint32_t idxA = threadIdx.x + blockIdx.x * blockDim.x;
  if (idxA >= countA) return;

  uint32_t nodeIDA = d_indicesA[idxA];
  cuBQL::box3f boxA = d_nodesA[nodeIDA].bounds;
  uint32_t localIntersections = 0;

  for (uint32_t idxB = 0; idxB < countB; ++idxB) {
    uint32_t nodeIDB = d_indicesB[idxB];
    cuBQL::box3f boxB = d_nodesB[nodeIDB].bounds;
    if (boxBoxesIntersect(boxA, boxB)) { 
      localIntersections++; 
    }
  }
  d_perThreadCounts[idxA] = localIntersections;
}

// Pass 2: Fill parallel compact allocations completely collision-free
__global__ void fillTargetNodesIntersections(
    const BvhNodeT* d_nodesA, const uint32_t* d_indicesA, uint32_t countA,
    const BvhNodeT* d_nodesB, const uint32_t* d_indicesB, uint32_t countB,
    const uint32_t* d_offsets,
    uint32_t* d_outPairsA, uint32_t* d_outPairsB) 
{
  uint32_t idxA = threadIdx.x + blockIdx.x * blockDim.x;
  if (idxA >= countA) return;

  // We still need nodeIDA to fetch the bounding box for intersection testing
  uint32_t nodeIDA = d_indicesA[idxA];
  cuBQL::box3f boxA = d_nodesA[nodeIDA].bounds;
  uint32_t writePos = d_offsets[idxA];

  for (uint32_t idxB = 0; idxB < countB; ++idxB) {
    uint32_t nodeIDB = d_indicesB[idxB];
    cuBQL::box3f boxB = d_nodesB[nodeIDB].bounds;
    
    if (boxBoxesIntersect(boxA, boxB)) { 
      // Write the index of the marked node array, NOT the BVH Node ID
      d_outPairsA[writePos] = idxA;
      d_outPairsB[writePos] = idxB;
      writePos++;
    }
  }
}

// Orchestrator function called by your main pipeline code
extern "C" uint32_t executeCrissCrossIntersection(
    const cuBQL::bvh3f& bvhA,
    const thrust::device_vector<uint32_t>& d_markedNodeIndicesA,
    uint32_t h_outMarkedCountA,
    const cuBQL::bvh3f& bvhB,
    const thrust::device_vector<uint32_t>& d_markedNodeIndicesB,
    uint32_t h_outMarkedCountB,
    thrust::device_vector<uint32_t>& d_outPairsA,
    thrust::device_vector<uint32_t>& d_outPairsB)
{
  if (h_outMarkedCountA == 0 || h_outMarkedCountB == 0) return 0;

  uint32_t threadsCross = 256;
  uint32_t blocksCross = (h_outMarkedCountA + threadsCross - 1) / threadsCross;

  // Pass 1: Count
  thrust::device_vector<uint32_t> d_perThreadCounts(h_outMarkedCountA, 0);
  countTargetNodesIntersections<<<blocksCross, threadsCross>>>(
      bvhA.nodes, thrust::raw_pointer_cast(d_markedNodeIndicesA.data()), h_outMarkedCountA,
      bvhB.nodes, thrust::raw_pointer_cast(d_markedNodeIndicesB.data()), h_outMarkedCountB,
      thrust::raw_pointer_cast(d_perThreadCounts.data()));
  cudaDeviceSynchronize();

  // Prefix Scan
  thrust::device_vector<uint32_t> d_offsets(h_outMarkedCountA, 0);
  thrust::exclusive_scan(d_perThreadCounts.begin(), d_perThreadCounts.end(), d_offsets.begin());
  
  uint32_t lastCount = d_perThreadCounts.back();
  uint32_t lastOffset = d_offsets.back();
  uint32_t totalIntersections = lastOffset + lastCount;

  if (totalIntersections > 0) {
    d_outPairsA.resize(totalIntersections);
    d_outPairsB.resize(totalIntersections);

    // Pass 2: Fill
    fillTargetNodesIntersections<<<blocksCross, threadsCross>>>(
        bvhA.nodes, thrust::raw_pointer_cast(d_markedNodeIndicesA.data()), h_outMarkedCountA,
        bvhB.nodes, thrust::raw_pointer_cast(d_markedNodeIndicesB.data()), h_outMarkedCountB,
        thrust::raw_pointer_cast(d_offsets.data()),
        thrust::raw_pointer_cast(d_outPairsA.data()),
        thrust::raw_pointer_cast(d_outPairsB.data()));
    cudaDeviceSynchronize();
  }

  return totalIntersections;
}