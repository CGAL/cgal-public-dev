#include "crossCheckFlexible.h"
#include <thrust/scan.h>

// Standard bounding box overlap check (unmodified)
__device__ inline bool boxBoxesIntersect(const cuBQL::box3f& a, const cuBQL::box3f& b) {
  return (a.lower[0] <= b.upper[0] && a.upper[0] >= b.lower[0]) &&
         (a.lower[1] <= b.upper[1] && a.upper[1] >= b.lower[1]) &&
         (a.lower[2] <= b.upper[2] && a.upper[2] >= b.lower[2]);
}

// Pass 1: Count intersections
__global__ void countTargetNodesIntersections(
    const BvhNodeT* d_nodesA, const uint32_t* d_indicesA, uint32_t countA,
    const BvhNodeT* d_nodesB, const uint32_t* d_indicesB, uint32_t countB,
    uint32_t* d_perThreadCounts, float tx, float ty, float tz) 
{
  uint32_t idxB = threadIdx.x + blockIdx.x * blockDim.x;
  if (idxB >= countB) return;

  uint32_t nodeIDB = d_indicesB[idxB];
  
  // Fetch boxB and shift it ONCE per thread before entering the loop
  cuBQL::box3f boxB = d_nodesB[nodeIDB].bounds;
  boxB.lower[0] += tx; boxB.upper[0] += tx;
  boxB.lower[1] += ty; boxB.upper[1] += ty;
  boxB.lower[2] += tz; boxB.upper[2] += tz;

  uint32_t localIntersections = 0;

  for (uint32_t idxA = 0; idxA < countA; ++idxA) {
    uint32_t nodeIDA = d_indicesA[idxA];
    cuBQL::box3f boxA = d_nodesA[nodeIDA].bounds;
    
    // Standard intersection check!
    if (boxBoxesIntersect(boxA, boxB)) { 
      localIntersections++; 
    }
  }
  d_perThreadCounts[idxB] = localIntersections;
}

// Pass 2: Fill parallel allocations
__global__ void fillTargetNodesIntersections(
    const BvhNodeT* d_nodesA, const uint32_t* d_indicesA, uint32_t countA,
    const BvhNodeT* d_nodesB, const uint32_t* d_indicesB, uint32_t countB,
    const uint32_t* d_offsets,
    uint32_t* d_outPairsA, uint32_t* d_outPairsB,
    float tx, float ty, float tz) 
{
  uint32_t idxB = threadIdx.x + blockIdx.x * blockDim.x;
  if (idxB >= countB) return;

  uint32_t nodeIDB = d_indicesB[idxB];

  // Fetch and shift boxB ONCE per thread
  cuBQL::box3f boxB = d_nodesB[nodeIDB].bounds;
  boxB.lower[0] += tx; boxB.upper[0] += tx;
  boxB.lower[1] += ty; boxB.upper[1] += ty;
  boxB.lower[2] += tz; boxB.upper[2] += tz;

  uint32_t writePos = d_offsets[idxB];

  for (uint32_t idxA = 0; idxA < countA; ++idxA) {
    uint32_t nodeIDA = d_indicesA[idxA];
    cuBQL::box3f boxA = d_nodesA[nodeIDA].bounds;
    
    if (boxBoxesIntersect(boxA, boxB)) { 
      d_outPairsA[writePos] = nodeIDA;
      d_outPairsB[writePos] = nodeIDB;
      writePos++;
    }
  }
}

// Orchestrator function
extern "C" uint32_t executeCrossIntersectionFlexible(
    const cuBQL::bvh3f& bvhA,
    const thrust::device_vector<uint32_t>& d_markedNodeIndicesA,
    uint32_t h_outMarkedCountA,
    const cuBQL::bvh3f& bvhB,
    const thrust::device_vector<uint32_t>& d_markedNodeIndicesB,
    uint32_t h_outMarkedCountB,
    thrust::device_vector<uint32_t>& d_outPairsA,
    thrust::device_vector<uint32_t>& d_outPairsB, 
    float tx, float ty, float tz)
{
  if (h_outMarkedCountA == 0 || h_outMarkedCountB == 0) return 0;

  uint32_t threadsCross = 256;
  uint32_t blocksCross = (h_outMarkedCountB + threadsCross - 1) / threadsCross;

  // Pass 1: Count
  thrust::device_vector<uint32_t> d_perThreadCounts(h_outMarkedCountB, 0);
  countTargetNodesIntersections<<<blocksCross, threadsCross>>>(
      bvhA.nodes, thrust::raw_pointer_cast(d_markedNodeIndicesA.data()), h_outMarkedCountA,
      bvhB.nodes, thrust::raw_pointer_cast(d_markedNodeIndicesB.data()), h_outMarkedCountB,
      thrust::raw_pointer_cast(d_perThreadCounts.data()),
      tx, ty, tz);
  cudaDeviceSynchronize();

  // Prefix Scan
  thrust::device_vector<uint32_t> d_offsets(h_outMarkedCountB, 0);
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
        thrust::raw_pointer_cast(d_outPairsB.data()),
        tx, ty, tz);
    cudaDeviceSynchronize();
  }

  return totalIntersections;
}