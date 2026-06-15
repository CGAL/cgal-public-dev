#include "volumeSanityCheck.h"
#include <thrust/host_vector.h>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <cmath>

void runVolumeSanityCheck(
    const cuBQL::bvh3f& bvhB,
    uint32_t h_outMarkedCountB,
    const thrust::device_vector<uint32_t>& d_markedNodeIndicesB,
    const thrust::device_vector<uint32_t>& d_nodeDescendantCountsB,
    const thrust::device_vector<uint32_t>& d_outOffsetsB,
    const thrust::device_vector<uint32_t>& d_outPrimsFlatB,
    const cuBQL::Triangle* hMeshB)
{
  if (h_outMarkedCountB == 0) return;

  std::cout << "\n--------------------------------------------------------------------\n";
  std::cout << "--- EXTRACTED TARGET THRESHOLD NODES (GEOMETRIC EVALUATION - MESH B) ---\n";
  std::cout << "--------------------------------------------------------------------\n";
  
  auto oldFlags = std::cout.flags();
  std::cout.setf(std::ios::fixed, std::ios::floatfield);
std::streamsize oldPrec = std::cout.precision(5);

  std::vector<BvhNodeT> h_nodesVectorB(bvhB.numNodes);
  cudaMemcpy(h_nodesVectorB.data(), bvhB.nodes, bvhB.numNodes * sizeof(BvhNodeT), cudaMemcpyDefault);
  
  thrust::host_vector<uint32_t> h_countsHostB = d_nodeDescendantCountsB;
  thrust::host_vector<uint32_t> h_markedIndicesHostB = d_markedNodeIndicesB;

  auto& rootNode = h_nodesVectorB[0];
  float rdx = rootNode.bounds.upper[0] - rootNode.bounds.lower[0];
  float rdy = rootNode.bounds.upper[1] - rootNode.bounds.lower[1];
  float rdz = rootNode.bounds.upper[2] - rootNode.bounds.lower[2];
  float globalSceneVolume = (rdx > 0.0f && rdy > 0.0f && rdz > 0.0f) ? (rdx * rdy * rdz) : 1.0f;

  thrust::host_vector<uint32_t> h_outOffsetsB = d_outOffsetsB;
  thrust::host_vector<uint32_t> h_outPrimsFlatB = d_outPrimsFlatB;

  for(size_t i = 0; i < h_outMarkedCountB; ++i) {
    uint32_t nodeID = h_markedIndicesHostB[i];
    uint32_t expectedPrimsCount = h_countsHostB[nodeID];

    auto& node = h_nodesVectorB[nodeID];
    float dx = node.bounds.upper[0] - node.bounds.lower[0];
    float dy = node.bounds.upper[1] - node.bounds.lower[1];
    float dz = node.bounds.upper[2] - node.bounds.lower[2];
    float nodeVolume = (dx > 0.0f && dy > 0.0f && dz > 0.0f) ? (dx * dy * dz) : 0.0f;
    float metadataVolumePercentage = (nodeVolume / globalSceneVolume) * 100.0f;

    uint32_t startOffset = h_outOffsetsB[i];
    uint32_t endOffset = (i + 1 < h_outMarkedCountB) ? h_outOffsetsB[i + 1] : (uint32_t)h_outPrimsFlatB.size();
    uint32_t harvestedPrimsCount = endOffset - startOffset;

    cuBQL::box3f primitiveHarvestBBox;
    primitiveHarvestBBox.lower[0] = primitiveHarvestBBox.lower[1] = primitiveHarvestBBox.lower[2] = INFINITY;
    primitiveHarvestBBox.upper[0] = primitiveHarvestBBox.upper[1] = primitiveHarvestBBox.upper[2] = -INFINITY;

    for (uint32_t pIdx = startOffset; pIdx < endOffset; ++pIdx) {
      uint32_t actualPrimID = h_outPrimsFlatB[pIdx];
      cuBQL::box3f primBox = hMeshB[actualPrimID].bounds();
      
      for(int axis = 0; axis < 3; ++axis) {
        primitiveHarvestBBox.lower[axis] = std::min(primitiveHarvestBBox.lower[axis], primBox.lower[axis]);
        primitiveHarvestBBox.upper[axis] = std::max(primitiveHarvestBBox.upper[axis], primBox.upper[axis]);
      }
    }

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

  std::cout.flags(oldFlags);
  std::cout.precision(oldPrec);
}