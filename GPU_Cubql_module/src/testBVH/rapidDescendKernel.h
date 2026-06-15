#pragma once
#include "cuBQL/bvh.h"
#include <thrust/device_vector.h>

// Struct tracking the layout state across parallel levels
struct FrontierElement
{
  uint32_t nodeID;     
  uint32_t writeCursor; 
};

// C-linkage interface wrapper exposing the pipeline to main/other modules
extern "C" void executeRapidDescentBFS(
    const cuBQL::bvh3f& bvhA,
    int numTrianglesA,
    const thrust::device_vector<uint32_t>& d_markedNodeIndicesA,
    const thrust::device_vector<uint32_t>& d_nodeDescendantCountsA, 
    uint32_t h_outMarkedCountA,
    thrust::device_vector<uint32_t>& d_outOffsets,                 
    thrust::device_vector<uint32_t>& d_outPrimsFlat
);