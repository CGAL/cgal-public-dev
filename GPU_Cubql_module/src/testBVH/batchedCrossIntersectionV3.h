#ifndef BATCHED_CROSS_INTERSECTION_V3_H
#define BATCHED_CROSS_INTERSECTION_V3_H

#include "cuBQL/bvh.h"
#include <thrust/device_vector.h>
#include <vector_types.h>
#include <iostream>
#include <vector>
#include "samples/common/loadOBJ.h"

#include <tbb/concurrent_vector.h>

#include "IntersectionTimeTracker.h"
#include "../src/CPU/CgalDefinitions.h"

// Simple functor definition required for Thrust compaction filtering

uint64_t executeBatchedCrossIntersectionLoopV3(Mesh & meshAcpu, Mesh & MeshBcpu,
    int batchMultiplier,
    int totalBatches,
    const thrust::device_vector<uint32_t>& d_outPairsA,
    const thrust::device_vector<uint32_t>& d_outPairsB,
    const thrust::device_vector<uint32_t>& d_reverseMapB,
    const thrust::device_vector<uint32_t>& d_markedNodeIndicesB, 
    const thrust::device_vector<uint32_t>& d_outOffsetsB,
    const thrust::device_vector<uint32_t>& d_outPrimsFlatB,
    const thrust::device_vector<uint32_t>& d_nodeDescendantCountsB, 
    uint32_t h_outMarkedCountB,
    const cuBQL::bvh3f& bvhA,
    const cuBQL::Triangle* dMeshA,
    const cuBQL::Triangle* dMeshB,
    const float2 *triAMetrics,
    const float2 *triBMetrics,
    tbb::concurrent_vector<int2> & finalExactPairs,
    IntersectionTimeTracker& tracker,cudaStream_t stream 
);

//std::vector<int2>& hGreenPairs,
//    std::vector<int2>& hYellowPairs,

#endif // BATCHED_CROSS_INTERSECTION_H