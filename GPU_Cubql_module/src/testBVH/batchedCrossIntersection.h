#ifndef BATCHED_CROSS_INTERSECTION_H
#define BATCHED_CROSS_INTERSECTION_H

#include "cuBQL/bvh.h"
#include <thrust/device_vector.h>
#include <vector_types.h>
#include <iostream>
#include <vector>
#include "samples/common/loadOBJ.h"

#include "IntersectionTimeTracker.h"

// Simple functor definition required for Thrust compaction filtering
struct IsTargetPairStatus {
    int target;
    __host__ __device__ bool operator()(const int status) const {
        return status == target;
    }
};

uint64_t executeBatchedCrossIntersectionLoop(
    int batchMultiplier,
    int totalBatches,
    const thrust::device_vector<uint32_t>& d_outPairsA,
    const thrust::device_vector<uint32_t>& d_outPairsB,
    const thrust::device_vector<uint32_t>& d_markedNodeIndicesA,
    const thrust::device_vector<uint32_t>& d_markedNodeIndicesB,
    const thrust::device_vector<uint32_t>& d_outOffsetsB,
    const thrust::device_vector<uint32_t>& d_outPrimsFlatB,
    const thrust::device_vector<uint32_t>& d_nodeDescendantCountsB, 
    uint32_t h_outMarkedCountB,
    const cuBQL::bvh3f& bvhA,
    const cuBQL::Triangle* dMeshA,
    const cuBQL::Triangle* dMeshB,
    std::vector<int2>& hGreenPairs,  // Output target for confirmed intersections
    std::vector<int2>& hYellowPairs, // Output target for coplanar / boundary elements
    IntersectionTimeTracker& outLoopTime
);

#endif // BATCHED_CROSS_INTERSECTION_H