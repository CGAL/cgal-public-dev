#ifndef BRUTE_BATCHED_CROSS_INTERSECTION_H
#define BRUTE_BATCHED_CROSS_INTERSECTION_H

#include <thrust/device_vector.h>
#include <vector_types.h>
#include <iostream>
#include <vector>
#include "samples/common/loadOBJ.h"

#include "IntersectionTimeTracker.h"



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
    const thrust::device_vector<uint32_t>& d_outOffsetsA,
    const thrust::device_vector<uint32_t>& d_outPrimsFlatA,
    const thrust::device_vector<uint32_t>& d_nodeDescendantCountsA, 
    uint32_t h_outMarkedCountA,
    uint32_t h_outMarkedCountB,
    const cuBQL::Triangle* dMeshA,
    const cuBQL::Triangle* dMeshB,
    std::vector<int2>& hGreenPairs,
    std::vector<int2>& hYellowPairs,
    IntersectionTimeTracker& tracker 
);

#endif // BATCHED_CROSS_INTERSECTION_H