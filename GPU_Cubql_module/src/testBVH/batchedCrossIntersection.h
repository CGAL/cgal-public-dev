#ifndef BATCHED_CROSS_INTERSECTION_H
#define BATCHED_CROSS_INTERSECTION_H

#include "cuBQL/bvh.h"
#include <thrust/device_vector.h>

#include "samples/common/loadOBJ.h"

/**
 * Executes the multi-stream batched bounding box cross-intersection loop.
 * Returns the final total number of AABB candidate pairs found.
 */
uint64_t executeBatchedCrossIntersectionLoop(
    int batchMultiplier,
    int totalBatches,
    const thrust::device_vector<uint32_t>& d_outPairsA,
    const thrust::device_vector<uint32_t>& d_outPairsB,
    const thrust::device_vector<uint32_t>& d_markedNodeIndicesA,
    const thrust::device_vector<uint32_t>& d_markedNodeIndicesB,
    const thrust::device_vector<uint32_t>& d_outOffsetsB,
    const thrust::device_vector<uint32_t>& d_outPrimsFlatB,
    const thrust::device_vector<uint32_t>& d_nodeDescendantCountsB, // <-- ADD THIS PARAMETER
    uint32_t h_outMarkedCountB,
    const cuBQL::bvh3f& bvhA,
    const cuBQL::Triangle* dMeshA,
    const cuBQL::Triangle* dMeshB,
    double& outLoopTime
);



#endif // BATCHED_CROSS_INTERSECTION_H