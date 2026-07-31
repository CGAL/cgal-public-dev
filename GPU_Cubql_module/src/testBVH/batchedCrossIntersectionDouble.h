#ifndef BATCHED_CROSS_INTERSECTION_DOUBLE_H
#define BATCHED_CROSS_INTERSECTION_DOUBLE_H

#include "cuBQL/bvh.h"
#include <thrust/device_vector.h>
#include <vector_types.h>
#include <iostream>
#include <vector>
#include <tbb/concurrent_vector.h>
#include "IntersectionTimeTracker.h"

#include "../custom_pipeline/GPUPredicatesCheckDouble.h"
#include "../src/CPU/CgalDefinitions.h"
#include "samples/common/loadOBJ.h"

/**
 * Batched Cross-Intersection Loop using direct precomputed cuBQL::box3f arrays
 * for coarse AABB overlap checks, and double3 vertex & uint3 index buffers for exact GPU predicates.
 */
uint64_t executeBatchedCrossIntersectionLoopDouble(
    Mesh & meshAcpu, Mesh & meshBcpu,
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
    const cuBQL::box3f* d_boxesA,
    const cuBQL::box3f* d_boxesB,
    const double3* d_vertsA,
    const uint3* d_indicesA,
    const double3* d_vertsB,
    const uint3* d_indicesB,
    tbb::concurrent_vector<int2> & finalExactPairs,
    IntersectionTimeTracker& tracker,
    cudaStream_t stream = 0
);

#endif // BATCHED_CROSS_INTERSECTION_DOUBLE_H