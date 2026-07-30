#ifndef BATCHED_CROSS_INTERSECTION_DOUBLE_H
#define BATCHED_CROSS_INTERSECTION_DOUBLE_H

#include "cuBQL/bvh.h"
#include <thrust/device_vector.h>
#include <vector_types.h>
#include <iostream>
#include <vector>
#include <tbb/concurrent_vector.h>

#include "IntersectionTimeTracker.h"
#include "../custom_pipeline/TriangleDouble.h"
#include "../custom_pipeline/GPUPredicatesCheckDouble.h"
#include "../src/CPU/CgalDefinitions.h"
#include "samples/common/loadOBJ.h"

/**
 * Batched Cross-Intersection Loop with GPU Double-Precision Filter Stage.
 *
 * If `dMeshADouble` and `dMeshBDouble` are provided (!= nullptr), candidates flagged
 * as PAIR_YELLOW during the float evaluation pass are re-evaluated on the GPU in 
 * 64-bit precision interval arithmetic before resorting to CPU CGAL exact filtering.
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
    const cuBQL::Triangle* dMeshA,
    const cuBQL::Triangle* dMeshB,
    const float2 *triAMetrics,
    const float2 *triBMetrics,
    tbb::concurrent_vector<int2> & finalExactPairs,
    IntersectionTimeTracker& tracker,
    cudaStream_t stream = 0,
    const TriangleDouble* dMeshADouble = nullptr, // Pass nullptr to skip GPU Double evaluation
    const TriangleDouble* dMeshBDouble = nullptr
);

#endif // BATCHED_CROSS_INTERSECTION_DOUBLE_H