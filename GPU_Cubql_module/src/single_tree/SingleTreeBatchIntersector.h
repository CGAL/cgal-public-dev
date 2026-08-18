#ifndef BATCHED_CROSS_INTERSECTION_DOUBLE_H
#define BATCHED_CROSS_INTERSECTION_DOUBLE_H
// Standard C++ & System Headers

#include <iostream>
#include <vector>

// CUDA & Thrust Headers
#include <vector_types.h>
#include <thrust/device_vector.h>

// Intel TBB Headers
#include <tbb/concurrent_vector.h>

// Vendor / Library Headers
#include <cuBQL/bvh.h>
#include <loadOBJ.h>

// Custom Project Modules (relative to src/single_tree/)
#include "common/IntersectionTimeTracker.h"
#include "GPU_exact_predicates/TwoSumDoublePredicate.h"
#include "CPU/CgalDefinitions.h"
/**
 * Batched Cross-Intersection Loop using direct precomputed cuBQL::box3f arrays
 * for coarse AABB overlap checks, and double3 vertex & uint3 index buffers for exact GPU predicates.
 */
uint64_t executeSingleTreeBatchedTraversalWithPredicates(
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
    int2*& outFinalExactPairs,       // Fast raw pointer return
    size_t& outFinalCount, 
    IntersectionTimeTracker& tracker, Point3 m_centerA, Point3 m_centerB, double3 m_rotA, double3 m_transA, double3 m_rotB, double3 m_transB,
    int exactPredicateComputeMode = 1,
    cudaStream_t stream = 0
);

#endif // BATCHED_CROSS_INTERSECTION_DOUBLE_H