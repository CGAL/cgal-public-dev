#ifndef BATCHED_CROSS_INTERSECTION_DOUBLE_OLD_H
#define BATCHED_CROSS_INTERSECTION_DOUBLE_OLD_H

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

uint64_t executeBatchedCrossIntersectionLoopDoubleOld(
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
    int2*& outFinalExactPairs,       
    size_t& outFinalCount, 
    IntersectionTimeTracker& tracker, Point3 m_centerA, Point3 m_centerB, double3 m_rotA, double3 m_transA, double3 m_rotB, double3 m_transB,
    cudaStream_t stream = 0
);

#endif // BATCHED_CROSS_INTERSECTION_DOUBLE_OLD_H