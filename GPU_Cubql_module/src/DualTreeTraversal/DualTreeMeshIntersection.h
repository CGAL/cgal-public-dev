#ifndef DUAL_TREE_MESH_INTERSECTION_H
#define DUAL_TREE_MESH_INTERSECTION_H

#include "cuBQL/bvh.h"
#include <thrust/device_vector.h>
#include <vector_types.h>
#include <iostream>
#include <vector>
#include "samples/common/loadOBJ.h"
#include "../testBVH/batchedCrossIntersection.h"



uint64_t executeDualTreeTraversal(
    int batchMultiplier,
    int totalBatches,
    uint32_t maxDescendantsA,      // Max total primitives across subtrees in A
    uint32_t maxDescendantsB,      // Max total primitives across subtrees in B
    float expectedIntersectionDensity, // Value between (0.0f, 1.0f] matching expected overlap complexity
    const thrust::device_vector<uint32_t>& d_outPairsA,
    const thrust::device_vector<uint32_t>& d_outPairsB,
    const thrust::device_vector<uint32_t>& d_markedNodeIndicesA,
    const thrust::device_vector<uint32_t>& d_markedNodeIndicesB,
    const thrust::device_vector<uint32_t>& d_nodeDescendantCountsA, 
    const thrust::device_vector<uint32_t>& d_nodeDescendantCountsB, 
    uint32_t h_outMarkedCountA,
    uint32_t h_outMarkedCountB,
    const cuBQL::bvh3f& bvhA,
    const cuBQL::bvh3f& bvhB,
    const cuBQL::Triangle* dMeshA,
    const cuBQL::Triangle* dMeshB,
    std::vector<int2>& hGreenPairs,  
    std::vector<int2>& hYellowPairs, 
    IntersectionTimeTracker& outLoopTime
);

#endif // BDUAL_TREE_MESH_INTERSECTION_H