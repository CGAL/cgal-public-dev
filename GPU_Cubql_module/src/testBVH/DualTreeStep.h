#ifndef DUAL_TREE_STEP_H
#define DUAL_TREE_STEP_H

#include "cuBQL/bvh.h"
#include <thrust/device_vector.h>
#include <vector_types.h>
#include <iostream>
#include <vector>
#include "samples/common/loadOBJ.h"

uint64_t executeDualTreeStep(
    int numSteps,
    uint32_t maxDescendantsA,      // Max total primitives across subtrees in A
    uint32_t maxDescendantsB,      // Max total primitives across subtrees in B
    thrust::device_vector<uint32_t>& d_outPairsA,                 // Removed const
    thrust::device_vector<uint32_t>& d_outPairsB,                 // Removed const
    thrust::device_vector<uint32_t>& d_markedNodeIndicesA,       // Removed const
    thrust::device_vector<uint32_t>& d_markedNodeIndicesB,       // Removed const
    thrust::device_vector<uint32_t>& d_nodeDescendantCountsA,    // Removed const
    thrust::device_vector<uint32_t>& d_nodeDescendantCountsB,    // Removed const
    uint32_t& h_outMarkedCountA,                                 // Changed to non-const reference
    uint32_t& h_outMarkedCountB,                                 // Changed to non-const reference
    const cuBQL::bvh3f& bvhA,
    const cuBQL::bvh3f& bvhB                                       // Matches implementation's void*
);

#endif // DUAL_TREE_STEP_H