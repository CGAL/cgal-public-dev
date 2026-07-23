#pragma once
#include "cuBQL/bvh.h"
#include <thrust/device_vector.h>

// Mirror the node type declaration
typedef typename cuBQL::BinaryBVH<float, 3>::Node BvhNodeT;

extern "C" uint32_t executeCrossIntersectionFlexible(
    const cuBQL::bvh3f& bvhA,
    const thrust::device_vector<uint32_t>& d_markedNodeIndicesA,
    uint32_t h_outMarkedCountA,
    const cuBQL::bvh3f& bvhB,
    const thrust::device_vector<uint32_t>& d_markedNodeIndicesB,
    uint32_t h_outMarkedCountB,
    thrust::device_vector<uint32_t>& d_outPairsA,
    thrust::device_vector<uint32_t>& d_outPairsB,
    float tx, float ty, float tz
);