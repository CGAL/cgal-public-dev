#pragma once
#include "cuBQL/bvh.h"
#include <thrust/device_vector.h>
#include <vector>
#include "samples/common/loadOBJ.h"

typedef typename cuBQL::BinaryBVH<float, 3>::Node BvhNodeT;

void runVolumeSanityCheck(
    const cuBQL::bvh3f& bvhB,
    uint32_t h_outMarkedCountB,
    const thrust::device_vector<uint32_t>& d_markedNodeIndicesB,
    const thrust::device_vector<uint32_t>& d_nodeDescendantCountsB,
    const thrust::device_vector<uint32_t>& d_outOffsetsB,
    const thrust::device_vector<uint32_t>& d_outPrimsFlatB,
    const cuBQL::Triangle* hMeshB
);