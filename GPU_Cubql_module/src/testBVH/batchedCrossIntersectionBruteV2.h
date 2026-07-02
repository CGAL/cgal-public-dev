#ifndef BRUTE_V2_BATCHED_CROSS_INTERSECTION_H
#define BRUTE_V2_BATCHED_CROSS_INTERSECTION_H

#include <thrust/device_vector.h>
#include <vector_types.h>
#include <iostream>
#include <vector>
#include "samples/common/loadOBJ.h"




template <typename T, int D>
uint64_t executePrimitiveBoxCrossIntersectionCheck(
    int batchMultiplier,
    const thrust::device_vector<uint32_t>& d_intersectPairsA,
    const thrust::device_vector<uint32_t>& d_intersectPairsB,
    const uint32_t* d_outSortedPrimIDsA, const uint32_t* d_outNodeOffsetsA, const cuBQL::box_t<T,D>* d_boxesA,
    const uint32_t* d_outSortedPrimIDsB, const uint32_t* d_outNodeOffsetsB, const cuBQL::box_t<T,D>* d_boxesB,
    cudaStream_t s);

#endif // BATCHED_V2_CROSS_INTERSECTION_H