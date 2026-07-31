#ifndef GPU_PREDICATES_CHECK_DOUBLE_ASSISTED_H
#define GPU_PREDICATES_CHECK_DOUBLE_ASSISTED_H

#include <cuda_runtime.h>
#include <vector_types.h>
#include "PairStatus.h"

void evaluateAndCompactPairsDoubleAssisted(
    const int2* dCandidatePairs,
    int* dPairStatuses,
    const double3* dVertsA,
    const uint3* dIndicesA,
    const double3* dVertsB,
    const uint3* dIndicesB,
    int totalBatchPairs,
    cudaStream_t stream = 0
);

#endif // GPU_PREDICATES_CHECK_DOUBLE_ASSISTED_H