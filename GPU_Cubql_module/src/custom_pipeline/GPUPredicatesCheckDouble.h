#ifndef GPU_PREDICATES_CHECK_DOUBLE_H
#define GPU_PREDICATES_CHECK_DOUBLE_H

#include <cuda_runtime.h>
#include <vector_types.h>
#include "PairStatus.h"

/**
 * Executes double-precision geometric interval evaluation on overlapping candidate pairs.
 * Fetches triangle vertices directly on GPU using vertex (double3) and index (uint3) buffers.
 * Pairs containing 0 in their determinant intervals return PAIR_YELLOW for CPU exact fallback.
 */
void evaluateAndCompactPairsDouble(
    int2* dCandidatePairs,
    int* dPairStatuses,
    const double3* dVertsA,
    const uint3* dIndicesA,
    const double3* dVertsB,
    const uint3* dIndicesB,
    int totalBatchPairs,
    cudaStream_t stream = 0);

#endif // GPU_PREDICATES_CHECK_DOUBLE_H