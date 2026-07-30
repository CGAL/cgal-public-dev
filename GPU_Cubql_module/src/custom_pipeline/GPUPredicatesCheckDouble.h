#ifndef GPU_PREDICATES_CHECK_DOUBLE_H
#define GPU_PREDICATES_CHECK_DOUBLE_H

#include <cuda_runtime.h>
#include <vector_types.h>
#include "PairStatus.h"
#include "TriangleDouble.h"

// Double-precision triangle primitive layout matching 64-bit coordinates


/**
 * Executes double-precision geometric interval evaluation on overlapping candidate pairs.
 * Evaluates pairs in raw 64-bit precision without downcast drift or metric buffers.
 * Pairs containing 0 in their determinant intervals return PAIR_YELLOW for CPU exact fallback.
 */
void evaluateAndCompactPairsDouble(
    int2* dCandidatePairs,
    int* dPairStatuses,
    const TriangleDouble* dA,
    const TriangleDouble* dB_batch,
    int totalBatchPairs,
    cudaStream_t stream = 0);

#endif // GPU_PREDICATES_CHECK_DOUBLE_H