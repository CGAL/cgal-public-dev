#ifndef GPU_PREDICATES_CHECK_DOUBLE_H
#define GPU_PREDICATES_CHECK_DOUBLE_H

#include <cuda_runtime.h>
#include <vector_types.h>
#include "PairStatus.h"

// Executes a double-single (two-float) precision evaluation on the yellow list.
// Outputs PAIR_GREEN, PAIR_NO, or PAIR_YELLOW (for the final exact CPU fallback).
void evaluateAndCompactPairsDouble(
    const int2* dYellowCandidatePairs,
    int* dPairStatuses,
    const double3* dVertsA,
    const uint3* dIndicesA,
    const double3* dVertsB,
    const uint3* dIndicesB,
    int numYellowPairs,
    cudaStream_t stream = 0
);

#endif // GPU_PREDICATES_CHECK_DOUBLE_H