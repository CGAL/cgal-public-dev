#ifndef GPU_PREDICATE_DOUBLE_H
#define GPU_PREDICATE_DOUBLE_H

#include <cuda_runtime.h>
#include <vector_types.h>
#include "PairStatus.h"

// Executes a hardware FP64 double-precision evaluation on the yellow candidate list.
// Outputs PAIR_GREEN (intersecting), PAIR_NO (separated), or PAIR_YELLOW (uncertain/exact CPU fallback).
void evaluateAndCompactPairsDoubleFull(
    const int2* dYellowCandidatePairs,
    int* dPairStatuses,
    const double3* dVertsA,
    const uint3* dIndicesA,
    const double3* dVertsB,
    const uint3* dIndicesB,
    int numYellowPairs,
    cudaStream_t stream = 0
);

#endif // GPU_PREDICATE_DOUBLE_H