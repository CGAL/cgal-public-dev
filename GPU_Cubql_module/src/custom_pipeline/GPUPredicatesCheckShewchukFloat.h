#ifndef GPU_PREDICATES_CHECK_SHEWCHUK_FLOAT_H
#define GPU_PREDICATES_CHECK_SHEWCHUK_FLOAT_H

#include <cuda_runtime.h>
#include <vector_types.h>
#include "PairStatus.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Evaluate a list of geometric pairs (triangles) using Shewchuk's exact
 * predicates in single precision. This is intended for small lists of
 * "yellow" pairs that could not be resolved by floating‑point filters.
 */
void evaluateAndCompactPairsShewchukFloat(
    const int2* dYellowCandidatePairs,
    int* dPairStatuses,
    const double3* dVertsA,
    const uint3* dIndicesA,
    const double3* dVertsB,
    const uint3* dIndicesB,
    int numYellowPairs,
    cudaStream_t stream = 0
);

#ifdef __cplusplus
}
#endif

#endif // GPU_PREDICATES_CHECK_SHEWCHUK_FLOAT_H