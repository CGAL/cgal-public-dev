#ifndef GPU_PREDICATES_FUSED_H
#define GPU_PREDICATES_FUSED_H

#include <cuda_runtime.h>
#include <vector_types.h>
#include "cuBQL/bvh.h"
#include "PairStatus.h"
#include "TriangleDouble.h"

/**
 * Fused single/double-precision geometric evaluation kernel wrapper.
 * Evaluates pairs in float interval arithmetic; if ambiguous (PAIR_YELLOW),
 * immediately escalates within the same thread to FP64 interval arithmetic.
 */
void evaluateAndCompactPairsFused(
    int2* dCandidatePairs,
    int* dPairStatuses,
    const cuBQL::Triangle* dA_float,
    const cuBQL::Triangle* dB_float,
    const float2* triAMetrics,
    const float2* triBMetrics,
    const TriangleDouble* dA_double,
    const TriangleDouble* dB_double,
    int totalBatchPairs,
    cudaStream_t stream = 0
);

#endif // GPU_PREDICATES_FUSED_H