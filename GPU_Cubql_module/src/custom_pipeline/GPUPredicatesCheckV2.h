#ifndef GPU_PREDICATES_CHECK_V2_H
#define GPU_PREDICATES_CHECK_V2_H

#include <vector>
#include "cuBQL/bvh.h" // Needed for cuBQL::Triangle reference
#include "samples/common/loadOBJ.h"
#include "kernelsV2.h"
#include "PairStatus.h"

// Struct to store final intersection pairs


// Enum to manage classification criteria statuses


/**
 * Executes geometric evaluation on overlapping candidate pairs, filters them using Thrust,
 * and appends the resulting green and yellow pairs directly into the host output vectors.
 * * This function also cleans up the provided dCandidatePairs, dEvaluatedPairs, and dPairStatuses device memory.
 */
void evaluateAndCompactPairsV2(
    int2* dCandidatePairs,
    int* dPairStatuses,
    const cuBQL::Triangle* dA,
    const cuBQL::Triangle* dB_batch,
    const float2 *triAMetrics,
    const float2 *triBMetrics,
    int totalBatchPairs,
    cudaStream_t stream);

#endif // GPU_PREDICATES_CHECK_H