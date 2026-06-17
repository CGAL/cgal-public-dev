#ifndef GPU_PREDICATES_CHECK_H
#define GPU_PREDICATES_CHECK_H

#include <vector>
#include "cuBQL/bvh.h" // Needed for cuBQL::Triangle reference
#include "samples/common/loadOBJ.h"
#include "kernelsV2.h"

// Struct to store final intersection pairs


// Enum to manage classification criteria statuses
enum PairStatus {
    PAIR_NO     = 0,
    PAIR_GREEN  = 1,
    PAIR_YELLOW = 2
};

/**
 * Executes geometric evaluation on overlapping candidate pairs, filters them using Thrust,
 * and appends the resulting green and yellow pairs directly into the host output vectors.
 * * This function also cleans up the provided dCandidatePairs, dEvaluatedPairs, and dPairStatuses device memory.
 */
void evaluateAndCompactPairs(
    int2* dCandidatePairs,
    int* dPairStatuses,
    const cuBQL::Triangle* dA,
    const cuBQL::Triangle* dB_batch,
    int totalBatchPairs,
    double& outEvaluateGeometricTime);

#endif // GPU_PREDICATES_CHECK_H