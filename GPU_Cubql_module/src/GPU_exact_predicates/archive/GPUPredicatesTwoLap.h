#ifndef GPU_PREDICATES_TWO_LAP_H
#define GPU_PREDICATES_TWO_LAP_H

#include <cuda_runtime.h>
#include <vector_types.h>
#include "PairStatus.h"

// Executes a two-lap exact coplanar filtering pass on candidate pairs:
// Lap 1: Checks 3D plane determinant = 0 in GPU registers.
// Lap 2: Projects coplanar pairs to 2D and executes orient2d_ds SAT.
// Non-zero ambiguous cases are routed directly to PAIR_YELLOW.
void evaluateTwoLapPairs(
    const int2* dCandidatePairs,
    int* dPairStatuses,
    const double3* dVertsA,
    const uint3* dIndicesA,
    const double3* dVertsB,
    const uint3* dIndicesB,
    int numPairs,
    cudaStream_t stream = 0
);

#endif // GPU_PREDICATES_TWO_LAP_H