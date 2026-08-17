#ifndef GPU_PREDICATES_CHECK_BIGINTEGER_V2_H
#define GPU_PREDICATES_CHECK_BIGINTEGER_V2_H

#include <cuda_runtime.h>
#include <vector_types.h>
#include "PairStatus.h"

#ifdef __cplusplus
extern "C" {
#endif

void evaluateGeometricPairsKernelBigIntV2(
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

#endif // GPU_PREDICATES_CHECK_BIGINTEGER_V2_H