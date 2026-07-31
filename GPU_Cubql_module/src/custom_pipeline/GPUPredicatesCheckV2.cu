#include "GPUPredicatesCheckV2.h"
#include "GPUPredicatesCommon.h"
#include <thrust/device_ptr.h>
#include <thrust/copy.h>
#include <thrust/execution_policy.h>
#include <thrust/device_vector.h>
#include <algorithm>

__device__ inline PairStatus classifyPair(
    const cuBQL::Triangle& A, 
    const cuBQL::Triangle& B,
    const float2 metricsA,
    const float2 metricsB)
{    
    if (!A.bounds().overlaps(B.bounds())) return PAIR_NO;

    const cuBQL::vec3f A0 = A.a, A1 = A.b, A2 = A.c;
    const cuBQL::vec3f B0 = B.a, B1 = B.b, B2 = B.c;

    float L = fmaxf(metricsA.x, metricsB.x);
    float E2 = fmaxf(metricsA.y, metricsB.y);
    const float float_machine_epsilon = 1.1920929e-7f; 
    float eps = 8.0f * L * E2 * float_machine_epsilon;

    int ob0 = orient3d_interval(A0, A1, A2, B0, eps);
    int ob1 = orient3d_interval(A0, A1, A2, B1, eps);
    int ob2 = orient3d_interval(A0, A1, A2, B2, eps);
    if (ob0 != 0 && ob1 != 0 && ob2 != 0) {
        if ((ob0 > 0 && ob1 > 0 && ob2 > 0) || (ob0 < 0 && ob1 < 0 && ob2 < 0)) return PAIR_NO;
    }

    int oa0 = orient3d_interval(B0, B1, B2, A0, eps);
    int oa1 = orient3d_interval(B0, B1, B2, A1, eps);
    int oa2 = orient3d_interval(B0, B1, B2, A2, eps);
    if (oa0 != 0 && oa1 != 0 && oa2 != 0) {
        if ((oa0 > 0 && oa1 > 0 && oa2 > 0) || (oa0 < 0 && oa1 < 0 && oa2 < 0)) return PAIR_NO;
    }

    int r;
    r = edgeTriInterval(A0, A1, B0, B1, B2, eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;
    r = edgeTriInterval(A1, A2, B0, B1, B2, eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;
    r = edgeTriInterval(A2, A0, B0, B1, B2, eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;

    r = edgeTriInterval(B0, B1, A0, A1, A2, eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;
    r = edgeTriInterval(B1, B2, A0, A1, A2, eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;
    r = edgeTriInterval(B2, B0, A0, A1, A2, eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;

    return PAIR_NO;
}

__global__ void evaluateGeometricPairsKernelV2(
    int *outStatuses, 
    const int2 *candidatePairs, 
    const cuBQL::Triangle *triA, 
    const cuBQL::Triangle *triB, 
    const float2 *triAMetrics, 
    const float2 *triBMetrics, 
    int numPairs) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= numPairs) return;

    int2 pair = candidatePairs[tid];
    
    PairStatus status = classifyPair(
        triA[pair.x], 
        triB[pair.y], 
        triAMetrics[pair.x], 
        triBMetrics[pair.y]
    );

    outStatuses[tid] = (int)status;
}

struct IsTargetPairStatus {
    int targetStatus;
    __host__ __device__ bool operator()(const int &status) const {
        return status == targetStatus;
    }
};

void evaluateAndCompactPairsV2(
    int2* dCandidatePairs,
    int* dPairStatuses,
    const cuBQL::Triangle* dA,
    const cuBQL::Triangle* dB_batch,
    const float2 *triAMetrics,
    const float2 *triBMetrics,
    int totalBatchPairs, 
    cudaStream_t stream)
{
    int threadsPerBlock = 256;
    int blocksPerGrid = (totalBatchPairs + threadsPerBlock - 1) / threadsPerBlock;
    
    evaluateGeometricPairsKernelV2<<<blocksPerGrid, threadsPerBlock, 0, stream>>>(
         dPairStatuses, dCandidatePairs, dA, dB_batch, triAMetrics, triBMetrics, totalBatchPairs
    );
}