#include "GPUPredicatesCheckV2.h"
#include <cuda_runtime.h>
#include <thrust/device_ptr.h>
#include <thrust/copy.h>
#include <thrust/execution_policy.h>
#include <thrust/device_vector.h>
#include <algorithm>

// Ensure proper macro safety mapping context for CUBQL allocations within this module
#ifndef CUBQL_CUDA_CALL
#define CUBQL_CUDA_CALL(X) { cudaError_t err = cuda ## X; if(err != cudaSuccess) { printf("CUDA Error at %s:%d\n", __FILE__, __LINE__); } }
#endif

// --------------------------------------------------------------------
// INTERVAL ARITHMETIC & GEOMETRIC CLASS CONTEXT (DEVICE INTERNAL)
// --------------------------------------------------------------------
struct GPUInterval {
    float lo;
    float hi;
};

__device__ inline GPUInterval interval_sub(float a, float b) { return { __fsub_rd(a, b), __fsub_ru(a, b) }; }
__device__ inline GPUInterval interval_add(GPUInterval a, GPUInterval b) { return { __fadd_rd(a.lo, b.lo), __fadd_ru(a.hi, b.hi) }; }
__device__ inline GPUInterval interval_sub(GPUInterval a, GPUInterval b) { return { __fsub_rd(a.lo, b.hi), __fsub_ru(a.hi, b.lo) }; }

__device__ inline GPUInterval interval_mul(GPUInterval a, GPUInterval b) {
    float lo1 = __fmul_rd(a.lo, b.lo); float lo2 = __fmul_rd(a.lo, b.hi);
    float lo3 = __fmul_rd(a.hi, b.lo); float lo4 = __fmul_rd(a.hi, b.hi);
    float hi1 = __fmul_ru(a.lo, b.lo); float hi2 = __fmul_ru(a.lo, b.hi);
    float hi3 = __fmul_ru(a.hi, b.lo); float hi4 = __fmul_ru(a.hi, b.hi);
    return { fminf(fminf(lo1, lo2), fminf(lo3, lo4)), fmaxf(fmaxf(hi1, hi2), fmaxf(hi3, hi4)) };
}

__device__ int orient3d_interval(
    const cuBQL::vec3f& p, const cuBQL::vec3f& q, 
    const cuBQL::vec3f& r, const cuBQL::vec3f& s, float eps) 
{
    GPUInterval pdx = interval_sub(p.x, s.x); GPUInterval pdy = interval_sub(p.y, s.y); GPUInterval pdz = interval_sub(p.z, s.z);
    GPUInterval qdx = interval_sub(q.x, s.x); GPUInterval qdy = interval_sub(q.y, s.y); GPUInterval qdz = interval_sub(q.z, s.z);
    GPUInterval rdx = interval_sub(r.x, s.x); GPUInterval rdy = interval_sub(r.y, s.y); GPUInterval rdz = interval_sub(r.z, s.z);

    GPUInterval det = interval_add(
        interval_add(
            interval_mul(pdx, interval_sub(interval_mul(qdy, rdz), interval_mul(qdz, rdy))),
            interval_mul(pdy, interval_sub(interval_mul(qdz, rdx), interval_mul(qdx, rdz)))
        ),
        interval_mul(pdz, interval_sub(interval_mul(qdx, rdy), interval_mul(qdy, rdx)))
    );

    if (det.lo > eps) return 1;
    if (det.hi < -eps) return -1;
    return 0; 
}

__device__ inline int edgeTriInterval(
    const cuBQL::vec3f &a, const cuBQL::vec3f &b,
    const cuBQL::vec3f &p, const cuBQL::vec3f &q, const cuBQL::vec3f &r, float eps) 
{
    int s0 = orient3d_interval(p, q, r, a, eps);
    int s1 = orient3d_interval(p, q, r, b, eps);

    if (s0 == 0 || s1 == 0) return PAIR_YELLOW;
    if ((s0 > 0 && s1 > 0) || (s0 < 0 && s1 < 0)) return PAIR_NO;

    int e0 = orient3d_interval(a, b, p, q, eps);
    int e1 = orient3d_interval(a, b, q, r, eps);
    int e2 = orient3d_interval(a, b, r, p, eps);

    if (e0 == 0 || e1 == 0 || e2 == 0) return PAIR_YELLOW;
    if ((e0 >= 0 && e1 >= 0 && e2 >= 0) || (e0 <= 0 && e1 <= 0 && e2 <= 0)) return PAIR_GREEN;

    return PAIR_NO;
}

__device__ inline PairStatus classifyPair(
    const cuBQL::Triangle& A, 
    const cuBQL::Triangle& B,
    const float2 metricsA, // Packed: .x = LA, .y = E2A
    const float2 metricsB) // Packed: .x = LB, .y = E2B
{    
    if (!A.bounds().overlaps(B.bounds())) return PAIR_NO;

    // Fetch coordinates directly for interval evaluations
    const cuBQL::vec3f A0 = A.a, A1 = A.b, A2 = A.c;
    const cuBQL::vec3f B0 = B.a, B1 = B.b, B2 = B.c;

    // O(1) Vector-optimized Scale Bounds Extraction
    //float L = fmaxf(metricsA.x, metricsB.x);
    float L = fmaxf(metricsA.x, metricsB.x);
    //float E2 = fmaxf(metricsA.y, metricsB.y); // max(EA, EB)^2 is equivalent to max(EA^2, EB^2)
    float E2 = fmaxf(metricsA.y, metricsB.y);
    // Compute identical tight volumetric error bound instantly
    const float float_machine_epsilon = 1.1920929e-7f; 
    //float eps = 48.0f * L * E2 * float_machine_epsilon;
    float eps = 8.0f * L * E2 * float_machine_epsilon;

    // --- Geometrical Predicate Evaluation Passes (Unchanged) ---
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
// --------------------------------------------------------------------
// CORE INTERSECTION KERNEL
// --------------------------------------------------------------------
__global__ void evaluateGeometricPairsKernel(
    int *outStatuses, 
    const int2 *candidatePairs, 
    const cuBQL::Triangle *triA, 
    const cuBQL::Triangle *triB, 
    const float2 *triAMetrics, // NEW: Companion metrics array 
    const float2 *triBMetrics, // NEW: Companion metrics array
    int numPairs) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= numPairs) return;

    int2 pair = candidatePairs[tid];
    
    // Low-overhead evaluation execution
    PairStatus status = classifyPair(
        triA[pair.x], 
        triB[pair.y], 
        triAMetrics[pair.x], 
        triBMetrics[pair.y]
    );

    outStatuses[tid] = (int)status;
}

// --------------------------------------------------------------------
// THRUST COMPACTION STRUCTS
// --------------------------------------------------------------------
struct IsTargetPairStatus {
    int targetStatus;
    __host__ __device__ bool operator()(const int &status) const {
        return status == targetStatus;
    }
};

// --------------------------------------------------------------------
// EXPORTED PIPELINE WRAPPER STAGE
// --------------------------------------------------------------------
void evaluateAndCompactPairsV2(
    int2* dCandidatePairs,
    int* dPairStatuses,
    const cuBQL::Triangle* dA,
    const cuBQL::Triangle* dB_batch,
    const float2 *triAMetrics,
    const float2 *triBMetrics,
    int totalBatchPairs,cudaStream_t stream 
    )
{
  
    int threadsPerBlock = 256;
    int blocksPerGrid = (totalBatchPairs + threadsPerBlock - 1) / threadsPerBlock;
    
    evaluateGeometricPairsKernel<<<blocksPerGrid, threadsPerBlock, 0, stream>>>(
         dPairStatuses, dCandidatePairs, dA, dB_batch, triAMetrics, triBMetrics, totalBatchPairs
    );
   

}