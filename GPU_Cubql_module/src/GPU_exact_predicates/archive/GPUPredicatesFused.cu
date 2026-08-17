#include "TriangleDouble.h"
#include <cuda_runtime.h>
#include <cstdio>
#include "samples/common/loadOBJ.h"

#ifndef CUBQL_CUDA_CALL
#define CUBQL_CUDA_CALL(X) { cudaError_t err = cuda ## X; if(err != cudaSuccess) { printf("CUDA Error at %s:%d\n", __FILE__, __LINE__); } }
#endif

// ====================================================================
// 1. FLOAT INTERVAL ARITHMETIC
// ====================================================================
struct GPUIntervalFloat {
    float lo;
    float hi;
};

__device__ inline GPUIntervalFloat interval_sub(float a, float b) { return { __fsub_rd(a, b), __fsub_ru(a, b) }; }
__device__ inline GPUIntervalFloat interval_add(GPUIntervalFloat a, GPUIntervalFloat b) { return { __fadd_rd(a.lo, b.lo), __fadd_ru(a.hi, b.hi) }; }
__device__ inline GPUIntervalFloat interval_sub(GPUIntervalFloat a, GPUIntervalFloat b) { return { __fsub_rd(a.lo, b.hi), __fsub_ru(a.hi, b.lo) }; }

__device__ inline GPUIntervalFloat interval_mul(GPUIntervalFloat a, GPUIntervalFloat b) {
    float lo1 = __fmul_rd(a.lo, b.lo); float lo2 = __fmul_rd(a.lo, b.hi);
    float lo3 = __fmul_rd(a.hi, b.lo); float lo4 = __fmul_rd(a.hi, b.hi);
    float hi1 = __fmul_ru(a.lo, b.lo); float hi2 = __fmul_ru(a.lo, b.hi);
    float hi3 = __fmul_ru(a.hi, b.lo); float hi4 = __fmul_ru(a.hi, b.hi);
    return { fminf(fminf(lo1, lo2), fminf(lo3, lo4)), fmaxf(fmaxf(hi1, hi2), fmaxf(hi3, hi4)) };
}

__device__ inline int orient3d_interval_float(
    const cuBQL::vec3f& p, const cuBQL::vec3f& q, 
    const cuBQL::vec3f& r, const cuBQL::vec3f& s, float eps) 
{
    GPUIntervalFloat pdx = interval_sub(p.x, s.x); GPUIntervalFloat pdy = interval_sub(p.y, s.y); GPUIntervalFloat pdz = interval_sub(p.z, s.z);
    GPUIntervalFloat qdx = interval_sub(q.x, s.x); GPUIntervalFloat qdy = interval_sub(q.y, s.y); GPUIntervalFloat qdz = interval_sub(q.z, s.z);
    GPUIntervalFloat rdx = interval_sub(r.x, s.x); GPUIntervalFloat rdy = interval_sub(r.y, s.y); GPUIntervalFloat rdz = interval_sub(r.z, s.z);

    GPUIntervalFloat det = interval_add(
        interval_add(
            interval_mul(pdx, interval_sub(interval_mul(qdy, rdz), interval_mul(qdz, rdy))),
            interval_mul(pdy, interval_sub(interval_mul(qdz, rdx), interval_mul(qdx, rdz)))
        ),
        interval_mul(pdz, interval_sub(interval_mul(qdx, rdy), interval_mul(qdy, rdx)))
    );

    if (det.lo > eps) return 1;
    if (det.hi < -eps) return -1;
    return 0; // FP32 uncertain
}

// ====================================================================
// 2. DOUBLE INTERVAL ARITHMETIC
// ====================================================================
struct GPUIntervalDouble {
    double lo;
    double hi;
};

__device__ inline GPUIntervalDouble interval_sub(double a, double b) { return { __dsub_rd(a, b), __dsub_ru(a, b) }; }
__device__ inline GPUIntervalDouble interval_add(GPUIntervalDouble a, GPUIntervalDouble b) { return { __dadd_rd(a.lo, b.lo), __dadd_ru(a.hi, b.hi) }; }
__device__ inline GPUIntervalDouble interval_sub(GPUIntervalDouble a, GPUIntervalDouble b) { return { __dsub_rd(a.lo, b.hi), __dsub_ru(a.hi, b.lo) }; }

__device__ inline GPUIntervalDouble interval_mul(GPUIntervalDouble a, GPUIntervalDouble b) {
    double lo1 = __dmul_rd(a.lo, b.lo); double lo2 = __dmul_rd(a.lo, b.hi);
    double lo3 = __dmul_rd(a.hi, b.lo); double lo4 = __dmul_rd(a.hi, b.hi);
    double hi1 = __dmul_ru(a.lo, b.lo); double hi2 = __dmul_ru(a.lo, b.hi);
    double hi3 = __dmul_ru(a.hi, b.lo); double hi4 = __dmul_ru(a.hi, b.hi);
    return { fmin(fmin(lo1, lo2), fmin(lo3, lo4)), fmax(fmax(hi1, hi2), fmax(hi3, hi4)) };
}

__device__ inline int orient3d_interval_double(
    const double3& p, const double3& q, 
    const double3& r, const double3& s) 
{
    GPUIntervalDouble pdx = interval_sub(p.x, s.x); GPUIntervalDouble pdy = interval_sub(p.y, s.y); GPUIntervalDouble pdz = interval_sub(p.z, s.z);
    GPUIntervalDouble qdx = interval_sub(q.x, s.x); GPUIntervalDouble qdy = interval_sub(q.y, s.y); GPUIntervalDouble qdz = interval_sub(q.z, s.z);
    GPUIntervalDouble rdx = interval_sub(r.x, s.x); GPUIntervalDouble rdy = interval_sub(r.y, s.y); GPUIntervalDouble rdz = interval_sub(r.z, s.z);

    GPUIntervalDouble det = interval_add(
        interval_add(
            interval_mul(pdx, interval_sub(interval_mul(qdy, rdz), interval_mul(qdz, rdy))),
            interval_mul(pdy, interval_sub(interval_mul(qdz, rdx), interval_mul(qdx, rdz)))
        ),
        interval_mul(pdz, interval_sub(interval_mul(qdx, rdy), interval_mul(qdy, rdx)))
    );

    if (det.lo > 0.0) return 1;
    if (det.hi < 0.0) return -1;
    return 0; // True hard degeneracy -> CPU CGAL fallback
}

// ====================================================================
// 3. PURE FP32 & FP64 TEST IMPLEMENTATIONS
// ====================================================================
__device__ inline int edgeTriIntervalFP32(
    const cuBQL::vec3f &a, const cuBQL::vec3f &b,
    const cuBQL::vec3f &p, const cuBQL::vec3f &q, const cuBQL::vec3f &r,
    float eps, bool &uncertain) 
{
    int s0 = orient3d_interval_float(p, q, r, a, eps);
    int s1 = orient3d_interval_float(p, q, r, b, eps);

    if (s0 == 0 || s1 == 0) uncertain = true;
    if ((s0 > 0 && s1 > 0) || (s0 < 0 && s1 < 0)) return PAIR_NO;

    int e0 = orient3d_interval_float(a, b, p, q, eps);
    int e1 = orient3d_interval_float(a, b, q, r, eps);
    int e2 = orient3d_interval_float(a, b, r, p, eps);

    if (e0 == 0 || e1 == 0 || e2 == 0) uncertain = true;

    if ((e0 >= 0 && e1 >= 0 && e2 >= 0) || (e0 <= 0 && e1 <= 0 && e2 <= 0)) {
        if (s0 != 0 && s1 != 0 && e0 != 0 && e1 != 0 && e2 != 0) {
            return PAIR_GREEN; // Unambiguous FP32 intersection!
        }
    }

    return PAIR_NO;
}

__device__ inline int edgeTriIntervalFP64(
    const double3 &a, const double3 &b,
    const double3 &p, const double3 &q, const double3 &r) 
{
    int s0 = orient3d_interval_double(p, q, r, a);
    int s1 = orient3d_interval_double(p, q, r, b);

    if (s0 == 0 || s1 == 0) return PAIR_YELLOW;
    if ((s0 > 0 && s1 > 0) || (s0 < 0 && s1 < 0)) return PAIR_NO;

    int e0 = orient3d_interval_double(a, b, p, q);
    int e1 = orient3d_interval_double(a, b, q, r);
    int e2 = orient3d_interval_double(a, b, r, p);

    if (e0 == 0 || e1 == 0 || e2 == 0) return PAIR_YELLOW;
    if ((e0 >= 0 && e1 >= 0 && e2 >= 0) || (e0 <= 0 && e1 <= 0 && e2 <= 0)) return PAIR_GREEN;

    return PAIR_NO;
}

// Full Double Fallback Evaluator
__device__ inline PairStatus classifyPairFP64(
    const TriangleDouble& A_d, 
    const TriangleDouble& B_d)
{
    int ob0 = orient3d_interval_double(A_d.a, A_d.b, A_d.c, B_d.a);
    int ob1 = orient3d_interval_double(A_d.a, A_d.b, A_d.c, B_d.b);
    int ob2 = orient3d_interval_double(A_d.a, A_d.b, A_d.c, B_d.c);
    if (ob0 != 0 && ob1 != 0 && ob2 != 0) {
        if ((ob0 > 0 && ob1 > 0 && ob2 > 0) || (ob0 < 0 && ob1 < 0 && ob2 < 0)) return PAIR_NO;
    }

    int oa0 = orient3d_interval_double(B_d.a, B_d.b, B_d.c, A_d.a);
    int oa1 = orient3d_interval_double(B_d.a, B_d.b, B_d.c, A_d.b);
    int oa2 = orient3d_interval_double(B_d.a, B_d.b, B_d.c, A_d.c);
    if (oa0 != 0 && oa1 != 0 && oa2 != 0) {
        if ((oa0 > 0 && oa1 > 0 && oa2 > 0) || (oa0 < 0 && oa1 < 0 && oa2 < 0)) return PAIR_NO;
    }

    int r;
    r = edgeTriIntervalFP64(A_d.a, A_d.b, B_d.a, B_d.b, B_d.c); if (r == PAIR_GREEN || r == PAIR_YELLOW) return (PairStatus)r;
    r = edgeTriIntervalFP64(A_d.b, A_d.c, B_d.a, B_d.b, B_d.c); if (r == PAIR_GREEN || r == PAIR_YELLOW) return (PairStatus)r;
    r = edgeTriIntervalFP64(A_d.c, A_d.a, B_d.a, B_d.b, B_d.c); if (r == PAIR_GREEN || r == PAIR_YELLOW) return (PairStatus)r;

    r = edgeTriIntervalFP64(B_d.a, B_d.b, A_d.a, A_d.b, A_d.c); if (r == PAIR_GREEN || r == PAIR_YELLOW) return (PairStatus)r;
    r = edgeTriIntervalFP64(B_d.b, B_d.c, A_d.a, A_d.b, A_d.c); if (r == PAIR_GREEN || r == PAIR_YELLOW) return (PairStatus)r;
    r = edgeTriIntervalFP64(B_d.c, B_d.a, A_d.a, A_d.b, A_d.c); if (r == PAIR_GREEN || r == PAIR_YELLOW) return (PairStatus)r;

    return PAIR_NO;
}

// ====================================================================
// 4. MACRO-CASCADED CLASSIFIER
// ====================================================================
__device__ inline PairStatus classifyPairCascaded(
    const cuBQL::Triangle& A, 
    const cuBQL::Triangle& B,
    const float2 metricsA,
    const float2 metricsB,
    const TriangleDouble *triA_double,
    const TriangleDouble *triB_double,
    int idxA, int idxB) 
{    
    if (!A.bounds().overlaps(B.bounds())) return PAIR_NO;

    const cuBQL::vec3f A0 = A.a, A1 = A.b, A2 = A.c;
    const cuBQL::vec3f B0 = B.a, B1 = B.b, B2 = B.c;

    float L = fmaxf(metricsA.x, metricsB.x);
    float E2 = fmaxf(metricsA.y, metricsB.y);
    const float float_machine_epsilon = 1.1920929e-7f; 
    float eps = 8.0f * L * E2 * float_machine_epsilon;

    bool fp32_uncertain = false;

    // --- Point-Plane Orientation Passes (A vs B) ---
    int ob0 = orient3d_interval_float(A0, A1, A2, B0, eps);
    int ob1 = orient3d_interval_float(A0, A1, A2, B1, eps);
    int ob2 = orient3d_interval_float(A0, A1, A2, B2, eps);
    if (ob0 != 0 && ob1 != 0 && ob2 != 0) {
        if ((ob0 > 0 && ob1 > 0 && ob2 > 0) || (ob0 < 0 && ob1 < 0 && ob2 < 0)) 
            return PAIR_NO; // Guaranteed separated in FP32
    } else {
        // Mark uncertain ONLY if we haven't already proven straddling
        if (!((ob0 > 0 || ob1 > 0 || ob2 > 0) && (ob0 < 0 || ob1 < 0 || ob2 < 0))) {
            fp32_uncertain = true;
        }
    }

    // --- Point-Plane Orientation Passes (B vs A) ---
    int oa0 = orient3d_interval_float(B0, B1, B2, A0, eps);
    int oa1 = orient3d_interval_float(B0, B1, B2, A1, eps);
    int oa2 = orient3d_interval_float(B0, B1, B2, A2, eps);
    if (oa0 != 0 && oa1 != 0 && oa2 != 0) {
        if ((oa0 > 0 && oa1 > 0 && oa2 > 0) || (oa0 < 0 && oa1 < 0 && oa2 < 0)) 
            return PAIR_NO; // Guaranteed separated in FP32
    } else {
        if (!((oa0 > 0 || oa1 > 0 || oa2 > 0) && (oa0 < 0 || oa1 < 0 || oa2 < 0))) {
            fp32_uncertain = true;
        }
    }

    // --- Edge-Triangle Intersection Passes ---
    int r;
    r = edgeTriIntervalFP32(A0, A1, B0, B1, B2, eps, fp32_uncertain); if (r == PAIR_GREEN) return PAIR_GREEN;
    r = edgeTriIntervalFP32(A1, A2, B0, B1, B2, eps, fp32_uncertain); if (r == PAIR_GREEN) return PAIR_GREEN;
    r = edgeTriIntervalFP32(A2, A0, B0, B1, B2, eps, fp32_uncertain); if (r == PAIR_GREEN) return PAIR_GREEN;

    r = edgeTriIntervalFP32(B0, B1, A0, A1, A2, eps, fp32_uncertain); if (r == PAIR_GREEN) return PAIR_GREEN;
    r = edgeTriIntervalFP32(B1, B2, A0, A1, A2, eps, fp32_uncertain); if (r == PAIR_GREEN) return PAIR_GREEN;
    r = edgeTriIntervalFP32(B2, B0, A0, A1, A2, eps, fp32_uncertain); if (r == PAIR_GREEN) return PAIR_GREEN;

    // If FP32 was completely decisive and showed no intersection:
    if (!fp32_uncertain) return PAIR_NO;

    // ====================================================================
    // LAZY FP64 FALLBACK (Only executes on ~2-5% of candidate pairs)
    // ====================================================================
    TriangleDouble A_d = triA_double[idxA];
    TriangleDouble B_d = triB_double[idxB];

    return classifyPairFP64(A_d, B_d);
}

// ====================================================================
// 5. CORE INTERSECTION KERNEL & EXPORTED HOST STAGE
// ====================================================================
__global__ void evaluateGeometricPairsKernelCascaded(
    int *outStatuses, 
    const int2 *candidatePairs, 
    const cuBQL::Triangle *triA, 
    const cuBQL::Triangle *triB, 
    const float2 *triAMetrics, 
    const float2 *triBMetrics, 
    const TriangleDouble *triA_double,
    const TriangleDouble *triB_double,
    int numPairs) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= numPairs) return;

    int2 pair = candidatePairs[tid];

    PairStatus status = classifyPairCascaded(
        triA[pair.x], 
        triB[pair.y], 
        triAMetrics[pair.x], 
        triBMetrics[pair.y],
        triA_double,
        triB_double,
        pair.x,
        pair.y
    );

    outStatuses[tid] = (int)status;
}

void evaluateAndCompactPairsV2(
    int2* dCandidatePairs,
    int* dPairStatuses,
    const cuBQL::Triangle* dA,
    const cuBQL::Triangle* dB_batch,
    const float2 *triAMetrics,
    const float2 *triBMetrics,
    const TriangleDouble *dA_double,
    const TriangleDouble *dB_double,
    int totalBatchPairs,
    cudaStream_t stream)
{
    if (totalBatchPairs <= 0) return;

    int threadsPerBlock = 256;
    int blocksPerGrid = (totalBatchPairs + threadsPerBlock - 1) / threadsPerBlock;

    evaluateGeometricPairsKernelCascaded<<<blocksPerGrid, threadsPerBlock, 0, stream>>>(
         dPairStatuses, 
         dCandidatePairs, 
         dA, 
         dB_batch, 
         triAMetrics, 
         triBMetrics, 
         dA_double, 
         dB_double, 
         totalBatchPairs
    );
}