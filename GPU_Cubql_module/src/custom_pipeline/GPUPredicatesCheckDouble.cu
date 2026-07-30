#include "GPUPredicatesCheckDouble.h"
#include <cuda_runtime.h>
#include <thrust/device_ptr.h>
#include <thrust/copy.h>
#include <thrust/execution_policy.h>
#include <thrust/device_vector.h>
#include <algorithm>
#include <cstdio>

#ifndef CUBQL_CUDA_CALL
#define CUBQL_CUDA_CALL(X) { cudaError_t err = cuda ## X; if(err != cudaSuccess) { printf("CUDA Error at %s:%d\n", __FILE__, __LINE__); } }
#endif

// --------------------------------------------------------------------
// DOUBLE INTERVAL ARITHMETIC (DEVICE INTERNAL)
// --------------------------------------------------------------------
struct GPUIntervalDouble {
    double lo;
    double hi;
};

__device__ inline GPUIntervalDouble interval_sub(double a, double b) { 
    return { __dsub_rd(a, b), __dsub_ru(a, b) }; 
}

__device__ inline GPUIntervalDouble interval_add(GPUIntervalDouble a, GPUIntervalDouble b) { 
    return { __dadd_rd(a.lo, b.lo), __dadd_ru(a.hi, b.hi) }; 
}

__device__ inline GPUIntervalDouble interval_sub(GPUIntervalDouble a, GPUIntervalDouble b) { 
    return { __dsub_rd(a.lo, b.hi), __dsub_ru(a.hi, b.lo) }; 
}

__device__ inline GPUIntervalDouble interval_mul(GPUIntervalDouble a, GPUIntervalDouble b) {
    double lo1 = __dmul_rd(a.lo, b.lo); double lo2 = __dmul_rd(a.lo, b.hi);
    double lo3 = __dmul_rd(a.hi, b.lo); double lo4 = __dmul_rd(a.hi, b.hi);
    double hi1 = __dmul_ru(a.lo, b.lo); double hi2 = __dmul_ru(a.lo, b.hi);
    double hi3 = __dmul_ru(a.hi, b.lo); double hi4 = __dmul_ru(a.hi, b.hi);
    return { fmin(fmin(lo1, lo2), fmin(lo3, lo4)), fmax(fmax(hi1, hi2), fmax(hi3, hi4)) };
}

// Double-precision Orient3D with strict interval zero-exclusion check
__device__ inline int orient3d_interval_double(
    const double3& p, const double3& q, 
    const double3& r, const double3& s) 
{
    GPUIntervalDouble pdx = interval_sub(p.x, s.x); 
    GPUIntervalDouble pdy = interval_sub(p.y, s.y); 
    GPUIntervalDouble pdz = interval_sub(p.z, s.z);

    GPUIntervalDouble qdx = interval_sub(q.x, s.x); 
    GPUIntervalDouble qdy = interval_sub(q.y, s.y); 
    GPUIntervalDouble qdz = interval_sub(q.z, s.z);

    GPUIntervalDouble rdx = interval_sub(r.x, s.x); 
    GPUIntervalDouble rdy = interval_sub(r.y, s.y); 
    GPUIntervalDouble rdz = interval_sub(r.z, s.z);

    GPUIntervalDouble det = interval_add(
        interval_add(
            interval_mul(pdx, interval_sub(interval_mul(qdy, rdz), interval_mul(qdz, rdy))),
            interval_mul(pdy, interval_sub(interval_mul(qdz, rdx), interval_mul(qdx, rdz)))
        ),
        interval_mul(pdz, interval_sub(interval_mul(qdx, rdy), interval_mul(qdy, rdx)))
    );

    // If 0 is strictly outside [lo, hi], exact sign is guaranteed
    if (det.lo > 0.0) return 1;
    if (det.hi < 0.0) return -1;
    
    // Determinant interval contains 0 -> defer to CPU CGAL exact predicate
    return 0; 
}

__device__ inline int edgeTriIntervalDouble(
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

__device__ inline PairStatus classifyPairDouble(
    const TriangleDouble& A, 
    const TriangleDouble& B) 
{    
    const double3 A0 = A.a, A1 = A.b, A2 = A.c;
    const double3 B0 = B.a, B1 = B.b, B2 = B.c;

    // Reject early if vertices are strictly separated by plane orientation
    int ob0 = orient3d_interval_double(A0, A1, A2, B0);
    int ob1 = orient3d_interval_double(A0, A1, A2, B1);
    int ob2 = orient3d_interval_double(A0, A1, A2, B2);
    if (ob0 != 0 && ob1 != 0 && ob2 != 0) {
        if ((ob0 > 0 && ob1 > 0 && ob2 > 0) || (ob0 < 0 && ob1 < 0 && ob2 < 0)) return PAIR_NO;
    }

    int oa0 = orient3d_interval_double(B0, B1, B2, A0);
    int oa1 = orient3d_interval_double(B0, B1, B2, A1);
    int oa2 = orient3d_interval_double(B0, B1, B2, A2);
    if (oa0 != 0 && oa1 != 0 && oa2 != 0) {
        if ((oa0 > 0 && oa1 > 0 && oa2 > 0) || (oa0 < 0 && oa1 < 0 && oa2 < 0)) return PAIR_NO;
    }

    // Edge-triangle intersection passes
    int r;
    r = edgeTriIntervalDouble(A0, A1, B0, B1, B2); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;
    r = edgeTriIntervalDouble(A1, A2, B0, B1, B2); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;
    r = edgeTriIntervalDouble(A2, A0, B0, B1, B2); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;

    r = edgeTriIntervalDouble(B0, B1, A0, A1, A2); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;
    r = edgeTriIntervalDouble(B1, B2, A0, A1, A2); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;
    r = edgeTriIntervalDouble(B2, B0, A0, A1, A2); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;

    return PAIR_NO;
}

// --------------------------------------------------------------------
// CORE INTERSECTION KERNEL
// --------------------------------------------------------------------
__global__ void evaluateGeometricPairsKernelDouble(
    int *outStatuses, 
    const int2 *candidatePairs, 
    const TriangleDouble *triA, 
    const TriangleDouble *triB, 
    int numPairs) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= numPairs) return;

    int2 pair = candidatePairs[tid];
    
    PairStatus status = classifyPairDouble(
        triA[pair.x], 
        triB[pair.y]
    );

    outStatuses[tid] = (int)status;
}

// --------------------------------------------------------------------
// EXPORTED PIPELINE WRAPPER STAGE
// --------------------------------------------------------------------
void evaluateAndCompactPairsDouble(
    int2* dCandidatePairs,
    int* dPairStatuses,
    const TriangleDouble* dA,
    const TriangleDouble* dB_batch,
    int totalBatchPairs,
    cudaStream_t stream)
{
    if (totalBatchPairs <= 0) return;

    int threadsPerBlock = 256;
    int blocksPerGrid = (totalBatchPairs + threadsPerBlock - 1) / threadsPerBlock;
    
    evaluateGeometricPairsKernelDouble<<<blocksPerGrid, threadsPerBlock, 0, stream>>>(
         dPairStatuses, dCandidatePairs, dA, dB_batch, totalBatchPairs
    );
}