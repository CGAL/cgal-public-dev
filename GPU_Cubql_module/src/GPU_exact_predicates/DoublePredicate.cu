#include "DoublePredicate.h"
#include "GPUPredicatesCommon.h"
#include <cmath>

// --------------------------------------------------------------------
// MEMORY HELPER (Fixes __ldg overload issue for 24-byte double3)
// --------------------------------------------------------------------
__device__ inline double3 load_double3(const double3* ptr) {
    return make_double3(__ldg(&(ptr->x)), __ldg(&(ptr->y)), __ldg(&(ptr->z)));
}

// --------------------------------------------------------------------
// FP64 DOUBLE PREDICATES WITH DYNAMIC ERROR BOUNDING
// --------------------------------------------------------------------
__device__ inline int orient3d_double(
    const double3& p, const double3& q, 
    const double3& r, const double3& s, double double_eps) 
{
    const double pdx = p.x - s.x; const double pdy = p.y - s.y; const double pdz = p.z - s.z;
    const double qdx = q.x - s.x; const double qdy = q.y - s.y; const double qdz = q.z - s.z;
    const double rdx = r.x - s.x; const double rdy = r.y - s.y; const double rdz = r.z - s.z;

    const double det = 
        pdx * (qdy * rdz - qdz * rdy) +
        pdy * (qdz * rdx - qdx * rdz) +
        pdz * (qdx * rdy - qdy * rdx);

    if (det > double_eps) return 1;
    if (det < -double_eps) return -1;
    return 0; // Fallback to Exact CPU Predicate (PAIR_YELLOW)
}

__device__ inline int edgeTri_double(
    const double3 &a, const double3 &b,
    const double3 &p, const double3 &q, const double3 &r, double double_eps) 
{
    int s0 = orient3d_double(p, q, r, a, double_eps);
    int s1 = orient3d_double(p, q, r, b, double_eps);

    if (s0 == 0 || s1 == 0) return PAIR_YELLOW;
    if ((s0 > 0 && s1 > 0) || (s0 < 0 && s1 < 0)) return PAIR_NO;

    int e0 = orient3d_double(a, b, p, q, double_eps);
    int e1 = orient3d_double(a, b, q, r, double_eps);
    int e2 = orient3d_double(a, b, r, p, double_eps);

    if (e0 == 0 || e1 == 0 || e2 == 0) return PAIR_YELLOW;
    if ((e0 >= 0 && e1 >= 0 && e2 >= 0) || (e0 <= 0 && e1 <= 0 && e2 <= 0)) return PAIR_GREEN;

    return PAIR_NO;
}

// --------------------------------------------------------------------
// CLASSIFICATION FUNCTION (Hardware FP64 Precision)
// --------------------------------------------------------------------
__device__ inline PairStatus classifyPairDoubleFull(
    const double3& Aa, const double3& Ab, const double3& Ac,
    const double3& Ba, const double3& Bb, const double3& Bc) 
{    
    // 1. Exact Double Anchor Shift (Eliminates World-Space Offset Cancellation)
    const double3 anchor = Aa;
    
    const double3 relA0 = { 0.0, 0.0, 0.0 };
    const double3 relA1 = { Ab.x - anchor.x, Ab.y - anchor.y, Ab.z - anchor.z };
    const double3 relA2 = { Ac.x - anchor.x, Ac.y - anchor.y, Ac.z - anchor.z };

    const double3 relB0 = { Ba.x - anchor.x, Ba.y - anchor.y, Ba.z - anchor.z };
    const double3 relB1 = { Bb.x - anchor.x, Bb.y - anchor.y, Bb.z - anchor.z };
    const double3 relB2 = { Bc.x - anchor.x, Bc.y - anchor.y, Bc.z - anchor.z };

    // 2. Compute Local Bounding Scale for Machine Epsilon Scaling
    double maxA1 = fmax(fabs(relA1.x), fmax(fabs(relA1.y), fabs(relA1.z)));
    double maxA2 = fmax(fabs(relA2.x), fmax(fabs(relA2.y), fabs(relA2.z)));
    double maxB0 = fmax(fabs(relB0.x), fmax(fabs(relB0.y), fabs(relB0.z)));
    double maxB1 = fmax(fabs(relB1.x), fmax(fabs(relB1.y), fabs(relB1.z)));
    double maxB2 = fmax(fabs(relB2.x), fmax(fabs(relB2.y), fabs(relB2.z)));

    double M_local = fmax(maxA1, fmax(maxA2, fmax(maxB0, fmax(maxB1, maxB2))));

    // 3. Dynamic FP64 Error Bound: Scale^3 * 64 * 2^-53 (64 * 1.110223e-16 = 7.105427e-15)
    const double double_eps = M_local * M_local * M_local * 7.105427357601002e-15;

    // 4. Early Rejection Checks
    int ob0 = orient3d_double(relA0, relA1, relA2, relB0, double_eps);
    int ob1 = orient3d_double(relA0, relA1, relA2, relB1, double_eps);
    int ob2 = orient3d_double(relA0, relA1, relA2, relB2, double_eps);
    if (ob0 != 0 && ob1 != 0 && ob2 != 0) {
        if ((ob0 > 0 && ob1 > 0 && ob2 > 0) || (ob0 < 0 && ob1 < 0 && ob2 < 0)) return PAIR_NO;
    }

    int oa0 = orient3d_double(relB0, relB1, relB2, relA0, double_eps);
    int oa1 = orient3d_double(relB0, relB1, relB2, relA1, double_eps);
    int oa2 = orient3d_double(relB0, relB1, relB2, relA2, double_eps);
    if (oa0 != 0 && oa1 != 0 && oa2 != 0) {
        if ((oa0 > 0 && oa1 > 0 && oa2 > 0) || (oa0 < 0 && oa1 < 0 && oa2 < 0)) return PAIR_NO;
    }

    // 5. Detailed Edge-Triangle Predicates
    int r;
    r = edgeTri_double(relA0, relA1, relB0, relB1, relB2, double_eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;
    r = edgeTri_double(relA1, relA2, relB0, relB1, relB2, double_eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;
    r = edgeTri_double(relA2, relA0, relB0, relB1, relB2, double_eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;

    r = edgeTri_double(relB0, relB1, relA0, relA1, relA2, double_eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;
    r = edgeTri_double(relB1, relB2, relA0, relA1, relA2, double_eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;
    r = edgeTri_double(relB2, relB0, relA0, relA1, relA2, double_eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;

    return PAIR_NO;
}

// --------------------------------------------------------------------
// CORE INTERSECTION KERNEL
// --------------------------------------------------------------------
__global__ void evaluateGeometricPairsKernelDoubleFull(
    int* dPairStatuses, 
    const int2* dYellowCandidatePairs, 
    const double3* dVertsA,
    const uint3* dIndicesA,
    const double3* dVertsB,
    const uint3* dIndicesB,
    int numYellowPairs) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= numYellowPairs) return;

    int2 pair = dYellowCandidatePairs[tid];
    
    uint3 idxA = dIndicesA[pair.x];
    uint3 idxB = dIndicesB[pair.y];

    // Load double3 using component-wise __ldg loads
    double3 Aa = load_double3(&dVertsA[idxA.x]);
    double3 Ab = load_double3(&dVertsA[idxA.y]);
    double3 Ac = load_double3(&dVertsA[idxA.z]);

    double3 Ba = load_double3(&dVertsB[idxB.x]);
    double3 Bb = load_double3(&dVertsB[idxB.y]);
    double3 Bc = load_double3(&dVertsB[idxB.z]);

    PairStatus status = classifyPairDoubleFull(Aa, Ab, Ac, Ba, Bb, Bc);

    dPairStatuses[tid] = (int)status; 
}

// --------------------------------------------------------------------
// PIPELINE EXPORTED WRAPPER
// --------------------------------------------------------------------
void evaluateAndCompactPairsDoubleFull(
    const int2* dYellowCandidatePairs,
    int* dPairStatuses,
    const double3* dVertsA,
    const uint3* dIndicesA,
    const double3* dVertsB,
    const uint3* dIndicesB,
    int numYellowPairs,
    cudaStream_t stream)
{
    if (numYellowPairs <= 0) return;

    int threadsPerBlock = 256;
    int blocksPerGrid = (numYellowPairs + threadsPerBlock - 1) / threadsPerBlock;
    
    evaluateGeometricPairsKernelDoubleFull<<<blocksPerGrid, threadsPerBlock, 0, stream>>>(
         dPairStatuses, dYellowCandidatePairs, 
         dVertsA, dIndicesA, 
         dVertsB, dIndicesB, 
         numYellowPairs
    );
}