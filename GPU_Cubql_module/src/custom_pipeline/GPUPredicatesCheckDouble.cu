#include "GPUPredicatesCheckDouble.h"
#include "GPUPredicatesCommon.h"
#include <cmath>

// --------------------------------------------------------------------
// MEMORY HELPER (Fixes __ldg overload error for 24-byte double3)
// --------------------------------------------------------------------
__device__ inline double3 load_double3(const double3* ptr) {
    return make_double3(__ldg(&(ptr->x)), __ldg(&(ptr->y)), __ldg(&(ptr->z)));
}

// --------------------------------------------------------------------
// DOUBLE-SINGLE (DS) SIMULATION ARITHMETIC (1:1 FP32 Throughput)
// --------------------------------------------------------------------
struct ds_float {
    float hi;
    float lo;
};

struct ds_vec3 {
    ds_float x, y, z;
};

// Converts native double (53-bit mantissa) into two exact FP32 components
__device__ inline ds_float double_to_ds(double val) {
    float hi = static_cast<float>(val);
    float lo = static_cast<float>(val - static_cast<double>(hi));
    return { hi, lo };
}

__device__ inline ds_vec3 to_ds_vec(double3 v) {
    return { double_to_ds(v.x), double_to_ds(v.y), double_to_ds(v.z) };
}

__device__ inline ds_float quick_two_sum(float a, float b) {
    float s = __fadd_rn(a, b);
    float e = __fsub_rn(b, __fsub_rn(s, a));
    return { s, e };
}

__device__ inline ds_float two_sum(float a, float b) {
    float s = __fadd_rn(a, b);
    float v = __fsub_rn(s, a);
    float e = __fadd_rn(__fsub_rn(a, __fsub_rn(s, v)), __fsub_rn(b, v));
    return { s, e };
}

__device__ inline ds_float two_prod(float a, float b) {
    float p = __fmul_rn(a, b);
    float e = __fmaf_rn(a, b, -p); 
    return { p, e };
}

__device__ inline ds_float ds_add(ds_float a, ds_float b) {
    ds_float s = two_sum(a.hi, b.hi);
    ds_float t = two_sum(a.lo, b.lo);
    s.lo = __fadd_rn(s.lo, t.hi);
    s = quick_two_sum(s.hi, s.lo);
    s.lo = __fadd_rn(s.lo, t.lo);
    return quick_two_sum(s.hi, s.lo);
}

__device__ inline ds_float ds_sub(ds_float a, ds_float b) {
    return ds_add(a, { -b.hi, -b.lo });
}

__device__ inline ds_float ds_mul(ds_float a, ds_float b) {
    ds_float p = two_prod(a.hi, b.hi);
    p.lo = __fadd_rn(p.lo, __fadd_rn(__fmul_rn(a.hi, b.lo), __fmul_rn(a.lo, b.hi)));
    return quick_two_sum(p.hi, p.lo);
}

// --------------------------------------------------------------------
// PRECISION DS PREDICATES WITH DYNAMIC ERROR BOUNDING
// --------------------------------------------------------------------
__device__ inline int orient3d_ds(
    const ds_vec3& p, const ds_vec3& q, 
    const ds_vec3& r, const ds_vec3& s, float ds_eps) 
{
    ds_float pdx = ds_sub(p.x, s.x); ds_float pdy = ds_sub(p.y, s.y); ds_float pdz = ds_sub(p.z, s.z);
    ds_float qdx = ds_sub(q.x, s.x); ds_float qdy = ds_sub(q.y, s.y); ds_float qdz = ds_sub(q.z, s.z);
    ds_float rdx = ds_sub(r.x, s.x); ds_float rdy = ds_sub(r.y, s.y); ds_float rdz = ds_sub(r.z, s.z);

    ds_float det = ds_add(
        ds_add(
            ds_mul(pdx, ds_sub(ds_mul(qdy, rdz), ds_mul(qdz, rdy))),
            ds_mul(pdy, ds_sub(ds_mul(qdz, rdx), ds_mul(qdx, rdz)))
        ),
        ds_mul(pdz, ds_sub(ds_mul(qdx, rdy), ds_mul(qdy, rdx)))
    );

    if (det.hi > ds_eps) return 1;
    if (det.hi < -ds_eps) return -1;
    return 0; // Fallback to Exact CPU Predicate (PAIR_YELLOW)
}

__device__ inline int edgeTri_ds(
    const ds_vec3 &a, const ds_vec3 &b,
    const ds_vec3 &p, const ds_vec3 &q, const ds_vec3 &r, float ds_eps) 
{
    int s0 = orient3d_ds(p, q, r, a, ds_eps);
    int s1 = orient3d_ds(p, q, r, b, ds_eps);

    if (s0 == 0 || s1 == 0) return PAIR_YELLOW;
    if ((s0 > 0 && s1 > 0) || (s0 < 0 && s1 < 0)) return PAIR_NO;

    int e0 = orient3d_ds(a, b, p, q, ds_eps);
    int e1 = orient3d_ds(a, b, q, r, ds_eps);
    int e2 = orient3d_ds(a, b, r, p, ds_eps);

    if (e0 == 0 || e1 == 0 || e2 == 0) return PAIR_YELLOW;
    if ((e0 >= 0 && e1 >= 0 && e2 >= 0) || (e0 <= 0 && e1 <= 0 && e2 <= 0)) return PAIR_GREEN;

    return PAIR_NO;
}

// --------------------------------------------------------------------
// CLASSIFICATION FUNCTION
// --------------------------------------------------------------------
__device__ inline PairStatus classifyPairDS(
    const double3& Aa, const double3& Ab, const double3& Ac,
    const double3& Ba, const double3& Bb, const double3& Bc) 
{    
    // 1. Exact Double Anchor Shift (Eliminates World-Space Offset Cancellation)
    const double3 anchor = Aa;
    
    double3 relA1 = { Ab.x - anchor.x, Ab.y - anchor.y, Ab.z - anchor.z };
    double3 relA2 = { Ac.x - anchor.x, Ac.y - anchor.y, Ac.z - anchor.z };

    double3 relB0 = { Ba.x - anchor.x, Ba.y - anchor.y, Ba.z - anchor.z };
    double3 relB1 = { Bb.x - anchor.x, Bb.y - anchor.y, Bb.z - anchor.z };
    double3 relB2 = { Bc.x - anchor.x, Bc.y - anchor.y, Bc.z - anchor.z };

    // 2. Compute Local Bounding Scale for Strict Machine Epsilon
    double maxA1 = fmax(fabs(relA1.x), fmax(fabs(relA1.y), fabs(relA1.z)));
    double maxA2 = fmax(fabs(relA2.x), fmax(fabs(relA2.y), fabs(relA2.z)));
    double maxB0 = fmax(fabs(relB0.x), fmax(fabs(relB0.y), fabs(relB0.z)));
    double maxB1 = fmax(fabs(relB1.x), fmax(fabs(relB1.y), fabs(relB1.z)));
    double maxB2 = fmax(fabs(relB2.x), fmax(fabs(relB2.y), fabs(relB2.z)));

    double M_local = fmax(maxA1, fmax(maxA2, fmax(maxB0, fmax(maxB1, maxB2))));

    // 3. Dynamic DS Error Bound: Scale^3 * 64 * 2^-48
    const float ds_eps = static_cast<float>(M_local * M_local * M_local * 2.27e-13);

    // 4. Promote Relative Coordinates to Double-Single Vectors
    ds_vec3 ds_A0 = { {0.0f, 0.0f}, {0.0f, 0.0f}, {0.0f, 0.0f} };
    ds_vec3 ds_A1 = to_ds_vec(relA1);
    ds_vec3 ds_A2 = to_ds_vec(relA2);

    ds_vec3 ds_B0 = to_ds_vec(relB0);
    ds_vec3 ds_B1 = to_ds_vec(relB1);
    ds_vec3 ds_B2 = to_ds_vec(relB2);

    // 5. Early Rejection Checks
    int ob0 = orient3d_ds(ds_A0, ds_A1, ds_A2, ds_B0, ds_eps);
    int ob1 = orient3d_ds(ds_A0, ds_A1, ds_A2, ds_B1, ds_eps);
    int ob2 = orient3d_ds(ds_A0, ds_A1, ds_A2, ds_B2, ds_eps);
    if (ob0 != 0 && ob1 != 0 && ob2 != 0) {
        if ((ob0 > 0 && ob1 > 0 && ob2 > 0) || (ob0 < 0 && ob1 < 0 && ob2 < 0)) return PAIR_NO;
    }

    int oa0 = orient3d_ds(ds_B0, ds_B1, ds_B2, ds_A0, ds_eps);
    int oa1 = orient3d_ds(ds_B0, ds_B1, ds_B2, ds_A1, ds_eps);
    int oa2 = orient3d_ds(ds_B0, ds_B1, ds_B2, ds_A2, ds_eps);
    if (oa0 != 0 && oa1 != 0 && oa2 != 0) {
        if ((oa0 > 0 && oa1 > 0 && oa2 > 0) || (oa0 < 0 && oa1 < 0 && oa2 < 0)) return PAIR_NO;
    }

    // 6. Detailed Edge-Triangle Predicates
    int r;
    r = edgeTri_ds(ds_A0, ds_A1, ds_B0, ds_B1, ds_B2, ds_eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;
    r = edgeTri_ds(ds_A1, ds_A2, ds_B0, ds_B1, ds_B2, ds_eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;
    r = edgeTri_ds(ds_A2, ds_A0, ds_B0, ds_B1, ds_B2, ds_eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;

    r = edgeTri_ds(ds_B0, ds_B1, ds_A0, ds_A1, ds_A2, ds_eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;
    r = edgeTri_ds(ds_B1, ds_B2, ds_A0, ds_A1, ds_A2, ds_eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;
    r = edgeTri_ds(ds_B2, ds_B0, ds_A0, ds_A1, ds_A2, ds_eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;

    return PAIR_NO;
}

// --------------------------------------------------------------------
// CORE INTERSECTION KERNEL
// --------------------------------------------------------------------
__global__ void evaluateGeometricPairsKernelDouble(
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

    PairStatus status = classifyPairDS(Aa, Ab, Ac, Ba, Bb, Bc);

    dPairStatuses[tid] = (int)status; 
}

// --------------------------------------------------------------------
// PIPELINE EXPORTED WRAPPER
// --------------------------------------------------------------------
void evaluateAndCompactPairsDouble(
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
    
    evaluateGeometricPairsKernelDouble<<<blocksPerGrid, threadsPerBlock, 0, stream>>>(
         dPairStatuses, dYellowCandidatePairs, 
         dVertsA, dIndicesA, 
         dVertsB, dIndicesB, 
         numYellowPairs
    );
}