#include "GPUPredicatesTwoLap.h"
#include "GPUPredicatesCommon.h"
#include <cmath>

// --------------------------------------------------------------------
// MEMORY HELPER
// --------------------------------------------------------------------
__device__ inline double3 load_double3(const double3* ptr) {
    return make_double3(__ldg(&(ptr->x)), __ldg(&(ptr->y)), __ldg(&(ptr->z)));
}

// --------------------------------------------------------------------
// DOUBLE-SINGLE (DS) REGISTER ARITHMETIC PRIMITIVES
// --------------------------------------------------------------------
struct ds_float {
    float hi;
    float lo;
};

struct ds_vec3 {
    ds_float x, y, z;
};

struct ds_vec2 {
    ds_float x, y;
};

__device__ inline ds_float double_to_ds(double val) {
    float hi = static_cast<float>(val);
    float lo = static_cast<float>(val - static_cast<double>(hi));
    return { hi, lo };
}

__device__ inline ds_vec3 to_ds_vec3(double3 v) {
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
// 2D AND 3D REGISTER ORIENTATION PREDICATES
// --------------------------------------------------------------------

// LAP 1 PRIMITIVE: Evaluates 3D 4-point determinant in registers
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
    
    // If within [-ds_eps, ds_eps], check if high part is strictly zero
    if (det.hi == 0.0f && det.lo == 0.0f) return 0; // True Exact Zero

    return 2; // Ambiguous non-zero in (0, ds_eps] range -> Route to PAIR_YELLOW
}

// LAP 2 PRIMITIVE: Evaluates 2D 3-point determinant in registers
__device__ inline int orient2d_ds(
    const ds_vec2& a, const ds_vec2& b, const ds_vec2& c, float ds_eps_2d) 
{
    ds_float acx = ds_sub(a.x, c.x); ds_float acy = ds_sub(a.y, c.y);
    ds_float bcx = ds_sub(b.x, c.x); ds_float bcy = ds_sub(b.y, c.y);

    ds_float det = ds_sub(ds_mul(acx, bcy), ds_mul(acy, bcx));

    if (det.hi > ds_eps_2d) return 1;
    if (det.hi < -ds_eps_2d) return -1;
    if (det.hi == 0.0f && det.lo == 0.0f) return 0;

    return 2; // Ambiguous sub-grid -> Route to PAIR_YELLOW
}

// Project 3D Double-Single vector to 2D by dropping dominant normal axis
__device__ inline ds_vec2 project2D(const ds_vec3& v, int drop_axis) {
    if (drop_axis == 0) return { v.y, v.z }; // Drop X
    if (drop_axis == 1) return { v.x, v.z }; // Drop Y
    return { v.x, v.y };                    // Drop Z
}

// --------------------------------------------------------------------
// LAP 2: COPLANAR 2D SEPARATION AXIS TEST (SAT)
// --------------------------------------------------------------------
__device__ inline PairStatus evaluateCoplanar2D(
    const ds_vec3& A0, const ds_vec3& A1, const ds_vec3& A2,
    const ds_vec3& B0, const ds_vec3& B1, const ds_vec3& B2,
    float ds_eps) 
{
    // 1. Calculate approximate normal of Triangle A to find dominant axis
    ds_float e1x = ds_sub(A1.x, A0.x); ds_float e1y = ds_sub(A1.y, A0.y); ds_float e1z = ds_sub(A1.z, A0.z);
    ds_float e2x = ds_sub(A2.x, A0.x); ds_float e2y = ds_sub(A2.y, A0.y); ds_float e2z = ds_sub(A2.z, A0.z);

    float nx = fabsf(ds_sub(ds_mul(e1y, e2z), ds_mul(e1z, e2y)).hi);
    float ny = fabsf(ds_sub(ds_mul(e1z, e2x), ds_mul(e1x, e2z)).hi);
    float nz = fabsf(ds_sub(ds_mul(e1x, e2y), ds_mul(e1y, e2x)).hi);

    int drop_axis = 2; // Default drop Z
    if (nx > ny && nx > nz) drop_axis = 0;
    else if (ny > nz) drop_axis = 1;

    // 2. Project points to 2D
    ds_vec2 a0 = project2D(A0, drop_axis);
    ds_vec2 a1 = project2D(A1, drop_axis);
    ds_vec2 a2 = project2D(A2, drop_axis);

    ds_vec2 b0 = project2D(B0, drop_axis);
    ds_vec2 b1 = project2D(B1, drop_axis);
    ds_vec2 b2 = project2D(B2, drop_axis);

    float ds_eps_2d = ds_eps; // Scaled 2D threshold

    // 3. Test edges of Triangle A against vertices of B using orient2d
    const ds_vec2 A_verts[3] = { a0, a1, a2 };
    const ds_vec2 B_verts[3] = { b0, b1, b2 };

    for (int i = 0; i < 3; ++i) {
        ds_vec2 eA = A_verts[i];
        ds_vec2 eB = A_verts[(i + 1) % 3];

        int o0 = orient2d_ds(eA, eB, b0, ds_eps_2d);
        int o1 = orient2d_ds(eA, eB, b1, ds_eps_2d);
        int o2 = orient2d_ds(eA, eB, b2, ds_eps_2d);

        if (o0 == 2 || o1 == 2 || o2 == 2) return PAIR_YELLOW; // Sub-eps float miss

        // If all vertices of B are strictly outside this edge, triangles are disjoint
        if ((o0 > 0 && o1 > 0 && o2 > 0) || (o0 < 0 && o1 < 0 && o2 < 0)) {
            return PAIR_NO;
        }
    }

    // 4. Test edges of Triangle B against vertices of A
    for (int i = 0; i < 3; ++i) {
        ds_vec2 eA = B_verts[i];
        ds_vec2 eB = B_verts[(i + 1) % 3];

        int o0 = orient2d_ds(eA, eB, a0, ds_eps_2d);
        int o1 = orient2d_ds(eA, eB, a1, ds_eps_2d);
        int o2 = orient2d_ds(eA, eB, a2, ds_eps_2d);

        if (o0 == 2 || o1 == 2 || o2 == 2) return PAIR_YELLOW;

        if ((o0 > 0 && o1 > 0 && o2 > 0) || (o0 < 0 && o1 < 0 && o2 < 0)) {
            return PAIR_NO;
        }
    }

    return PAIR_GREEN; // Exact 2D intersection confirmed
}

// --------------------------------------------------------------------
// TWO-LAP CLASSIFICATION PIPELINE
// --------------------------------------------------------------------
__device__ inline PairStatus classifyPairTwoLap(
    const double3& Aa, const double3& Ab, const double3& Ac,
    const double3& Ba, const double3& Bb, const double3& Bc) 
{    
    // Step 1: Anchor shift & local scale bounds
    const double3 anchor = Aa;
    
    double3 relA1 = { Ab.x - anchor.x, Ab.y - anchor.y, Ab.z - anchor.z };
    double3 relA2 = { Ac.x - anchor.x, Ac.y - anchor.y, Ac.z - anchor.z };

    double3 relB0 = { Ba.x - anchor.x, Ba.y - anchor.y, Ba.z - anchor.z };
    double3 relB1 = { Bb.x - anchor.x, Bb.y - anchor.y, Bb.z - anchor.z };
    double3 relB2 = { Bc.x - anchor.x, Bc.y - anchor.y, Bc.z - anchor.z };

    double maxA1 = fmax(fabs(relA1.x), fmax(fabs(relA1.y), fabs(relA1.z)));
    double maxA2 = fmax(fabs(relA2.x), fmax(fabs(relA2.y), fabs(relA2.z)));
    double maxB0 = fmax(fabs(relB0.x), fmax(fabs(relB0.y), fabs(relB0.z)));
    double maxB1 = fmax(fabs(relB1.x), fmax(fabs(relB1.y), fabs(relB1.z)));
    double maxB2 = fmax(fabs(relB2.x), fmax(fabs(relB2.y), fabs(relB2.z)));

    double M_local = fmax(maxA1, fmax(maxA2, fmax(maxB0, fmax(maxB1, maxB2))));
    const float ds_eps = static_cast<float>(M_local * M_local * M_local * 2.27e-13);

    ds_vec3 ds_A0 = { {0.0f, 0.0f}, {0.0f, 0.0f}, {0.0f, 0.0f} };
    ds_vec3 ds_A1 = to_ds_vec3(relA1);
    ds_vec3 ds_A2 = to_ds_vec3(relA2);

    ds_vec3 ds_B0 = to_ds_vec3(relB0);
    ds_vec3 ds_B1 = to_ds_vec3(relB1);
    ds_vec3 ds_B2 = to_ds_vec3(relB2);

    // ================================================================
    // LAP 1: FIND ALL EXACT 0-DETERMINANT CASES
    // ================================================================
    int ob0 = orient3d_ds(ds_A0, ds_A1, ds_A2, ds_B0, ds_eps);
    int ob1 = orient3d_ds(ds_A0, ds_A1, ds_A2, ds_B1, ds_eps);
    int ob2 = orient3d_ds(ds_A0, ds_A1, ds_A2, ds_B2, ds_eps);

    // If any orientation is in the ambiguous (0, eps] noise band, push to PAIR_YELLOW
    if (ob0 == 2 || ob1 == 2 || ob2 == 2) {
        return PAIR_YELLOW;
    }

    // QUALIFY LAP 1: Are all 3 vertices strictly coplanar (determinants == 0)?
    if (ob0 == 0 && ob1 == 0 && ob2 == 0) {
        // ============================================================
        // LAP 2: EXECUTE 2D ORIENT2D SAT PROJECTION FOR COPLANAR PAIR
        // ============================================================
        return evaluateCoplanar2D(ds_A0, ds_A1, ds_A2, ds_B0, ds_B1, ds_B2, ds_eps);
    }

    // Standard 3D Early Rejection (all vertices on one side of plane)
    if ((ob0 > 0 && ob1 > 0 && ob2 > 0) || (ob0 < 0 && ob1 < 0 && ob2 < 0)) {
        return PAIR_NO;
    }

    // Check converse orientations (Triangle A against Plane B)
    int oa0 = orient3d_ds(ds_B0, ds_B1, ds_B2, ds_A0, ds_eps);
    int oa1 = orient3d_ds(ds_B0, ds_B1, ds_B2, ds_A1, ds_eps);
    int oa2 = orient3d_ds(ds_B0, ds_B1, ds_B2, ds_A2, ds_eps);

    if (oa0 == 2 || oa1 == 2 || oa2 == 2) {
        return PAIR_YELLOW;
    }

    if (oa0 == 0 && oa1 == 0 && oa2 == 0) {
        return evaluateCoplanar2D(ds_A0, ds_A1, ds_A2, ds_B0, ds_B1, ds_B2, ds_eps);
    }

    if ((oa0 > 0 && oa1 > 0 && oa2 > 0) || (oa0 < 0 && oa1 < 0 && oa2 < 0)) {
        return PAIR_NO;
    }

    // Remaining non-coplanar 3D straddling cases -> Send to CPU exact solver
    return PAIR_YELLOW;
}

// --------------------------------------------------------------------
// CUDA KERNEL & PIPELINE EXPORT
// --------------------------------------------------------------------
__global__ void evaluateTwoLapPairsKernel(
    int* dPairStatuses, 
    const int2* dCandidatePairs, 
    const double3* dVertsA,
    const uint3* dIndicesA,
    const double3* dVertsB,
    const uint3* dIndicesB,
    int numPairs) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= numPairs) return;

    int2 pair = dCandidatePairs[tid];
    
    uint3 idxA = dIndicesA[pair.x];
    uint3 idxB = dIndicesB[pair.y];

    double3 Aa = load_double3(&dVertsA[idxA.x]);
    double3 Ab = load_double3(&dVertsA[idxA.y]);
    double3 Ac = load_double3(&dVertsA[idxA.z]);

    double3 Ba = load_double3(&dVertsB[idxB.x]);
    double3 Bb = load_double3(&dVertsB[idxB.y]);
    double3 Bc = load_double3(&dVertsB[idxB.z]);

    PairStatus status = classifyPairTwoLap(Aa, Ab, Ac, Ba, Bb, Bc);

    dPairStatuses[tid] = (int)status; 
}

void evaluateTwoLapPairs(
    const int2* dCandidatePairs,
    int* dPairStatuses,
    const double3* dVertsA,
    const uint3* dIndicesA,
    const double3* dVertsB,
    const uint3* dIndicesB,
    int numPairs,
    cudaStream_t stream)
{
    if (numPairs <= 0) return;

    int threadsPerBlock = 256;
    int blocksPerGrid = (numPairs + threadsPerBlock - 1) / threadsPerBlock;
    
    evaluateTwoLapPairsKernel<<<blocksPerGrid, threadsPerBlock, 0, stream>>>(
         dPairStatuses, dCandidatePairs, 
         dVertsA, dIndicesA, 
         dVertsB, dIndicesB, 
         numPairs
    );
}