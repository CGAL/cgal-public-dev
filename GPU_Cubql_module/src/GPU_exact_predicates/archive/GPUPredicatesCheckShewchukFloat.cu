#include "GPUPredicatesCheckShewchukFloat.h"
#include "GPUPredicatesCommon.h"
#include <cuda_runtime.h>
#include <cmath>
#include <stdint.h>

// --------------------------------------------------------------------
// HELPERS
// --------------------------------------------------------------------
__device__ inline float3 load_float3_from_double3(const double3* ptr) {
    return make_float3((float)__ldg(&(ptr->x)),
                       (float)__ldg(&(ptr->y)),
                       (float)__ldg(&(ptr->z)));
}

// --------------------------------------------------------------------
// SHEWCHUK EXPANSION ARITHMETIC (float version, max length 64)
// --------------------------------------------------------------------
#define MAX_EXP_LEN 64

struct Expansion {
    float parts[MAX_EXP_LEN];
    int len;
};

// TwoSum: (a+b) exactly as an expansion of length ≤ 2
__device__ inline Expansion two_sum(float a, float b) {
    Expansion e;
    float s = a + b;
    float v = s - a;
    float e0 = (a - (s - v)) + (b - v);
    e.parts[0] = s;
    e.parts[1] = e0;
    e.len = (e0 != 0.0f) ? 2 : 1;
    return e;
}

// TwoProd: (a*b) exactly as an expansion of length ≤ 2
__device__ inline Expansion two_prod(float a, float b) {
    Expansion e;
    float p = a * b;
    float e0 = __fmaf_rn(a, b, -p); // exact product error
    e.parts[0] = p;
    e.parts[1] = e0;
    e.len = (e0 != 0.0f) ? 2 : 1;
    return e;
}

// Grow expansion by one term (add a scalar to an expansion)
__device__ inline Expansion grow_expansion(const Expansion& e, float x) {
    Expansion out;
    out.len = 0;
    float q = x;
    for (int i = 0; i < e.len; ++i) {
        Expansion sum = two_sum(e.parts[i], q);
        out.parts[out.len++] = sum.parts[1];
        q = sum.parts[0];
    }
    out.parts[out.len++] = q;
    // Compress to remove zeros and sort by absolute value descending
    // (simplified compression; for exact sign we need sorted parts)
    int new_len = 0;
    for (int i = 0; i < out.len; ++i) {
        if (out.parts[i] != 0.0f) {
            out.parts[new_len++] = out.parts[i];
        }
    }
    out.len = new_len;
    // Sort descending by absolute value (simple bubble sort)
    for (int i = 0; i < out.len; ++i) {
        for (int j = i+1; j < out.len; ++j) {
            if (fabsf(out.parts[i]) < fabsf(out.parts[j])) {
                float tmp = out.parts[i];
                out.parts[i] = out.parts[j];
                out.parts[j] = tmp;
            }
        }
    }
    return out;
}

// Add two expansions
__device__ inline Expansion add_expansions(const Expansion& a, const Expansion& b) {
    Expansion out = { {0.0f}, 0 };
    const Expansion* bigger = (a.len >= b.len) ? &a : &b;
    const Expansion* smaller = (a.len >= b.len) ? &b : &a;
    out = grow_expansion(*bigger, smaller->parts[0]);
    for (int i = 1; i < smaller->len; ++i) {
        out = grow_expansion(out, smaller->parts[i]);
    }
    return out;
}

// Multiply two expansions (naive: O(n*m) but n,m ≤ 2)
__device__ inline Expansion multiply_expansions(const Expansion& a, const Expansion& b) {
    Expansion out = { {0.0f}, 0 };
    for (int i = 0; i < a.len; ++i) {
        for (int j = 0; j < b.len; ++j) {
            Expansion prod = two_prod(a.parts[i], b.parts[j]);
            out = add_expansions(out, prod);
        }
    }
    return out;
}

// Negate expansion
__device__ inline Expansion negate_expansion(const Expansion& e) {
    Expansion neg;
    neg.len = e.len;
    for (int i = 0; i < e.len; ++i) neg.parts[i] = -e.parts[i];
    return neg;
}

// --------------------------------------------------------------------
// EXACT ORIENTATION PREDICATES (float version)
// --------------------------------------------------------------------

__device__ inline int orient2d_exact_float(
    float ax, float ay,
    float bx, float by,
    float cx, float cy)
{
    Expansion acx = two_sum(ax, -cx);
    Expansion acy = two_sum(ay, -cy);
    Expansion bcx = two_sum(bx, -cx);
    Expansion bcy = two_sum(by, -cy);

    Expansion term1 = multiply_expansions(acx, bcy);
    Expansion term2 = multiply_expansions(acy, bcx);
    Expansion neg_term2 = negate_expansion(term2);
    Expansion det = add_expansions(term1, neg_term2);

    // Sum expansion parts (sorted descending by absolute value)
    float sum = 0.0f;
    for (int i = 0; i < det.len; ++i) sum += det.parts[i];
    if (sum > 0.0f) return 1;
    if (sum < 0.0f) return -1;
    return 0;
}

__device__ inline int orient3d_exact_float(
    float ax, float ay, float az,
    float bx, float by, float bz,
    float cx, float cy, float cz,
    float dx, float dy, float dz)
{
    Expansion adx = two_sum(ax, -dx);
    Expansion ady = two_sum(ay, -dy);
    Expansion adz = two_sum(az, -dz);
    Expansion bdx = two_sum(bx, -dx);
    Expansion bdy = two_sum(by, -dy);
    Expansion bdz = two_sum(bz, -dz);
    Expansion cdx = two_sum(cx, -dx);
    Expansion cdy = two_sum(cy, -dy);
    Expansion cdz = two_sum(cz, -dz);

    Expansion m0 = add_expansions(
        multiply_expansions(bdy, cdz),
        multiply_expansions(bdz, negate_expansion(cdy))
    );
    Expansion m1 = add_expansions(
        multiply_expansions(bdz, cdx),
        multiply_expansions(bdx, negate_expansion(cdz))
    );
    Expansion m2 = add_expansions(
        multiply_expansions(bdx, cdy),
        multiply_expansions(bdy, negate_expansion(cdx))
    );

    Expansion term0 = multiply_expansions(adz, m0);
    Expansion term1 = multiply_expansions(ady, m1);
    Expansion term2 = multiply_expansions(adx, m2);

    Expansion neg_term1 = negate_expansion(term1);
    Expansion det = add_expansions(add_expansions(term0, neg_term1), term2);

    float sum = 0.0f;
    for (int i = 0; i < det.len; ++i) sum += det.parts[i];
    if (sum > 0.0f) return 1;
    if (sum < 0.0f) return -1;
    return 0;
}

// --------------------------------------------------------------------
// EXACT TRIANGLE–TRIANGLE INTERSECTION (float predicates)
// --------------------------------------------------------------------

__device__ inline bool all_outside_edge_exact_float(
    float Eu[3], float Ev[3],
    float Ou[3], float Ov[3],
    int i, int j, int k)
{
    int o_ref = orient2d_exact_float(Eu[i], Ev[i], Eu[j], Ev[j], Eu[k], Ev[k]);
    if (o_ref == 0) return false;

    int o0 = orient2d_exact_float(Eu[i], Ev[i], Eu[j], Ev[j], Ou[0], Ov[0]);
    int o1 = orient2d_exact_float(Eu[i], Ev[i], Eu[j], Ev[j], Ou[1], Ov[1]);
    int o2 = orient2d_exact_float(Eu[i], Ev[i], Eu[j], Ev[j], Ou[2], Ov[2]);

    bool out0 = (o_ref > 0) ? (o0 < 0) : (o0 > 0);
    bool out1 = (o_ref > 0) ? (o1 < 0) : (o1 > 0);
    bool out2 = (o_ref > 0) ? (o2 < 0) : (o2 > 0);
    return out0 && out1 && out2;
}

__device__ inline PairStatus coplanar_tri_tri_exact_float(
    float Aa[3], float Ab[3], float Ac[3],
    float Ba[3], float Bb[3], float Bc[3])
{
    float nx = (Ab[1]-Aa[1])*(Ac[2]-Aa[2]) - (Ab[2]-Aa[2])*(Ac[1]-Aa[1]);
    float ny = (Ab[2]-Aa[2])*(Ac[0]-Aa[0]) - (Ab[0]-Aa[0])*(Ac[2]-Aa[2]);
    float nz = (Ab[0]-Aa[0])*(Ac[1]-Aa[1]) - (Ab[1]-Aa[1])*(Ac[0]-Aa[0]);
    float abs_nx = fabsf(nx), abs_ny = fabsf(ny), abs_nz = fabsf(nz);

    float Au[3], Av[3], Bu[3], Bv[3];
    if (abs_nx >= abs_ny && abs_nx >= abs_nz) {
        Au[0]=Aa[1]; Av[0]=Aa[2]; Au[1]=Ab[1]; Av[1]=Ab[2]; Au[2]=Ac[1]; Av[2]=Ac[2];
        Bu[0]=Ba[1]; Bv[0]=Ba[2]; Bu[1]=Bb[1]; Bv[1]=Bb[2]; Bu[2]=Bc[1]; Bv[2]=Bc[2];
    } else if (abs_ny >= abs_nx && abs_ny >= abs_nz) {
        Au[0]=Aa[0]; Av[0]=Aa[2]; Au[1]=Ab[0]; Av[1]=Ab[2]; Au[2]=Ac[0]; Av[2]=Ac[2];
        Bu[0]=Ba[0]; Bv[0]=Ba[2]; Bu[1]=Bb[0]; Bv[1]=Bb[2]; Bu[2]=Bc[0]; Bv[2]=Bc[2];
    } else {
        Au[0]=Aa[0]; Av[0]=Aa[1]; Au[1]=Ab[0]; Av[1]=Ab[1]; Au[2]=Ac[0]; Av[2]=Ac[1];
        Bu[0]=Ba[0]; Bv[0]=Ba[1]; Bu[1]=Bb[0]; Bv[1]=Bb[1]; Bu[2]=Bc[0]; Bv[2]=Bc[1];
    }

    for (int i=0; i<3; ++i) {
        int j = (i+1)%3, k = (i+2)%3;
        if (all_outside_edge_exact_float(Au, Av, Bu, Bv, i, j, k)) return PAIR_NO;
    }
    for (int i=0; i<3; ++i) {
        int j = (i+1)%3, k = (i+2)%3;
        if (all_outside_edge_exact_float(Bu, Bv, Au, Av, i, j, k)) return PAIR_NO;
    }
    return PAIR_GREEN;
}

__device__ inline bool point_in_or_on_tri_exact_float(
    float px, float py,
    float Tu[3], float Tv[3])
{
    for (int i=0; i<3; ++i) {
        int j=(i+1)%3, k=(i+2)%3;
        int o_ref = orient2d_exact_float(Tu[i], Tv[i], Tu[j], Tv[j], Tu[k], Tv[k]);
        if (o_ref == 0) continue;
        int o_pt = orient2d_exact_float(Tu[i], Tv[i], Tu[j], Tv[j], px, py);
        bool outside = (o_ref > 0) ? (o_pt < 0) : (o_pt > 0);
        if (outside) return false;
    }
    return true;
}

__device__ inline bool on_segment_collinear_exact_float(
    float px, float py, float qx, float qy, float rx, float ry)
{
    float minx = (px < qx) ? px : qx, maxx = (px > qx) ? px : qx;
    float miny = (py < qy) ? py : qy, maxy = (py > qy) ? py : qy;
    return (rx >= minx && rx <= maxx && ry >= miny && ry <= maxy);
}

__device__ inline bool segments_intersect_2d_exact_float(
    float ax, float ay, float bx, float by,
    float cx, float cy, float dx, float dy)
{
    int o1 = orient2d_exact_float(ax, ay, bx, by, cx, cy);
    int o2 = orient2d_exact_float(ax, ay, bx, by, dx, dy);
    int o3 = orient2d_exact_float(cx, cy, dx, dy, ax, ay);
    int o4 = orient2d_exact_float(cx, cy, dx, dy, bx, by);

    if (o1 != 0 && o2 != 0 && o3 != 0 && o4 != 0 && o1 != o2 && o3 != o4)
        return true;

    if (o1 == 0 && on_segment_collinear_exact_float(ax, ay, bx, by, cx, cy)) return true;
    if (o2 == 0 && on_segment_collinear_exact_float(ax, ay, bx, by, dx, dy)) return true;
    if (o3 == 0 && on_segment_collinear_exact_float(cx, cy, dx, dy, ax, ay)) return true;
    if (o4 == 0 && on_segment_collinear_exact_float(cx, cy, dx, dy, bx, by)) return true;
    return false;
}

__device__ inline PairStatus edgeTri_coplanar_exact_float(
    float a[3], float b[3],
    float p[3], float q[3], float r[3])
{
    float nx = (q[1]-p[1])*(r[2]-p[2]) - (q[2]-p[2])*(r[1]-p[1]);
    float ny = (q[2]-p[2])*(r[0]-p[0]) - (q[0]-p[0])*(r[2]-p[2]);
    float nz = (q[0]-p[0])*(r[1]-p[1]) - (q[1]-p[1])*(r[0]-p[0]);
    float abs_nx = fabsf(nx), abs_ny = fabsf(ny), abs_nz = fabsf(nz);

    float Tu[3], Tv[3], au, av, bu, bv;
    if (abs_nx >= abs_ny && abs_nx >= abs_nz) {
        Tu[0]=p[1]; Tv[0]=p[2]; Tu[1]=q[1]; Tv[1]=q[2]; Tu[2]=r[1]; Tv[2]=r[2];
        au=a[1]; av=a[2]; bu=b[1]; bv=b[2];
    } else if (abs_ny >= abs_nx && abs_ny >= abs_nz) {
        Tu[0]=p[0]; Tv[0]=p[2]; Tu[1]=q[0]; Tv[1]=q[2]; Tu[2]=r[0]; Tv[2]=r[2];
        au=a[0]; av=a[2]; bu=b[0]; bv=b[2];
    } else {
        Tu[0]=p[0]; Tv[0]=p[1]; Tu[1]=q[0]; Tv[1]=q[1]; Tu[2]=r[0]; Tv[2]=r[1];
        au=a[0]; av=a[1]; bu=b[0]; bv=b[1];
    }

    if (point_in_or_on_tri_exact_float(au, av, Tu, Tv)) return PAIR_GREEN;
    if (point_in_or_on_tri_exact_float(bu, bv, Tu, Tv)) return PAIR_GREEN;

    for (int i=0; i<3; ++i) {
        int j=(i+1)%3;
        if (segments_intersect_2d_exact_float(au, av, bu, bv, Tu[i], Tv[i], Tu[j], Tv[j]))
            return PAIR_GREEN;
    }
    return PAIR_NO;
}

__device__ inline int edgeTri_exact_float(
    float a[3], float b[3],
    float p[3], float q[3], float r[3])
{
    int s0 = orient3d_exact_float(p[0],p[1],p[2], q[0],q[1],q[2], r[0],r[1],r[2], a[0],a[1],a[2]);
    int s1 = orient3d_exact_float(p[0],p[1],p[2], q[0],q[1],q[2], r[0],r[1],r[2], b[0],b[1],b[2]);

    if ((s0 > 0 && s1 > 0) || (s0 < 0 && s1 < 0)) return PAIR_NO;
    if (s0 == 0 && s1 == 0) {
        return edgeTri_coplanar_exact_float(a, b, p, q, r);
    }

    int e0 = orient3d_exact_float(a[0],a[1],a[2], b[0],b[1],b[2], p[0],p[1],p[2], q[0],q[1],q[2]);
    int e1 = orient3d_exact_float(a[0],a[1],a[2], b[0],b[1],b[2], q[0],q[1],q[2], r[0],r[1],r[2]);
    int e2 = orient3d_exact_float(a[0],a[1],a[2], b[0],b[1],b[2], r[0],r[1],r[2], p[0],p[1],p[2]);

    if ((e0 >= 0 && e1 >= 0 && e2 >= 0) || (e0 <= 0 && e1 <= 0 && e2 <= 0))
        return PAIR_GREEN;
    return PAIR_NO;
}

// --------------------------------------------------------------------
// MAIN CLASSIFICATION (float exact)
// --------------------------------------------------------------------
__device__ inline PairStatus classifyPairShewchukFloat(
    const double3& Aa, const double3& Ab, const double3& Ac,
    const double3& Ba, const double3& Bb, const double3& Bc)
{
    // Convert to float arrays
    float A[3][3] = {
        {(float)Aa.x, (float)Aa.y, (float)Aa.z},
        {(float)Ab.x, (float)Ab.y, (float)Ab.z},
        {(float)Ac.x, (float)Ac.y, (float)Ac.z}
    };
    float B[3][3] = {
        {(float)Ba.x, (float)Ba.y, (float)Ba.z},
        {(float)Bb.x, (float)Bb.y, (float)Bb.z},
        {(float)Bc.x, (float)Bc.y, (float)Bc.z}
    };

    // Early rejection tests
    int ob0 = orient3d_exact_float(A[0][0],A[0][1],A[0][2], A[1][0],A[1][1],A[1][2], A[2][0],A[2][1],A[2][2], B[0][0],B[0][1],B[0][2]);
    int ob1 = orient3d_exact_float(A[0][0],A[0][1],A[0][2], A[1][0],A[1][1],A[1][2], A[2][0],A[2][1],A[2][2], B[1][0],B[1][1],B[1][2]);
    int ob2 = orient3d_exact_float(A[0][0],A[0][1],A[0][2], A[1][0],A[1][1],A[1][2], A[2][0],A[2][1],A[2][2], B[2][0],B[2][1],B[2][2]);
    if ((ob0 > 0 && ob1 > 0 && ob2 > 0) || (ob0 < 0 && ob1 < 0 && ob2 < 0))
        return PAIR_NO;

    int oa0 = orient3d_exact_float(B[0][0],B[0][1],B[0][2], B[1][0],B[1][1],B[1][2], B[2][0],B[2][1],B[2][2], A[0][0],A[0][1],A[0][2]);
    int oa1 = orient3d_exact_float(B[0][0],B[0][1],B[0][2], B[1][0],B[1][1],B[1][2], B[2][0],B[2][1],B[2][2], A[1][0],A[1][1],A[1][2]);
    int oa2 = orient3d_exact_float(B[0][0],B[0][1],B[0][2], B[1][0],B[1][1],B[1][2], B[2][0],B[2][1],B[2][2], A[2][0],A[2][1],A[2][2]);
    if ((oa0 > 0 && oa1 > 0 && oa2 > 0) || (oa0 < 0 && oa1 < 0 && oa2 < 0))
        return PAIR_NO;

    // Coplanar check
    if (ob0 == 0 && ob1 == 0 && ob2 == 0 && oa0 == 0 && oa1 == 0 && oa2 == 0) {
        return coplanar_tri_tri_exact_float(A[0], A[1], A[2], B[0], B[1], B[2]);
    }

    // Edge-triangle tests
    int r;
    r = edgeTri_exact_float(A[0], A[1], B[0], B[1], B[2]); if (r == PAIR_GREEN) return PAIR_GREEN;
    r = edgeTri_exact_float(A[1], A[2], B[0], B[1], B[2]); if (r == PAIR_GREEN) return PAIR_GREEN;
    r = edgeTri_exact_float(A[2], A[0], B[0], B[1], B[2]); if (r == PAIR_GREEN) return PAIR_GREEN;

    r = edgeTri_exact_float(B[0], B[1], A[0], A[1], A[2]); if (r == PAIR_GREEN) return PAIR_GREEN;
    r = edgeTri_exact_float(B[1], B[2], A[0], A[1], A[2]); if (r == PAIR_GREEN) return PAIR_GREEN;
    r = edgeTri_exact_float(B[2], B[0], A[0], A[1], A[2]); if (r == PAIR_GREEN) return PAIR_GREEN;

    return PAIR_NO;
}

// --------------------------------------------------------------------
// KERNEL
// --------------------------------------------------------------------
__global__ void evaluatePairsShewchukFloatKernel(
    int* dPairStatuses,
    const int2* dYellowCandidatePairs,
    const double3* dVertsA,
    const uint3* dIndicesA,
    const double3* dVertsB,
    const uint3* dIndicesB,
    int numPairs)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= numPairs) return;

    int2 pair = dYellowCandidatePairs[tid];
    uint3 idxA = dIndicesA[pair.x];
    uint3 idxB = dIndicesB[pair.y];

    // Load double3 and convert to float3
    float3 Aa = load_float3_from_double3(&dVertsA[idxA.x]);
    float3 Ab = load_float3_from_double3(&dVertsA[idxA.y]);
    float3 Ac = load_float3_from_double3(&dVertsA[idxA.z]);
    float3 Ba = load_float3_from_double3(&dVertsB[idxB.x]);
    float3 Bb = load_float3_from_double3(&dVertsB[idxB.y]);
    float3 Bc = load_float3_from_double3(&dVertsB[idxB.z]);

    // Need to pass double3 to classifier; we'll convert inside
    double3 dAa = make_double3(Aa.x, Aa.y, Aa.z);
    double3 dAb = make_double3(Ab.x, Ab.y, Ab.z);
    double3 dAc = make_double3(Ac.x, Ac.y, Ac.z);
    double3 dBa = make_double3(Ba.x, Ba.y, Ba.z);
    double3 dBb = make_double3(Bb.x, Bb.y, Bb.z);
    double3 dBc = make_double3(Bc.x, Bc.y, Bc.z);

    PairStatus status = classifyPairShewchukFloat(dAa, dAb, dAc, dBa, dBb, dBc);
    dPairStatuses[tid] = (int)status;
}

// --------------------------------------------------------------------
// HOST LAUNCHER
// --------------------------------------------------------------------
extern "C" void evaluateAndCompactPairsShewchukFloat(
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
    int blockSize = 256;
    int gridSize = (numYellowPairs + blockSize - 1) / blockSize;
    evaluatePairsShewchukFloatKernel<<<gridSize, blockSize, 0, stream>>>(
        dPairStatuses, dYellowCandidatePairs,
        dVertsA, dIndicesA,
        dVertsB, dIndicesB,
        numYellowPairs
    );
}