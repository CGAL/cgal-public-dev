#include "GPUPredicatesCheckBigInteger.h"
#include "GPUPredicatesCommon.h"
#include <cuda_runtime.h>
#include <stdint.h>
#include <cmath>

// --------------------------------------------------------------------
// MEMORY HELPER 
// --------------------------------------------------------------------
__device__ inline double3 load_double3(const double3* ptr) {
    return make_double3(__ldg(&(ptr->x)), __ldg(&(ptr->y)), __ldg(&(ptr->z)));
}

// --------------------------------------------------------------------
// 128-BIT EXACT INTEGER PREDICATES
// --------------------------------------------------------------------
struct int3_128 {
    int64_t x, y, z;
};

__device__ inline int orient2d_exact(
    int64_t ax, int64_t ay, 
    int64_t bx, int64_t by, 
    int64_t cx, int64_t cy) 
{
    int64_t acx = ax - cx;
    int64_t acy = ay - cy;
    int64_t bcx = bx - cx;
    int64_t bcy = by - cy;

    // Widened to __int128_t: at higher quantization bit-depths, acx*bcy can
    // exceed int64_t range. Doing this in plain int64_t was silently
    // overflowing/wrapping at high precision, which is why the caller used
    // to be capped at a tiny 21-bit quantization budget.
    __int128_t det = (__int128_t)acx * bcy - (__int128_t)acy * bcx;
    if (det > 0) return 1;
    if (det < 0) return -1;
    return 0;
}

__device__ inline int orient3d_exact(
    const int3_128& a, const int3_128& b, 
    const int3_128& c, const int3_128& d) 
{
    int64_t adx = a.x - d.x; int64_t ady = a.y - d.y; int64_t adz = a.z - d.z;
    int64_t bdx = b.x - d.x; int64_t bdy = b.y - d.y; int64_t bdz = b.z - d.z;
    int64_t cdx = c.x - d.x; int64_t cdy = c.y - d.y; int64_t cdz = c.z - d.z;

    __int128_t m0 = (__int128_t)bdx * cdy - (__int128_t)bdy * cdx;
    __int128_t m1 = (__int128_t)bdx * cdz - (__int128_t)bdz * cdx;
    __int128_t m2 = (__int128_t)bdy * cdz - (__int128_t)bdz * cdy;

    __int128_t det = (__int128_t)adz * m0 
                   - (__int128_t)ady * m1 
                   + (__int128_t)adx * m2;

    if (det > 0) return 1;
    if (det < 0) return -1;
    return 0; 
}

__device__ inline __int128_t abs128(__int128_t v) {
    return (v < 0) ? -v : v;
}

// --------------------------------------------------------------------
// FULL 2D SAT COPLANAR INTERSECTION (winding-independent)
// --------------------------------------------------------------------
// For each edge (i,j) of a triangle, the "inside" side is derived from
// that triangle's own opposite vertex k, rather than assumed to be a
// fixed sign. This makes the test correct regardless of which axis got
// dropped during projection, or the resulting winding order.
__device__ inline bool all_outside_edge(
    const int64_t Eu[3], const int64_t Ev[3],  // edge-owner triangle (u,v) coords
    const int64_t Ou[3], const int64_t Ov[3],  // other triangle's (u,v) coords being tested
    int i, int j, int k)
{
    int o_ref = orient2d_exact(Eu[i], Ev[i], Eu[j], Ev[j], Eu[k], Ev[k]);
    if (o_ref == 0) return false; // degenerate edge basis, can't separate on this axis

    int o0 = orient2d_exact(Eu[i], Ev[i], Eu[j], Ev[j], Ou[0], Ov[0]);
    int o1 = orient2d_exact(Eu[i], Ev[i], Eu[j], Ev[j], Ou[1], Ov[1]);
    int o2 = orient2d_exact(Eu[i], Ev[i], Eu[j], Ev[j], Ou[2], Ov[2]);

    bool out0 = (o_ref > 0) ? (o0 < 0) : (o0 > 0);
    bool out1 = (o_ref > 0) ? (o1 < 0) : (o1 > 0);
    bool out2 = (o_ref > 0) ? (o2 < 0) : (o2 > 0);

    return out0 && out1 && out2;
}

__device__ inline PairStatus coplanar_tri_tri_exact(
    const int3_128& Aa, const int3_128& Ab, const int3_128& Ac,
    const int3_128& Ba, const int3_128& Bb, const int3_128& Bc) 
{
    // Widened to __int128_t for the same overflow reason as orient2d_exact:
    // these are only used for a magnitude comparison to pick the dominant
    // projection axis, but at higher quantization precision the raw
    // int64_t products can exceed int64_t range.
    __int128_t nx = (__int128_t)(Ab.y - Aa.y) * (Ac.z - Aa.z) - (__int128_t)(Ab.z - Aa.z) * (Ac.y - Aa.y);
    __int128_t ny = (__int128_t)(Ab.z - Aa.z) * (Ac.x - Aa.x) - (__int128_t)(Ab.x - Aa.x) * (Ac.z - Aa.z);
    __int128_t nz = (__int128_t)(Ab.x - Aa.x) * (Ac.y - Aa.y) - (__int128_t)(Ab.y - Aa.y) * (Ac.x - Aa.x);

    __int128_t abs_nx = abs128(nx), abs_ny = abs128(ny), abs_nz = abs128(nz);
    
    int64_t Au[3], Av[3], Bu[3], Bv[3];

    if (abs_nx >= abs_ny && abs_nx >= abs_nz) { 
        Au[0] = Aa.y; Av[0] = Aa.z; Au[1] = Ab.y; Av[1] = Ab.z; Au[2] = Ac.y; Av[2] = Ac.z;
        Bu[0] = Ba.y; Bv[0] = Ba.z; Bu[1] = Bb.y; Bv[1] = Bb.z; Bu[2] = Bc.y; Bv[2] = Bc.z;
    } else if (abs_ny >= abs_nx && abs_ny >= abs_nz) { 
        Au[0] = Aa.x; Av[0] = Aa.z; Au[1] = Ab.x; Av[1] = Ab.z; Au[2] = Ac.x; Av[2] = Ac.z;
        Bu[0] = Ba.x; Bv[0] = Ba.z; Bu[1] = Bb.x; Bv[1] = Bb.z; Bu[2] = Bc.x; Bv[2] = Bc.z;
    } else { 
        Au[0] = Aa.x; Av[0] = Aa.y; Au[1] = Ab.x; Av[1] = Ab.y; Au[2] = Ac.x; Av[2] = Ac.y;
        Bu[0] = Ba.x; Bv[0] = Ba.y; Bu[1] = Bb.x; Bv[1] = Bb.y; Bu[2] = Bc.x; Bv[2] = Bc.y;
    }

    // SAT: is B completely outside any edge of A? (sign derived from A's own opposite vertex)
    for (int i = 0; i < 3; i++) {
        int j = (i + 1) % 3;
        int k = (i + 2) % 3;
        if (all_outside_edge(Au, Av, Bu, Bv, i, j, k)) return PAIR_NO;
    }

    // SAT: is A completely outside any edge of B?
    for (int i = 0; i < 3; i++) {
        int j = (i + 1) % 3;
        int k = (i + 2) % 3;
        if (all_outside_edge(Bu, Bv, Au, Av, i, j, k)) return PAIR_NO;
    }

    return PAIR_GREEN;
}

// --------------------------------------------------------------------
// COPLANAR SEGMENT-VS-TRIANGLE (exact 2D)
// --------------------------------------------------------------------
// Used only for the degenerate case where a candidate edge lies exactly
// in the triangle's plane. In that situation orient3d_exact(a,b,p,q) etc.
// are always zero (four coplanar points), so the 3D containment test in
// edgeTri_exact can't distinguish "crosses the triangle" from "misses it
// entirely" -- it would trivially report GREEN every time. This does the
// equivalent test properly, in 2D, after projecting onto the triangle's
// dominant axis.

__device__ inline bool point_in_or_on_tri_2d(
    int64_t px, int64_t py,
    const int64_t Tu[3], const int64_t Tv[3])
{
    // Winding-independent, like all_outside_edge: derive the "inside" sign
    // for each edge from the triangle's own opposite vertex.
    for (int i = 0; i < 3; i++) {
        int j = (i + 1) % 3;
        int k = (i + 2) % 3;
        int o_ref = orient2d_exact(Tu[i], Tv[i], Tu[j], Tv[j], Tu[k], Tv[k]);
        if (o_ref == 0) continue; // degenerate triangle guard
        int o_pt = orient2d_exact(Tu[i], Tv[i], Tu[j], Tv[j], px, py);
        bool outside = (o_ref > 0) ? (o_pt < 0) : (o_pt > 0);
        if (outside) return false;
    }
    return true;
}

__device__ inline bool on_segment_collinear(
    int64_t px, int64_t py, int64_t qx, int64_t qy,
    int64_t rx, int64_t ry)
{
    // Precondition: p, q, r are already known to be collinear.
    // Checks whether r lies within the bounding box of segment [p,q].
    int64_t minx = (px < qx) ? px : qx, maxx = (px > qx) ? px : qx;
    int64_t miny = (py < qy) ? py : qy, maxy = (py > qy) ? py : qy;
    return (rx >= minx && rx <= maxx && ry >= miny && ry <= maxy);
}

__device__ inline bool segments_intersect_2d(
    int64_t ax, int64_t ay, int64_t bx, int64_t by,
    int64_t cx, int64_t cy, int64_t dx, int64_t dy)
{
    int o1 = orient2d_exact(ax, ay, bx, by, cx, cy);
    int o2 = orient2d_exact(ax, ay, bx, by, dx, dy);
    int o3 = orient2d_exact(cx, cy, dx, dy, ax, ay);
    int o4 = orient2d_exact(cx, cy, dx, dy, bx, by);

    // Proper crossing: c,d strictly straddle line ab and vice versa
    if (o1 != 0 && o2 != 0 && o3 != 0 && o4 != 0 && o1 != o2 && o3 != o4)
        return true;

    // Collinear touching / overlapping cases
    if (o1 == 0 && on_segment_collinear(ax, ay, bx, by, cx, cy)) return true;
    if (o2 == 0 && on_segment_collinear(ax, ay, bx, by, dx, dy)) return true;
    if (o3 == 0 && on_segment_collinear(cx, cy, dx, dy, ax, ay)) return true;
    if (o4 == 0 && on_segment_collinear(cx, cy, dx, dy, bx, by)) return true;

    return false;
}

__device__ inline PairStatus edgeTri_coplanar_exact(
    const int3_128 &a, const int3_128 &b,
    const int3_128 &p, const int3_128 &q, const int3_128 &r)
{
    __int128_t nx = (__int128_t)(q.y - p.y) * (r.z - p.z) - (__int128_t)(q.z - p.z) * (r.y - p.y);
    __int128_t ny = (__int128_t)(q.z - p.z) * (r.x - p.x) - (__int128_t)(q.x - p.x) * (r.z - p.z);
    __int128_t nz = (__int128_t)(q.x - p.x) * (r.y - p.y) - (__int128_t)(q.y - p.y) * (r.x - p.x);

    __int128_t abs_nx = abs128(nx), abs_ny = abs128(ny), abs_nz = abs128(nz);

    int64_t Tu[3], Tv[3], au, av, bu, bv;

    if (abs_nx >= abs_ny && abs_nx >= abs_nz) {
        Tu[0] = p.y; Tv[0] = p.z; Tu[1] = q.y; Tv[1] = q.z; Tu[2] = r.y; Tv[2] = r.z;
        au = a.y; av = a.z; bu = b.y; bv = b.z;
    } else if (abs_ny >= abs_nx && abs_ny >= abs_nz) {
        Tu[0] = p.x; Tv[0] = p.z; Tu[1] = q.x; Tv[1] = q.z; Tu[2] = r.x; Tv[2] = r.z;
        au = a.x; av = a.z; bu = b.x; bv = b.z;
    } else {
        Tu[0] = p.x; Tv[0] = p.y; Tu[1] = q.x; Tv[1] = q.y; Tu[2] = r.x; Tv[2] = r.y;
        au = a.x; av = a.y; bu = b.x; bv = b.y;
    }

    // Either endpoint inside (or on the boundary of) the triangle?
    if (point_in_or_on_tri_2d(au, av, Tu, Tv)) return PAIR_GREEN;
    if (point_in_or_on_tri_2d(bu, bv, Tu, Tv)) return PAIR_GREEN;

    // Otherwise, does the segment cross any triangle edge?
    for (int i = 0; i < 3; i++) {
        int j = (i + 1) % 3;
        if (segments_intersect_2d(au, av, bu, bv, Tu[i], Tv[i], Tu[j], Tv[j]))
            return PAIR_GREEN;
    }

    return PAIR_NO;
}

// --------------------------------------------------------------------
// 3D EDGE-TRIANGLE EXACT PREDICATE
// --------------------------------------------------------------------
__device__ inline int edgeTri_exact(
    const int3_128 &a, const int3_128 &b,
    const int3_128 &p, const int3_128 &q, const int3_128 &r) 
{
    int s0 = orient3d_exact(p, q, r, a);
    int s1 = orient3d_exact(p, q, r, b);

    // If both segment endpoints are on the same (non-zero) side of the plane
    if ((s0 > 0 && s1 > 0) || (s0 < 0 && s1 < 0)) return PAIR_NO;

    // Special case: the edge lies exactly in the triangle's plane.
    // orient3d_exact(a,b,p,q)-style tests below are always zero here
    // (four coplanar points), so they can't tell containment apart from
    // a clean miss -- fall through to the exact 2D test instead.
    if (s0 == 0 && s1 == 0) {
        return edgeTri_coplanar_exact(a, b, p, q, r);
    }

    int e0 = orient3d_exact(a, b, p, q);
    int e1 = orient3d_exact(a, b, q, r);
    int e2 = orient3d_exact(a, b, r, p);

    // If segment pierces the triangle interior
    if ((e0 >= 0 && e1 >= 0 && e2 >= 0) || (e0 <= 0 && e1 <= 0 && e2 <= 0)) return PAIR_GREEN;

    return PAIR_NO;
}

// --------------------------------------------------------------------
// EXACT INTEGER CLASSIFICATION
// --------------------------------------------------------------------
__device__ inline PairStatus classifyPairBigInt(
    const double3& Aa, const double3& Ab, const double3& Ac,
    const double3& Ba, const double3& Bb, const double3& Bc) 
{    
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

    if (M_local == 0.0) return PAIR_YELLOW; // Degenerate geometry

    // 1. Dynamic Quantization Scale
    //
    // The candidates reaching this function are exactly the near-touching,
    // borderline pairs the coarse float filter couldn't resolve -- i.e.
    // precision matters most right here. The old 21-bit budget gave a grid
    // spacing of M_local / 2^21 (~1 part in 2 million): any two points
    // closer than that snapped to the SAME integer, turning genuinely
    // separate-but-close geometry into "touching" and inflating the
    // intersection count. That's why GPU/orange was finding more hits than
    // the CPU reference, never fewer.
    //
    // 37 bits is the largest budget that stays safely inside __int128_t
    // for the worst-case term in orient3d_exact: quantized coordinate
    // differences are < 2^(B+1), so a 2x2 minor is < 2^(2B+3), and the
    // final adz*m0-style triple product summed over 3 terms is bounded by
    // roughly 2^(3B+6). At B=37 that's ~2^117, versus __int128_t's ~2^127
    // ceiling -- about a 1000x safety margin. (orient2d_exact and the
    // dominant-axis normal computation are now widened to __int128_t too,
    // so they're no longer the binding constraint the way they were when
    // this was capped at 21 bits.)
    const double MAX_SAFE_VAL = (double)((1LL << 37) - 1); 
    const double dynamic_scale = MAX_SAFE_VAL / M_local;

    // Round to nearest instead of truncating toward zero, to avoid a
    // consistent directional bias in the quantized coordinates.
    int3_128 qA0 = {0, 0, 0};
    int3_128 qA1 = { llround(relA1.x * dynamic_scale), llround(relA1.y * dynamic_scale), llround(relA1.z * dynamic_scale) };
    int3_128 qA2 = { llround(relA2.x * dynamic_scale), llround(relA2.y * dynamic_scale), llround(relA2.z * dynamic_scale) };

    int3_128 qB0 = { llround(relB0.x * dynamic_scale), llround(relB0.y * dynamic_scale), llround(relB0.z * dynamic_scale) };
    int3_128 qB1 = { llround(relB1.x * dynamic_scale), llround(relB1.y * dynamic_scale), llround(relB1.z * dynamic_scale) };
    int3_128 qB2 = { llround(relB2.x * dynamic_scale), llround(relB2.y * dynamic_scale), llround(relB2.z * dynamic_scale) };

    // 2. Disproportional Scale Guard (SYMMETRIC)
    // M_local is the max extent over BOTH triangles, so whichever triangle
    // is much smaller than the other can collapse to a single point after
    // quantization. The original code only checked B for this; a small A
    // next to a large B slipped through as a degenerate triangle at the
    // anchor. Check both directions and offload to CPU exact (yellow) if
    // either collapses.
    bool aDegenerate =
        (qA1.x == 0 && qA1.y == 0 && qA1.z == 0 &&
         qA2.x == 0 && qA2.y == 0 && qA2.z == 0);

    bool bDegenerate =
        (qB0.x == qB1.x && qB0.y == qB1.y && qB0.z == qB1.z &&
         qB0.x == qB2.x && qB0.y == qB2.y && qB0.z == qB2.z);

    if (aDegenerate || bDegenerate) {
        return PAIR_YELLOW;
    }

    // 3. Early Rejection Checks
    int ob0 = orient3d_exact(qA0, qA1, qA2, qB0);
    int ob1 = orient3d_exact(qA0, qA1, qA2, qB1);
    int ob2 = orient3d_exact(qA0, qA1, qA2, qB2);

    if (ob0 == 0 && ob1 == 0 && ob2 == 0) {
        return coplanar_tri_tri_exact(qA0, qA1, qA2, qB0, qB1, qB2);
    }
    if ((ob0 > 0 && ob1 > 0 && ob2 > 0) || (ob0 < 0 && ob1 < 0 && ob2 < 0)) return PAIR_NO;

    int oa0 = orient3d_exact(qB0, qB1, qB2, qA0);
    int oa1 = orient3d_exact(qB0, qB1, qB2, qA1);
    int oa2 = orient3d_exact(qB0, qB1, qB2, qA2);
    
    if ((oa0 > 0 && oa1 > 0 && oa2 > 0) || (oa0 < 0 && oa1 < 0 && oa2 < 0)) return PAIR_NO;

    // 4. Restored Detailed Edge-Triangle Predicates
    int r;
    r = edgeTri_exact(qA0, qA1, qB0, qB1, qB2); if (r == PAIR_GREEN) return PAIR_GREEN;
    r = edgeTri_exact(qA1, qA2, qB0, qB1, qB2); if (r == PAIR_GREEN) return PAIR_GREEN;
    r = edgeTri_exact(qA2, qA0, qB0, qB1, qB2); if (r == PAIR_GREEN) return PAIR_GREEN;

    r = edgeTri_exact(qB0, qB1, qA0, qA1, qA2); if (r == PAIR_GREEN) return PAIR_GREEN;
    r = edgeTri_exact(qB1, qB2, qA0, qA1, qA2); if (r == PAIR_GREEN) return PAIR_GREEN;
    r = edgeTri_exact(qB2, qB0, qA0, qA1, qA2); if (r == PAIR_GREEN) return PAIR_GREEN;

    return PAIR_NO;
}

// --------------------------------------------------------------------
// KERNEL IMPLEMENTATION
// --------------------------------------------------------------------
__global__ void evaluateGeometricPairsKernelBigInt_Kernel(
    const int2* dYellowCandidatePairs, 
    int* dPairStatuses, 
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

    double3 Aa = load_double3(&dVertsA[idxA.x]);
    double3 Ab = load_double3(&dVertsA[idxA.y]);
    double3 Ac = load_double3(&dVertsA[idxA.z]);

    double3 Ba = load_double3(&dVertsB[idxB.x]);
    double3 Bb = load_double3(&dVertsB[idxB.y]);
    double3 Bc = load_double3(&dVertsB[idxB.z]);

    PairStatus status = classifyPairBigInt(Aa, Ab, Ac, Ba, Bb, Bc);
    dPairStatuses[tid] = (int)status; 
}

// --------------------------------------------------------------------
// HOST LAUNCHER (C LINKAGE)
// --------------------------------------------------------------------
extern "C" void evaluateGeometricPairsKernelBigInt(
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

    evaluateGeometricPairsKernelBigInt_Kernel<<<gridSize, blockSize, 0, stream>>>(
        dYellowCandidatePairs, 
        dPairStatuses, 
        dVertsA, 
        dIndicesA, 
        dVertsB, 
        dIndicesB, 
        numYellowPairs
    );
}