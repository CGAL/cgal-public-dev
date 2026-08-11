#include "GPUPredicatesCheckBigIntegerV2.h"
#include "GPUPredicatesCommon.h"
#include <cuda_runtime.h>
#include <stdint.h>
#include <cmath>

// ======================================================================
// FIXED-WIDTH EXACT INTEGER ARITHMETIC (384-bit signed-magnitude)
//
// Coordinates are exact bit-for-bit reconstructions of the input doubles
// (mantissa + shared exponent, no rounding). That means a single
// coordinate can need up to ~113 bits (53-bit mantissa + up to 60 bits of
// exponent spread, see SPREAD_CAP below) -- still fits in __int128_t.
// But orient3d's triple product multiplies three such coordinates
// together, which can need up to ~340 bits. Hence this wide type, used
// only for the internal predicate arithmetic. WLIMBS=6 (384 bits) gives
// generous headroom above that ~340-bit worst case.
// ======================================================================
#define WLIMBS_V2 6

struct wideU_V2 { uint64_t w[WLIMBS_V2]; };

__device__ inline wideU_V2 wideU_zero_V2() {
    wideU_V2 r;
    #pragma unroll
    for (int i = 0; i < WLIMBS_V2; i++) r.w[i] = 0;
    return r;
}

__device__ inline bool wideU_is_zero_V2(const wideU_V2& a) {
    uint64_t acc = 0;
    #pragma unroll
    for (int i = 0; i < WLIMBS_V2; i++) acc |= a.w[i];
    return acc == 0;
}

__device__ inline int wideU_cmp_V2(const wideU_V2& a, const wideU_V2& b) {
    for (int i = WLIMBS_V2 - 1; i >= 0; i--) {
        if (a.w[i] != b.w[i]) return (a.w[i] < b.w[i]) ? -1 : 1;
    }
    return 0;
}

__device__ inline wideU_V2 wideU_add_V2(const wideU_V2& a, const wideU_V2& b) {
    wideU_V2 r;
    unsigned __int128 carry = 0;
    #pragma unroll
    for (int i = 0; i < WLIMBS_V2; i++) {
        unsigned __int128 s = (unsigned __int128)a.w[i] + b.w[i] + carry;
        r.w[i] = (uint64_t)s;
        carry = s >> 64;
    }
    return r;
}

__device__ inline wideU_V2 wideU_sub_V2(const wideU_V2& a, const wideU_V2& b) {
    // Precondition: a >= b (magnitude)
    wideU_V2 r;
    __int128 borrow = 0;
    #pragma unroll
    for (int i = 0; i < WLIMBS_V2; i++) {
        __int128 d = (__int128)a.w[i] - (__int128)b.w[i] - borrow;
        if (d < 0) { d += ((__int128)1 << 64); borrow = 1; } else { borrow = 0; }
        r.w[i] = (uint64_t)d;
    }
    return r;
}

__device__ inline wideU_V2 wideU_from_u128_V2(unsigned __int128 v) {
    wideU_V2 r = wideU_zero_V2();
    r.w[0] = (uint64_t)v;
    r.w[1] = (uint64_t)(v >> 64);
    return r;
}

// Schoolbook multiply, WLIMBS x WLIMBS -> keep low WLIMBS limbs.
// Correct (no truncation) as long as the true product fits in WLIMBS*64
// bits, which is guaranteed upstream by SPREAD_CAP.
__device__ inline wideU_V2 wideU_mul_V2(const wideU_V2& a, const wideU_V2& b) {
    uint64_t acc[WLIMBS_V2 * 2];
    #pragma unroll
    for (int i = 0; i < WLIMBS_V2 * 2; i++) acc[i] = 0;

    #pragma unroll
    for (int i = 0; i < WLIMBS_V2; i++) {
        if (a.w[i] == 0) continue;
        unsigned __int128 carry = 0;
        #pragma unroll
        for (int j = 0; j < WLIMBS_V2; j++) {
            unsigned __int128 cur = (unsigned __int128)a.w[i] * b.w[j] + acc[i + j] + carry;
            acc[i + j] = (uint64_t)cur;
            carry = cur >> 64;
        }
        int k = i + WLIMBS_V2;
        while (carry != 0 && k < WLIMBS_V2 * 2) {
            unsigned __int128 cur = (unsigned __int128)acc[k] + carry;
            acc[k] = (uint64_t)cur;
            carry = cur >> 64;
            k++;
        }
    }
    wideU_V2 r;
    #pragma unroll
    for (int i = 0; i < WLIMBS_V2; i++) r.w[i] = acc[i];
    return r;
}

struct wideS_V2 { wideU_V2 mag; bool neg; };

__device__ inline wideS_V2 wideS_from_i128_V2(__int128_t v) {
    wideS_V2 r;
    bool neg = v < 0;
    unsigned __int128 uv = neg ? (unsigned __int128)(-v) : (unsigned __int128)v;
    r.mag = wideU_from_u128_V2(uv);
    r.neg = neg && !wideU_is_zero_V2(r.mag);
    return r;
}

__device__ inline wideS_V2 wideS_mul_V2(const wideS_V2& a, const wideS_V2& b) {
    wideS_V2 r;
    r.mag = wideU_mul_V2(a.mag, b.mag);
    r.neg = (a.neg != b.neg) && !wideU_is_zero_V2(r.mag);
    return r;
}

__device__ inline wideS_V2 wideS_add_V2(const wideS_V2& a, const wideS_V2& b) {
    wideS_V2 r;
    if (a.neg == b.neg) {
        r.mag = wideU_add_V2(a.mag, b.mag);
        r.neg = a.neg && !wideU_is_zero_V2(r.mag);
    } else {
        int c = wideU_cmp_V2(a.mag, b.mag);
        if (c == 0) { r.mag = wideU_zero_V2(); r.neg = false; }
        else if (c > 0) { r.mag = wideU_sub_V2(a.mag, b.mag); r.neg = a.neg; }
        else { r.mag = wideU_sub_V2(b.mag, a.mag); r.neg = b.neg; }
    }
    return r;
}

__device__ inline wideS_V2 wideS_neg_V2(wideS_V2 a) {
    if (!wideU_is_zero_V2(a.mag)) a.neg = !a.neg;
    return a;
}

__device__ inline int wideS_sign_V2(const wideS_V2& a) {
    if (wideU_is_zero_V2(a.mag)) return 0;
    return a.neg ? -1 : 1;
}

__device__ inline wideS_V2 wide_mul_i128_V2(__int128_t a, __int128_t b) {
    return wideS_mul_V2(wideS_from_i128_V2(a), wideS_from_i128_V2(b));
}

// a1*b2 - a2*b1, exact, wide result
__device__ inline wideS_V2 wide_cross_term_V2(__int128_t a1, __int128_t a2, __int128_t b1, __int128_t b2) {
    return wideS_add_V2(wide_mul_i128_V2(a1, b2), wideS_neg_V2(wide_mul_i128_V2(a2, b1)));
}

__device__ inline wideS_V2 wide_mul_scalar_wide_V2(__int128_t s, const wideS_V2& w) {
    return wideS_mul_V2(wideS_from_i128_V2(s), w);
}

// --------------------------------------------------------------------
// EXACT DOUBLE -> INTEGER DECOMPOSITION (lossless)
//
// x == mant * 2^exp2 exactly, for any finite double x. mant fits in
// ~54 signed bits. This is a pure bit-reinterpretation via frexp, not a
// rounding/quantization step -- no precision is lost here.
// --------------------------------------------------------------------
__device__ inline void exact_decompose_V2(double x, __int128_t& mant_out, int& exp_out, bool& nonzero_out) {
    if (x == 0.0) { mant_out = 0; exp_out = 0; nonzero_out = false; return; }
    int e;
    double m = frexp(x, &e);                    // x = m * 2^e, 0.5 <= |m| < 1
    double scaled = m * 9007199254740992.0;      // m * 2^53 -- exact (binary point shift only)
    mant_out = (__int128_t)(int64_t)scaled;      // exact, |mant| < 2^53
    exp_out = e - 53;                            // x == mant_out * 2^exp_out, exactly
    nonzero_out = true;
}

// --------------------------------------------------------------------
// MEMORY HELPER
// --------------------------------------------------------------------
__device__ inline double3 load_double3_V2(const double3* ptr) {
    return make_double3(__ldg(&(ptr->x)), __ldg(&(ptr->y)), __ldg(&(ptr->z)));
}

// --------------------------------------------------------------------
// EXACT COORDINATE TYPE AND PREDICATES
// --------------------------------------------------------------------
struct int3_wide_V2 { __int128_t x, y, z; };

__device__ inline int orient2d_exact_V2(
    __int128_t ax, __int128_t ay,
    __int128_t bx, __int128_t by,
    __int128_t cx, __int128_t cy)
{
    __int128_t acx = ax - cx, acy = ay - cy;
    __int128_t bcx = bx - cx, bcy = by - cy;
    wideS_V2 det = wideS_add_V2(wide_mul_i128_V2(acx, bcy), wideS_neg_V2(wide_mul_i128_V2(acy, bcx)));
    return wideS_sign_V2(det);
}

__device__ inline int orient3d_exact_V2(
    const int3_wide_V2& a, const int3_wide_V2& b,
    const int3_wide_V2& c, const int3_wide_V2& d)
{
    __int128_t adx = a.x - d.x, ady = a.y - d.y, adz = a.z - d.z;
    __int128_t bdx = b.x - d.x, bdy = b.y - d.y, bdz = b.z - d.z;
    __int128_t cdx = c.x - d.x, cdy = c.y - d.y, cdz = c.z - d.z;

    wideS_V2 m0 = wide_cross_term_V2(bdx, bdy, cdx, cdy);
    wideS_V2 m1 = wide_cross_term_V2(bdx, bdz, cdx, cdz);
    wideS_V2 m2 = wide_cross_term_V2(bdy, bdz, cdy, cdz);

    wideS_V2 t0 = wide_mul_scalar_wide_V2(adz, m0);
    wideS_V2 t1 = wide_mul_scalar_wide_V2(ady, m1);
    wideS_V2 t2 = wide_mul_scalar_wide_V2(adx, m2);

    wideS_V2 det = wideS_add_V2(wideS_add_V2(t0, wideS_neg_V2(t1)), t2);
    return wideS_sign_V2(det);
}

// --------------------------------------------------------------------
// FULL 2D SAT COPLANAR INTERSECTION (winding-independent)
// --------------------------------------------------------------------
__device__ inline bool all_outside_edge_V2(
    const __int128_t Eu[3], const __int128_t Ev[3],
    const __int128_t Ou[3], const __int128_t Ov[3],
    int i, int j, int k)
{
    int o_ref = orient2d_exact_V2(Eu[i], Ev[i], Eu[j], Ev[j], Eu[k], Ev[k]);
    if (o_ref == 0) return false;

    int o0 = orient2d_exact_V2(Eu[i], Ev[i], Eu[j], Ev[j], Ou[0], Ov[0]);
    int o1 = orient2d_exact_V2(Eu[i], Ev[i], Eu[j], Ev[j], Ou[1], Ov[1]);
    int o2 = orient2d_exact_V2(Eu[i], Ev[i], Eu[j], Ev[j], Ou[2], Ov[2]);

    bool out0 = (o_ref > 0) ? (o0 < 0) : (o0 > 0);
    bool out1 = (o_ref > 0) ? (o1 < 0) : (o1 > 0);
    bool out2 = (o_ref > 0) ? (o2 < 0) : (o2 > 0);

    return out0 && out1 && out2;
}

__device__ inline PairStatus coplanar_tri_tri_exact_V2(
    const int3_wide_V2& Aa, const int3_wide_V2& Ab, const int3_wide_V2& Ac,
    const int3_wide_V2& Ba, const int3_wide_V2& Bb, const int3_wide_V2& Bc)
{
    __int128_t d1x = Ab.x - Aa.x, d1y = Ab.y - Aa.y, d1z = Ab.z - Aa.z;
    __int128_t d2x = Ac.x - Aa.x, d2y = Ac.y - Aa.y, d2z = Ac.z - Aa.z;

    wideS_V2 nx = wide_cross_term_V2(d1y, d1z, d2y, d2z);
    wideS_V2 ny = wide_cross_term_V2(d1z, d1x, d2z, d2x);
    wideS_V2 nz = wide_cross_term_V2(d1x, d1y, d2x, d2y);

    bool useX = (wideU_cmp_V2(nx.mag, ny.mag) >= 0) && (wideU_cmp_V2(nx.mag, nz.mag) >= 0);
    bool useY = (!useX) && (wideU_cmp_V2(ny.mag, nz.mag) >= 0);

    __int128_t Au[3], Av[3], Bu[3], Bv[3];

    if (useX) {
        Au[0] = Aa.y; Av[0] = Aa.z; Au[1] = Ab.y; Av[1] = Ab.z; Au[2] = Ac.y; Av[2] = Ac.z;
        Bu[0] = Ba.y; Bv[0] = Ba.z; Bu[1] = Bb.y; Bv[1] = Bb.z; Bu[2] = Bc.y; Bv[2] = Bc.z;
    } else if (useY) {
        Au[0] = Aa.x; Av[0] = Aa.z; Au[1] = Ab.x; Av[1] = Ab.z; Au[2] = Ac.x; Av[2] = Ac.z;
        Bu[0] = Ba.x; Bv[0] = Ba.z; Bu[1] = Bb.x; Bv[1] = Bb.z; Bu[2] = Bc.x; Bv[2] = Bc.z;
    } else {
        Au[0] = Aa.x; Av[0] = Aa.y; Au[1] = Ab.x; Av[1] = Ab.y; Au[2] = Ac.x; Av[2] = Ac.y;
        Bu[0] = Ba.x; Bv[0] = Ba.y; Bu[1] = Bb.x; Bv[1] = Bb.y; Bu[2] = Bc.x; Bv[2] = Bc.y;
    }

    for (int i = 0; i < 3; i++) {
        int j = (i + 1) % 3, k = (i + 2) % 3;
        if (all_outside_edge_V2(Au, Av, Bu, Bv, i, j, k)) return PAIR_NO;
    }
    for (int i = 0; i < 3; i++) {
        int j = (i + 1) % 3, k = (i + 2) % 3;
        if (all_outside_edge_V2(Bu, Bv, Au, Av, i, j, k)) return PAIR_NO;
    }
    return PAIR_GREEN;
}

// --------------------------------------------------------------------
// COPLANAR SEGMENT-VS-TRIANGLE (exact 2D)
// --------------------------------------------------------------------
__device__ inline bool point_in_or_on_tri_2d_V2(
    __int128_t px, __int128_t py,
    const __int128_t Tu[3], const __int128_t Tv[3])
{
    for (int i = 0; i < 3; i++) {
        int j = (i + 1) % 3, k = (i + 2) % 3;
        int o_ref = orient2d_exact_V2(Tu[i], Tv[i], Tu[j], Tv[j], Tu[k], Tv[k]);
        if (o_ref == 0) continue;
        int o_pt = orient2d_exact_V2(Tu[i], Tv[i], Tu[j], Tv[j], px, py);
        bool outside = (o_ref > 0) ? (o_pt < 0) : (o_pt > 0);
        if (outside) return false;
    }
    return true;
}

__device__ inline bool on_segment_collinear_V2(
    __int128_t px, __int128_t py, __int128_t qx, __int128_t qy,
    __int128_t rx, __int128_t ry)
{
    __int128_t minx = (px < qx) ? px : qx, maxx = (px > qx) ? px : qx;
    __int128_t miny = (py < qy) ? py : qy, maxy = (py > qy) ? py : qy;
    return (rx >= minx && rx <= maxx && ry >= miny && ry <= maxy);
}

__device__ inline bool segments_intersect_2d_V2(
    __int128_t ax, __int128_t ay, __int128_t bx, __int128_t by,
    __int128_t cx, __int128_t cy, __int128_t dx, __int128_t dy)
{
    int o1 = orient2d_exact_V2(ax, ay, bx, by, cx, cy);
    int o2 = orient2d_exact_V2(ax, ay, bx, by, dx, dy);
    int o3 = orient2d_exact_V2(cx, cy, dx, dy, ax, ay);
    int o4 = orient2d_exact_V2(cx, cy, dx, dy, bx, by);

    if (o1 != 0 && o2 != 0 && o3 != 0 && o4 != 0 && o1 != o2 && o3 != o4)
        return true;

    if (o1 == 0 && on_segment_collinear_V2(ax, ay, bx, by, cx, cy)) return true;
    if (o2 == 0 && on_segment_collinear_V2(ax, ay, bx, by, dx, dy)) return true;
    if (o3 == 0 && on_segment_collinear_V2(cx, cy, dx, dy, ax, ay)) return true;
    if (o4 == 0 && on_segment_collinear_V2(cx, cy, dx, dy, bx, by)) return true;

    return false;
}

__device__ inline PairStatus edgeTri_coplanar_exact_V2(
    const int3_wide_V2& a, const int3_wide_V2& b,
    const int3_wide_V2& p, const int3_wide_V2& q, const int3_wide_V2& r)
{
    __int128_t d1x = q.x - p.x, d1y = q.y - p.y, d1z = q.z - p.z;
    __int128_t d2x = r.x - p.x, d2y = r.y - p.y, d2z = r.z - p.z;

    wideS_V2 nx = wide_cross_term_V2(d1y, d1z, d2y, d2z);
    wideS_V2 ny = wide_cross_term_V2(d1z, d1x, d2z, d2x);
    wideS_V2 nz = wide_cross_term_V2(d1x, d1y, d2x, d2y);

    bool useX = (wideU_cmp_V2(nx.mag, ny.mag) >= 0) && (wideU_cmp_V2(nx.mag, nz.mag) >= 0);
    bool useY = (!useX) && (wideU_cmp_V2(ny.mag, nz.mag) >= 0);

    __int128_t Tu[3], Tv[3], au, av, bu, bv;

    if (useX) {
        Tu[0] = p.y; Tv[0] = p.z; Tu[1] = q.y; Tv[1] = q.z; Tu[2] = r.y; Tv[2] = r.z;
        au = a.y; av = a.z; bu = b.y; bv = b.z;
    } else if (useY) {
        Tu[0] = p.x; Tv[0] = p.z; Tu[1] = q.x; Tv[1] = q.z; Tu[2] = r.x; Tv[2] = r.z;
        au = a.x; av = a.z; bu = b.x; bv = b.z;
    } else {
        Tu[0] = p.x; Tv[0] = p.y; Tu[1] = q.x; Tv[1] = q.y; Tu[2] = r.x; Tv[2] = r.y;
        au = a.x; av = a.y; bu = b.x; bv = b.y;
    }

    if (point_in_or_on_tri_2d_V2(au, av, Tu, Tv)) return PAIR_GREEN;
    if (point_in_or_on_tri_2d_V2(bu, bv, Tu, Tv)) return PAIR_GREEN;

    for (int i = 0; i < 3; i++) {
        int j = (i + 1) % 3;
        if (segments_intersect_2d_V2(au, av, bu, bv, Tu[i], Tv[i], Tu[j], Tv[j]))
            return PAIR_GREEN;
    }
    return PAIR_NO;
}

// --------------------------------------------------------------------
// 3D EDGE-TRIANGLE EXACT PREDICATE
// --------------------------------------------------------------------
__device__ inline int edgeTri_exact_V2(
    const int3_wide_V2& a, const int3_wide_V2& b,
    const int3_wide_V2& p, const int3_wide_V2& q, const int3_wide_V2& r)
{
    int s0 = orient3d_exact_V2(p, q, r, a);
    int s1 = orient3d_exact_V2(p, q, r, b);

    if ((s0 > 0 && s1 > 0) || (s0 < 0 && s1 < 0)) return PAIR_NO;

    if (s0 == 0 && s1 == 0) {
        return edgeTri_coplanar_exact_V2(a, b, p, q, r);
    }

    int e0 = orient3d_exact_V2(a, b, p, q);
    int e1 = orient3d_exact_V2(a, b, q, r);
    int e2 = orient3d_exact_V2(a, b, r, p);

    if ((e0 >= 0 && e1 >= 0 && e2 >= 0) || (e0 <= 0 && e1 <= 0 && e2 <= 0)) return PAIR_GREEN;

    return PAIR_NO;
}

// --------------------------------------------------------------------
// EXACT INTEGER CLASSIFICATION V2 (fully lossless quantization)
// --------------------------------------------------------------------
__device__ inline PairStatus classifyPairBigIntV2(
    const double3& Aa, const double3& Ab, const double3& Ac,
    const double3& Ba, const double3& Bb, const double3& Bc)
{
    const double3 anchor = Aa;

    double3 relA1 = { Ab.x - anchor.x, Ab.y - anchor.y, Ab.z - anchor.z };
    double3 relA2 = { Ac.x - anchor.x, Ac.y - anchor.y, Ac.z - anchor.z };
    double3 relB0 = { Ba.x - anchor.x, Ba.y - anchor.y, Ba.z - anchor.z };
    double3 relB1 = { Bb.x - anchor.x, Bb.y - anchor.y, Bb.z - anchor.z };
    double3 relB2 = { Bc.x - anchor.x, Bc.y - anchor.y, Bc.z - anchor.z };

    // --- EXACT DECOMPOSITION of every relative coordinate ---
    double coordVals[15] = {
        relA1.x, relA1.y, relA1.z,
        relA2.x, relA2.y, relA2.z,
        relB0.x, relB0.y, relB0.z,
        relB1.x, relB1.y, relB1.z,
        relB2.x, relB2.y, relB2.z
    };

    __int128_t mant[15];
    int exp2[15];
    bool nz[15];
    #pragma unroll
    for (int i = 0; i < 15; i++) exact_decompose_V2(coordVals[i], mant[i], exp2[i], nz[i]);

    int e_min = 0, e_max = 0;
    bool any = false;
    for (int i = 0; i < 15; i++) {
        if (!nz[i]) continue;
        if (!any) { e_min = e_max = exp2[i]; any = true; }
        else {
            if (exp2[i] < e_min) e_min = exp2[i];
            if (exp2[i] > e_max) e_max = exp2[i];
        }
    }
    if (!any) return PAIR_YELLOW; // every relative coordinate is exactly zero

    // Deterministic exactness/overflow guard -- NOT a geometric epsilon.
    // Coordinates end up with |value| < 2^(53 + spread). orient3d's triple
    // product then needs roughly 3*(54+spread) bits, which must stay under
    // the WLIMBS_V2*64 = 384-bit budget. spread <= 60 keeps that at ~342
    // bits with margin to spare. Only pairs with genuinely huge dynamic
    // range between their coordinates (~10^18x) ever hit this.
    const int SPREAD_CAP_V2 = 60;
    if ((e_max - e_min) > SPREAD_CAP_V2) return PAIR_YELLOW;

    __int128_t q[15];
    #pragma unroll
    for (int i = 0; i < 15; i++) {
        q[i] = nz[i] ? (mant[i] * (((__int128_t)1) << (exp2[i] - e_min))) : (__int128_t)0;
    }

    int3_wide_V2 qA0 = { 0, 0, 0 };
    int3_wide_V2 qA1 = { q[0], q[1], q[2] };
    int3_wide_V2 qA2 = { q[3], q[4], q[5] };
    int3_wide_V2 qB0 = { q[6], q[7], q[8] };
    int3_wide_V2 qB1 = { q[9], q[10], q[11] };
    int3_wide_V2 qB2 = { q[12], q[13], q[14] };

    // --- CHECK 1: 0D POINT & 1D LINE COLLAPSE GUARD (now exact facts, not rounding artifacts) ---
    bool aDegenerate =
        (qA1.x == 0 && qA1.y == 0 && qA1.z == 0 &&
         qA2.x == 0 && qA2.y == 0 && qA2.z == 0);

    bool bDegenerate =
        (qB0.x == qB1.x && qB0.y == qB1.y && qB0.z == qB1.z &&
         qB0.x == qB2.x && qB0.y == qB2.y && qB0.z == qB2.z);

    if (aDegenerate || bDegenerate) return PAIR_YELLOW;

    wideS_V2 nxA = wide_cross_term_V2(qA1.y, qA1.z, qA2.y, qA2.z);
    wideS_V2 nyA = wide_cross_term_V2(qA1.z, qA1.x, qA2.z, qA2.x);
    wideS_V2 nzA = wide_cross_term_V2(qA1.x, qA1.y, qA2.x, qA2.y);

    __int128_t b10x = qB1.x - qB0.x, b10y = qB1.y - qB0.y, b10z = qB1.z - qB0.z;
    __int128_t b20x = qB2.x - qB0.x, b20y = qB2.y - qB0.y, b20z = qB2.z - qB0.z;

    wideS_V2 nxB = wide_cross_term_V2(b10y, b10z, b20y, b20z);
    wideS_V2 nyB = wide_cross_term_V2(b10z, b10x, b20z, b20x);
    wideS_V2 nzB = wide_cross_term_V2(b10x, b10y, b20x, b20y);

    if ((wideU_is_zero_V2(nxA.mag) && wideU_is_zero_V2(nyA.mag) && wideU_is_zero_V2(nzA.mag)) ||
        (wideU_is_zero_V2(nxB.mag) && wideU_is_zero_V2(nyB.mag) && wideU_is_zero_V2(nzB.mag))) {
        return PAIR_YELLOW; // triangle collapsed to a line
    }

    // --- CHECK 2: PREDICATES ---
    // ob0==ob1==ob2==0 below is now a genuine exact-coplanarity fact
    // (the quantization is lossless), so there's no spurious-coplanarity
    // artifact to filter -- the old double-precision distance re-check is
    // gone.
    int ob0 = orient3d_exact_V2(qA0, qA1, qA2, qB0);
    int ob1 = orient3d_exact_V2(qA0, qA1, qA2, qB1);
    int ob2 = orient3d_exact_V2(qA0, qA1, qA2, qB2);

    if (ob0 == 0 && ob1 == 0 && ob2 == 0) {
        return coplanar_tri_tri_exact_V2(qA0, qA1, qA2, qB0, qB1, qB2);
    }
    if ((ob0 > 0 && ob1 > 0 && ob2 > 0) || (ob0 < 0 && ob1 < 0 && ob2 < 0)) return PAIR_NO;

    int oa0 = orient3d_exact_V2(qB0, qB1, qB2, qA0);
    int oa1 = orient3d_exact_V2(qB0, qB1, qB2, qA1);
    int oa2 = orient3d_exact_V2(qB0, qB1, qB2, qA2);

    if (oa0 == 0 && oa1 == 0 && oa2 == 0) {
        return coplanar_tri_tri_exact_V2(qA0, qA1, qA2, qB0, qB1, qB2);
    }
    if ((oa0 > 0 && oa1 > 0 && oa2 > 0) || (oa0 < 0 && oa1 < 0 && oa2 < 0)) return PAIR_NO;

    // --- CHECK 3: EDGE-TRIANGLE PREDICATES ---
    int r;
    r = edgeTri_exact_V2(qA0, qA1, qB0, qB1, qB2); if (r == PAIR_GREEN) return PAIR_GREEN;
    r = edgeTri_exact_V2(qA1, qA2, qB0, qB1, qB2); if (r == PAIR_GREEN) return PAIR_GREEN;
    r = edgeTri_exact_V2(qA2, qA0, qB0, qB1, qB2); if (r == PAIR_GREEN) return PAIR_GREEN;

    r = edgeTri_exact_V2(qB0, qB1, qA0, qA1, qA2); if (r == PAIR_GREEN) return PAIR_GREEN;
    r = edgeTri_exact_V2(qB1, qB2, qA0, qA1, qA2); if (r == PAIR_GREEN) return PAIR_GREEN;
    r = edgeTri_exact_V2(qB2, qB0, qA0, qA1, qA2); if (r == PAIR_GREEN) return PAIR_GREEN;

    return PAIR_NO;
}

// --------------------------------------------------------------------
// KERNEL IMPLEMENTATION
// --------------------------------------------------------------------
__global__ void evaluateGeometricPairsKernelBigInt_KernelV2(
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

    double3 Aa = load_double3_V2(&dVertsA[idxA.x]);
    double3 Ab = load_double3_V2(&dVertsA[idxA.y]);
    double3 Ac = load_double3_V2(&dVertsA[idxA.z]);

    double3 Ba = load_double3_V2(&dVertsB[idxB.x]);
    double3 Bb = load_double3_V2(&dVertsB[idxB.y]);
    double3 Bc = load_double3_V2(&dVertsB[idxB.z]);

    PairStatus status = classifyPairBigIntV2(Aa, Ab, Ac, Ba, Bb, Bc);
    dPairStatuses[tid] = (int)status;
}

// --------------------------------------------------------------------
// HOST LAUNCHER (C LINKAGE)
// --------------------------------------------------------------------
extern "C" void evaluateGeometricPairsKernelBigIntV2(
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

    evaluateGeometricPairsKernelBigInt_KernelV2<<<gridSize, blockSize, 0, stream>>>(
        dYellowCandidatePairs,
        dPairStatuses,
        dVertsA,
        dIndicesA,
        dVertsB,
        dIndicesB,
        numYellowPairs
    );
}