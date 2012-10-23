// ============================================================================
//
// Copyright (c) 2001-2008 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// this file is not part of any library ;-)
//
// ----------------------------------------------------------------------------
//
// Library       : CUDA MP
//
// File          : 
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

 
#ifndef _ZMOD_DEV_CU_
#define _ZMOD_DEV_CU_

//! contant memory space for moduli set
__constant__ unsigned dev_const_mem[CUMP_DEVICE_MODULI_SZ];

#if CUMP_USE_32BIT_MODULI_SET
typedef double2 fp_limb2;
#else
typedef float2 fp_limb2;
#endif

#define CUMP_USE_SAD_TRICK 1  // whether to use 'sum of absolute different'
                              // trick to compute FFT butterfly
#define CUMP_USE_FAST_TRUNC 1 // whether to use fast 23-bit float2int
                              // truncation in reduce_m operation

#define CUMP_D2I_TRUNC (double)(3ll << 51)

#if CUMP_USE_32BIT_MODULI_SET
#define MODARITHM_INIT(var) do { \
    mx100 = 0, invm = CUMP_D2I_TRUNC + (var >> 20);  \
    invk = __hiloint2double(mods[3], mods[2]); \
    } while(0);
#else
#define MODARITHM_INIT(var) do { \
    mx100 = UMUL(100, m); \
    invk = __int_as_float(mods[2]); \
    invm = __int_as_float(mods[3]); \
    } while(0);
#endif

__device__ __forceinline__ unsigned sub_m(unsigned a, unsigned b, volatile unsigned m) {

#if CUMP_USE_32BIT_MODULI_SET
    unsigned c = a - b;
//     VADDx(c, c, m, c, "min") // == umin(r, r + m);
//     c = umin(c, c + m); // replace by c < 0 ??
    if((int)c < 0)
        c += m;
    return c;
#else
    unsigned c = a - b;
    if((int)c < 0)
        c += m;
    return c;
#endif
}

__device__ __forceinline__ unsigned add_m(unsigned a, unsigned b, volatile unsigned m) {

#if CUMP_USE_32BIT_MODULI_SET
    // VSUBx(c, c, m, c, "min") // == umin(r, r + m);
//     unsigned mb = m - b;
//     if(a >= mb)
//         return a - mb;
//     else
//         return m - mb + a; // a + b ??
    unsigned c = a + b;
    if(c >= m)
        c -= m;
    return c;
#else
    unsigned c = a + b;
    if(c >= m)
        c -= m;
    return c;
#endif
}

#if !CUMP_USE_32BIT_MODULI_SET
__device__ __forceinline__ unsigned mul_no_reduce_m(unsigned a, unsigned b,
    volatile unsigned m, volatile fp_limb invk) {

    unsigned mulhi, l;
    UMUL24HI(mulhi, a, b)
    float prodf = __uint2float_rn(mulhi); // _rz sometimes does not work...

    // _rz does not work here either
#if CUDART_VERSION < 2000
    prodf = prodf * invk;
#else
    prodf = __fmul_rn(prodf, invk);
#endif
    l = __float2uint_rz(prodf);

    unsigned resid = __umul24(a, b) - __umul24(l, m);
    return resid;
}
#endif

//! when \c CUMP_USE_32BIT_MODULI_SET is on, the pair (invk; aux) represents
//! a single double value: 1 / m, in 24-bit mode \c aux is not used
__device__ __forceinline__ unsigned mul_m(unsigned a, unsigned b, volatile unsigned m,
    volatile fp_limb invk) {

#if CUMP_USE_32BIT_MODULI_SET

//     unsigned hi = __umulhi(a, b);
//     double rf = __uint2double_rn(hi*4) * invk + CUMP_D2I_TRUNC;
//     unsigned r = (unsigned)__double2loint(rf);
//     r = a * b - r * m;
//     if((int)r < 0)
//         r += m;
// 
//     if((int)r >= m)
//         r -= m;

    unsigned hi = __umulhi(a*2, b*2); // 3 flops

    // 2 double instructions
    double rf = __uint2double_rn(hi) * invk + CUMP_D2I_TRUNC;
    unsigned r = (unsigned)__double2loint(rf);
    r = a * b - r * m; // 2 flops
    if((int)r < 0) // smth more
        r += m;
    //VADDx(r, r, m, r, "min") // == umin(r, r + m);

    return r;
#else
    unsigned r = mul_no_reduce_m(a, b, m, invk);
    if((int)r < 0)
        r = r + __umul24(m, 0x1000002);

    r = umin(r, r - m);
    return r;
#endif
}

//! computes: a - b * c mod m
__device__ __forceinline__ unsigned submul_m(unsigned a, unsigned b, unsigned c,
    volatile unsigned m, volatile fp_limb invk) {

#if CUMP_USE_32BIT_MODULI_SET

    unsigned r = mul_m(b, c, m, invk);
    r = sub_m(a, r, m);

    return r;
#else

    unsigned mulhi, l;
    UMUL24HI(mulhi, b, c)
    float prodf = __uint2float_rn(mulhi); // _rz sometimes does not work...

    if(0) { // UseExpTrunc
        limb _2e23 = 0x4b000000;
        prodf = prodf * invk + __int_as_float(_2e23);
        l = __float_as_int(prodf);
    } else {
#if CUDART_VERSION < 2000
        prodf = prodf * invk;
#else
        prodf = __fmul_rn(prodf, invk);
#endif
        l = __float2uint_rz(prodf);
    }

    unsigned r = a - __umul24(b, c) + __umul24(l, m);
    if((int)r < 0)
        r = r + __umul24(m, 0x1000002);

    r = umin(r, r - m);
    r = umin(r, r - m);
    return r;

#endif
}


#if CUMP_USE_32BIT_MODULI_SET
#define mul_small_m mul_m
#else
//! the same as \c mul_m but \c b is assumed to be small, s.t.
//! \c l = a*b*2^16 / m fits 23-bits
__device__ __forceinline__ unsigned mul_small_m(unsigned a, unsigned b,
    volatile unsigned m, volatile fp_limb invk) {

    unsigned mulhi, l;
    UMUL24HI(mulhi, a, b)
    float prodf = __uint2float_rn(mulhi); // _rz sometimes does not work...

    // _rz does not work here either
    limb _2e23 = 0x4b000000;
    prodf = prodf * invk + __int_as_float(_2e23);
    l = __float_as_int(prodf);

    // NOTE: for some reason fast float2int trick results in wider range
    // residue..
    unsigned resid = __umul24(a, b) - __umul24(l, m);

    if((int)resid < 0)
       resid = resid + __umul24(m, 0x1000002);

    resid = umin(resid, resid - m);
    return resid;
}
#endif

//! compute a + b*c mod m
__device__  __forceinline__ unsigned addmul_m(unsigned a, unsigned b, unsigned c,
    volatile unsigned m, volatile fp_limb invk) {

#if CUMP_USE_32BIT_MODULI_SET

//     typedef unsigned long long u64;
//     return (unsigned)((((u64)b * (u64)c) + a)%m);

    unsigned r = mul_m(b, c, m, invk);
    r = add_m(a, r, m);
    return r;
#else
    unsigned mulhi, l;
    UMUL24HI(mulhi, b, c);
    float prodf = __uint2float_rn(mulhi); // _rz sometimes does not work...
#if CUDART_VERSION >= 2000
    prodf = __fmul_rn(prodf, invk);
#else
    prodf = prodf * invk;
#endif
    l = __float2uint_rz(prodf);

    unsigned r = a + __umul24(b, c) - __umul24(l, m);
    if((int)r < 0)
        r = r + __umul24(m, 0x1000002);

    r = umin(r, r - m);
    r = umin(r, r - m);
    return r;
#endif
}

#if CUMP_USE_32BIT_MODULI_SET
__device__ __forceinline__ unsigned sub_mul_reduce_mX(unsigned a1, unsigned b1,
        unsigned a2, unsigned b2, volatile unsigned m, volatile fp_limb invk,
        volatile fp_limb invm, volatile unsigned mx100) {

//! here \c CUMP_D2I_TRUNC is loaded to invm which reduces register count

    unsigned h1 = __umulhi(a1*2, b1*2),
             h2 = __umulhi(a2*2, b2*2);
    double rf1 = __uint2double_rn(h1) * invk + invm,
           rf2 = __uint2double_rn(h2) * invk + invm;
    unsigned r1 = (unsigned)__double2loint(rf1),
             r2 = (unsigned)__double2loint(rf2);

    r1 = a1 * b1 - r1 * m;
//     VADDx(r1, r1, m, r1, "min") // == umin(r, r + m);

    if((int)r1 < 0)
        r1 += m;
    r2 = a2 * b2 - r2 * m;
//     VADDx(r2, r2, m, r2, "min") // == umin(r, r + m);
    if((int)r2 < 0)
        r2 += m;

    r1 -= r2;
//     VSUBx(r1, r1, m, r1, "min") // == umin(r, r - m);

    if((int)r1 < 0)
        r1 += m; // umin(r, r - m);

    return r1;

}
#else
#define sub_mul_reduce_mX sub_mul_reduce_m
#endif

//! computes (a1 * b1 - a2 * b2) mod m
//! \c mx100 = 100 * m, \c _2e23 = (float)(1 << 23)
__device__ __forceinline__ unsigned sub_mul_reduce_m(unsigned a1, unsigned b1,
        unsigned a2, unsigned b2, volatile unsigned m, volatile fp_limb invk,
        volatile fp_limb invm, volatile unsigned mx100) {

#if CUMP_USE_32BIT_MODULI_SET

    unsigned h1 = __umulhi(a1*2, b1*2),
             h2 = __umulhi(a2*2, b2*2);
    double rf1 = __uint2double_rn(h1) * invk + CUMP_D2I_TRUNC,
           rf2 = __uint2double_rn(h2) * invk + CUMP_D2I_TRUNC;
    unsigned r1 = (unsigned)__double2loint(rf1),
             r2 = (unsigned)__double2loint(rf2);

    r1 = a1 * b1 - r1 * m;
//     r1 = min(r1, r1 + m);
    if((int)r1 < 0)
        r1 += m;
    r2 = a2 * b2 - r2 * m;
//     r2 = min(r2, r2 + m);
    if((int)r2 < 0)
        r2 += m;

    r1 -= r2;
    if((int)r1 < 0)
        r1 += m; // umin(r, r - m);

    return r1;
#else

    unsigned m1hi, m2hi, l1, l2;
    UMUL24HI(m1hi, a1, b1); UMUL24HI(m2hi, a2, b2);
    float prodf1 = __uint2float_rz(m1hi),
        prodf2 = __uint2float_rz(m2hi);

    prodf1 = __fmul_rz(prodf1, invk);
    prodf2 = __fmul_rz(prodf2, invk);
    l1 = __float2uint_rz(prodf1), l2 = __float2uint_rz(prodf2);

    unsigned a = mx100 + __umul24(a1, b1) - __umul24(l1, m) -
        __umul24(a2, b2) + __umul24(l2, m); // 4 fma instructions

    float af = __uint2float_rn(a); //  af = __fmul_rn(af, invm);

    // fast float2int trunc: add 2**23 to normalize mantissa
    volatile limb _2e23 = 0x4b000000;
    af = af * invm + __int_as_float(_2e23); //  mad.rn
    unsigned l = __float_as_int(af);

    a -= __umul24(l, m); // fused into signed mad24

    if((int)a < 0)
        a += m;
    return a;
#endif // CUMP_USE_32BIT_MODULI_SET
}

//! computes (a1 * b1 + a2 * b2) mod m
//! \c mx100 = 100 * m, \c _2e23 = (float)(1 << 23)
__device__ __forceinline__ unsigned add_mul_reduce_m(unsigned a1, unsigned b1,
        unsigned a2, unsigned b2, volatile unsigned m, volatile fp_limb invk,
        volatile fp_limb invm, volatile unsigned mx100) {

#if CUMP_USE_32BIT_MODULI_SET

#if 1

    unsigned h1 = __umulhi(a1*2, b1*2),
            h2 = __umulhi(a2*2, b2*2);
    double rf1 = __uint2double_rn(h1) * invk + CUMP_D2I_TRUNC,
           rf2 = __uint2double_rn(h2) * invk + CUMP_D2I_TRUNC;
    unsigned r1 = (unsigned)__double2loint(rf1),
             r2 = (unsigned)__double2loint(rf2);

    r1 = a1 * b1 - r1 * m;
//     VADDx(r1, r1, m, r1, "min") // == umin(r, r + m);

    if((int)r1 < 0)
        r1 += m;
    r2 = a2 * b2 - r2 * m;
//     VADDx(r2, r2, m, r2, "min") // == umin(r, r + m);
    if((int)r2 < 0)
        r2 += m;

    r1 += r2;
//     VSUBx(r1, r1, m, r1, "min") // == umin(r, r - m);
    if(r1 >= m)
        r1 -= m; // umin(r, r - m);

    return r1;

#else
    unsigned r1, r2;
    double f1 = __uint2double_rn(a1) * __uint2double_rn(b1);
    f1 = f1 * invk + CUMP_D2I_TRUNC;
    double f2 = __uint2double_rn(a2) * __uint2double_rn(b2);
    f2 = f2 * invk + CUMP_D2I_TRUNC;

    r1 = a1 * b1 - (unsigned)__double2loint(f1) * m;
    // interesting but results differ from vadd.. 
//     VADDx(r1, r1, m, r1, "min") // == umin(r, r + m);
    if((int)r1 < 0)
        r1 += m;

    r2 = a2 * b2 - (unsigned)__double2loint(f2) * m;
    if((int)r2 < 0)
        r2 += m;

    r1 += r2;
    if(r1 >= m)
        r1 -= m; // umin(r, r - m);
    return r1;
#endif

#else

    unsigned m1hi, m2hi, l1, l2;
    UMUL24HI(m1hi, a1, b1); UMUL24HI(m2hi, a2, b2);
    float prodf1 = __uint2float_rz(m1hi),
        prodf2 = __uint2float_rz(m2hi);

    prodf1 = __fmul_rz(prodf1, invk);
    prodf2 = __fmul_rz(prodf2, invk);
    l1 = __float2uint_rz(prodf1), l2 = __float2uint_rz(prodf2);

//! NOTE NOTE NOTE: this needs to be checked out if it works correctly
//! because this was never tested before
    unsigned a = mx100 + __umul24(a1, b1) - __umul24(l1, m) +
        __umul24(a2, b2) - __umul24(l2, m); // 4 fma instructions

    float af = __uint2float_rn(a); //  af = __fmul_rn(af, invm);

    // fast float2int trunc: add 2**23 to normalize mantissa
    volatile limb _2e23 = 0x4b000000;
    af = af * invm + __int_as_float(_2e23); //  mad.rn
    unsigned l = __float_as_int(af);

    a -= __umul24(l, m); // fused into signed mad24

    if((int)a < 0)
        a += m;
    return a;
#endif // CUMP_USE_32BIT_MODULI_SET
}

// fma has different accuracy than that of CPU -> errors
__device__ __forceinline__ limb fma_m(limb a, limb b, limb c, limb m, fp_limb invk) {
    // 6 flops

    limb mulhi, l;
    UMUL24HI(mulhi, a, b)
    //! NOTE NOTE was previously rn / rn !!
    float prodf = __uint2float_rz(mulhi);

#if CUMP_USE_SAD_TRICK && (CUDART_VERSION >= 2000)
    //prodf = __fmul_rn(prodf, __int_as_float(invk)) + 4.0f;
    prodf = (prodf * invk) + 2.0f;
#else

#if CUDART_VERSION < 2000
    prodf = prodf * invk;
#else
    prodf = __fmul_rn(prodf, invk);
#endif

#endif // CUMP_USE_SAD_TRICK
    l = __float2uint_rz(prodf);

    // must be fused into 2 mad24s!!!
    limb resid = c + __umul24(a, b);
    resid = resid - __umul24(l, m);
    
    // result is [-2m + e; 2m + e] if c is potitive
    return resid;
}

//! the same as \c fma_m but performs reduction immediately after
//! assumes that \c [a,b,c] are valid residues !!
__device__ __forceinline__ limb fma_with_reduce_m(limb a, limb b, limb c, limb m,
        fp_limb invk) {

    // result is in [-2m + e; 2m + e] if c is potitive
    limb resid = fma_m(a, b, c, m, invk);

    if((int)resid < 0)
        resid = resid + __umul24(m, 0x1000002);

    resid = umin(resid, resid - m);
    resid = umin(resid, resid - m);
    return resid;
}

// 6 flops ?
__device__ __forceinline__ void reduce_m(limb& a, limb m, fp_limb invm) {

    // in principle one can modify the method, s.t. a has the leading
    // bit set => then conversion to floating-point should be easy
    a += __umul24(100, m);
    // NOTE FIXME ATTENTION: does not work for round-to-zero mode
    float af = __uint2float_rn(a);

    // NOTE NOTE: this is an experimental version of reduce_m
    // use it carefully
#if CUMP_USE_FAST_TRUNC
    float addf = (float)(1 << 23);
    af = af * invm + addf;
    unsigned l = __float_as_int(af);
#else
    af = __fmul_rz(af, invm);
    unsigned l = __float2uint_rd(af); // at most 8 bits long
#endif
 
    a -= __umul24(l, m); // fused into mad24
    if((int)a < 0)
        a += m;
}

//! reduces 48-bit number (hi << 16) + lo modulo m
//! \c lo and \c hi must be 32 bits and overlap in the middle
__device__ __forceinline__ limb large_reduce_m(limb hi, limb lo, limb m, fp_limb invk) {

    float prodf = __uint2float_rz(hi);
    prodf = __fmul_rz(prodf, invk);
    unsigned l = __float2uint_rz(prodf);
    unsigned resid = lo - __umul24(l, m);
    //TODO: you forgot post-reduction step
    if((int)resid < 0)
        resid = resid + __umul24(m, 0x1000002);

    resid = umin(resid, resid - m); // this is seemingly faster
    return resid;
}

//! computes x = c + a * b; y = c - a * b
//! with \c a and \c b are residues modulo \c m
__device__ __forceinline__ void fma_bfy2_m(limb& x, limb &y, limb a, limb b,
        limb c, limb m, limb invk) {

    x = fma_m(a, b, c, m, invk);
#if CUMP_USE_SAD_TRICK && (CUDART_VERSION >= 2000)
    y = __sad(x, c, c);
    //y = c * 2 - x;
#else
    y = c * 2 - x;
#endif
}

//! computes x = c + a * b; y = c - a * b
//! with \c a reduced modulo \c m, \c b must be correct residue !
__device__ __forceinline__ void fma_reduce_bfy2_m(limb& x, limb &y, limb a, limb b,
        limb c, limb m, limb invk, limb invm) {

    reduce_m(a, m, invm);
    fma_bfy2_m(x, y, a, b, c, m, invk);
}

//! multiplies 24-bit numbers, result in [hi; lo]
__device__ __forceinline__ void mul_24_wide(limb& hi, limb& lo, limb a, limb b) {

    lo = __umul24(a, b);
    UMUL24HI(hi, a, b)
}

//! multiplies 24-bit numbers and adds 32-bit number to the result, carry
//! is not saved. Can be acquired as (lo < c)
__device__ __forceinline__ void mad_24_wide(limb& hi, limb& lo, limb a, limb b, limb c) {

    lo = __umul24(a, b) + c;
    UMUL24HI(hi, a, b);
}

//! multiplies 48-bit number [hi; lo] by 24-bit number \c a , returns 72-bit
//! product as [r2; r1; r0] having 16, 32 and 24 bits respectively
//! NOTE: r0's MSB 8 are not masked!!!
__device__ __forceinline__ void mul_48_24_wide(limb& r2, limb& r1, limb& r0,
        limb hi, limb lo, limb a) {

    limb h1;
    mul_24_wide(h1, r0, lo, a);
    h1 >>= 8; // keep only 24 MSB bits
    //r0 &= 0xffffff; // mask out 8 MSB bits that overlap with h1

    // [h3; l3] = m1m2_hi * c.x + h2
    mad_24_wide(r2, r1, hi, a, h1); // [r2;r1] = hi * a + h1
    r2 >>= 16, r2 += (r1 < h1);
}

//! same as before but result is: r2(8), r1(32), r0(32)
__device__ __forceinline__ void mul_48_32_wide(limb& r2, limb& r1, limb& r0,
        limb hi, limb lo, limb a) {

    limb h1;
    mul_24_wide(h1, r0, lo, a);
    h1 >>= 8; // keep only 24 MSB bits

    // [h3; l3] = m1m2_hi * c.x + h2
    mad_24_wide(r2, r1, hi, a, h1); // [r2;r1] = hi * a + h1
    if(r1 < h1)
        r2 += 0x10000;

    r0 = __umul24(r0, 0x1000001) + (r1 << 24);
    r1 = __umul24(r2 & 0xff0000, 0x1000100) + (r1 >> 8);
    r2 >>= 24;
}


//! computes modular inverse using montgomery shift algorithm
//! \c mu = -m^-1 mod beta; \c beta = 2^24
__device__ __forceinline__ unsigned montgomery_inverse(unsigned x, unsigned m, unsigned mu) {

    limb v = x, u = m, s = 1, r = 0;

    volatile int k = 0;
    while(1) {

        int tmprs = r;
        if(v & 1) {
            int safeuv = v;

            if(((int)(safeuv ^ u)) < 0)
                v += u;
            else
                v -= u;

            if(((int)(v ^ safeuv)) < 0) {
                u = safeuv, tmprs = s;
            }
            s += r;
        }
        v = (int)v >> 1, k++;
        r = tmprs << 1;

        if(v == 0)
            break;
    } // warp serialization is a problem here ??

    if(r != 0)
        r = m - r;

    if((int)r < 0)
        r += m;

#if CUMP_USE_32BIT_MODULI_SET
   if(k >= 32) { // r = MonPro(r, 1)

        s = r * mu;
        u = s * m, v = __umulhi(s, m);
        u = u + r;
        r = v + (u < r); // addc ?
        k -= 32;
    }
    if(k == 0)
        return r;

    // r = MonPro(r, 1 << (32 - k))
    limb hi2;
    s = r << (32 - k);
    hi2 = s * mu;
    u = hi2 * m, v = __umulhi(hi2, m);

    hi2 = r >> k;
    u = u + s;
    v = v + hi2 + (u < s);

#else

    limb hi2;
    if(k > 24) {
        s = __umul24(r, mu);
        mul_24_wide(v, u, s, m);
        u = __umul24(u, 0x1000001) + r;
        r = (v >> 8) + (u >> 24);
        k -= 24;
    }
    s = r << (24 - k);
    hi2 = __umul24(s, mu);

    mul_24_wide(v, u, hi2, m);
    hi2 = r >> k;

    u &= 0xffffff;
    u = __umul24(s, 0x1000001) + u;
    v = (v >> 8) + hi2 + (u >> 24);
#endif // CUMP_USE_32BIT_MODULI_SET

    return v;
}

#endif // _ZMOD_DEV_CU_

