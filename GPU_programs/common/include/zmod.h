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

#ifndef ZMOD_H
#define ZMOD_H

#include <include/macros.h>
#include <math.h>
// #include <include/ieee754_fp.h>
// #include <CGAL/FPU.h>

#define D2I_TRUNC (double)(3ll << 51)
// strict operand check for modular arithmetic
#define STRICT_MOD_CHECK 0

namespace CGAL {

namespace internal {

#if STRICT_MOD_CHECK
inline bool residue_check(unsigned x, unsigned m) {
    return ((int)x >= 0 && x < m);
}
#else
inline bool residue_check(unsigned x, unsigned m) {
    return true;
}
#endif

typedef unsigned umod_t;
typedef int      smod_t;

static double g_inv_mod;
static double g_inv_mod_e16;

inline umod_t incr_mod(umod_t a, umod_t m)
{ a++; if ( a==m )  a = 0;  return a; }

inline umod_t decr_mod(umod_t a, umod_t m)
{ if ( a==0 )  a = m - 1; else a--;  return a; }

inline umod_t neg_mod(umod_t b, umod_t m)
{ if ( 0==b )  return 0;  else  return m - b; }


inline umod_t my_umul24(umod_t a, umod_t b) {

    umod_t mask = (1ull << 24) - 1;
    umod_t res = (a & mask) * (b & mask); 
    return (res & 0xffffffff);
}

inline umod_t my_min(umod_t a, umod_t b) {
    return (a < b ? a : b);
}

inline umod_t my_umul24hi(umod_t a, umod_t b) {

    umod_t mask = (1ull << 24) - 1;
    uint64 res = (uint64)(a & mask) * (uint64)(b & mask);
    return (umod_t)(res >> 16); // drop 16 LSBs
}

// returns |a - b| + c
inline smod_t my_sad(smod_t a, smod_t b, smod_t c) {

    int diff = a - b;
    if(diff < 0)
        diff = -diff;
    return diff + c;
}

inline void mul_24_wide_host(unsigned& hi, unsigned& lo,
        unsigned a, unsigned b) {

    lo = my_umul24(a, b);
    hi = my_umul24hi(a, b);
}

//! multiplies 24-bit numbers and adds 32-bit number to the result, carry
//! is not saved. Can be acquired as (lo < c)
inline void mad_24_wide_host(unsigned& hi, unsigned& lo,
        unsigned a, unsigned b, unsigned c) {

    lo = my_umul24(a, b) + c;
    hi = my_umul24hi(a, b);
}

inline umod_t sub_mod(umod_t a, umod_t b, umod_t m) {

#if STRICT_MOD_CHECK
    if(a >= m || b >= m) {
        printf("sub_mod error: out of bounds: a: %d; b: %d; m: %d\n", a, b, m);
        exit(1);
    }
#endif
    umod_t c;
    if ( a>=b )  c =  a - b;
    else         c = m - b + a;
    return c;

//     a -= b;
//     return my_min(a, a + m);
/*    umod_t c = a - b;
    if(a < b)
        c += m;
    return c;*/
}

inline umod_t add_mod(umod_t a, umod_t b, umod_t m) {

#if STRICT_MOD_CHECK
    if(a >= m || b >= m) {
        printf("add_mod error: out of bounds: a: %d; b: %d; m: %d\n", a, b, m);
        exit(1);
    }
#endif
//     a += b;
//     return my_min(a, a - m);
    if(0 == b)
        return a;

    umod_t mb = m - b;
    if(a >= mb)
        return a - mb;
    else
        return m - mb + a;
}

inline umod_t mul_mod(umod_t a, umod_t b, umod_t m) {

#if STRICT_MOD_CHECK
    if(a >= m || b >= m) {
        printf("mul_mod error: out of bounds: a: %u; b: %u; m: %u\n", a, b, m);
        exit(1);
    }
#endif

#ifdef CUMP_USE_32BIT_MODULI_SET

//     modifyFPUStateX86(__FPU_CW_ROUND_MASK__, __FPU_CW_ROUND_CHOP__);
//      CGAL::FPU_set_cw(CGAL_FE_TONEAREST);
    unsigned x = (unsigned)a * (unsigned)b, r;
// some magic crap:  work only when mode is set to nearest here
//     double f = ((double)a * (double)b * g_inv_mod + D2I_TRUNC);
//     r = x - (unsigned&)f * m;

//     g_inv_mod = 1.0 / m;
//     printf("g_inv_mod: %f\n", 1.0/g_inv_mod);

    unsigned f = (unsigned)((double)a * (double)b * g_inv_mod + 0.5);
    r = x - f * m;
    
    if((int)r < 0)
        r += m;

//     unsigned c1 = ((uint64)a * (uint64)b) % m;
//     if(c1 != r) {
//        printf("mul_mod: m: %x a1: %x b1: %x\n",
//                 m, a, b);
//         throw "fuck";
//     }

    return (unsigned)r;
    
#else

#if 1
    unsigned x = (unsigned)a * (unsigned)b, r;
// some magic crap:  work only when mode is set to nearest here
//     double f = ((double)a * (double)b * g_inv_mod + D2I_TRUNC);
//     r = x - (unsigned&)f * m;

//     g_inv_mod = 1.0 / m;
    unsigned f = (unsigned)((double)a * (double)b * g_inv_mod + 0.5);
    r = x - f * m;
    
    if((int)r < 0)
        r += m;

    return (unsigned)r;

#else
    unsigned mullo = my_umul24(a, b), mulhi = my_umul24hi(a, b), l;
    float prodf = float(mulhi);

    float invk = (float)65536.0f / m;
    prodf *= invk;
    l = (unsigned)floorf(prodf);
    unsigned resid = mullo - my_umul24(l, m);

    resid = my_min(resid, resid + m*2);
    resid = my_min(resid, resid - m);
    return resid;
#endif
#endif // CUMP_USE_32BIT_MODULI_SET
}

inline unsigned addmul_mod(unsigned a, unsigned b, unsigned c, unsigned m) {

    b = mul_mod(b, c, m);
    return add_mod(a, b, m);
}

inline void reduce_mod(unsigned& a, unsigned m) {

    unsigned au = a + 100*m;
    float inv = 1.0f / m;

    float af = (float)au;
    af = af * inv;
    unsigned l = (unsigned)floorf(af); // at most 8 bits long
    unsigned resid = au - l * m;

    if((int)resid < 0)
        resid += m;

    a = resid;
    if(a - m < a)
        a -= m;
}

inline umod_t add_mul_reduce_mod(unsigned a1, unsigned b1,
        unsigned a2, unsigned b2, unsigned m) {

    unsigned c1 = mul_mod(a1, b1, m),
            c2 = mul_mod(a2, b2, m);
    return add_mod(c1, c2, m);

/*
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
        r1 -= m; // umin(r, r - m);*/
}


inline umod_t sub_mul_reduce_mod(unsigned a1, unsigned b1,
        unsigned a2, unsigned b2, unsigned m) {

#if CUMP_USE_32BIT_MODULI_SET

    unsigned m1 = a1 * b1, m2 = a2 * b2;
    double mf1 = (double)a1 * (double)b1, mf2 = (double)a2 * (double)b2;

//    g_inv_mod = 1.0 / m;

#if 1
    double f1 = mf1 * g_inv_mod + 0.5,
           f2 = mf2 * g_inv_mod + 0.5;

//     CGAL::FPU_set_cw(CGAL_FE_TONEAREST);    
//  printf("m1 = %x m2 = %x mf1: %f mf2: %f g_inv_mod: %f\n", m1, m2, mf1, mf2,
//              g_inv_mod);

    m1 -= (unsigned)f1 * m;
    if((int)m1 < 0)
        m1 += m;

    m2 -= (unsigned)f2 * m;
    if((int)m2 < 0)
        m2 += m;

    unsigned r = sub_mod(m1, m2, m);
//       unsigned c1 = ((uint64)a1 * (uint64)b1) % m,
//                 c2 = ((uint64)a2 * (uint64)b2) % m;
//       unsigned r2 = sub_mod(c1, c2, m);
// 
//     if(!residue_check(r, m) || r != r2) {
//         printf("sub_mul_reduce_mod: m: %x a1: %x b1: %x a2: %x b2: %x\n",
//                 m, a1, b1, a2, b2);
//         throw "fuck";
//     }
    return r;

#else

    unsigned r1, r2;

    double f12 = floor((mf1 - mf2) * g_inv_mod);
    r2 = m * (unsigned)f12;

    r2 = m1 - m2 - r2;
    // or we need to work with extended precision here..
    if((int)r2 < 0)
        r2 += m;

    if(r2 >= m)
        r2 -= m;
    return r2;
#endif

#else
    unsigned m1hi = my_umul24hi(a1, b1), m2hi = my_umul24hi(a2, b2), l1, l2;
    float prodf1 = float(m1hi), prodf2 = float(m2hi);

    float invk = (float)65536.0f / m;
    prodf1 *= invk, prodf2 *= invk;

    l1 = (unsigned)floorf(prodf1);
    l2 = (unsigned)floorf(prodf2);

    unsigned a = my_umul24(a1, b1) - my_umul24(l1, m) -
        my_umul24(a2, b2) + my_umul24(l2, m);

    reduce_mod(a, m);
    return a;
#endif  // CUMP_USE_32BIT_MODULI_SET
}

// computes a * b + c partially reduced
inline unsigned fma_mod(unsigned a, unsigned b, unsigned c, unsigned m) {

#if STRICT_MOD_CHECK
    if(a >= m || b >= m) {
        printf("fma_mod error: a: %x; b: %x; c: %x; m: %x\n",
                a, b, c, m);
        exit(1);
    }
#endif

    unsigned mullo = my_umul24(a, b), mulhi = my_umul24hi(a, b), l;
    unsigned k = (m >> 9);
    // TODO remove fftmod variable ? use only m instead
    float invk = 128.0f / k;
    float prodf = (float)mulhi;

    prodf = prodf * invk + 2.0f;
    l = (unsigned)floorf(prodf);

    // possibly fused into 2 mad24s
    unsigned resid = mullo - my_umul24(l, m) + c;
    return resid;
}

//! computes (a mod m) * b + c partially reduced
//! \c b must be already reduced !!
// look if there is smth to benefit from..
inline unsigned fma_reduce_mod(unsigned a, unsigned b, unsigned c,
        unsigned m) {

    reduce_mod(a, m);
    return fma_mod(a, b, c, m);
}

inline void fma_bfy2_mod(unsigned& x, unsigned &y, unsigned a, unsigned b,
        unsigned c, unsigned m) {

    x = fma_mod(a, b, c, m);
    y = 2 * c - x;
}

inline void fma_reduce_bfy2_mod(unsigned& x, unsigned &y, unsigned a,
        unsigned b, unsigned c, unsigned m) {

    x = fma_reduce_mod(a, b, c, m);

    if((int)c < (int)x)
        printf("Less: c: %d x: %d\n", c, x);
    //y = 2 * c - x;
    y = my_sad(x, c, c);
}

inline unsigned large_reduce_mod(unsigned hi, unsigned lo, unsigned m) {

    float prodf = float(hi);
    unsigned k = (m >> 9);
    float invk = 128.0f / k;

    prodf = prodf * invk;
    unsigned l = (unsigned)floorf(prodf);

    unsigned resid = lo - my_umul24(l, m);

    //TODO: you forgot post-reduction step
    if((int)resid < 0)
        resid = resid + m * 2;

    if((int)resid >= m)
        resid -= m;

    return resid;
}

//! multiplies 48-bit number [hi; lo] by 24-bit number \c a , returns 72-bit
//! product as [r2; r1; r0] having 16, 32 and 24 bits respectively
//! NOTE: r0's 8 MSB bits are not masked !!!!
inline void mul_48_24_wide_host(unsigned& r2, unsigned& r1, unsigned& r0,
        unsigned hi, unsigned lo, unsigned a) {

    unsigned h1;
    mul_24_wide_host(h1, r0, lo, a);
    h1 >>= 8; // keep only 24 MSB bits
    //r0 &= 0xffffff; // mask out 8 MSB bits that overlap with h1

    // [h3; l3] = m1m2_hi * c.x + h2
    mad_24_wide_host(r2, r1, hi, a, h1); // [r2;r1] = hi * a + h1
    r2 >>= 16, r2 += (r1 < h1);
}

//! same as before but result is: r2(8), r1(32), r0(32)
inline void mul_48_32_wide_host(unsigned& r2, unsigned& r1, unsigned& r0,
        unsigned hi, unsigned lo, unsigned a) {

    unsigned h1;
    mul_24_wide_host(h1, r0, lo, a);
    h1 >>= 8; // keep only 24 MSB bits

    // [h3; l3] = m1m2_hi * c.x + h2
    mad_24_wide_host(r2, r1, hi, a, h1); // [r2;r1] = hi * a + h1
    if(r1 < h1)
        r2 += 0x10000;

    r0 = (r0 & 0xffffff) + (r1 << 24);
    r1 = (r1 >> 8) + ((r2 & 0xff0000) << 8);
    r2 >>= 24;
}


inline umod_t pow_mod(umod_t a, unsigned e, umod_t m) {

    if(e == 0)  
        return 1;
    umod_t z = a;
    umod_t y = 1;
    while(1) {
        if(e & 1)  
            y = mul_mod(y, z, m);  // y *= z;
        e >>= 1;
        if(e == 0)  
            break;
        z = mul_mod(z, z, m);  
    }
    return  y;
}

// Return u3 and set u1,v1 so that
//   gcd(u,v) == u3 == u*u1 + v*u2
// Type must be a signed type.
template <typename Type>
Type egcd(Type u, Type v, Type &tu1, Type &tu2) {
    Type u1 = 1,  u2 = 0;
    Type v1 = 0,  v3 = v;
    Type u3 = u,  v2 = 1;
    while ( v3!=0 )
    {
        Type q = u3 / v3;

        Type t1 = u1 - v1 * q;
        u1 = v1;  v1 = t1;

        Type t3 = u3 - v3 * q;
        u3 = v3;  v3 = t3;

        Type t2 = u2 - v2 * q;
        u2 = v2;  v2 = t2;
    }
    tu1 = u1;  tu2 = u2;
    return u3;
}

struct zmod // arithmetic modulo word-sized prime
{
    umod_t x;

public:
    static umod_t modulus;        // 0 <= x < modulus
    
    zmod()  { ; }

    explicit zmod(const smod_t i) { 

//         g_inv_mod = 1.0 / modulus;
//         umod_t y = (umod_t)modulus *
//           (umod_t)((double)i * g_inv_mod + 0.5);
//         x = i - y;
//         if((smod_t)x < 0)
//             x += modulus;
        x = i;
    }

    zmod(const zmod &m) : x(m.x)  { ; }
    ~zmod()  { ; }

    static void set_mod(umod_t m) {
        modulus = m;
        g_inv_mod = 1.0 / m;
//         g_inv_mod_e16 = 65536.0 / m;
    }

    zmod pow(unsigned e)  const
    {
        zmod s;
        s.x = pow_mod(x, e, zmod::modulus);
        return  s;
    }

    inline zmod & operator = (const zmod &h)  { x = h.x;  return *this; }
    inline zmod & operator = (umod_t i)  { (*this) = zmod(i);  return *this; }
    
    friend inline zmod & operator ++ (zmod &z)
    { ++z.x; if ( z.x==zmod::modulus ) z.x=0;  return z; }

    friend inline zmod & operator -- (zmod &z)
    { if ( z.x==0 ) z.x=zmod::modulus;  z.x--;  return z; };

    friend inline zmod & operator += (zmod &z, const zmod &h)
    { z.x = add_mod(z.x, h.x, zmod::modulus);  return z; }

    friend inline zmod & operator -= (zmod &z, const zmod &h)
    { z.x = sub_mod(z.x, h.x, zmod::modulus);  return z; }

    friend inline zmod & operator *= (zmod &z, const zmod &h)
    { z.x = mul_mod(z.x, h.x, zmod::modulus);  return z; }

    friend inline zmod operator + (const zmod &h1, const zmod &h2)
    { zmod z(h1); z += h2; return z; }

    friend inline zmod operator - (const zmod &h1, const zmod &h2)
    { zmod z(h1); z -= h2; return z; }

    friend inline zmod operator * (const zmod &h1, const zmod &h2)
    { zmod z(h1); z *= h2; return z; }

    zmod & negate()  { x = neg_mod(x, zmod::modulus);
        return *this; }

    friend inline zmod operator - (const zmod &h)
    { zmod n(h);  n.negate();  return n; }

    friend inline zmod operator + (const zmod &h)  { return h; }

    friend inline bool operator == (const zmod &h1, const zmod &h2)
    { return  h1.x == h2.x; }

    friend inline bool operator != (const zmod &h1, const zmod &h2)
    { return  h1.x != h2.x; }

    friend inline bool operator < (const zmod &h1, const zmod &h2)
    { return  h1.x < h2.x; }

    friend inline bool operator <= (const zmod &h1, const zmod &h2)
    { return  h1.x <= h2.x; }

    friend inline bool operator > (const zmod &h1, const zmod &h2)
    { return  h1.x > h2.x; }

    friend inline bool operator >= (const zmod &h1, const zmod &h2)
    { return  h1.x >= h2.x; }
};

std::ostream& operator <<(std::ostream& os, const zmod& x);

} // namespace internal

} // namespace CGAL

#endif // ZMOD_H
