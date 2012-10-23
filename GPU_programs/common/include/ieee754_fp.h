/*
 * Copyright 1993-2008 NVIDIA Corporation.  All rights reserved.
 *
 * NOTICE TO USER:   
 *
 * This source code is subject to NVIDIA ownership rights under U.S. and 
 * international Copyright laws.  Users and possessors of this source code 
 * are hereby granted a nonexclusive, royalty-free license to use this code 
 * in individual and commercial software.
 *
 * NVIDIA MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE 
 * CODE FOR ANY PURPOSE.  IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR 
 * IMPLIED WARRANTY OF ANY KIND.  NVIDIA DISCLAIMS ALL WARRANTIES WITH 
 * REGARD TO THIS SOURCE CODE, INCLUDING ALL IMPLIED WARRANTIES OF 
 * MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.
 * IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL, 
 * OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS 
 * OF USE, DATA OR PROFITS,  WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE 
 * OR OTHER TORTIOUS ACTION,  ARISING OUT OF OR IN CONNECTION WITH THE USE 
 * OR PERFORMANCE OF THIS SOURCE CODE.  
 *
 * U.S. Government End Users.   This source code is a "commercial item" as 
 * that term is defined at  48 C.F.R. 2.101 (OCT 1995), consisting  of 
 * "commercial computer  software"  and "commercial computer software 
 * documentation" as such terms are  used in 48 C.F.R. 12.212 (SEPT 1995) 
 * and is provided to the U.S. Government only as a commercial end item.  
 * Consistent with 48 C.F.R.12.212 and 48 C.F.R. 227.7202-1 through 
 * 227.7202-4 (JUNE 1995), all U.S. Government End Users acquire the 
 * source code with only those rights set forth herein. 
 *
 * Any use of this source code in individual and commercial software must 
 * include, in the user documentation and internal comments to the code,
 * the above Disclaimer and U.S. Government End Users Notice.
 */

#if !defined(__2111111111111111__)
#define __2111111111111111__

#define __internal_clamp11(val, max, min, nan)                                         \
       if (val >= max) return max;                                                   \
       if (val <= min) return min

#define __device_func__(x) x


#define __FPU_CW_EXCEPTION_MASK__   (0x003f)
#define __FPU_CW_INVALID__          (0x0001)
#define __FPU_CW_DENORMAL__         (0x0002)
#define __FPU_CW_ZERODIVIDE__       (0x0004)
#define __FPU_CW_OVERFLOW__         (0x0008)
#define __FPU_CW_UNDERFLOW__        (0x0010)
#define __FPU_CW_INEXACT__          (0x0020)

#define __FPU_CW_PREC_MASK__        (0x0300)
#define __FPU_CW_PREC_SINGLE__      (0x0000)
#define __FPU_CW_PREC_DOUBLE__      (0x0200)
#define __FPU_CW_PREC_EXTENDED__    (0x0300)

#define __FPU_CW_ROUND_MASK__       (0x0c00)
#define __FPU_CW_ROUND_NEAR__       (0x0000)
#define __FPU_CW_ROUND_DOWN__       (0x0400)
#define __FPU_CW_ROUND_UP__         (0x0800)
#define __FPU_CW_ROUND_CHOP__       (0x0c00)

#define __FPU_CW_MASK_ALL__         (0x1f3f)


#define __SSE_CW_FLUSHZERO__        (0x8000)
    
#define __SSE_CW_ROUND_MASK__       (0x6000)
#define __SSE_CW_ROUND_NEAR__       (0x0000)
#define __SSE_CW_ROUND_DOWN__       (0x2000)
#define __SSE_CW_ROUND_UP__         (0x4000)
#define __SSE_CW_ROUND_CHOP__       (0x6000)

#define __SSE_CW_EXCEPTION_MASK__   (0x1f80)
#define __SSE_CW_PRECISION__        (0x1000)
#define __SSE_CW_UNDERFLOW__        (0x0800)
#define __SSE_CW_OVERFLOW__         (0x0400)
#define __SSE_CW_DIVIDEZERO__       (0x0200)
#define __SSE_CW_DENORMAL__         (0x0100)
#define __SSE_CW_INVALID__          (0x0080)
// not on all SSE machines
// #define __SSE_CW_DENORMALZERO__     (0x0040)

#define __SSE_CW_MASK_ALL__         (0xffc0)

#define __MOD_FPU_CW_DEFAULT__ \
    (__FPU_CW_EXCEPTION_MASK__ + __FPU_CW_PREC_DOUBLE__ + __FPU_CW_ROUND_CHOP__)

#define __MOD_SSE_CW_DEFAULT__ \
    (__SSE_CW_EXCEPTION_MASK__ + __SSE_CW_ROUND_CHOP__ + __SSE_CW_FLUSHZERO__)

#ifdef USE_SSE
inline unsigned int getSSEStateX86(void);
inline void setSSEModDefault(unsigned int control);
inline void modifySSEStateX86(unsigned int control, unsigned int mask);
#endif // USE_SSE

inline void setRoundingMode(unsigned int round);

inline unsigned int getFPUStateX86(void);
inline void setFPUStateX86(unsigned int control);

inline void setFPUModDefault(void);
inline void assertFPUModDefault(void);

inline void modifyFPUStateX86(const unsigned int control, const unsigned int mask);

inline int FastFtol(const float a);

// assume for now that we are running on an x86
// #ifdef __i386__

#ifdef USE_SSE

inline
unsigned int
getSSEStateX86
    (void)
{
    return _mm_getcsr();
}

inline
void
setSSEStateX86
    (unsigned int control)
{
    _mm_setcsr(control);
}


inline
void
modifySSEStateX86
    (unsigned int control, unsigned int mask)
{
    unsigned int oldControl = getFPUStateX86();
    unsigned int newControl = ((oldControl & (~mask)) | (control & mask));
    setFPUStateX86(newControl);
}


inline
void
setSSEModDefault
    (void)
{
    modifySSEStateX86(__MOD_SSE_CW_DEFAULT__, __SSE_CW_MASK_ALL__);
}    
#endif // USE_SSE


inline
void
setRoundingMode
    (unsigned round)
{
    unsigned mask = 0x3;

    unsigned fpuControl = getFPUStateX86();
    fpuControl &= ~(mask << 10);
    fpuControl |= round << 10;
    setFPUStateX86(fpuControl);

#ifdef USE_SSE
    Uint32 sseControl = getSSEStateX86();
    sseControl &= ~(mask << 13);
    sseControl |= round << 13;
    setSSEStateX86(sseControl);
#endif // USE_SSE
}


inline
unsigned int
getFPUStateX86
    (void)
{
    unsigned int control = 0;
#if defined(_MSVC)
    __asm fnstcw control;
#elif defined(__GNUG__)
    __asm__ __volatile__ ("fnstcw %0" : "=m" (control));
#endif
    return control;
}

/* set status */
inline
void
setFPUStateX86
    (unsigned int control)
{
#if defined(_MSVC)
    __asm fldcw control;
#elif defined(__GNUG__)
    __asm__ __volatile__ ("fldcw %0" : : "m" (control));
#endif
}


inline
void
modifyFPUStateX86
    (const unsigned int control, const unsigned int mask)
{
    unsigned int oldControl = getFPUStateX86();
    unsigned int newControl = ((oldControl & (~mask)) | (control & mask));
    setFPUStateX86(newControl);
}


inline
void
setFPUModDefault
    (void)
{
    modifyFPUStateX86(__MOD_FPU_CW_DEFAULT__, __FPU_CW_MASK_ALL__);
    assertFPUModDefault();
}


inline
void
assertFPUModDefault
    (void)
{
//     assert((getFPUStateX86() & (__FPU_CW_MASK_ALL__)) == 
//             __MOD_FPU_CW_DEFAULT__);
}

// taken from http://www.stereopsis.com/FPU.html
// this assumes the CPU is in double precision mode
inline
int
FastFtol(const float a)
{
    static int    b;
    
#if defined(_MSVC)
    __asm fld a
    __asm fistp b
#elif defined(__GNUG__)
    // use AT&T inline assembly style, document that
    // we use memory as output (=m) and input (m)
    __asm__ __volatile__ (
        "flds %1        \n\t"
        "fistpl %0      \n\t"
        : "=m" (b)
        : "m" (a));
#endif
    return b;
}

/*******************************************************************************
*                                                                              *
* HOST IMPLEMENTATIONS FOR FUNCTIONS WITH BUILTIN NVOPENCC OPREATIONS          *
*                                                                              *
*******************************************************************************/

#ifndef DO_NOT_USE_THIS_STUFF

#include <math_constants.h>

namespace XXX {

__device_func__(int __mulhi(int a, int b))
{
  long long int c = (long long int)a * (long long int)b;

  return (int)(c >> 32);
}

__device_func__(unsigned int __umulhi(unsigned int a, unsigned int b))
{
  unsigned long long int c = (unsigned long long int)a * (unsigned long long int)b;

  return (unsigned int)(c >> 32);
}

__device_func__(unsigned long long int __umul64hi(unsigned long long int a, unsigned long long int b))
{
  unsigned int           a_lo = (unsigned int)a;
  unsigned long long int a_hi = a >> 32;
  unsigned int           b_lo = (unsigned int)b;
  unsigned long long int b_hi = b >> 32;
  unsigned long long int m1 = a_lo * b_hi;
  unsigned long long int m2 = a_hi * b_lo;
  unsigned int           carry;

  carry = (0ULL + __umulhi(a_lo, b_lo) + (unsigned int)m1 + (unsigned int)m2) >> 32;

  return a_hi * b_hi + (m1 >> 32) + (m2 >> 32) + carry;
}

__device_func__(long long int __mul64hi(long long int a, long long int b))
{
  long long int res;
  res = __umul64hi(a, b);
  if (a < 0LL) res = res - b;
  if (b < 0LL) res = res - a;
  return res;
}

__device_func__(float __saturatef(float a))
{
  return a >= 1.0f ? 1.0f : a <= 0.0f ? 0.0f : a;
}

__device_func__(int __mul24(int a, int b))
{
#if !defined(__MULTI_CORE__)
  a &= 0xffffff;
  a = (a & 0x800000) != 0 ? a | ~0xffffff : a;
  b &= 0xffffff;
  b = (b & 0x800000) != 0 ? b | ~0xffffff : b;
#endif /* !__MULTI_CORE__ */

  return a * b;
}

__device_func__(unsigned int __umul24(unsigned int a, unsigned int b))
{
#if !defined(__MULTI_CORE__)
  a &= 0xffffff;
  b &= 0xffffff;
#endif /* !__MULTI_CORE__ */

  return a * b;
}

__device_func__(unsigned int __umul24hi(unsigned int a, unsigned int b))
{
#if !defined(__MULTI_CORE__)
  a &= 0xffffff;
  b &= 0xffffff;
#endif /* !__MULTI_CORE__ */

   unsigned long long int res =
            (unsigned long long int)(a) *
            (unsigned long long int)(b);
   return (unsigned int)(res >> 16); // drop 16 LSBs
}

__device_func__(float __int_as_float(int a))
{
  volatile union {int a; float b;} u;

  u.a = a;

  return u.b;
}

__device_func__(int __float_as_int(float a))
{
  volatile union {float a; int b;} u;

  u.a = a;

  return u.b;
}

__device_func__(long long int __internal_float2ll_kernel(float a, long long int max, long long int min, long long int nan, enum cudaRoundMode rndMode))
{
  unsigned long long int res, t = 0ULL;
  int shift;
  unsigned int ia;

  __internal_clamp11(a, max, min, nan);
  ia = __float_as_int(a);
  shift = 189 - ((ia >> 23) & 0xff);
  res = (unsigned long long int)(((ia << 8) | 0x80000000) >> 1) << 32;
  if (shift >= 64) {
    t = res;
    res = 0;
  } else if (shift) {
    t = res << (64 - shift);
    res = res >> shift;
  }
  if (rndMode == cudaRoundNearest && (long long int)t < 0LL) {
    res += t == 0x8000000000000000ULL ? res & 1ULL : 1ULL;
  }
  else if (rndMode == cudaRoundMinInf && t != 0ULL && ia > 0x80000000) {
    res++;
  }
  else if (rndMode == cudaRoundPosInf && t != 0ULL && (int)ia > 0) {
    res++;
  }
  if ((int)ia < 0) res = (unsigned long long int)-(long long int)res;
  return (long long int)res;
}

__device_func__(int __internal_float2int(float a, enum cudaRoundMode rndMode))
{
  return (int)__internal_float2ll_kernel(a, 2147483647LL, -2147483648LL, 0LL, rndMode);
}

__device_func__(int __float2int_rz(float a))
{
#if defined(__MULTI_CORE__)
  return (int)a;
#else /* __MULTI_CORE__ */
  return __internal_float2int(a, cudaRoundZero);
#endif /* __MULTI_CORE__ */
}

__device_func__(int __float2int_ru(float a))
{
  return __internal_float2int(a, cudaRoundPosInf);
}

__device_func__(int __float2int_rd(float a))
{
  return __internal_float2int(a, cudaRoundMinInf);
}

__device_func__(int __float2int_rn(float a))
{
  return __internal_float2int(a, cudaRoundNearest);
}

__device_func__(long long int __internal_float2ll(float a, enum cudaRoundMode rndMode))
{
  return __internal_float2ll_kernel(a, 9223372036854775807LL, -9223372036854775807LL -1LL, -9223372036854775807LL -1LL, rndMode);
}

__device_func__(long long int __float2ll_rz(float a))
{
#if defined(__MULTI_CORE__)
  return (long long int)a;
#else /* __MULTI_CORE__ */
  return __internal_float2ll(a, cudaRoundZero);
#endif /* __MULTI_CORE__ */
}

__device_func__(long long int __float2ll_ru(float a))
{
  return __internal_float2ll(a, cudaRoundPosInf);
}

__device_func__(long long int __float2ll_rd(float a))
{
  return __internal_float2ll(a, cudaRoundMinInf);
}

__device_func__(long long int __float2ll_rn(float a))
{
  return __internal_float2ll(a, cudaRoundNearest);
}

__device_func__(unsigned long long int __internal_float2ull_kernel(float a, unsigned long long int max, unsigned long long int nan, enum cudaRoundMode rndMode))
{
  unsigned long long int res, t = 0ULL;
  int shift;
  unsigned int ia;

  __internal_clamp11(a, max, 0LL, nan);
  ia = __float_as_int(a);
    // this is an exponent
  shift = 190 - ((ia >> 23) & 0xff);
    // this is mantissa with set hidden bit
  res = (unsigned long long int)((ia << 8) | 0x80000000) << 32;
  if (shift >= 64) {
    t = res >> (int)(shift > 64);
    res = 0;
  } else if (shift) {
    t = res << (64 - shift); // t is the LSB of res
    res = res >> shift;
  }
  if (rndMode == cudaRoundNearest && (long long int)t < 0LL) {
    res += t == 0x8000000000000000ULL ? res & 1ULL : 1ULL;
  }
  else if (rndMode == cudaRoundPosInf && t != 0ULL) {
    res++;
  }
  return res;
}

__device_func__(unsigned int __internal_float2uint(float a, enum cudaRoundMode rndMode))
{
  return (unsigned int)__internal_float2ull_kernel(a, 4294967295U, 0U, rndMode);
}

__device_func__(unsigned int __float2uint_rz(float a))
{
#if defined(__MULTI_CORE__)
  return (unsigned int)a;
#else /* __MULTI_CORE__ */
  return __internal_float2uint(a, cudaRoundZero);
#endif /* __MULTI_CORE__ */
}

__device_func__(unsigned int __float2uint_ru(float a))
{
  return __internal_float2uint(a, cudaRoundPosInf);
}

__device_func__(unsigned int __float2uint_rd(float a))
{
  return __internal_float2uint(a, cudaRoundMinInf);
}

__device_func__(unsigned int __float2uint_rn(float a))
{
  return __internal_float2uint(a, cudaRoundNearest);
}

__device_func__(unsigned long long int __internal_float2ull(float a, enum cudaRoundMode rndMode))
{
  return __internal_float2ull_kernel(a, 18446744073709551615ULL, 9223372036854775808ULL, rndMode);
}

__device_func__(unsigned long long int __float2ull_rz(float a))
{
#if defined(__MULTI_CORE__)
  return (unsigned long long int)a;
#else /* __MULTI_CORE__ */
  return __internal_float2ull(a, cudaRoundZero);
#endif /* __MULTI_CORE__ */
}

__device_func__(unsigned long long int __float2ull_ru(float a))
{
  return __internal_float2ull(a, cudaRoundPosInf);
}

__device_func__(unsigned long long int __float2ull_rd(float a))
{
  return __internal_float2ull(a, cudaRoundMinInf);
}

__device_func__(unsigned long long int __float2ull_rn(float a))
{
  return __internal_float2ull(a, cudaRoundNearest);
}

__device_func__(int __internal_normalize64(unsigned long long int *a))
{
  int lz = 0;

  if ((*a & 0xffffffff00000000ULL) == 0ULL) {
    *a <<= 32;
    lz += 32;
  }
  if ((*a & 0xffff000000000000ULL) == 0ULL) {
    *a <<= 16;
    lz += 16;
  }
  if ((*a & 0xff00000000000000ULL) == 0ULL) {
    *a <<= 8;
    lz += 8;
  }
  if ((*a & 0xf000000000000000ULL) == 0ULL) {
    *a <<= 4;
    lz += 4;
  }
  if ((*a & 0xC000000000000000ULL) == 0ULL) {
    *a <<= 2;
    lz += 2;
  }
  if ((*a & 0x8000000000000000ULL) == 0ULL) {
    *a <<= 1;
    lz += 1;
  }
  return lz;
}

__device_func__(int __internal_normalize(unsigned int *a))
{
  unsigned long long int t = (unsigned long long int)*a;

  int lz = __internal_normalize64(&t);

  // shifts out all zeros to the right 
  // returns # of leading zeros  
  *a = (unsigned int)(t >> 32);

  return lz - 32;
}

__device_func__(float __internal_int2float_kernel(int a, enum cudaRoundMode rndMode))
{
  volatile union {
    float f;
    unsigned int i;
  } res;
  int shift;
  unsigned int t;
  res.i = a;
  if (a == 0) return res.f;
  if (a < 0) res.i = (unsigned int)-a;
  shift = __internal_normalize((unsigned int*)&res.i);
  t = res.i << 24;
  res.i = (res.i >> 8);
  res.i += (127 + 30 - shift) << 23;
  if (a < 0) res.i |= 0x80000000;
  if ((rndMode == cudaRoundNearest) && (t >= 0x80000000)) {
    res.i += (t == 0x80000000) ? (res.i & 1) : (t >> 31);
  }
  else if ((rndMode == cudaRoundMinInf) && t && (a < 0)) {
    res.i++;
  }
  else if ((rndMode == cudaRoundPosInf) && t && (a > 0)) {
    res.i++;
  }
  return res.f;
}

__device_func__(float __int2float_rz(int a))
{
  return __internal_int2float_kernel(a, cudaRoundZero);
}

__device_func__(float __int2float_ru(int a))
{
  return __internal_int2float_kernel(a, cudaRoundPosInf);
}

__device_func__(float __int2float_rd(int a))
{
  return __internal_int2float_kernel(a, cudaRoundMinInf);
}

__device_func__(float __int2float_rn(int a))
{
#if defined(__MULTI_CORE__)
  return (float)a;
#else /* __MULTI_CORE__ */
  return __internal_int2float_kernel(a, cudaRoundNearest);
#endif /* __MULTI_CORE__ */
}

__device_func__(float __internal_uint2float_kernel(unsigned int a, enum cudaRoundMode rndMode))
{
  volatile union {
    float f;
    unsigned int i;
  } res;
  int shift;
  unsigned int t;
  res.i = a;
  if (a == 0) return res.f;
  shift = __internal_normalize((unsigned int*)&res.i);

// shift = 0;
//     printf("shift = %d\n", shift);
  t = res.i << 24;
  res.i = (res.i >> 8);
  res.i += (127 + 30 - shift) << 23;
  if ((rndMode == cudaRoundNearest) && (t >= 0x80000000)) {
    res.i += (t == 0x80000000) ? (res.i & 1) : (t >> 31);
  }
  else if ((rndMode == cudaRoundPosInf) && t) {
    res.i++;
  }
  return res.f;
}

__device_func__(float __uint2float_rz(unsigned int a))
{
  return __internal_uint2float_kernel(a, cudaRoundZero);
}

__device_func__(float __uint2float_ru(unsigned int a))
{
  return __internal_uint2float_kernel(a, cudaRoundPosInf);
}

__device_func__(float __uint2float_rd(unsigned int a))
{
  return __internal_uint2float_kernel(a, cudaRoundMinInf);
}

__device_func__(float __uint2float_rn(unsigned int a))
{
#if defined(__MULTI_CORE__)
  return (float)a;
#else /* __MULTI_CORE__ */
  return __internal_uint2float_kernel(a, cudaRoundNearest);
#endif /* __MULTI_CORE__ */
}

__device_func__(float __ll2float_rn(long long int a))
{
  return (float)a;
}      

__device_func__(float __ull2float_rn(unsigned long long int a))
{
#if defined(__MULTI_CORE__)
  return (float)a;
#else /* __MULTI_CORE__ */
  unsigned long long int temp;
  unsigned int res, t;
  int shift;
  if (a == 0ULL) return 0.0f;
  temp = a;
  shift = __internal_normalize64(&temp);
  temp = (temp >> 8) | ((temp & 0xffULL) ? 1ULL : 0ULL);
  res = (unsigned int)(temp >> 32);
  t = (unsigned int)temp;
  res += (127 + 62 - shift) << 23; /* add in exponent */
  res += t == 0x80000000 ? res & 1 : t >> 31;
  return __int_as_float(res);
#endif /* __MULTI_CORE__ */
}      

__device_func__(float __internal_fmul_kernel(float a, float b, int rndNearest))
{
  unsigned long long product;
  volatile union {
    float f;
    unsigned int i;
  } xx, yy;
  unsigned expo_x, expo_y;
    
  xx.f = a;
  yy.f = b;

  expo_y = 0xFF;
  expo_x = expo_y & (xx.i >> 23);
  expo_x = expo_x - 1;
  expo_y = expo_y & (yy.i >> 23);
  expo_y = expo_y - 1;
    
  if ((expo_x <= 0xFD) && 
      (expo_y <= 0xFD)) {
multiply:
    expo_x = expo_x + expo_y;
    expo_y = xx.i ^ yy.i;
    xx.i = xx.i & 0x00ffffff;
    yy.i = yy.i << 8;
    xx.i = xx.i | 0x00800000;
    yy.i = yy.i | 0x80000000;
    /* compute product */
    product = ((unsigned long long)xx.i) * yy.i;
    expo_x = expo_x - 127 + 2;
    expo_y = expo_y & 0x80000000;
    xx.i = (unsigned int)(product >> 32);
    yy.i = (unsigned int)(product & 0xffffffff);
    /* normalize mantissa */
    if (xx.i < 0x00800000) {
      xx.i = (xx.i << 1) | (yy.i >> 31);
      yy.i = (yy.i << 1);
      expo_x--;
    }
    if (expo_x <= 0xFD) {
      xx.i = xx.i | expo_y;          /* OR in sign bit */
      xx.i = xx.i + (expo_x << 23);  /* add in exponent */
      /* round result to nearest or even */
      if (yy.i < 0x80000000) return xx.f;
      xx.i += (((yy.i == 0x80000000) ? (xx.i & 1) : (yy.i >> 31)) 
               && rndNearest);
      return xx.f;
    } else if ((int)expo_x >= 254) {
      /* overflow: return infinity */
      xx.i = (expo_y | 0x7F800000) - (!rndNearest);
      return xx.f;
    } else {
      /* zero, denormal, or smallest normal */
      expo_x = ((unsigned int)-((int)expo_x));
      if (expo_x > 25) {
        /* massive underflow: return 0 */
        xx.i = expo_y;
        return xx.f;
      } else {
        yy.i = (xx.i << (32 - expo_x)) | ((yy.i) ? 1 : 0);
        xx.i = expo_y + (xx.i >> expo_x);
        xx.i += (((yy.i == 0x80000000) ? (xx.i & 1) : (yy.i >> 31)) 
                 && rndNearest);
        return xx.f;
      }
    }
  } else {
    product = xx.i ^ yy.i;
    product = product & 0x80000000;
    if (!(xx.i & 0x7fffffff)) {
      if (expo_y != 254) {
        xx.i = (unsigned int)product;
        return xx.f;
      }
      expo_y = yy.i << 1;
      if (expo_y == 0xFF000000) {
        xx.i = expo_y | 0x00C00000;
      } else {
        xx.i = yy.i | 0x00400000;
      }
      return xx.f;
    }
    if (!(yy.i & 0x7fffffff)) {
      if (expo_x != 254) {
        xx.i = (unsigned int)product;
        return xx.f;
      }
      expo_x = xx.i << 1;
      if (expo_x == 0xFF000000) {
        xx.i = expo_x | 0x00C00000;
      } else {
        xx.i = xx.i | 0x00400000;
      }
      return xx.f;
    }
    if ((expo_y != 254) && (expo_x != 254)) {
      expo_y++;
      expo_x++;
      if (expo_x == 0) {
        expo_y |= xx.i & 0x80000000;
        /*
         * If both operands are denormals, we only need to normalize 
         * one of them as the result will be either a denormal or zero.
         */
        xx.i = xx.i << 8;
        while (!(xx.i & 0x80000000)) {
          xx.i <<= 1;
          expo_x--;
        }
        xx.i = (xx.i >> 8) | (expo_y & 0x80000000);
        expo_y &= ~0x80000000;
        expo_y--;
        goto multiply;
      }
      if (expo_y == 0) {
        expo_x |= yy.i & 0x80000000;
        yy.i = yy.i << 8;
        while (!(yy.i & 0x80000000)) {
          yy.i <<= 1;
          expo_y--;
        }
        yy.i = (yy.i >> 8) | (expo_x & 0x80000000);
        expo_x &= ~0x80000000;
        expo_x--;
        goto multiply;
      }
    }
    expo_x = xx.i << 1;
    expo_y = yy.i << 1;
    /* if x is NaN, return x */
    if (expo_x > 0xFF000000) {
      /* cvt any SNaNs to QNaNs */
      xx.i = xx.i | 0x00400000;
      return xx.f;
    }
    /* if y is NaN, return y */
    if (expo_y > 0xFF000000) {
      /* cvt any SNaNs to QNaNs */
      xx.i = yy.i | 0x00400000;
      return xx.f;
    } 
    xx.i = (unsigned int)product | 0x7f800000;
    return xx.f;
  }
}

__device_func__(float __internal_fadd_kernel(float a, float b, int rndNearest))
{
  volatile union {
    float f;
    unsigned int i;
  } xx, yy;
  unsigned int expo_x;
  unsigned int expo_y;
  unsigned int temp;

  xx.f = a;
  yy.f = b;

  /* make bigger operand the augend */
  expo_y = yy.i << 1;
  if (expo_y > (xx.i << 1)) {
    expo_y = xx.i;
    xx.i   = yy.i;
    yy.i   = expo_y;
  }
    
  temp = 0xff;
  expo_x = temp & (xx.i >> 23);
  expo_x = expo_x - 1;
  expo_y = temp & (yy.i >> 23);
  expo_y = expo_y - 1;
    
  if ((expo_x <= 0xFD) && 
      (expo_y <= 0xFD)) {
        
add:
    expo_y = expo_x - expo_y;
    if (expo_y > 25) {
      expo_y = 31;
    }
    temp = xx.i ^ yy.i;
    xx.i = xx.i & ~0x7f000000;
    xx.i = xx.i |  0x00800000;
    yy.i = yy.i & ~0xff000000;
    yy.i = yy.i |  0x00800000;
        
    if ((int)temp < 0) {
      /* signs differ, effective subtraction */
      temp = 32 - expo_y;
      temp = (expo_y) ? (yy.i << temp) : 0;
      temp = (unsigned int)(-((int)temp));
      xx.i = xx.i - (yy.i >> expo_y) - (temp ? 1 : 0);
      if (xx.i & 0x00800000) {
        if (expo_x <= 0xFD) {
          xx.i = xx.i & ~0x00800000; /* lop off integer bit */
          xx.i = (xx.i + (expo_x << 23)) + 0x00800000;
          if (temp < 0x80000000) return xx.f;
          xx.i += (((temp == 0x80000000) ? (xx.i & 1) : (temp >> 31))
                   && rndNearest);
          return xx.f;
        }
      } else {
        if ((temp | (xx.i << 1)) == 0) {
          /* operands cancelled, resulting in a clean zero */
          xx.i = 0;
          return xx.f;
        }
        /* normalize result */
        yy.i = xx.i & 0x80000000;
        do {
          xx.i = (xx.i << 1) | (temp >> 31);
          temp <<= 1;
          expo_x--;
        } while (!(xx.i & 0x00800000));
        xx.i = xx.i | yy.i;
      }
    } else {
      /* signs are the same, effective addition */
      temp = 32 - expo_y;
      temp = (expo_y) ? (yy.i << temp) : 0;
      xx.i = xx.i + (yy.i >> expo_y);
      if (!(xx.i & 0x01000000)) {
        if (expo_x <= 0xFD) {
          expo_y = xx.i & 1;
          xx.i = xx.i + (expo_x << 23);
          if (temp < 0x80000000) return xx.f;
          xx.i += (((temp == 0x80000000) ? expo_y : (temp >> 31)) 
                   && rndNearest);
          return xx.f;
        }
      } else {
        /* normalize result */
        temp = (xx.i << 31) | (temp >> 1);
        /* not ANSI compliant: xx.i = (((int)xx.i)>>1) & ~0x40000000 */
        xx.i = ((xx.i & 0x80000000) | (xx.i >> 1)) & ~0x40000000;
        expo_x++;
      }
    }
    if (expo_x <= 0xFD) {
      expo_y = xx.i & 1;
      xx.i += (((temp == 0x80000000) ? expo_y : (temp >> 31)) 
               && rndNearest);
      xx.i = xx.i + (expo_x << 23);
      return xx.f;
    }
    if ((int)expo_x >= 254) {
      /* overflow: return infinity */
        xx.i = ((xx.i & 0x80000000) | 0x7f800000) - (!rndNearest);
        return xx.f;
    }
    /* underflow: denormal, or smallest normal */
    expo_y = expo_x + 32;
    yy.i = xx.i &  0x80000000;
    xx.i = xx.i & ~0xff000000;
        
    expo_x = (unsigned int)(-((int)expo_x));
    temp = xx.i << expo_y | ((temp) ? 1 : 0);
    xx.i = yy.i | (xx.i >> expo_x);
    xx.i += (((temp == 0x80000000) ? (xx.i & 1) : (temp >> 31)) 
             && rndNearest);
    return xx.f;
  } else {
    /* handle special cases separately */
    if (!(yy.i << 1)) {
      if (xx.i == 0x80000000) {
          xx.i = yy.i;
      }
      if ((xx.i << 1) > 0xff000000) {
          xx.i |= 0x00400000;
      }
      return xx.f;
    }
    if ((expo_y != 254) && (expo_x != 254)) {
      /* remove sign bits */
      if (expo_x == (unsigned int) -1) {
        temp = xx.i & 0x80000000;
        xx.i = xx.i << 8;
        while (!(xx.i & 0x80000000)) {
          xx.i <<= 1;
          expo_x--;
        }
        expo_x++;
        xx.i = (xx.i >> 8) | temp;
      }
      if (expo_y == (unsigned int) -1) {
        temp = yy.i & 0x80000000;
        yy.i = yy.i << 8;
        while (!(yy.i & 0x80000000)) {
          yy.i <<= 1;
          expo_y--;
        }
        expo_y++;
        yy.i = (yy.i >> 8) | temp;
      }
      goto add;
    }
    expo_x = xx.i << 1;
    expo_y = yy.i << 1;
    /* if x is NaN, return x */
    if (expo_x > 0xff000000) {
      /* cvt any SNaNs to QNaNs */
      xx.i = xx.i | 0x00400000;
      return xx.f;
    }
    /* if y is NaN, return y */
    if (expo_y > 0xff000000) {
      /* cvt any SNaNs to QNaNs */
      xx.i = yy.i | 0x00400000;
      return xx.f;
    }
    if ((expo_x == 0xff000000) && (expo_y == 0xff000000)) {
      /*
       * subtraction of infinities with the same sign, and addition of
       * infinities of unlike sign is undefined: return NaN INDEFINITE
       */
      expo_x = xx.i ^ yy.i;
      xx.i = xx.i | ((expo_x) ? 0xffc00000 : 0);
      return xx.f;
    }
    /* handle infinities */
    if (expo_y == 0xff000000) {
      xx.i = yy.i;
    }
    return xx.f;
  }
}

__device_func__(float __fadd_rz(float a, float b))
{
  return __internal_fadd_kernel(a, b, 0);
}

__device_func__(float __fmul_rz(float a, float b))
{
  return __internal_fmul_kernel(a, b, 0);
}

__device_func__(float __fadd_rn(float a, float b))
{
  return __internal_fadd_kernel(a, b, 1);
}

__device_func__(float __fmul_rn(float a, float b))
{
  return __internal_fmul_kernel(a, b, 1);
}

__device_func__(void __brkpt(int c))
{
  /* TODO */
}

/*******************************************************************************
*                                                                              *
* DEVICE IMPLEMENTATIONS FOR FUNCTIONS WITH BUILTIN NVOPENCC OPERATIONS        *
*                                                                              *
*******************************************************************************/

__device_func__(float __sinf(float a))
{
  return sinf(a);
}

__device_func__(float __cosf(float a))
{
  return cosf(a);
}

__device_func__(float __log2f(float a))
{
  return log2f(a);
}

/*******************************************************************************
*                                                                              *
* SHARED HOST AND DEVICE IMPLEMENTATIONS                                       *
*                                                                              *
*******************************************************************************/


__device_func__(float __tanf(float a))
{
#if defined(__MULTI_CORE__)
  return tanf(a);
#else /* __MULTI_CORE__ */
  return __sinf(a) / __cosf(a);
#endif /* __MULTI_CORE__ */
}

__device_func__(void __sincosf(float a, float *sptr, float *cptr))
{
#if defined(__MULTI_CORE__)
  sincosf(a, sptr, cptr);
#else /* __MULTI_CORE__ */
  *sptr = __sinf(a);
  *cptr = __cosf(a);
#endif /* __MULTI_CORE__ */
}

__device_func__(int __clz(int a))
{
  return (a)?(158-(__float_as_int(__uint2float_rz((unsigned int)a))>>23)):32;
}

__device_func__(int __ffs(int a))
{
  return 32 - __clz (a & -a);
}

__device_func__(int __popc(unsigned int a))
{
  a = a - ((a >> 1) & 0x55555555);
  a = (a & 0x33333333) + ((a >> 2) & 0x33333333);
  a = (a + (a >> 4)) & 0x0f0f0f0f;
  a = ((__umul24(a, 0x808080) << 1) + a) >> 24;
  return a;
}

__device_func__(int __clzll(long long int a))
{
  int ahi = ((int)((unsigned long long)a >> 32));
  int alo = ((int)((unsigned long long)a & 0xffffffffULL));
  int res;
  if (ahi) {
      res = 0;
  } else {
      res = 32;
      ahi = alo;
  }
  res = res + __clz(ahi);
  return res;
}

__device_func__(int __ffsll(long long int a))
{
  return 64 - __clzll (a & -a);
}

__device_func__(int __popcll(unsigned long long int a))
{
  unsigned int ahi = ((unsigned int)(a >> 32));
  unsigned int alo = ((unsigned int)(a & 0xffffffffULL));
  alo = alo - ((alo >> 1) & 0x55555555);
  alo = (alo & 0x33333333) + ((alo >> 2) & 0x33333333);
  ahi = ahi - ((ahi >> 1) & 0x55555555);
  ahi = (ahi & 0x33333333) + ((ahi >> 2) & 0x33333333);
  alo = alo + ahi;
  alo = (alo & 0x0f0f0f0f) + ((alo >> 4) & 0x0f0f0f0f);
  alo = ((__umul24(alo, 0x808080) << 1) + alo) >> 24;
  return alo;
}

__device_func__(int __double2int_rz(double a))
{
  return __float2int_rz((float)a);
}

__device_func__(unsigned int __double2uint_rz(double a))
{
  return __float2uint_rz((float)a);
}

__device_func__(long long int __double2ll_rz(double a))
{
  return __float2ll_rz((float)a);
}

__device_func__(unsigned long long int __double2ull_rz(double a))
{
  return __float2ull_rz((float)a);
}


} // namespace XXX

#endif // DO_NOT_USE_THIS_STUFF

#endif /* !__DEVICE_FUNCTIONS_H__ */
