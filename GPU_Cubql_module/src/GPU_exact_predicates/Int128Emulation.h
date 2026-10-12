#pragma once
// ----------------------------------------------------------------------------
// Portable 128-bit integer emulation for host AND device code.
//
// Used when the host compiler has no native __int128 (MSVC). nvcc cannot
// provide __int128 on top of MSVC either, so this implements it by hand.
// Semantics follow two's-complement __int128: arithmetic wraps modulo 2^128,
// signed/unsigned mixes are done in unsigned arithmetic (like the builtin),
// right shift of a signed value is arithmetic.
// ----------------------------------------------------------------------------
#include <cuda_runtime.h>
#include <cstdint>
#include <type_traits>

#define I128_HD __host__ __device__ inline

struct int128_t;

// high 64 bits of the 128-bit product a*b
I128_HD uint64_t i128_umulhi(uint64_t a, uint64_t b) {
#if defined(__CUDA_ARCH__)
  return __umul64hi(a, b);
#else
  const uint64_t a0 = (uint32_t)a, a1 = a >> 32;
  const uint64_t b0 = (uint32_t)b, b1 = b >> 32;
  const uint64_t p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0, p11 = a1 * b1;
  const uint64_t mid = (p00 >> 32) + (uint32_t)p01 + (uint32_t)p10;
  return p11 + (p01 >> 32) + (p10 >> 32) + (mid >> 32);
#endif
}

struct uint128_t {
  uint64_t lo;
  uint64_t hi;

  uint128_t() = default;
  // plain constructors: signed sources are sign-extended (as for __int128)
  I128_HD uint128_t(int v)                : lo((uint64_t)(int64_t)v), hi(v < 0 ? ~0ull : 0ull) {}
  I128_HD uint128_t(long v)               : lo((uint64_t)(int64_t)v), hi(v < 0 ? ~0ull : 0ull) {}
  I128_HD uint128_t(long long v)          : lo((uint64_t)v),          hi(v < 0 ? ~0ull : 0ull) {}
  I128_HD uint128_t(unsigned v)           : lo(v), hi(0) {}
  I128_HD uint128_t(unsigned long v)      : lo(v), hi(0) {}
  I128_HD uint128_t(unsigned long long v) : lo(v), hi(0) {}
  I128_HD uint128_t(const int128_t& v);   // implicit, defined below
  I128_HD explicit uint128_t(double d) {  // truncating, d >= 0 expected
    if (d < 0.0) d = -d;
    const double two64 = 18446744073709551616.0;
    const double h = d / two64;
    hi = (uint64_t)h;
    lo = (uint64_t)(d - (double)hi * two64);
  }
  I128_HD explicit uint128_t(float f) : uint128_t((double)f) {}

  // explicit conversions to any arithmetic type
  template <class T, class = typename std::enable_if<std::is_arithmetic<T>::value>::type>
  I128_HD explicit operator T() const {
    if constexpr (std::is_same<T, bool>::value) return (lo | hi) != 0;
    else if constexpr (std::is_floating_point<T>::value)
      return (T)((double)hi * 18446744073709551616.0 + (double)lo);
    else return (T)lo;
  }

  // ---- bitwise / arithmetic -------------------------------------------------
  I128_HD uint128_t operator~() const { uint128_t r; r.lo = ~lo; r.hi = ~hi; return r; }
  I128_HD uint128_t operator-() const { uint128_t r; r.lo = 0 - lo; r.hi = 0 - hi - (lo != 0 ? 1 : 0); return r; }
  I128_HD uint128_t operator+() const { return *this; }
  I128_HD bool operator!() const { return (lo | hi) == 0; }

  friend I128_HD uint128_t operator+(const uint128_t& a, const uint128_t& b) {
    uint128_t r; r.lo = a.lo + b.lo; r.hi = a.hi + b.hi + (r.lo < a.lo ? 1 : 0); return r;
  }
  friend I128_HD uint128_t operator-(const uint128_t& a, const uint128_t& b) {
    uint128_t r; r.lo = a.lo - b.lo; r.hi = a.hi - b.hi - (a.lo < b.lo ? 1 : 0); return r;
  }
  friend I128_HD uint128_t operator*(const uint128_t& a, const uint128_t& b) {
    uint128_t r;
    r.lo = a.lo * b.lo;
    r.hi = i128_umulhi(a.lo, b.lo) + a.lo * b.hi + a.hi * b.lo;
    return r;
  }
  friend I128_HD uint128_t operator&(const uint128_t& a, const uint128_t& b) { uint128_t r; r.lo = a.lo & b.lo; r.hi = a.hi & b.hi; return r; }
  friend I128_HD uint128_t operator|(const uint128_t& a, const uint128_t& b) { uint128_t r; r.lo = a.lo | b.lo; r.hi = a.hi | b.hi; return r; }
  friend I128_HD uint128_t operator^(const uint128_t& a, const uint128_t& b) { uint128_t r; r.lo = a.lo ^ b.lo; r.hi = a.hi ^ b.hi; return r; }

  // shifts: count may be any integer-like type (including int128_t/uint128_t)
  template <class S> friend I128_HD uint128_t operator<<(const uint128_t& a, S cnt) {
    const int n = (int)cnt; uint128_t r;
    if (n <= 0)        { r = a; }
    else if (n >= 128) { r.lo = 0; r.hi = 0; }
    else if (n >= 64)  { r.lo = 0; r.hi = a.lo << (n - 64); }
    else               { r.lo = a.lo << n; r.hi = (a.hi << n) | (a.lo >> (64 - n)); }
    return r;
  }
  template <class S> friend I128_HD uint128_t operator>>(const uint128_t& a, S cnt) {
    const int n = (int)cnt; uint128_t r;
    if (n <= 0)        { r = a; }
    else if (n >= 128) { r.lo = 0; r.hi = 0; }
    else if (n >= 64)  { r.hi = 0; r.lo = a.hi >> (n - 64); }
    else               { r.hi = a.hi >> n; r.lo = (a.lo >> n) | (a.hi << (64 - n)); }
    return r;
  }

  // ---- comparisons ----------------------------------------------------------
  friend I128_HD bool operator==(const uint128_t& a, const uint128_t& b) { return a.lo == b.lo && a.hi == b.hi; }
  friend I128_HD bool operator!=(const uint128_t& a, const uint128_t& b) { return !(a == b); }
  friend I128_HD bool operator<(const uint128_t& a, const uint128_t& b) { return a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo); }
  friend I128_HD bool operator>(const uint128_t& a, const uint128_t& b) { return b < a; }
  friend I128_HD bool operator<=(const uint128_t& a, const uint128_t& b) { return !(b < a); }
  friend I128_HD bool operator>=(const uint128_t& a, const uint128_t& b) { return !(a < b); }

  // ---- division (shift/subtract; slow but exact, not used on hot paths) ----
  static I128_HD void divmod(const uint128_t& n, const uint128_t& d, uint128_t& q, uint128_t& r) {
    q.lo = q.hi = 0; r.lo = r.hi = 0;
    if (d.lo == 0 && d.hi == 0) return;  // division by zero: return 0
    for (int i = 127; i >= 0; --i) {
      r = r << 1;
      const uint64_t bit = (i >= 64) ? ((n.hi >> (i - 64)) & 1ull) : ((n.lo >> i) & 1ull);
      r.lo |= bit;
      if (!(r < d)) {
        r = r - d;
        if (i >= 64) q.hi |= (1ull << (i - 64)); else q.lo |= (1ull << i);
      }
    }
  }
  friend I128_HD uint128_t operator/(const uint128_t& a, const uint128_t& b) { uint128_t q, r; divmod(a, b, q, r); return q; }
  friend I128_HD uint128_t operator%(const uint128_t& a, const uint128_t& b) { uint128_t q, r; divmod(a, b, q, r); return r; }

  // ---- compound assignment --------------------------------------------------
  I128_HD uint128_t& operator+=(const uint128_t& o) { *this = *this + o; return *this; }
  I128_HD uint128_t& operator-=(const uint128_t& o) { *this = *this - o; return *this; }
  I128_HD uint128_t& operator*=(const uint128_t& o) { *this = *this * o; return *this; }
  I128_HD uint128_t& operator/=(const uint128_t& o) { *this = *this / o; return *this; }
  I128_HD uint128_t& operator%=(const uint128_t& o) { *this = *this % o; return *this; }
  I128_HD uint128_t& operator&=(const uint128_t& o) { *this = *this & o; return *this; }
  I128_HD uint128_t& operator|=(const uint128_t& o) { *this = *this | o; return *this; }
  I128_HD uint128_t& operator^=(const uint128_t& o) { *this = *this ^ o; return *this; }
  template <class S> I128_HD uint128_t& operator<<=(S n) { *this = *this << n; return *this; }
  template <class S> I128_HD uint128_t& operator>>=(S n) { *this = *this >> n; return *this; }
  I128_HD uint128_t& operator++() { *this = *this + uint128_t(1); return *this; }
  I128_HD uint128_t& operator--() { *this = *this - uint128_t(1); return *this; }
  I128_HD uint128_t operator++(int) { uint128_t t = *this; ++*this; return t; }
  I128_HD uint128_t operator--(int) { uint128_t t = *this; --*this; return t; }
};

struct int128_t {
  uint64_t lo;
  int64_t  hi;

  int128_t() = default;
  I128_HD int128_t(int v)                : lo((uint64_t)(int64_t)v), hi(v < 0 ? -1 : 0) {}
  I128_HD int128_t(long v)               : lo((uint64_t)(int64_t)v), hi(v < 0 ? -1 : 0) {}
  I128_HD int128_t(long long v)          : lo((uint64_t)v),          hi(v < 0 ? -1 : 0) {}
  I128_HD int128_t(unsigned v)           : lo(v), hi(0) {}
  I128_HD int128_t(unsigned long v)      : lo(v), hi(0) {}
  I128_HD int128_t(unsigned long long v) : lo(v), hi(0) {}
  // unsigned -> signed reinterpretation must be requested explicitly
  I128_HD explicit int128_t(const uint128_t& v) : lo(v.lo), hi((int64_t)v.hi) {}
  I128_HD explicit int128_t(double d) {  // truncates toward zero
    const bool neg = d < 0.0;
    uint128_t m(neg ? -d : d);
    if (neg) m = -m;
    lo = m.lo; hi = (int64_t)m.hi;
  }
  I128_HD explicit int128_t(float f) : int128_t((double)f) {}

  I128_HD uint128_t bits() const { uint128_t r; r.lo = lo; r.hi = (uint64_t)hi; return r; }
  I128_HD static int128_t fromBits(const uint128_t& u) { int128_t r; r.lo = u.lo; r.hi = (int64_t)u.hi; return r; }
  I128_HD bool isNeg() const { return hi < 0; }

  template <class T, class = typename std::enable_if<std::is_arithmetic<T>::value>::type>
  I128_HD explicit operator T() const {
    if constexpr (std::is_same<T, bool>::value) return (lo | (uint64_t)hi) != 0;
    else if constexpr (std::is_floating_point<T>::value) {
      if (hi < 0) { int128_t n = -*this; return (T)(-((double)(n.bits().hi) * 18446744073709551616.0 + (double)n.lo)); }
      return (T)((double)hi * 18446744073709551616.0 + (double)lo);
    }
    else return (T)lo;
  }

  I128_HD int128_t operator~() const { return fromBits(~bits()); }
  I128_HD int128_t operator-() const { return fromBits(-bits()); }
  I128_HD int128_t operator+() const { return *this; }
  I128_HD bool operator!() const { return (lo | (uint64_t)hi) == 0; }

  friend I128_HD int128_t operator+(const int128_t& a, const int128_t& b) { return fromBits(a.bits() + b.bits()); }
  friend I128_HD int128_t operator-(const int128_t& a, const int128_t& b) { return fromBits(a.bits() - b.bits()); }
  friend I128_HD int128_t operator*(const int128_t& a, const int128_t& b) { return fromBits(a.bits() * b.bits()); }
  friend I128_HD int128_t operator&(const int128_t& a, const int128_t& b) { return fromBits(a.bits() & b.bits()); }
  friend I128_HD int128_t operator|(const int128_t& a, const int128_t& b) { return fromBits(a.bits() | b.bits()); }
  friend I128_HD int128_t operator^(const int128_t& a, const int128_t& b) { return fromBits(a.bits() ^ b.bits()); }

  template <class S> friend I128_HD int128_t operator<<(const int128_t& a, S cnt) { return fromBits(a.bits() << cnt); }
  template <class S> friend I128_HD int128_t operator>>(const int128_t& a, S cnt) {  // arithmetic shift
    const int n = (int)cnt;
    if (n <= 0) return a;
    int128_t r;
    if (n >= 128)     { r.hi = a.hi < 0 ? -1 : 0; r.lo = (uint64_t)r.hi; }
    else if (n >= 64) { r.hi = a.hi < 0 ? -1 : 0; r.lo = (uint64_t)(a.hi >> (n - 64)); }
    else              { r.hi = a.hi >> n; r.lo = (a.lo >> n) | ((uint64_t)a.hi << (64 - n)); }
    return r;
  }

  friend I128_HD bool operator==(const int128_t& a, const int128_t& b) { return a.lo == b.lo && a.hi == b.hi; }
  friend I128_HD bool operator!=(const int128_t& a, const int128_t& b) { return !(a == b); }
  friend I128_HD bool operator<(const int128_t& a, const int128_t& b) { return a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo); }
  friend I128_HD bool operator>(const int128_t& a, const int128_t& b) { return b < a; }
  friend I128_HD bool operator<=(const int128_t& a, const int128_t& b) { return !(b < a); }
  friend I128_HD bool operator>=(const int128_t& a, const int128_t& b) { return !(a < b); }

  friend I128_HD int128_t operator/(const int128_t& a, const int128_t& b) {
    const bool na = a.isNeg(), nb = b.isNeg();
    uint128_t q = (na ? (-a).bits() : a.bits()) / (nb ? (-b).bits() : b.bits());
    return fromBits((na != nb) ? -q : q);
  }
  friend I128_HD int128_t operator%(const int128_t& a, const int128_t& b) {  // sign follows dividend
    const bool na = a.isNeg(), nb = b.isNeg();
    uint128_t r = (na ? (-a).bits() : a.bits()) % (nb ? (-b).bits() : b.bits());
    return fromBits(na ? -r : r);
  }

  I128_HD int128_t& operator+=(const int128_t& o) { *this = *this + o; return *this; }
  I128_HD int128_t& operator-=(const int128_t& o) { *this = *this - o; return *this; }
  I128_HD int128_t& operator*=(const int128_t& o) { *this = *this * o; return *this; }
  I128_HD int128_t& operator/=(const int128_t& o) { *this = *this / o; return *this; }
  I128_HD int128_t& operator%=(const int128_t& o) { *this = *this % o; return *this; }
  I128_HD int128_t& operator&=(const int128_t& o) { *this = *this & o; return *this; }
  I128_HD int128_t& operator|=(const int128_t& o) { *this = *this | o; return *this; }
  I128_HD int128_t& operator^=(const int128_t& o) { *this = *this ^ o; return *this; }
  template <class S> I128_HD int128_t& operator<<=(S n) { *this = *this << n; return *this; }
  template <class S> I128_HD int128_t& operator>>=(S n) { *this = *this >> n; return *this; }
  I128_HD int128_t& operator++() { *this = *this + int128_t(1); return *this; }
  I128_HD int128_t& operator--() { *this = *this - int128_t(1); return *this; }
  I128_HD int128_t operator++(int) { int128_t t = *this; ++*this; return t; }
  I128_HD int128_t operator--(int) { int128_t t = *this; --*this; return t; }
};

I128_HD uint128_t::uint128_t(const int128_t& v) : lo(v.lo), hi((uint64_t)v.hi) {}
