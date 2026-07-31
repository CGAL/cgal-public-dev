#pragma once

#include <cuda_runtime.h>
#include <cmath>
#include "PairStatus.h"
#include "cuBQL/bvh.h" // Needed for cuBQL::Triangle reference
#include "samples/common/loadOBJ.h"

#ifndef CUBQL_CUDA_CALL
#define CUBQL_CUDA_CALL(X) { cudaError_t err = cuda ## X; if(err != cudaSuccess) { printf("CUDA Error at %s:%d\n", __FILE__, __LINE__); } }
#endif

// --------------------------------------------------------------------
// INTERVAL ARITHMETIC & GEOMETRIC CLASS CONTEXT (DEVICE INTERNAL)
// --------------------------------------------------------------------
struct GPUInterval {
    float lo;
    float hi;
};

__device__ inline GPUInterval interval_sub(float a, float b) { 
    return { __fsub_rd(a, b), __fsub_ru(a, b) }; 
}

__device__ inline GPUInterval interval_add(GPUInterval a, GPUInterval b) { 
    return { __fadd_rd(a.lo, b.lo), __fadd_ru(a.hi, b.hi) }; 
}

__device__ inline GPUInterval interval_sub(GPUInterval a, GPUInterval b) { 
    return { __fsub_rd(a.lo, b.hi), __fsub_ru(a.hi, b.lo) }; 
}

__device__ inline GPUInterval interval_mul(GPUInterval a, GPUInterval b) {
    float lo1 = __fmul_rd(a.lo, b.lo); float lo2 = __fmul_rd(a.lo, b.hi);
    float lo3 = __fmul_rd(a.hi, b.lo); float lo4 = __fmul_rd(a.hi, b.hi);
    float hi1 = __fmul_ru(a.lo, b.lo); float hi2 = __fmul_ru(a.lo, b.hi);
    float hi3 = __fmul_ru(a.hi, b.lo); float hi4 = __fmul_ru(a.hi, b.hi);
    return { fminf(fminf(lo1, lo2), fminf(lo3, lo4)), fmaxf(fmaxf(hi1, hi2), fmaxf(hi3, hi4)) };
}

__device__ inline int orient3d_interval(
    const cuBQL::vec3f& p, const cuBQL::vec3f& q, 
    const cuBQL::vec3f& r, const cuBQL::vec3f& s, float eps) 
{
    GPUInterval pdx = interval_sub(p.x, s.x); GPUInterval pdy = interval_sub(p.y, s.y); GPUInterval pdz = interval_sub(p.z, s.z);
    GPUInterval qdx = interval_sub(q.x, s.x); GPUInterval qdy = interval_sub(q.y, s.y); GPUInterval qdz = interval_sub(q.z, s.z);
    GPUInterval rdx = interval_sub(r.x, s.x); GPUInterval rdy = interval_sub(r.y, s.y); GPUInterval rdz = interval_sub(r.z, s.z);

    GPUInterval det = interval_add(
        interval_add(
            interval_mul(pdx, interval_sub(interval_mul(qdy, rdz), interval_mul(qdz, rdy))),
            interval_mul(pdy, interval_sub(interval_mul(qdz, rdx), interval_mul(qdx, rdz)))
        ),
        interval_mul(pdz, interval_sub(interval_mul(qdx, rdy), interval_mul(qdy, rdx)))
    );

    if (det.lo > eps) return 1;
    if (det.hi < -eps) return -1;
    return 0; 
}

__device__ inline int edgeTriInterval(
    const cuBQL::vec3f &a, const cuBQL::vec3f &b,
    const cuBQL::vec3f &p, const cuBQL::vec3f &q, const cuBQL::vec3f &r, float eps) 
{
    int s0 = orient3d_interval(p, q, r, a, eps);
    int s1 = orient3d_interval(p, q, r, b, eps);

    if (s0 == 0 || s1 == 0) return PAIR_YELLOW;
    if ((s0 > 0 && s1 > 0) || (s0 < 0 && s1 < 0)) return PAIR_NO;

    int e0 = orient3d_interval(a, b, p, q, eps);
    int e1 = orient3d_interval(a, b, q, r, eps);
    int e2 = orient3d_interval(a, b, r, p, eps);

    if (e0 == 0 || e1 == 0 || e2 == 0) return PAIR_YELLOW;
    if ((e0 >= 0 && e1 >= 0 && e2 >= 0) || (e0 <= 0 && e1 <= 0 && e2 <= 0)) return PAIR_GREEN;

    return PAIR_NO;
}