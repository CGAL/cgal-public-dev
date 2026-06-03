#define CUBQL_GPU_BUILDER_IMPLEMENTATION 1
#include "cuBQL/bvh.h"
#include "cuBQL/traversal/fixedBoxQuery.h"
#include <thrust/device_ptr.h>
#include <thrust/scan.h>
#include <thrust/reduce.h>
#include <thrust/device_vector.h>
#include <vector>
#include <iostream>
#include "samples/common/loadOBJ.h"

struct IntersectionPair {
    int idA;
    int idB;
};

struct GPUTimingBreakdown {
    double uploadTime = 0.0;
    double executionTime = 0.0; // BVH Build + Query Passes
    double downloadTime = 0.0;
};

// Vector math utilities from your pristine kernel configuration
__device__ inline cuBQL::vec3f crossProduct(const cuBQL::vec3f &u, const cuBQL::vec3f &v) {
    return cuBQL::vec3f(u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x);
}

__device__ inline float dotProduct(const cuBQL::vec3f &u, const cuBQL::vec3f &v) {
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

// ====================================================================
// HARDWARE INTERVAL ARITHMETIC PRIMITIVES (IEEE-754 Directed Rounding)
// ====================================================================
struct GPUInterval {
    float lo;
    float hi;
};

enum PairStatus
{
    PAIR_NO     = 0,
    PAIR_GREEN  = 1,
    PAIR_YELLOW = 2
};

__device__ inline GPUInterval interval_sub(float a, float b) {
    return { __fsub_rd(a, b), __fsub_ru(a, b) };
}

__device__ inline GPUInterval interval_sub(GPUInterval a, GPUInterval b) {
    return { __fsub_rd(a.lo, b.hi), __fsub_ru(a.hi, b.lo) };
}

__device__ inline GPUInterval interval_mul(
    GPUInterval a,
    GPUInterval b)
{
    float lo1 = __fmul_rd(a.lo,b.lo);
    float lo2 = __fmul_rd(a.lo,b.hi);
    float lo3 = __fmul_rd(a.hi,b.lo);
    float lo4 = __fmul_rd(a.hi,b.hi);

    float hi1 = __fmul_ru(a.lo,b.lo);
    float hi2 = __fmul_ru(a.lo,b.hi);
    float hi3 = __fmul_ru(a.hi,b.lo);
    float hi4 = __fmul_ru(a.hi,b.hi);

    return {
        fminf(fminf(lo1,lo2),fminf(lo3,lo4)),
        fmaxf(fmaxf(hi1,hi2),fmaxf(hi3,hi4))
    };
}

__device__ inline GPUInterval interval_add(GPUInterval a, GPUInterval b) {
    return { __fadd_rd(a.lo, b.lo), __fadd_ru(a.hi, b.hi) };
}

// ====================================================================
// ROBUST INTERVAL ORIENT3D GEOMETRIC PREDICATE
// ====================================================================
// ====================================================================
// ROBUST FIXED INTERVAL ORIENT3D GEOMETRIC PREDICATE
// ====================================================================
// __device__ int orient3d_interval(
//     const cuBQL::vec3f& p, const cuBQL::vec3f& q, 
//     const cuBQL::vec3f& r, const cuBQL::vec3f& s,
//     const float eps = 1e-6f) // Added conservative tolerance parameter
// {
//     GPUInterval pdx = interval_sub(p.x, s.x);
//     GPUInterval pdy = interval_sub(p.y, s.y);
//     GPUInterval pdz = interval_sub(p.z, s.z);

//     GPUInterval qdx = interval_sub(q.x, s.x);
//     GPUInterval qdy = interval_sub(q.y, s.y);
//     GPUInterval qdz = interval_sub(q.z, s.z);

//     GPUInterval rdx = interval_sub(r.x, s.x);
//     GPUInterval rdy = interval_sub(r.y, s.y);
//     GPUInterval rdz = interval_sub(r.z, s.z);

//     // Formally aligned with standard 3x3 matrix determinant rules
//     GPUInterval det = interval_add(
//         interval_add(
//             interval_mul(pdx, interval_sub(interval_mul(qdy, rdz), interval_mul(qdz, rdy))),
//             interval_mul(pdy, interval_sub(interval_mul(qdz, rdx), interval_mul(qdx, rdz)))
//         ),
//         interval_mul(pdz, interval_sub(interval_mul(qdx, rdy), interval_mul(qdy, rdx)))
//     );

//     // More conservative evaluation: 
//     // The bounds must clear the safety epsilon threshold, not just 0.0f
//     if (det.lo > eps) return 1;
//     if (det.hi < -eps) return -1;
    
//     // Spans across the [-eps, eps] window -> safely flag as uncertain/degenerate
//     return 0; 
// }

__device__ int orient3d_interval(
    const cuBQL::vec3f& p, const cuBQL::vec3f& q, 
    const cuBQL::vec3f& r, const cuBQL::vec3f& s) // Removed hardcoded default eps parameter
{
    GPUInterval pdx = interval_sub(p.x, s.x);
    GPUInterval pdy = interval_sub(p.y, s.y);
    GPUInterval pdz = interval_sub(p.z, s.z);

    GPUInterval qdx = interval_sub(q.x, s.x);
    GPUInterval qdy = interval_sub(q.y, s.y);
    GPUInterval qdz = interval_sub(q.z, s.z);

    GPUInterval rdx = interval_sub(r.x, s.x);
    GPUInterval rdy = interval_sub(r.y, s.y);
    GPUInterval rdz = interval_sub(r.z, s.z);

    // Formally aligned with standard 3x3 matrix determinant rules
    GPUInterval det = interval_add(
        interval_add(
            interval_mul(pdx, interval_sub(interval_mul(qdy, rdz), interval_mul(qdz, rdy))),
            interval_mul(pdy, interval_sub(interval_mul(qdz, rdx), interval_mul(qdx, rdz)))
        ),
        interval_mul(pdz, interval_sub(interval_mul(qdx, rdy), interval_mul(qdy, rdx)))
    );

    // -----------------------------------------------------------------
    // DYNAMIC ERROR ESTIMATION (Accounts for Double->Float Truncation)
    // -----------------------------------------------------------------
    // 1. Find the maximum absolute coordinate component among all 4 points
    float max_p = fmaxf(fabsf(p.x), fmaxf(fabsf(p.y), fabsf(p.z)));
    float max_q = fmaxf(fabsf(q.x), fmaxf(fabsf(q.y), fabsf(q.z)));
    float max_r = fmaxf(fabsf(r.x), fmaxf(fabsf(r.y), fabsf(r.z)));
    float max_s = fmaxf(fabsf(s.x), fmaxf(fabsf(s.y), fabsf(s.z)));
    float max_coord = fmaxf(fmaxf(max_p, max_q), fmaxf(max_r, max_s));

    // 2. Compute 32-bit float grid resolution at this coordinate magnitude.
    // (1.1920929e-7f is FLT_EPSILON)
    float float_grid_resolution = max_coord * 1.1920929e-7f;

    // 3. Scale the coordinate error to a 3D determinant volume bound.
    // We add a safety multiplier of 8.0f for multi-step matrix operation drift.
    float dynamic_eps = 8.0f * (float_grid_resolution * float_grid_resolution * float_grid_resolution);
    
    // 4. Enforce an absolute minimum floor baseline for meshes close to the origin
    if (dynamic_eps < 1e-7f) {
        dynamic_eps = 1e-7f;
    }
    // -----------------------------------------------------------------

    // Safe conservative containment evaluation against the dynamic error threshold
    if (det.lo > dynamic_eps) return 1;
    if (det.hi < -dynamic_eps) return -1;
    
    // Spans across the uncertainty window -> safely flag as yellow for CPU exact re-evaluation
    return 0; 
}

// Highly-tight conservative Triangle-Triangle intersection filter running on candidate lists
__device__ bool triangleIntersectsSAT(const cuBQL::Triangle &A, const cuBQL::Triangle &B) {
    // Ground Truth AABB Pre-Filter check
    if (!A.bounds().overlaps(B.bounds())) return false;


    //! Debugging
    return true;

}


__device__ inline int edgeTriInterval(
    const cuBQL::vec3f &a,
    const cuBQL::vec3f &b,
    const cuBQL::vec3f &p,
    const cuBQL::vec3f &q,
    const cuBQL::vec3f &r)
{
    // Step 1: plane crossing
    int s0 = orient3d_interval(p,q,r,a);
    int s1 = orient3d_interval(p,q,r,b);

    if (s0 == 0 || s1 == 0)
        return PAIR_YELLOW;

    if ((s0 > 0 && s1 > 0) || (s0 < 0 && s1 < 0))
        return PAIR_NO;

    // Step 2: inside-triangle test using edge orientation
    int e0 = orient3d_interval(a,b,p,q);
    int e1 = orient3d_interval(a,b,q,r);
    int e2 = orient3d_interval(a,b,r,p);

    if (e0 == 0 || e1 == 0 || e2 == 0)
        return PAIR_YELLOW;

    if ((e0 >= 0 && e1 >= 0 && e2 >= 0) ||
        (e0 <= 0 && e1 <= 0 && e2 <= 0))
        return PAIR_GREEN;

    return PAIR_NO;
}

// ====================================================================
// ROBUST GEOMETRIC TRIANGLE-TRIANGLE INTERSECTION CLASSIFIER
// ====================================================================
// ====================================================================
// ROBUST GEOMETRIC TRIANGLE-TRIANGLE INTERSECTION CLASSIFIER
// ====================================================================

// Ensure your types are defined. Assuming cuBQL::vec3f is a 3-float struct
// and PairStatus is an enum or typedef.
__device__ inline PairStatus classifyPair(const cuBQL::Triangle& A, const cuBQL::Triangle& B)
{
    // 1. AABB overlap check
    if (!A.bounds().overlaps(B.bounds()))
        return PAIR_NO;

    const cuBQL::vec3f A0 = A.a, A1 = A.b, A2 = A.c;
    const cuBQL::vec3f B0 = B.a, B1 = B.b, B2 = B.c;

    // ----------------------------------------------------
    // 2. FAST PLANE SEPARATION (FILTER ONLY)
    // ----------------------------------------------------
    int ob0 = orient3d_interval(A0, A1, A2, B0);
    int ob1 = orient3d_interval(A0, A1, A2, B1);
    int ob2 = orient3d_interval(A0, A1, A2, B2);

    if (ob0 != 0 && ob1 != 0 && ob2 != 0) {
        if ((ob0 > 0 && ob1 > 0 && ob2 > 0) ||
            (ob0 < 0 && ob1 < 0 && ob2 < 0))
            return PAIR_NO;
    }

    int oa0 = orient3d_interval(B0, B1, B2, A0);
    int oa1 = orient3d_interval(B0, B1, B2, A1);
    int oa2 = orient3d_interval(B0, B1, B2, A2);

    if (oa0 != 0 && oa1 != 0 && oa2 != 0) {
        if ((oa0 > 0 && oa1 > 0 && oa2 > 0) ||
            (oa0 < 0 && oa1 < 0 && oa2 < 0))
            return PAIR_NO;
    }

    // ----------------------------------------------------
    // 3. EDGE TEST (ACTUAL INTERSECTION)
    // ----------------------------------------------------
    int r;

    r = edgeTriInterval(A0, A1, B0, B1, B2);
    if (r == PAIR_GREEN) return PAIR_GREEN;
    if (r == PAIR_YELLOW) return PAIR_YELLOW;

    r = edgeTriInterval(A1, A2, B0, B1, B2);
    if (r == PAIR_GREEN) return PAIR_GREEN;
    if (r == PAIR_YELLOW) return PAIR_YELLOW;

    r = edgeTriInterval(A2, A0, B0, B1, B2);
    if (r == PAIR_GREEN) return PAIR_GREEN;
    if (r == PAIR_YELLOW) return PAIR_YELLOW;

    r = edgeTriInterval(B0, B1, A0, A1, A2);
    if (r == PAIR_GREEN) return PAIR_GREEN;
    if (r == PAIR_YELLOW) return PAIR_YELLOW;

    r = edgeTriInterval(B1, B2, A0, A1, A2);
    if (r == PAIR_GREEN) return PAIR_GREEN;
    if (r == PAIR_YELLOW) return PAIR_YELLOW;

    r = edgeTriInterval(B2, B0, A0, A1, A2);
    if (r == PAIR_GREEN) return PAIR_GREEN;
    if (r == PAIR_YELLOW) return PAIR_YELLOW;

    return PAIR_NO;
}
// You can safely leave this as a basic bounds overlap wrapper or let classifyPair handle it.
//__device__ bool triangleIntersectsSAT(const cuBQL::Triangle &A, const cuBQL::Triangle &B) {
 //   return A.bounds().overlaps(B.bounds());
//}


// ====================================================================
// STEP 1: PARALLEL PAIR GENERATION KERNELS
// ====================================================================
__global__ void generateBoxes(cuBQL::box3f *boxes, const cuBQL::Triangle *tris, int N) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < N) boxes[i] = tris[i].bounds();
}

__global__ void countPairsKernel(int *pairCounts, cuBQL::bvh3f bvhA, const cuBQL::Triangle *triB, int NB) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= NB) return;

    cuBQL::Triangle b = triB[tid];
    cuBQL::box3f query = b.bounds();
    int count = 0;

    cuBQL::fixedBoxQuery::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        for (uint32_t i = 0; i < num; i++) {
            count++;
        }
        return CUBQL_CONTINUE_TRAVERSAL;
    }, bvhA, query);

    pairCounts[tid] = count;
}

__global__ void fillPairsKernel(int2 *candidatePairs, const int *offsets, cuBQL::bvh3f bvhA, const cuBQL::Triangle *triB, int NB) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= NB) return;

    // Fixed: Pull stable global unique write position based on current exclusive scan offset tracker
    int wPos = offsets[tid];
    cuBQL::Triangle b = triB[tid];
    cuBQL::box3f query = b.bounds();

    cuBQL::fixedBoxQuery::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        for (uint32_t i = 0; i < num; i++) {
            candidatePairs[wPos++] = make_int2((int)ids[i], tid);
        }
        return CUBQL_CONTINUE_TRAVERSAL;
    }, bvhA, query);
}

// ====================================================================
// STEP 2: DIVERGENCE-FREE INTERSECTION FILTER KERNEL
// ====================================================================
__global__ void evaluatePairsKernel(
    int *statusArray, const int2 *candidatePairs, 
    const cuBQL::Triangle *triA, const cuBQL::Triangle *triB, int numPairs) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= numPairs) return;

    int2 pair = candidatePairs[tid];
    cuBQL::Triangle ta = triA[pair.x];
    cuBQL::Triangle tb = triB[pair.y];

    statusArray[tid] = (int)classifyPair(ta,tb);
}

struct IsGreen
{
    __device__ bool operator()(const int x)
    {
        return x == PAIR_GREEN;
    }
};

struct IsYellow
{
    __device__ bool operator()(const int x)
    {
        return x == PAIR_YELLOW;
    }
};

extern "C" void runPartitionedMeshIntersection(
    const cuBQL::Triangle* hA, int NA, 
    const cuBQL::Triangle* hB, int NB, 
    std::vector<IntersectionPair>& hGreenPairs,
    std::vector<IntersectionPair>& hYellowPairs,
    GPUTimingBreakdown& outTimings) 
{
    double t0, t_exec;

    t0 = cuBQL::getCurrentTime();
    cuBQL::Triangle *dA, *dB;
    CUBQL_CUDA_CALL(Malloc(&dA, NA * sizeof(cuBQL::Triangle)));
    CUBQL_CUDA_CALL(Malloc(&dB, NB * sizeof(cuBQL::Triangle)));
    CUBQL_CUDA_CALL(Memcpy(dA, hA, NA * sizeof(cuBQL::Triangle), cudaMemcpyDefault));
    CUBQL_CUDA_CALL(Memcpy(dB, hB, NB * sizeof(cuBQL::Triangle), cudaMemcpyDefault));
    outTimings.uploadTime = cuBQL::getCurrentTime() - t0;

    t_exec = cuBQL::getCurrentTime();
    cuBQL::box3f *dBoxes;
    CUBQL_CUDA_CALL(Malloc(&dBoxes, NA * sizeof(cuBQL::box3f)));
    generateBoxes<<<cuBQL::divRoundUp(NA, 256), 256>>>(dBoxes, dA, NA);
    
    cuBQL::bvh3f bvh;
    cuBQL::gpuBuilder(bvh, dBoxes, NA, cuBQL::BuildConfig());
    CUBQL_CUDA_CALL(Free(dBoxes));

    int *dPairCounts, *dOffsets;
    CUBQL_CUDA_CALL(Malloc(&dPairCounts, NB * sizeof(int)));
    CUBQL_CUDA_CALL(Malloc(&dOffsets, NB * sizeof(int)));

    countPairsKernel<<<cuBQL::divRoundUp(NB, 128), 128>>>(dPairCounts, bvh, dB, NB);
    
    thrust::device_ptr<int> dev_counts(dPairCounts);
    thrust::device_ptr<int> dev_offsets(dOffsets);
    thrust::exclusive_scan(dev_counts, dev_counts + NB, dev_offsets);

    int totalPairs = 0;
    CUBQL_CUDA_CALL(Memcpy(&totalPairs, dOffsets + NB - 1, sizeof(int), cudaMemcpyDeviceToHost));
    int lastCount = 0;
    CUBQL_CUDA_CALL(Memcpy(&lastCount, dPairCounts + NB - 1, sizeof(int), cudaMemcpyDeviceToHost));
    totalPairs += lastCount;

    if (totalPairs == 0) {
        outTimings.executionTime = cuBQL::getCurrentTime() - t_exec;
        CUBQL_CUDA_CALL(Free(dA)); CUBQL_CUDA_CALL(Free(dB)); 
        CUBQL_CUDA_CALL(Free(dPairCounts)); CUBQL_CUDA_CALL(Free(dOffsets));
        cuBQL::cuda::free(bvh);
        return;
    }

    int2 *dCandidatePairs;
    CUBQL_CUDA_CALL(Malloc(&dCandidatePairs, totalPairs * sizeof(int2)));
    fillPairsKernel<<<cuBQL::divRoundUp(NB, 128), 128>>>(dCandidatePairs, dOffsets, bvh, dB, NB);

    int *dStatusArray;
    CUBQL_CUDA_CALL(Malloc(&dStatusArray, totalPairs * sizeof(int)));
    evaluatePairsKernel<<<cuBQL::divRoundUp(totalPairs, 256), 256>>>(dStatusArray, dCandidatePairs, dA, dB, totalPairs);

    thrust::device_ptr<int2> dev_pairs(dCandidatePairs);
    thrust::device_ptr<int> dev_status(dStatusArray);

    thrust::device_vector<int2> greenResult(totalPairs);
    thrust::device_vector<int2> yellowResult(totalPairs);

    auto greenEnd = thrust::copy_if(dev_pairs, dev_pairs + totalPairs, dev_status, greenResult.begin(), IsGreen());
    auto yellowEnd = thrust::copy_if(dev_pairs, dev_pairs + totalPairs, dev_status, yellowResult.begin(), IsYellow());

    int numGreen = greenEnd - greenResult.begin();
    int numYellow = yellowEnd - yellowResult.begin();
    outTimings.executionTime = cuBQL::getCurrentTime() - t_exec;

    t0 = cuBQL::getCurrentTime();
    // Fixed: Reinterpret casting raw pointer base positions ensures safe alignment mappings
    if (numGreen > 0) {
        hGreenPairs.resize(numGreen);
        CUBQL_CUDA_CALL(Memcpy(reinterpret_cast<void*>(hGreenPairs.data()), thrust::raw_pointer_cast(greenResult.data()), numGreen * sizeof(IntersectionPair), cudaMemcpyDefault));
    }
    if (numYellow > 0) {
        hYellowPairs.resize(numYellow);
        CUBQL_CUDA_CALL(Memcpy(reinterpret_cast<void*>(hYellowPairs.data()), thrust::raw_pointer_cast(yellowResult.data()), numYellow * sizeof(IntersectionPair), cudaMemcpyDefault));
    }
    outTimings.downloadTime = cuBQL::getCurrentTime() - t0;

    CUBQL_CUDA_CALL(Free(dA)); CUBQL_CUDA_CALL(Free(dB)); 
    CUBQL_CUDA_CALL(Free(dPairCounts)); CUBQL_CUDA_CALL(Free(dOffsets));
    CUBQL_CUDA_CALL(Free(dCandidatePairs)); CUBQL_CUDA_CALL(Free(dStatusArray));
    cuBQL::cuda::free(bvh);
}