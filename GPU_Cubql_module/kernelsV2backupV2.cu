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

// Layout for the single-flight monolithic download channel
struct TaggedIntersectionPair {
    int idA;
    int idB;
    int status; // 1 = GREEN, 2 = YELLOW
};

struct GPUTimingBreakdown {
    double uploadTime = 0.0;
    double executionTime = 0.0; // BVH Build + Query Passes
    double downloadTime = 0.0;
};

// Vector math utilities
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

__device__ inline GPUInterval interval_mul(GPUInterval a, GPUInterval b)
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

__device__ int orient3d_interval(
    const cuBQL::vec3f& p, const cuBQL::vec3f& q, 
    const cuBQL::vec3f& r, const cuBQL::vec3f& s,
    float eps)
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
    const cuBQL::vec3f &p, const cuBQL::vec3f &q, const cuBQL::vec3f &r,
    float eps)
{
    int s0 = orient3d_interval(p, q, r, a, eps);
    int s1 = orient3d_interval(p, q, r, b, eps);

    if (s0 == 0 || s1 == 0) return PAIR_YELLOW;
    if ((s0 > 0 && s1 > 0) || (s0 < 0 && s1 < 0)) return PAIR_NO;

    int e0 = orient3d_interval(a, b, p, q, eps);
    int e1 = orient3d_interval(a, b, q, r, eps);
    int e2 = orient3d_interval(a, b, r, p, eps);

    if (e0 == 0 || e1 == 0 || e2 == 0) return PAIR_YELLOW;

    if ((e0 >= 0 && e1 >= 0 && e2 >= 0) || (e0 <= 0 && e1 <= 0 && e2 <= 0))
        return PAIR_GREEN;

    return PAIR_NO;
}

__device__ inline PairStatus classifyPair(const cuBQL::Triangle& A, const cuBQL::Triangle& B)
{
    if (!A.bounds().overlaps(B.bounds()))
        return PAIR_NO;

    const cuBQL::vec3f A0 = A.a, A1 = A.b, A2 = A.c;
    const cuBQL::vec3f B0 = B.a, B1 = B.b, B2 = B.c;

    float max_A = fmaxf(fmaxf(fabsf(A0.x), fabsf(A0.y)), fmaxf(fabsf(A0.z), fmaxf(fabsf(A1.x), fmaxf(fabsf(A1.y), fabsf(A1.z)))));
    max_A = fmaxf(max_A, fmaxf(fabsf(A2.x), fmaxf(fabsf(A2.y), fabsf(A2.z))));
    
    float max_B = fmaxf(fmaxf(fabsf(B0.x), fabsf(B0.y)), fmaxf(fabsf(B0.z), fmaxf(fabsf(B1.x), fmaxf(fabsf(B1.y), fabsf(B1.z)))));
    max_B = fmaxf(max_B, fmaxf(fabsf(B2.x), fmaxf(fabsf(B2.y), fabsf(B2.z))));
    
    float max_coord = fmaxf(max_A, max_B);
    float float_grid_resolution = max_coord * 1.1920929e-7f; 
    float eps = 8.0f * (float_grid_resolution * float_grid_resolution * float_grid_resolution);
    
    if (eps < 1e-7f) eps = 1e-7f;

    int ob0 = orient3d_interval(A0, A1, A2, B0, eps);
    int ob1 = orient3d_interval(A0, A1, A2, B1, eps);
    int ob2 = orient3d_interval(A0, A1, A2, B2, eps);

    if (ob0 != 0 && ob1 != 0 && ob2 != 0) {
        if ((ob0 > 0 && ob1 > 0 && ob2 > 0) || (ob0 < 0 && ob1 < 0 && ob2 < 0))
            return PAIR_NO;
    }

    int oa0 = orient3d_interval(B0, B1, B2, A0, eps);
    int oa1 = orient3d_interval(B0, B1, B2, A1, eps);
    int oa2 = orient3d_interval(B0, B1, B2, A2, eps);

    if (oa0 != 0 && oa1 != 0 && oa2 != 0) {
        if ((oa0 > 0 && oa1 > 0 && oa2 > 0) || (oa0 < 0 && oa1 < 0 && oa2 < 0))
            return PAIR_NO;
    }

    int r;
    r = edgeTriInterval(A0, A1, B0, B1, B2, eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;
    r = edgeTriInterval(A1, A2, B0, B1, B2, eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;
    r = edgeTriInterval(A2, A0, B0, B1, B2, eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;

    r = edgeTriInterval(B0, B1, A0, A1, A2, eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;
    r = edgeTriInterval(B1, B2, A0, A1, A2, eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;
    r = edgeTriInterval(B2, B0, A0, A1, A2, eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;

    return PAIR_NO;
}

__global__ void generateBoxes(cuBQL::box3f *boxes, const cuBQL::Triangle *tris, int N) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < N) boxes[i] = tris[i].bounds();
}

__global__ void countPairsKernel(
    int *greenCounts, int *yellowCounts, 
    cuBQL::bvh3f bvhA, const cuBQL::Triangle *triA, 
    const cuBQL::Triangle *triB, int NB) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= NB) return;

    cuBQL::Triangle b = triB[tid];
    cuBQL::box3f query = b.bounds();
    int gCount = 0;
    int yCount = 0;

    cuBQL::fixedBoxQuery::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        for (uint32_t i = 0; i < num; i++) {
            PairStatus status = classifyPair(triA[ids[i]], b);
            if (status == PAIR_GREEN) gCount++;
            else if (status == PAIR_YELLOW) yCount++;
        }
        return CUBQL_CONTINUE_TRAVERSAL;
    }, bvhA, query);

    greenCounts[tid] = gCount;
    yellowCounts[tid] = yCount;
}

__global__ void fillPairsKernel(
    IntersectionPair *greenPairs, const int *greenOffsets,
    IntersectionPair *yellowPairs, const int *yellowOffsets,
    cuBQL::bvh3f bvhA, const cuBQL::Triangle *triA, 
    const cuBQL::Triangle *triB, int NB) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= NB) return;

    int wgPos = greenOffsets[tid];
    int wyPos = yellowOffsets[tid];
    cuBQL::Triangle b = triB[tid];
    cuBQL::box3f query = b.bounds();

    cuBQL::fixedBoxQuery::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        for (uint32_t i = 0; i < num; i++) {
            PairStatus status = classifyPair(triA[ids[i]], b);
            if (status == PAIR_GREEN) {
                greenPairs[wgPos++] = { (int)ids[i], tid };
            } else if (status == PAIR_YELLOW) {
                yellowPairs[wyPos++] = { (int)ids[i], tid };
            }
        }
        return CUBQL_CONTINUE_TRAVERSAL;
    }, bvhA, query);
}

// ====================================================================
// NEW: ULTRA-FAST INTERNAL MEMORY CONSOLIDATION KERNELS (D2D)
// ====================================================================
__global__ void mergeGreenToUnifiedKernel(TaggedIntersectionPair* dest, const IntersectionPair* src, int N) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < N) {
        dest[idx] = { src[idx].idA, src[idx].idB, 1 }; // Status 1 = Green
    }
}

__global__ void mergeYellowToUnifiedKernel(TaggedIntersectionPair* dest, const IntersectionPair* src, int N) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < N) {
        dest[idx] = { src[idx].idA, src[idx].idB, 2 }; // Status 2 = Yellow
    }
}

extern "C" void runPartitionedMeshIntersection(
    const cuBQL::Triangle* hA, int NA, 
    const cuBQL::Triangle* hB, int NB, 
    std::vector<TaggedIntersectionPair>& hUnifiedPairs, // Unified Host Array Receiver
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

    int *dGreenCounts, *dGreenOffsets;
    int *dYellowCounts, *dYellowOffsets;
    CUBQL_CUDA_CALL(Malloc(&dGreenCounts, NB * sizeof(int)));
    CUBQL_CUDA_CALL(Malloc(&dGreenOffsets, NB * sizeof(int)));
    CUBQL_CUDA_CALL(Malloc(&dYellowCounts, NB * sizeof(int)));
    CUBQL_CUDA_CALL(Malloc(&dYellowOffsets, NB * sizeof(int)));

    countPairsKernel<<<cuBQL::divRoundUp(NB, 128), 128>>>(dGreenCounts, dYellowCounts, bvh, dA, dB, NB);
    
    thrust::device_ptr<int> dev_gCounts(dGreenCounts);
    thrust::device_ptr<int> dev_gOffsets(dGreenOffsets);
    thrust::exclusive_scan(dev_gCounts, dev_gCounts + NB, dev_gOffsets);

    thrust::device_ptr<int> dev_yCounts(dYellowCounts);
    thrust::device_ptr<int> dev_yOffsets(dYellowOffsets);
    thrust::exclusive_scan(dev_yCounts, dev_yCounts + NB, dev_yOffsets);

    int totalGreen = 0, lastGreen = 0;
    int totalYellow = 0, lastYellow = 0;
    CUBQL_CUDA_CALL(Memcpy(&totalGreen, dGreenOffsets + NB - 1, sizeof(int), cudaMemcpyDeviceToHost));
    CUBQL_CUDA_CALL(Memcpy(&lastGreen, dGreenCounts + NB - 1, sizeof(int), cudaMemcpyDeviceToHost));
    totalGreen += lastGreen;

    CUBQL_CUDA_CALL(Memcpy(&totalYellow, dYellowOffsets + NB - 1, sizeof(int), cudaMemcpyDeviceToHost));
    CUBQL_CUDA_CALL(Memcpy(&lastYellow, dYellowCounts + NB - 1, sizeof(int), cudaMemcpyDeviceToHost));
    totalYellow += lastYellow;

    IntersectionPair *dGreenPairs = nullptr;
    IntersectionPair *dYellowPairs = nullptr;

    if (totalGreen > 0) CUBQL_CUDA_CALL(Malloc(&dGreenPairs, totalGreen * sizeof(IntersectionPair)));
    if (totalYellow > 0) CUBQL_CUDA_CALL(Malloc(&dYellowPairs, totalYellow * sizeof(IntersectionPair)));

    if (totalGreen > 0 || totalYellow > 0) {
        fillPairsKernel<<<cuBQL::divRoundUp(NB, 128), 128>>>(
            dGreenPairs, dGreenOffsets, 
            dYellowPairs, dYellowOffsets, 
            bvh, dA, dB, NB
        );
    }
    CUBQL_CUDA_SYNC_CHECK();
    outTimings.executionTime = cuBQL::getCurrentTime() - t_exec;

    // ====================================================================
    // THE FIX: CONSOLIDATE ON THE GPU AND STREAM IN A SINGLE PCIe FLIGHT
    // ====================================================================
    t0 = cuBQL::getCurrentTime();
    
    int totalElements = totalGreen + totalYellow;
    if (totalElements > 0) {
        hUnifiedPairs.resize(totalElements);

        TaggedIntersectionPair* dUnifiedStaging = nullptr;
        CUBQL_CUDA_CALL(Malloc(&dUnifiedStaging, totalElements * sizeof(TaggedIntersectionPair)));

        int block = 256;
        if (totalGreen > 0) {
            int gridGreen = cuBQL::divRoundUp(totalGreen, block);
            mergeGreenToUnifiedKernel<<<gridGreen, block>>>(dUnifiedStaging, dGreenPairs, totalGreen);
        }
        if (totalYellow > 0) {
            int gridYellow = cuBQL::divRoundUp(totalYellow, block);
            mergeYellowToUnifiedKernel<<<gridYellow, block>>>(dUnifiedStaging + totalGreen, dYellowPairs, totalYellow);
        }
        CUBQL_CUDA_SYNC_CHECK();

        // One contiguous PCIe transaction burst
        CUBQL_CUDA_CALL(Memcpy(hUnifiedPairs.data(), dUnifiedStaging, totalElements * sizeof(TaggedIntersectionPair), cudaMemcpyDefault));
        
        CUBQL_CUDA_CALL(Free(dUnifiedStaging));
    }
    CUBQL_CUDA_SYNC_CHECK();
    outTimings.downloadTime = cuBQL::getCurrentTime() - t0;

    // Cleanup
    CUBQL_CUDA_CALL(Free(dA)); CUBQL_CUDA_CALL(Free(dB)); 
    CUBQL_CUDA_CALL(Free(dGreenCounts)); CUBQL_CUDA_CALL(Free(dGreenOffsets));
    CUBQL_CUDA_CALL(Free(dYellowCounts)); CUBQL_CUDA_CALL(Free(dYellowOffsets));
    if (dGreenPairs) CUBQL_CUDA_CALL(Free(dGreenPairs));
    if (dYellowPairs) CUBQL_CUDA_CALL(Free(dYellowPairs));
    cuBQL::cuda::free(bvh);
}