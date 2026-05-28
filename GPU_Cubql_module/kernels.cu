#define CUBQL_GPU_BUILDER_IMPLEMENTATION 1
#include "cuBQL/bvh.h"
#include "cuBQL/traversal/fixedBoxQuery.h"
#include "samples/common/loadOBJ.h"
#include <thrust/device_ptr.h>
#include <thrust/scan.h>
#include <thrust/reduce.h>
#include <vector>
#include <iostream>

struct IntersectionPair {
    int idA;
    int idB;
};

// Struct to hold individual GPU phase timings (in seconds)
struct GPUTimingBreakdown {
    double uploadTime = 0.0;
    double executionTime = 0.0; // BVH Build + Query Passes
    double downloadTime = 0.0;
};

__global__ void generateBoxes(cuBQL::box3f *boxes, const cuBQL::Triangle *tris, int N) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < N) boxes[i] = tris[i].bounds();
}

// Vector math utilities
__device__ inline cuBQL::vec3f crossProduct(const cuBQL::vec3f &u, const cuBQL::vec3f &v) {
    return cuBQL::vec3f(u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x);
}

__device__ inline float dotProduct(const cuBQL::vec3f &u, const cuBQL::vec3f &v) {
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

// // Fast, conservative Triangle-Triangle intersection filter
// __device__ bool triangleIntersects(const cuBQL::Triangle &A, const cuBQL::Triangle &B) {
//     // 1. First Pass: Quick AABB overlap rejection
//     if (!A.bounds().overlaps(B.bounds())) return false;

//     // 2. Compute Plane Equation for Triangle A
//     cuBQL::vec3f e1A = A.b - A.a;
//     cuBQL::vec3f e2A = A.c - A.a;
//     cuBQL::vec3f nA = crossProduct(e1A, e2A);
//     float dA = -dotProduct(nA, A.a);

//     // Compute signed distances of B's vertices to A's plane
//     float dB0 = dotProduct(nA, B.a) + dA;
//     float dB1 = dotProduct(nA, B.b) + dA;
//     float dB2 = dotProduct(nA, B.c) + dA;

//     // If all vertices of B are on the same side of A's plane, no intersection
//     if ((dB0 > 1e-5f && dB1 > 1e-5f && dB2 > 1e-5f) || 
//         (dB0 < -1e-5f && dB1 < -1e-5f && dB2 < -1e-5f)) return false;

//     // 3. Compute Plane Equation for Triangle B
//     cuBQL::vec3f e1B = B.b - B.a;
//     cuBQL::vec3f e2B = B.c - B.a;
//     cuBQL::vec3f nB = crossProduct(e1B, e2B);
//     float dB_plane = -dotProduct(nB, B.a);

//     // Compute signed distances of A's vertices to B's plane
//     float dA0 = dotProduct(nB, A.a) + dB_plane;
//     float dA1 = dotProduct(nB, A.b) + dB_plane;
//     float dA2 = dotProduct(nB, A.c) + dB_plane;

//     // If all vertices of A are on the same side of B's plane, no intersection
//     if ((dA0 > 1e-5f && dA1 > 1e-5f && dA2 > 1e-5f) || 
//         (dA0 < -1e-5f && dA1 < -1e-5f && dA2 < -1e-5f)) return false;

//     // 4. Coplanar or complex edge-overlap case
//     return true;
// }

// Fast, highly-tight conservative Triangle-Triangle intersection filter
__device__ bool triangleIntersects(const cuBQL::Triangle &A, const cuBQL::Triangle &B) {
    // 1. First Pass: Quick AABB overlap rejection (Keep this, it's incredibly fast)
    if (!A.bounds().overlaps(B.bounds())) return false;

    // 2. Compute Plane Equation for Triangle A
    cuBQL::vec3f e0A = A.b - A.a;
    cuBQL::vec3f e1A = A.c - A.b; // Using cyclic edges: b-a, c-b, a-c
    cuBQL::vec3f e2A = A.a - A.c;
    cuBQL::vec3f nA = crossProduct(e0A, -e2A); // Normal of A
    float dA = -dotProduct(nA, A.a);

    // Compute signed distances of B's vertices to A's plane
    float dB0 = dotProduct(nA, B.a) + dA;
    float dB1 = dotProduct(nA, B.b) + dA;
    float dB2 = dotProduct(nA, B.c) + dA;

    // Reject if all vertices of B are on the same side of A's plane
    // SIMD-friendly sign check (robust against epsilon issues)
    if ((dB0 > 1e-5f && dB1 > 1e-5f && dB2 > 1e-5f) || 
        (dB0 < -1e-5f && dB1 < -1e-5f && dB2 < -1e-5f)) return false;

    // 3. Compute Plane Equation for Triangle B
    cuBQL::vec3f e0B = B.b - B.a;
    cuBQL::vec3f e1B = B.c - B.b;
    cuBQL::vec3f e2B = B.a - B.c;
    cuBQL::vec3f nB = crossProduct(e0B, -e2B); // Normal of B
    float dB_plane = -dotProduct(nB, B.a);

    // Compute signed distances of A's vertices to B's plane
    float dA0 = dotProduct(nB, A.a) + dB_plane;
    float dA1 = dotProduct(nB, A.b) + dB_plane;
    float dA2 = dotProduct(nB, A.c) + dB_plane;

    // Reject if all vertices of A are on the same side of B's plane
    if ((dA0 > 1e-5f && dA1 > 1e-5f && dA2 > 1e-5f) || 
        (dA0 < -1e-5f && dA1 < -1e-5f && dA2 < -1e-5f)) return false;

    // =========================================================================
    // 4. NEW: Tighten the filter using Edge Cross Products (SAT)
    // =========================================================================
    // There are 9 combinations of edges (3 from A x 3 from B).
    // We project both triangles onto each axis. If their projections don't overlap, 
    // they don't intersect.
    
    cuBQL::vec3f edgesA[3] = {e0A, e1A, e2A};
    cuBQL::vec3f edgesB[3] = {e0B, e1B, e2B};
    cuBQL::vec3f vertsA[3] = {A.a, A.b, A.c};
    cuBQL::vec3f vertsB[3] = {B.a, B.b, B.c};

    #pragma unroll
    for (int i = 0; i < 3; ++i) {
        #pragma unroll
        for (int j = 0; j < 3; ++j) {
            cuBQL::vec3f axis = crossProduct(edgesA[i], edgesB[j]);
            
            // Skip degenerate axes (if edges are parallel)
            if (dotProduct(axis, axis) < 1e-6f) continue;

            // Project Triangle A
            float pA0 = dotProduct(axis, vertsA[0]);
            float pA1 = dotProduct(axis, vertsA[1]);
            float pA2 = dotProduct(axis, vertsA[2]);
            float minA = fminf(pA0, fminf(pA1, pA2));
            float maxA = fmaxf(pA0, fmaxf(pA1, pA2));

            // Project Triangle B
            float pB0 = dotProduct(axis, vertsB[0]);
            float pB1 = dotProduct(axis, vertsB[1]);
            float pB2 = dotProduct(axis, vertsB[2]);
            float minB = fminf(pB0, fminf(pB1, pB2));
            float maxB = fmaxf(pB0, fmaxf(pB1, pB2));

            // Test for separation along this axis
            if (maxA < minB - 1e-5f || maxB < minA - 1e-5f) {
                return false; // Found a separating axis! Definitely no intersection.
            }
        }
    }

    // If it survives all 2 face-planes and 9 edge-cross-product axes, 
    // the triangles are either intersecting or coplanar-overlapping.
    return true;
}

// Device router optimized away at compile time via template instantiation
template<bool Improved>
__device__ inline bool intersects(const cuBQL::Triangle &a, const cuBQL::Triangle &b) {
    if (Improved) {
        return triangleIntersects(a, b);
    } else {
        return a.bounds().overlaps(b.bounds());
    }
}

template<bool Improved>
__global__ void countKernel(int *counts, const cuBQL::Triangle *triA, cuBQL::bvh3f bvhA, const cuBQL::Triangle *triB, int N) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= N) return;
    cuBQL::Triangle b = triB[tid];
    cuBQL::box3f query = b.bounds();
    int count = 0;
    cuBQL::fixedBoxQuery::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        for (int i = 0; i < num; i++)
            if (intersects<Improved>(triA[ids[i]], b)) count++;
        return CUBQL_CONTINUE_TRAVERSAL;
    }, bvhA, query);
    counts[tid] = count;
}

template<bool Improved>
__global__ void fillKernel(IntersectionPair *pairs, const int *offsets, const cuBQL::Triangle *triA, cuBQL::bvh3f bvhA, const cuBQL::Triangle *triB, int N) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= N) return;
    int writePos = offsets[tid];
    cuBQL::Triangle b = triB[tid];
    cuBQL::box3f query = b.bounds();
    cuBQL::fixedBoxQuery::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        for (int i = 0; i < num; i++) {
            if (intersects<Improved>(triA[ids[i]], b)) {
                pairs[writePos++] = { (int)ids[i], tid };
            }
        }
        return CUBQL_CONTINUE_TRAVERSAL;
    }, bvhA, query);
}

// Internal pipeline runner tracking the GPUTimingBreakdown metrics
template<bool Improved>
long long executeIntersectionPipeline(
    const cuBQL::Triangle* hA, int NA, 
    const cuBQL::Triangle* hB, int NB, 
    std::vector<IntersectionPair>& hPairs,
    GPUTimingBreakdown& outTimings) 
{
    double t0, t_exec;

    // ==========================================
    // PHASE 1: GPU ALLOCATION & UPLOAD
    // ==========================================
    t0 = cuBQL::getCurrentTime();
    cuBQL::Triangle *dA, *dB;
    CUBQL_CUDA_CALL(Malloc(&dA, NA * sizeof(cuBQL::Triangle)));
    CUBQL_CUDA_CALL(Malloc(&dB, NB * sizeof(cuBQL::Triangle)));
    CUBQL_CUDA_CALL(Memcpy(dA, hA, NA * sizeof(cuBQL::Triangle), cudaMemcpyDefault));
    CUBQL_CUDA_CALL(Memcpy(dB, hB, NB * sizeof(cuBQL::Triangle), cudaMemcpyDefault));
    CUBQL_CUDA_SYNC_CHECK();
    
    outTimings.uploadTime = cuBQL::getCurrentTime() - t0;

    // ==========================================
    // PHASE 2: BVH BUILD & KERNEL SEARCH EXECUTION
    // ==========================================
    t_exec = cuBQL::getCurrentTime(); // Track total execution time start
    
    // BVH BUILD
    cuBQL::box3f *dBoxes;
    CUBQL_CUDA_CALL(Malloc(&dBoxes, NA * sizeof(cuBQL::box3f)));
    generateBoxes<<<cuBQL::divRoundUp(NA, 256), 256>>>(dBoxes, dA, NA);
    
    cuBQL::bvh3f bvh;
    cuBQL::gpuBuilder(bvh, dBoxes, NA, cuBQL::BuildConfig());
    CUBQL_CUDA_SYNC_CHECK();
    CUBQL_CUDA_CALL(Free(dBoxes));

    // QUERY PASS 1 (COUNTING)
    int *dCounts, *dOffsets;
    CUBQL_CUDA_CALL(Malloc(&dCounts, NB * sizeof(int)));
    CUBQL_CUDA_CALL(Malloc(&dOffsets, NB * sizeof(int)));
    
    // FIXED: Changed <<<< to <<< below
    countKernel<Improved><<<cuBQL::divRoundUp(NB, 128), 128>>>(dCounts, dA, bvh, dB, NB);
    CUBQL_CUDA_SYNC_CHECK();

    // THRUST REDUCE & SCAN
    thrust::device_ptr<int> dev_counts(dCounts);
    thrust::device_ptr<int> dev_offsets(dOffsets);
    long long totalIntersections = thrust::reduce(dev_counts, dev_counts + NB, 0LL, thrust::plus<long long>());
    thrust::exclusive_scan(dev_counts, dev_counts + NB, dev_offsets);
    CUBQL_CUDA_SYNC_CHECK();

    // QUERY PASS 2 (FILLING)
    IntersectionPair *dPairs = nullptr;
    if (totalIntersections > 0) {
        CUBQL_CUDA_CALL(Malloc(&dPairs, totalIntersections * sizeof(IntersectionPair)));
        // FIXED: Changed <<<< to <<< below
        fillKernel<Improved><<<cuBQL::divRoundUp(NB, 128), 128>>>(dPairs, dOffsets, dA, bvh, dB, NB);
        CUBQL_CUDA_SYNC_CHECK();
    }
    
    outTimings.executionTime = cuBQL::getCurrentTime() - t_exec;

    // ==========================================
    // PHASE 3: DOWNLOAD FROM GPU
    // ==========================================
    t0 = cuBQL::getCurrentTime();
    
    if (totalIntersections > 0) {
        hPairs.resize(totalIntersections);
        CUBQL_CUDA_CALL(Memcpy(hPairs.data(), dPairs, totalIntersections * sizeof(IntersectionPair), cudaMemcpyDefault));
    }
    CUBQL_CUDA_SYNC_CHECK();
    
    outTimings.downloadTime = cuBQL::getCurrentTime() - t0;

    // Cleanup memory
    CUBQL_CUDA_CALL(Free(dA)); CUBQL_CUDA_CALL(Free(dB));
    CUBQL_CUDA_CALL(Free(dCounts)); CUBQL_CUDA_CALL(Free(dOffsets));
    if (dPairs) CUBQL_CUDA_CALL(Free(dPairs));
    cuBQL::cuda::free(bvh);

    return totalIntersections;
}

// Host entry point featuring both the strategy flag and the timing profile references
extern "C" long long runMeshIntersection(
    const cuBQL::Triangle* hA, int NA, 
    const cuBQL::Triangle* hB, int NB, 
    std::vector<IntersectionPair>& hPairs,
    GPUTimingBreakdown& outTimings,
    bool useImprovedIntersect = true) 
{
    if (useImprovedIntersect) {
        return executeIntersectionPipeline<true>(hA, NA, hB, NB, hPairs, outTimings);
    } else {
        return executeIntersectionPipeline<false>(hA, NA, hB, NB, hPairs, outTimings);
    }
}