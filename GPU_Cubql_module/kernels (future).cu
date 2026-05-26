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

// Fast, conservative Triangle-Triangle intersection filter
// Adjusted to use .a, .b, and .c naming layouts
__device__ bool triangleIntersects(const cuBQL::Triangle &A, const cuBQL::Triangle &B) {
    // 1. First Pass: Quick AABB overlap rejection
    if (!A.bounds().overlaps(B.bounds())) return false;

    // 2. Compute Plane Equation for Triangle A
    cuBQL::vec3f e1A = A.b - A.a;
    cuBQL::vec3f e2A = A.c - A.a;
    cuBQL::vec3f nA = crossProduct(e1A, e2A);
    float dA = -dotProduct(nA, A.a);

    // Compute signed distances of B's vertices to A's plane
    float dB0 = dotProduct(nA, B.a) + dA;
    float dB1 = dotProduct(nA, B.b) + dA;
    float dB2 = dotProduct(nA, B.c) + dA;

    // If all vertices of B are on the same side of A's plane, no intersection
    if ((dB0 > 1e-5f && dB1 > 1e-5f && dB2 > 1e-5f) || 
        (dB0 < -1e-5f && dB1 < -1e-5f && dB2 < -1e-5f)) return false;

    // 3. Compute Plane Equation for Triangle B
    cuBQL::vec3f e1B = B.b - B.a;
    cuBQL::vec3f e2B = B.c - B.a;
    cuBQL::vec3f nB = crossProduct(e1B, e2B);
    float dB_plane = -dotProduct(nB, B.a);

    // Compute signed distances of A's vertices to B's plane
    float dA0 = dotProduct(nB, A.a) + dB_plane;
    float dA1 = dotProduct(nB, A.b) + dB_plane;
    float dA2 = dotProduct(nB, A.c) + dB_plane;

    // If all vertices of A are on the same side of B's plane, no intersection
    if ((dA0 > 1e-5f && dA1 > 1e-5f && dA2 > 1e-5f) || 
        (dA0 < -1e-5f && dA1 < -1e-5f && dA2 < -1e-5f)) return false;

    // 4. Coplanar or complex edge-overlap case
    // We safely return true to let CGAL robustly determine actual intersection.
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

// Internal dispatcher that launches the correctly compiled kernel versions
template<bool Improved>
long long executeIntersectionPipeline(
    const cuBQL::Triangle* hA, int NA, 
    const cuBQL::Triangle* hB, int NB, 
    std::vector<IntersectionPair>& hPairs) 
{
    double t0;

    // GPU ALLOCATION & UPLOAD
    t0 = cuBQL::getCurrentTime();
    cuBQL::Triangle *dA, *dB;
    CUBQL_CUDA_CALL(Malloc(&dA, NA * sizeof(cuBQL::Triangle)));
    CUBQL_CUDA_CALL(Malloc(&dB, NB * sizeof(cuBQL::Triangle)));
    CUBQL_CUDA_CALL(Memcpy(dA, hA, NA * sizeof(cuBQL::Triangle), cudaMemcpyDefault));
    CUBQL_CUDA_CALL(Memcpy(dB, hB, NB * sizeof(cuBQL::Triangle), cudaMemcpyDefault));
    CUBQL_CUDA_SYNC_CHECK();
    std::cout << "[Step 2] GPU Alloc & Upload:    " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - t0) << "s\n";

    // BVH BUILD
    t0 = cuBQL::getCurrentTime();
    cuBQL::box3f *dBoxes;
    CUBQL_CUDA_CALL(Malloc(&dBoxes, NA * sizeof(cuBQL::box3f)));
    generateBoxes<<<cuBQL::divRoundUp(NA, 256), 256>>>(dBoxes, dA, NA);
    
    cuBQL::bvh3f bvh;
    cuBQL::gpuBuilder(bvh, dBoxes, NA, cuBQL::BuildConfig());
    CUBQL_CUDA_SYNC_CHECK();
    CUBQL_CUDA_CALL(Free(dBoxes));
    std::cout << "[Step 3] BVH Generation:        " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - t0) << "s\n";

    // QUERY PASS 1 (COUNTING)
    t0 = cuBQL::getCurrentTime();
    int *dCounts, *dOffsets;
    CUBQL_CUDA_CALL(Malloc(&dCounts, NB * sizeof(int)));
    CUBQL_CUDA_CALL(Malloc(&dOffsets, NB * sizeof(int)));
    
    countKernel<Improved><<<cuBQL::divRoundUp(NB, 128), 128>>>(dCounts, dA, bvh, dB, NB);
    CUBQL_CUDA_SYNC_CHECK();
    std::cout << "[Step 4] Count Kernel (Pass 1): " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - t0) << "s\n";

    // THRUST REDUCE & SCAN
    t0 = cuBQL::getCurrentTime();
    thrust::device_ptr<int> dev_counts(dCounts);
    thrust::device_ptr<int> dev_offsets(dOffsets);
    long long totalIntersections = thrust::reduce(dev_counts, dev_counts + NB, 0LL, thrust::plus<long long>());
    thrust::exclusive_scan(dev_counts, dev_counts + NB, dev_offsets);
    CUBQL_CUDA_SYNC_CHECK();
    std::cout << "[Step 5] Thrust Reduce & Scan:  " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - t0) << "s\n";

    // QUERY PASS 2 (FILLING)
    t0 = cuBQL::getCurrentTime();
    IntersectionPair *dPairs = nullptr;
    if (totalIntersections > 0) {
        CUBQL_CUDA_CALL(Malloc(&dPairs, totalIntersections * sizeof(IntersectionPair)));
        fillKernel<Improved><<<cuBQL::divRoundUp(NB, 128), 128>>>(dPairs, dOffsets, dA, bvh, dB, NB);
        CUBQL_CUDA_SYNC_CHECK();
    }
    std::cout << "[Step 6] Fill Kernel (Pass 2):  " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - t0) << "s\n";

    // DOWNLOAD FROM GPU
    t0 = cuBQL::getCurrentTime();
    if (totalIntersections > 0) {
        hPairs.resize(totalIntersections);
        CUBQL_CUDA_CALL(Memcpy(hPairs.data(), dPairs, totalIntersections * sizeof(IntersectionPair), cudaMemcpyDefault));
    }
    std::cout << "[Step 7] Download Results:      " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - t0) << "s\n";

    CUBQL_CUDA_CALL(Free(dA)); CUBQL_CUDA_CALL(Free(dB));
    CUBQL_CUDA_CALL(Free(dCounts)); CUBQL_CUDA_CALL(Free(dOffsets));
    if (dPairs) CUBQL_CUDA_CALL(Free(dPairs));
    cuBQL::cuda::free(bvh);

    return totalIntersections;
}

// Host entry point called from main.cpp with a runtime default parameter
extern "C" long long runMeshIntersection(
    const cuBQL::Triangle* hA, int NA, 
    const cuBQL::Triangle* hB, int NB, 
    std::vector<IntersectionPair>& hPairs,
    bool useImprovedIntersect = false) 
{
    if (useImprovedIntersect) {
        return executeIntersectionPipeline<true>(hA, NA, hB, NB, hPairs);
    } else {
        return executeIntersectionPipeline<false>(hA, NA, hB, NB, hPairs);
    }
}