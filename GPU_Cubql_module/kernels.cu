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

__device__ bool intersects(const cuBQL::Triangle &a, const cuBQL::Triangle &b) {
    return a.bounds().overlaps(b.bounds());
}

__global__ void countKernel(int *counts, const cuBQL::Triangle *triA, cuBQL::bvh3f bvhA, const cuBQL::Triangle *triB, int N) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= N) return;
    cuBQL::Triangle b = triB[tid];
    cuBQL::box3f query = b.bounds();
    int count = 0;
    cuBQL::fixedBoxQuery::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        for (int i = 0; i < num; i++)
            if (intersects(triA[ids[i]], b)) count++;
        return CUBQL_CONTINUE_TRAVERSAL;
    }, bvhA, query);
    counts[tid] = count;
}

__global__ void fillKernel(IntersectionPair *pairs, const int *offsets, const cuBQL::Triangle *triA, cuBQL::bvh3f bvhA, const cuBQL::Triangle *triB, int N) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= N) return;
    int writePos = offsets[tid];
    cuBQL::Triangle b = triB[tid];
    cuBQL::box3f query = b.bounds();
    cuBQL::fixedBoxQuery::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        for (int i = 0; i < num; i++) {
            if (intersects(triA[ids[i]], b)) {
                pairs[writePos++] = { (int)ids[i], tid };
            }
        }
        return CUBQL_CONTINUE_TRAVERSAL;
    }, bvhA, query);
}

// Struct to hold individual GPU phase timings (in seconds)
struct GPUTimingBreakdown {
    double uploadTime = 0.0;
    double executionTime = 0.0; // BVH Build + Query Passes
    double downloadTime = 0.0;
};

// Update the external C signature to return our new timing struct
extern "C" long long runMeshIntersection(
    const cuBQL::Triangle* hA, int NA, 
    const cuBQL::Triangle* hB, int NB, 
    std::vector<IntersectionPair>& hPairs,
    GPUTimingBreakdown& outTimings) // <--- Passing the struct by reference
{
    double t0;

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
    t0 = cuBQL::getCurrentTime(); // Reset timer for execution
    
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
    
    countKernel<<<cuBQL::divRoundUp(NB, 128), 128>>>(dCounts, dA, bvh, dB, NB);
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
        fillKernel<<<cuBQL::divRoundUp(NB, 128), 128>>>(dPairs, dOffsets, dA, bvh, dB, NB);
        CUBQL_CUDA_SYNC_CHECK();
    }
    
    outTimings.executionTime = cuBQL::getCurrentTime() - t0;

    // ==========================================
    // PHASE 3: DOWNLOAD FROM GPU
    // ==========================================
    t0 = cuBQL::getCurrentTime(); // Reset timer for download
    
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