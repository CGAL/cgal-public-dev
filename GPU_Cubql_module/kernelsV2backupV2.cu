#define CUBQL_GPU_BUILDER_IMPLEMENTATION 1
#include "cuBQL/bvh.h"
#include "cuBQL/traversal/fixedBoxQuery.h"
#include <thrust/device_ptr.h>
#include <thrust/scan.h>
#include <thrust/copy.h>
#include <thrust/execution_policy.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <vector>
#include <iostream>
#include "samples/common/loadOBJ.h"

struct IntersectionPair {
    int idA;
    int idB;
};

struct GPUTimingBreakdown {
    double uploadTime = 0.0;
    double executionTime = 0.0; 
    double downloadTime = 0.0;
    double bvhBuildTime = 0.0;  
    double queryTime = 0.0;     
};

struct GPUInterval {
    float lo;
    float hi;
};

enum PairStatus {
    PAIR_NO     = 0,
    PAIR_GREEN  = 1,
    PAIR_YELLOW = 2
};

__device__ inline GPUInterval interval_sub(float a, float b) { return { __fsub_rd(a, b), __fsub_ru(a, b) }; }
__device__ inline GPUInterval interval_add(GPUInterval a, GPUInterval b) { return { __fadd_rd(a.lo, b.lo), __fadd_ru(a.hi, b.hi) }; }
__device__ inline GPUInterval interval_sub(GPUInterval a, GPUInterval b) { return { __fsub_rd(a.lo, b.hi), __fsub_ru(a.hi, b.lo) }; }

__device__ inline GPUInterval interval_mul(GPUInterval a, GPUInterval b) {
    float lo1 = __fmul_rd(a.lo, b.lo); float lo2 = __fmul_rd(a.lo, b.hi);
    float lo3 = __fmul_rd(a.hi, b.lo); float lo4 = __fmul_rd(a.hi, b.hi);
    float hi1 = __fmul_ru(a.lo, b.lo); float hi2 = __fmul_ru(a.lo, b.hi);
    float hi3 = __fmul_ru(a.hi, b.lo); float hi4 = __fmul_ru(a.hi, b.hi);
    return { fminf(fminf(lo1, lo2), fminf(lo3, lo4)), fmaxf(fmaxf(hi1, hi2), fmaxf(hi3, hi4)) };
}

__device__ int orient3d_interval(
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

__device__ inline PairStatus classifyPair(const cuBQL::Triangle& A, const cuBQL::Triangle& B) {
    if (!A.bounds().overlaps(B.bounds())) return PAIR_NO;

    const cuBQL::vec3f A0 = A.a, A1 = A.b, A2 = A.c;
    const cuBQL::vec3f B0 = B.a, B1 = B.b, B2 = B.c;

    float max_A = fmaxf(fmaxf(fabsf(A0.x), fabsf(A0.y)), fmaxf(fabsf(A0.z), fabsf(A1.x)));
    max_A = fmaxf(max_A, fmaxf(fabsf(A1.y), fabsf(A1.z)));
    max_A = fmaxf(max_A, fmaxf(fmaxf(fabsf(A2.x), fabsf(A2.y)), fabsf(A2.z)));
    
    float max_B = fmaxf(fmaxf(fabsf(B0.x), fabsf(B0.y)), fmaxf(fabsf(B0.z), fabsf(B1.x)));
    max_B = fmaxf(max_B, fmaxf(fabsf(B1.y), fabsf(B1.z)));
    max_B = fmaxf(max_B, fmaxf(fmaxf(fabsf(B2.x), fabsf(B2.y)), fabsf(B2.z)));
    
    float max_coord = fmaxf(max_A, max_B);
    float float_grid_resolution = max_coord * 1.1920929e-7f; 
    float eps = 8.0f * (float_grid_resolution * float_grid_resolution * float_grid_resolution);
    if (eps < 1e-7f) eps = 1e-7f;

    int ob0 = orient3d_interval(A0, A1, A2, B0, eps);
    int ob1 = orient3d_interval(A0, A1, A2, B1, eps);
    int ob2 = orient3d_interval(A0, A1, A2, B2, eps);
    if (ob0 != 0 && ob1 != 0 && ob2 != 0) {
        if ((ob0 > 0 && ob1 > 0 && ob2 > 0) || (ob0 < 0 && ob1 < 0 && ob2 < 0)) return PAIR_NO;
    }

    int oa0 = orient3d_interval(B0, B1, B2, A0, eps);
    int oa1 = orient3d_interval(B0, B1, B2, A1, eps);
    int oa2 = orient3d_interval(B0, B1, B2, A2, eps);
    if (oa0 != 0 && oa1 != 0 && oa2 != 0) {
        if ((oa0 > 0 && oa1 > 0 && oa2 > 0) || (oa0 < 0 && oa1 < 0 && oa2 < 0)) return PAIR_NO;
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

// --------------------------------------------------------------------
// KERNELS: Accelerated Vector Clean Structure Passes
// --------------------------------------------------------------------
__global__ void generateBoxes(cuBQL::box3f *boxes, const cuBQL::Triangle *tris, int N) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < N) boxes[i] = tris[i].bounds();
}

__global__ void countAABBOverlapsKernel(int *pairCounts, cuBQL::bvh3f bvhA, const cuBQL::Triangle *triA, const cuBQL::Triangle *triB, int NB) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= NB) return;

    cuBQL::Triangle b = triB[tid];
    cuBQL::box3f query = b.bounds();
    int count = 0;

    cuBQL::fixedBoxQuery::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        for (uint32_t i = 0; i < num; i++) {
            if (triA[ids[i]].bounds().overlaps(query)) {
                count++;
            }
        }
        return CUBQL_CONTINUE_TRAVERSAL;
    }, bvhA, query);

    pairCounts[tid] = count;
}

__global__ void fillAABBOverlapsKernel(int2 *candidatePairs, const int *offsets, cuBQL::bvh3f bvhA, const cuBQL::Triangle *triA, const cuBQL::Triangle *triB, int NB) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= NB) return;

    int wPos = offsets[tid];
    cuBQL::Triangle b = triB[tid];
    cuBQL::box3f query = b.bounds();

    cuBQL::fixedBoxQuery::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        for (uint32_t i = 0; i < num; i++) {
            if (triA[ids[i]].bounds().overlaps(query)) {
                candidatePairs[wPos++] = make_int2((int)ids[i], tid);
            }
        }
        return CUBQL_CONTINUE_TRAVERSAL;
    }, bvhA, query);
}

__global__ void evaluateGeometricPairsKernel(
    IntersectionPair *outPairs, int *outStatuses, const int2 *candidatePairs, 
    const cuBQL::Triangle *triA, const cuBQL::Triangle *triB, int numPairs) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= numPairs) return;

    int2 pair = candidatePairs[tid];
    PairStatus status = classifyPair(triA[pair.x], triB[pair.y]);

    outPairs[tid] = { pair.x, pair.y };
    outStatuses[tid] = (int)status;
}

struct IsTargetPairStatus {
    int targetStatus;
    __host__ __device__ bool operator()(const int &status) const {
        return status == targetStatus;
    }
};

// Expands a 10-bit integer into 30 bits by inserting 2 zeros between each bit
__device__ inline uint32_t dev_expandBits3D(uint32_t v) {
    v = (v * 0x00010001u) & 0xFF0000FFu;
    v = (v * 0x00000101u) & 0x0F00F00Fu;
    v = (v * 0x00000011u) & 0xC30C30C3u;
    v = (v * 0x00000005u) & 0x49249249u;
    return v;
}

__global__ void computeMeshBMortonCodesKernel(
    uint32_t* dMortonCodes, 
    int* dIndices, 
    const cuBQL::Triangle* dB, 
    int NB, 
    cuBQL::box3f sceneBounds) 
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= NB) return;

    cuBQL::vec3f center = (dB[i].a + dB[i].b + dB[i].c) * (1.0f / 3.0f);
    
    float sizeX = sceneBounds.upper.x - sceneBounds.lower.x;
    float sizeY = sceneBounds.upper.y - sceneBounds.lower.y;
    float sizeZ = sceneBounds.upper.z - sceneBounds.lower.z;

    float normX = (sizeX > 1e-6f) ? (center.x - sceneBounds.lower.x) / sizeX : 0.0f;
    float normY = (sizeY > 1e-6f) ? (center.y - sceneBounds.lower.y) / sizeY : 0.0f;
    float normZ = (sizeZ > 1e-6f) ? (center.z - sceneBounds.lower.z) / sizeZ : 0.0f;

    uint32_t x = dev_expandBits3D(static_cast<uint32_t>(fminf(1023.0f, fmaxf(0.0f, normX * 1023.0f))));
    uint32_t y = dev_expandBits3D(static_cast<uint32_t>(fminf(1023.0f, fmaxf(0.0f, normY * 1023.0f))));
    uint32_t z = dev_expandBits3D(static_cast<uint32_t>(fminf(1023.0f, fmaxf(0.0f, normZ * 1023.0f))));
    
    dMortonCodes[i] = (x << 2) | (y << 1) | z;
    dIndices[i] = i; 
}

// OPTIMIZATION: Realizes full physical out-of-place reordering on the GPU
__global__ void permuteMeshBKernel(
    cuBQL::Triangle* dDest, 
    const cuBQL::Triangle* dSrc, 
    const int* dSortedIndices, 
    int NB) 
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= NB) return;
    dDest[i] = dSrc[dSortedIndices[i]];
}

// --------------------------------------------------------------------
// MAIN ENTRY PIPELINE
// --------------------------------------------------------------------
extern "C" void runPartitionedMeshIntersection(
    const cuBQL::Triangle* hA, int NA, 
    const cuBQL::Triangle* hB, int NB, 
    std::vector<IntersectionPair>& hGreenPairs,
    std::vector<IntersectionPair>& hYellowPairs,
    GPUTimingBreakdown& outTimings,
    int pipelineMode,
    int batchSize = 256000) 
{
    int currentDevice = 0;
    if (cudaGetDevice(&currentDevice) == cudaSuccess) {
        cudaSetDevice(currentDevice);
    }

    // 1. Upload both meshes to GPU
    double t0 = cuBQL::getCurrentTime();
    cuBQL::Triangle *dA, *dB_full;
    CUBQL_CUDA_CALL(Malloc(&dA, NA * sizeof(cuBQL::Triangle)));
    CUBQL_CUDA_CALL(Memcpy(dA, hA, NA * sizeof(cuBQL::Triangle), cudaMemcpyDefault));

    CUBQL_CUDA_CALL(Malloc(&dB_full, NB * sizeof(cuBQL::Triangle)));
    CUBQL_CUDA_CALL(Memcpy(dB_full, hB, NB * sizeof(cuBQL::Triangle), cudaMemcpyDefault));
    outTimings.uploadTime = cuBQL::getCurrentTime() - t0;

    double t_exec = cuBQL::getCurrentTime();
    
    // =================================================================
    // OPTIMIZATION: ONE-TIME HARDWARE PHYSICAL PERMUTATION
    // =================================================================
    cuBQL::box3f boundsB;
    for (int i = 0; i < NB; ++i) {
        boundsB.extend(hB[i].bounds());
    }

    uint32_t* dMortonCodes;
    int* dSortedIndices;
    CUBQL_CUDA_CALL(Malloc(&dMortonCodes, NB * sizeof(uint32_t)));
    CUBQL_CUDA_CALL(Malloc(&dSortedIndices, NB * sizeof(int)));

    computeMeshBMortonCodesKernel<<<cuBQL::divRoundUp(NB, 256), 256>>>(
        dMortonCodes, dSortedIndices, dB_full, NB, boundsB
    );
    cudaDeviceSynchronize();

    thrust::device_ptr<uint32_t> keys_ptr(dMortonCodes);
    thrust::device_ptr<int> values_ptr(dSortedIndices);
    thrust::sort_by_key(thrust::device, keys_ptr, keys_ptr + NB, values_ptr);
    CUBQL_CUDA_CALL(Free(dMortonCodes));

    // Physically rewrite Mesh B into perfectly aligned VRAM space
    cuBQL::Triangle* dB_permuted;
    CUBQL_CUDA_CALL(Malloc(&dB_permuted, NB * sizeof(cuBQL::Triangle)));
    permuteMeshBKernel<<<cuBQL::divRoundUp(NB, 256), 256>>>(dB_permuted, dB_full, dSortedIndices, NB);
    cudaDeviceSynchronize();

    CUBQL_CUDA_CALL(Free(dB_full)); 
    dB_full = dB_permuted; // Ready for continuous high-speed layout reading

    // ==========================================
    // PHASE A: BVH Build
    // ==========================================
    double t_bvh_start = cuBQL::getCurrentTime();
    cuBQL::box3f *dBoxes;
    CUBQL_CUDA_CALL(Malloc(&dBoxes, NA * sizeof(cuBQL::box3f)));
    generateBoxes<<<cuBQL::divRoundUp(NA, 256), 256>>>(dBoxes, dA, NA);
    cudaDeviceSynchronize(); 
    
    cuBQL::bvh3f bvh;
    cuBQL::gpuBuilder(bvh, dBoxes, NA, cuBQL::BuildConfig());
    CUBQL_CUDA_CALL(Free(dBoxes));
    cudaDeviceSynchronize(); 
    outTimings.bvhBuildTime = cuBQL::getCurrentTime() - t_bvh_start;

    hGreenPairs.clear();
    hYellowPairs.clear();
    
    // ==========================================
    // OPTIMIZATION: HOIST ENGINE SCRATCH WORKSPACES
    // ==========================================
    int *dPairCounts, *dOffsets;
    CUBQL_CUDA_CALL(Malloc(&dPairCounts, batchSize * sizeof(int)));
    CUBQL_CUDA_CALL(Malloc(&dOffsets, batchSize * sizeof(int)));

    // Set a predictable max pair bounds per batch (e.g., 4x batch scale) to prevent internal allocations
    int maxExpectedPairsPerBatch = batchSize * 4; 
    int2 *dCandidatePairs;
    IntersectionPair *dEvaluatedPairs;
    int *dPairStatuses;
    CUBQL_CUDA_CALL(Malloc(&dCandidatePairs, maxExpectedPairsPerBatch * sizeof(int2)));
    CUBQL_CUDA_CALL(Malloc(&dEvaluatedPairs, maxExpectedPairsPerBatch * sizeof(IntersectionPair)));
    CUBQL_CUDA_CALL(Malloc(&dPairStatuses, maxExpectedPairsPerBatch * sizeof(int)));

    thrust::device_vector<IntersectionPair> dev_green(maxExpectedPairsPerBatch);
    thrust::device_vector<IntersectionPair> dev_yellow(maxExpectedPairsPerBatch);

    // ==========================================
    // PHASE B: Query & Evaluation Loop
    // ==========================================
    double t_query_start = cuBQL::getCurrentTime();
    std::vector<int> hLocalIndices(batchSize);

    for (int bStart = 0; bStart < NB; bStart += batchSize) {
        int currentBatchSize = std::min(batchSize, NB - bStart);

        // OPTIMIZATION: Seamless pointer math layout. Zero gather overhead!
        const cuBQL::Triangle* dB_batch = dB_full + bStart;

        countAABBOverlapsKernel<<<cuBQL::divRoundUp(currentBatchSize, 128), 128>>>(
            dPairCounts, bvh, dA, dB_batch, currentBatchSize
        );
        
        thrust::device_ptr<int> dev_counts(dPairCounts);
        thrust::device_ptr<int> dev_offsets(dOffsets);
        thrust::exclusive_scan(thrust::device, dev_counts, dev_counts + currentBatchSize, dev_offsets);
        
        int totalBatchPairs = 0, lastCount = 0;
        CUBQL_CUDA_CALL(Memcpy(&totalBatchPairs, dOffsets + currentBatchSize - 1, sizeof(int), cudaMemcpyDeviceToHost));
        CUBQL_CUDA_CALL(Memcpy(&lastCount, dPairCounts + currentBatchSize - 1, sizeof(int), cudaMemcpyDeviceToHost));
        totalBatchPairs += lastCount;

        if (totalBatchPairs == 0) continue;

        // Dynamic resize protective buffer if allocations spike beyond capacity
        if (totalBatchPairs > maxExpectedPairsPerBatch) {
            CUBQL_CUDA_CALL(Free(dCandidatePairs)); CUBQL_CUDA_CALL(Free(dEvaluatedPairs)); CUBQL_CUDA_CALL(Free(dPairStatuses));
            maxExpectedPairsPerBatch = totalBatchPairs * 1.2f;
            CUBQL_CUDA_CALL(Malloc(&dCandidatePairs, maxExpectedPairsPerBatch * sizeof(int2)));
            CUBQL_CUDA_CALL(Malloc(&dEvaluatedPairs, maxExpectedPairsPerBatch * sizeof(IntersectionPair)));
            CUBQL_CUDA_CALL(Malloc(&dPairStatuses, maxExpectedPairsPerBatch * sizeof(int)));
            dev_green.resize(maxExpectedPairsPerBatch);
            dev_yellow.resize(maxExpectedPairsPerBatch);
        }

        fillAABBOverlapsKernel<<<cuBQL::divRoundUp(currentBatchSize, 128), 128>>>(
            dCandidatePairs, dOffsets, bvh, dA, dB_batch, currentBatchSize
        );

        evaluateGeometricPairsKernel<<<cuBQL::divRoundUp(totalBatchPairs, 256), 256>>>(
            dEvaluatedPairs, dPairStatuses, dCandidatePairs, dA, dB_batch, totalBatchPairs
        );

        thrust::device_ptr<IntersectionPair> dev_evaluated(dEvaluatedPairs);
        thrust::device_ptr<int> dev_statuses(dPairStatuses);
        thrust::device_ptr<IntersectionPair> dev_green_out(thrust::raw_pointer_cast(dev_green.data()));
        thrust::device_ptr<IntersectionPair> dev_yellow_out(thrust::raw_pointer_cast(dev_yellow.data()));

        auto green_end = thrust::copy_if(thrust::device, dev_evaluated, dev_evaluated + totalBatchPairs, dev_statuses, dev_green_out, IsTargetPairStatus{(int)PAIR_GREEN});
        auto yellow_end = thrust::copy_if(thrust::device, dev_evaluated, dev_evaluated + totalBatchPairs, dev_statuses, dev_yellow_out, IsTargetPairStatus{(int)PAIR_YELLOW});

        int totalGreen = green_end - dev_green_out;
        int totalYellow = yellow_end - dev_yellow_out;

        CUBQL_CUDA_CALL(Memcpy(hLocalIndices.data(), dSortedIndices + bStart, currentBatchSize * sizeof(int), cudaMemcpyDefault));

        if (totalGreen > 0) {
            size_t oldSize = hGreenPairs.size();
            hGreenPairs.resize(oldSize + totalGreen);
            std::vector<IntersectionPair> hTmpGreen(totalGreen);
            CUBQL_CUDA_CALL(Memcpy(hTmpGreen.data(), thrust::raw_pointer_cast(dev_green.data()), totalGreen * sizeof(IntersectionPair), cudaMemcpyDefault));
            for (auto& pair : hTmpGreen) {
                pair.idB = hLocalIndices[pair.idB];
            }
            std::copy(hTmpGreen.begin(), hTmpGreen.end(), hGreenPairs.begin() + oldSize);
        }

        if (totalYellow > 0) {
            size_t oldSize = hYellowPairs.size();
            hYellowPairs.resize(oldSize + totalYellow);
            std::vector<IntersectionPair> hTmpYellow(totalYellow);
            CUBQL_CUDA_CALL(Memcpy(hTmpYellow.data(), thrust::raw_pointer_cast(dev_yellow.data()), totalYellow * sizeof(IntersectionPair), cudaMemcpyDefault));
            for (auto& pair : hTmpYellow) {
                pair.idB = hLocalIndices[pair.idB]; 
            }
            std::copy(hTmpYellow.begin(), hTmpYellow.end(), hYellowPairs.begin() + oldSize);
        }
    }
    cudaDeviceSynchronize(); 

    outTimings.queryTime = cuBQL::getCurrentTime() - t_query_start;
    outTimings.executionTime = cuBQL::getCurrentTime() - t_exec;
    outTimings.downloadTime = 0.0; 

    // Global Structural Clean-up
    CUBQL_CUDA_CALL(Free(dA));
    CUBQL_CUDA_CALL(Free(dB_full));
    CUBQL_CUDA_CALL(Free(dSortedIndices));
    CUBQL_CUDA_CALL(Free(dPairCounts));
    CUBQL_CUDA_CALL(Free(dOffsets));
    CUBQL_CUDA_CALL(Free(dCandidatePairs));
    CUBQL_CUDA_CALL(Free(dEvaluatedPairs));
    CUBQL_CUDA_CALL(Free(dPairStatuses));
    cuBQL::cuda::free(bvh);
}