#define CUBQL_GPU_BUILDER_IMPLEMENTATION 1
#include "cuBQL/bvh.h"
#include "cuBQL/traversal/fixedBoxQuery.h"
#include "cuBQL/queries/triangleData/boxInsideOutsideIntersects.h"
#include <thrust/device_ptr.h>
#include <thrust/scan.h>
#include <thrust/copy.h>
#include <thrust/execution_policy.h>
#include <thrust/device_vector.h>
#include <vector>
#include <iostream>

using cuBQL::Triangle;
using cuBQL::vec3i;
using cuBQL::vec3f;
using cuBQL::box3f;
using cuBQL::bvh3f;
using cuBQL::divRoundUp;

struct VoxelContentMetrics {
    int voxelIndex;
    int countA;
    int countB;
    int firstPairOffset; 
};

struct VoxelTimingBreakdown {
    double uploadTime = 0.0;
    double bvhBuildTimeA = 0.0;
    double bvhBuildTimeB = 0.0;
    double countingPassTime = 0.0;
    double prefixScanTime = 0.0;
    double fillingPassTime = 0.0;
    double downloadTime = 0.0;
    double totalExecutionTime = 0.0;
};

// --------------------------------------------------------------------
// KERNELS
// --------------------------------------------------------------------

__global__ void fillBoundsVoxelizer(box3f *d_bounds, int numTriangles, const Triangle *d_tris) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= numTriangles) return;
    d_bounds[tid] = d_tris[tid].bounds();
}

// Pass 1: Optimized 1D Flat Voxel Counting Engine
__global__ void countVoxelPrimitivesKernel(
    int *voxelCountsA, 
    int *voxelCountsB, 
    vec3i dims, 
    box3f worldBounds, 
    bvh3f bvhA, 
    bvh3f bvhB,
    int numCells) 
{
    int voxelLinearIdx = threadIdx.x + blockIdx.x * blockDim.x;
    if (voxelLinearIdx >= numCells) return;

    // Unflatten the 1D index back to structured grid discrete coords
    int ix = voxelLinearIdx % dims.x;
    int iy = (voxelLinearIdx / dims.x) % dims.y;
    int iz = voxelLinearIdx / (dims.x * dims.y);

    vec3f f0 = vec3f(ix, iy, iz) / vec3f(dims);
    vec3f f1 = vec3f(ix + 1, iy + 1, iz + 1) / vec3f(dims);
    box3f queryBox { worldBounds.lerp(f0), worldBounds.lerp(f1) };

    int countA = 0;
    cuBQL::fixedBoxQuery::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        countA += num;
        return CUBQL_CONTINUE_TRAVERSAL;
    }, bvhA, queryBox);

    int countB = 0;
    cuBQL::fixedBoxQuery::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        countB += num;
        return CUBQL_CONTINUE_TRAVERSAL;
    }, bvhB, queryBox);

    voxelCountsA[voxelLinearIdx] = countA;
    voxelCountsB[voxelLinearIdx] = countB;
}

// Pass 2: Optimized 1D Flat Linear Array Packing
__global__ void fillVoxelPrimitivesKernel(
    int *d_flatListA, 
    int *d_flatListB,
    const int *offsetsA, 
    const int *offsetsB,
    vec3i dims, 
    box3f worldBounds, 
    bvh3f bvhA, 
    bvh3f bvhB,
    int numCells) 
{
    int voxelLinearIdx = threadIdx.x + blockIdx.x * blockDim.x;
    if (voxelLinearIdx >= numCells) return;

    int ix = voxelLinearIdx % dims.x;
    int iy = (voxelLinearIdx / dims.x) % dims.y;
    int iz = voxelLinearIdx / (dims.x * dims.y);

    vec3f f0 = vec3f(ix, iy, iz) / vec3f(dims);
    vec3f f1 = vec3f(ix + 1, iy + 1, iz + 1) / vec3f(dims);
    box3f queryBox { worldBounds.lerp(f0), worldBounds.lerp(f1) };

    int wPosA = offsetsA[voxelLinearIdx];
    cuBQL::fixedBoxQuery::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        for (uint32_t i = 0; i < num; i++) {
            d_flatListA[wPosA++] = (int)ids[i];
        }
        return CUBQL_CONTINUE_TRAVERSAL;
    }, bvhA, queryBox);

    int wPosB = offsetsB[voxelLinearIdx];
    cuBQL::fixedBoxQuery::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        for (uint32_t i = 0; i < num; i++) {
            d_flatListB[wPosB++] = (int)ids[i];
        }
        return CUBQL_CONTINUE_TRAVERSAL;
    }, bvhB, queryBox);
}

bvh3f buildLocalVoxelizerBVH(int numTriangles, const Triangle *dT_raw) {
    box3f *d_boxes;
    CUBQL_CUDA_CALL(Malloc(&d_boxes, numTriangles * sizeof(box3f)));
    fillBoundsVoxelizer<<<divRoundUp(numTriangles, 256), 256>>>(d_boxes, numTriangles, dT_raw);
    cudaDeviceSynchronize();

    bvh3f bvh;
    ::cuBQL::gpuBuilder(bvh, d_boxes, numTriangles);
    cudaDeviceSynchronize();
    CUBQL_CUDA_CALL(Free(d_boxes));
    return bvh;
}

// --------------------------------------------------------------------
// EXTERNAL INTERACTION INTERFACE
// --------------------------------------------------------------------
extern "C" void runDualPassVoxelizer(
    const Triangle* hA, int NA,
    const Triangle* hB, int NB,
    vec3i dims, // Now forced to strict N x N x N externally or internally
    std::vector<VoxelContentMetrics>& hVoxelMetrics,
    std::vector<int>& hFlatTrianglesA,
    std::vector<int>& hFlatTrianglesB,
    VoxelTimingBreakdown& outTimings) 
{
    double t_start = cuBQL::getCurrentTime();
    
    // FORCE STRICT N x N x N GRID REGARDLESS OF BOUNDING BOX ASPECT RATIO
    int numCells = dims.x * dims.y * dims.z;

    // 1. Calculate Combined Simulation World Spatial Bounds
    box3f worldBounds;
    for (int i = 0; i < NA; i++) worldBounds.extend(hA[i].bounds());
    for (int i = 0; i < NB; i++) worldBounds.extend(hB[i].bounds());

    // 2. Allocate and upload structural meshes
    double t0 = cuBQL::getCurrentTime();
    Triangle *dA, *dB;
    CUBQL_CUDA_CALL(Malloc(&dA, NA * sizeof(Triangle)));
    CUBQL_CUDA_CALL(Malloc(&dB, NB * sizeof(Triangle)));
    CUBQL_CUDA_CALL(Memcpy(dA, hA, NA * sizeof(Triangle), cudaMemcpyDefault));
    CUBQL_CUDA_CALL(Memcpy(dB, hB, NB * sizeof(Triangle), cudaMemcpyDefault));
    outTimings.uploadTime = cuBQL::getCurrentTime() - t0;

    // 3. Build Static Acceleration Structures
    double t_bvhA = cuBQL::getCurrentTime();
    bvh3f bvhA = buildLocalVoxelizerBVH(NA, dA);
    outTimings.bvhBuildTimeA = cuBQL::getCurrentTime() - t_bvhA;

    double t_bvhB = cuBQL::getCurrentTime();
    bvh3f bvhB = buildLocalVoxelizerBVH(NB, dB);
    outTimings.bvhBuildTimeB = cuBQL::getCurrentTime() - t_bvhB;

    // Allocate spatial matrix trackers
    int *dVoxelCountsA, *dVoxelCountsB;
    int *dVoxelOffsetsA, *dVoxelOffsetsB;
    CUBQL_CUDA_CALL(Malloc(&dVoxelCountsA, numCells * sizeof(int)));
    CUBQL_CUDA_CALL(Malloc(&dVoxelCountsB, numCells * sizeof(int)));
    CUBQL_CUDA_CALL(Malloc(&dVoxelOffsetsA, numCells * sizeof(int)));
    CUBQL_CUDA_CALL(Malloc(&dVoxelOffsetsB, numCells * sizeof(int)));

    // Setup highly parallel 1D execution blocks (256 threads per block)
    int blockSize = 256;
    int gridBlocks = divRoundUp(numCells, blockSize);

    // ================================================================
    // PASS 1: Execute The Counting Stage
    // ================================================================
    double t_count = cuBQL::getCurrentTime();
    countVoxelPrimitivesKernel<<<gridBlocks, blockSize>>>(
        dVoxelCountsA, dVoxelCountsB, dims, worldBounds, bvhA, bvhB, numCells
    );
    cudaDeviceSynchronize();
    outTimings.countingPassTime = cuBQL::getCurrentTime() - t_count;

    // ================================================================
    // BRIDGE: High-Speed Prefix Scan To Allocate Output Space
    // ================================================================
    double t_scan = cuBQL::getCurrentTime();
    thrust::device_ptr<int> dev_countsA(dVoxelCountsA);
    thrust::device_ptr<int> dev_countsB(dVoxelCountsB);
    thrust::device_ptr<int> dev_offsetsA(dVoxelOffsetsA);
    thrust::device_ptr<int> dev_offsetsB(dVoxelOffsetsB);

    thrust::exclusive_scan(thrust::device, dev_countsA, dev_countsA + numCells, dev_offsetsA);
    thrust::exclusive_scan(thrust::device, dev_countsB, dev_countsB + numCells, dev_offsetsB);

    // Read terminal index offsets to determine absolute packing requirements
    int totalElementsA = 0, lastCountA = 0;
    int totalElementsB = 0, lastCountB = 0;
    CUBQL_CUDA_CALL(Memcpy(&totalElementsA, dVoxelOffsetsA + numCells - 1, sizeof(int), cudaMemcpyDeviceToHost));
    CUBQL_CUDA_CALL(Memcpy(&lastCountA, dVoxelCountsA + numCells - 1, sizeof(int), cudaMemcpyDeviceToHost));
    totalElementsA += lastCountA;

    CUBQL_CUDA_CALL(Memcpy(&totalElementsB, dVoxelOffsetsB + numCells - 1, sizeof(int), cudaMemcpyDeviceToHost));
    CUBQL_CUDA_CALL(Memcpy(&lastCountB, dVoxelCountsB + numCells - 1, sizeof(int), cudaMemcpyDeviceToHost));
    totalElementsB += lastCountB;
    outTimings.prefixScanTime = cuBQL::getCurrentTime() - t_scan;

    // Allocate perfectly sized contiguous linear array spaces
    int *d_flatListA = nullptr;
    int *d_flatListB = nullptr;
    if (totalElementsA > 0) CUBQL_CUDA_CALL(Malloc(&d_flatListA, totalElementsA * sizeof(int)));
    if (totalElementsB > 0) CUBQL_CUDA_CALL(Malloc(&d_flatListB, totalElementsB * sizeof(int)));

    // ================================================================
    // PASS 2: Execute The Packing Stage
    // ================================================================
    double t_fill = cuBQL::getCurrentTime();
    fillVoxelPrimitivesKernel<<<gridBlocks, blockSize>>>(
        d_flatListA, d_flatListB, dVoxelOffsetsA, dVoxelOffsetsB, dims, worldBounds, bvhA, bvhB, numCells
    );
    cudaDeviceSynchronize();
    outTimings.fillingPassTime = cuBQL::getCurrentTime() - t_fill;

    // ================================================================
    // DOWNLOAD: Recover Structured Metadata Map back to the Host Vector
    // ================================================================
    double t_down = cuBQL::getCurrentTime();
    std::vector<int> h_countsA(numCells);
    std::vector<int> h_countsB(numCells);
    std::vector<int> h_offsetsA(numCells);

    CUBQL_CUDA_CALL(Memcpy(h_countsA.data(), dVoxelCountsA, numCells * sizeof(int), cudaMemcpyDeviceToHost));
    CUBQL_CUDA_CALL(Memcpy(h_countsB.data(), dVoxelCountsB, numCells * sizeof(int), cudaMemcpyDeviceToHost));
    CUBQL_CUDA_CALL(Memcpy(h_offsetsA.data(), dVoxelOffsetsA, numCells * sizeof(int), cudaMemcpyDeviceToHost));

    hVoxelMetrics.resize(numCells);
    for (int i = 0; i < numCells; i++) {
        hVoxelMetrics[i] = { i, h_countsA[i], h_countsB[i], h_offsetsA[i] };
    }

    if (totalElementsA > 0) {
        hFlatTrianglesA.resize(totalElementsA);
        CUBQL_CUDA_CALL(Memcpy(hFlatTrianglesA.data(), d_flatListA, totalElementsA * sizeof(int), cudaMemcpyDeviceToHost));
    }
    if (totalElementsB > 0) {
        hFlatTrianglesB.resize(totalElementsB);
        CUBQL_CUDA_CALL(Memcpy(hFlatTrianglesB.data(), d_flatListB, totalElementsB * sizeof(int), cudaMemcpyDeviceToHost));
    }
    outTimings.downloadTime = cuBQL::getCurrentTime() - t_down;

    // Cleanup Device allocations
    cuBQL::cuda::free(bvhA);
    cuBQL::cuda::free(bvhB);
    CUBQL_CUDA_CALL(Free(dA));
    CUBQL_CUDA_CALL(Free(dB));
    CUBQL_CUDA_CALL(Free(dVoxelCountsA));
    CUBQL_CUDA_CALL(Free(dVoxelCountsB));
    CUBQL_CUDA_CALL(Free(dVoxelOffsetsA));
    CUBQL_CUDA_CALL(Free(dVoxelOffsetsB));
    if (d_flatListA) CUBQL_CUDA_CALL(Free(d_flatListA));
    if (d_flatListB) CUBQL_CUDA_CALL(Free(d_flatListB));

    outTimings.totalExecutionTime = cuBQL::getCurrentTime() - t_start;
}