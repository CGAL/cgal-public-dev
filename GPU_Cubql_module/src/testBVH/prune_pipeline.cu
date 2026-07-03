#include <cuda_runtime.h>
#include <thrust/scan.h>
#include <thrust/execution_policy.h>

#include "samples/common/loadOBJ.h"
#include "cuBQL/bvh.h"
#include "cuBQL/builder/cuda.h"

// 1. Mark active cells using atomic flags to prevent write contention if pairs are duplicated
__global__ void markActiveCellsKernel_Opt(const uint32_t* d_intersectPairs, int numIntersections, uint32_t* d_flags) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < numIntersections) {
        uint32_t cellId = d_intersectPairs[idx];
        d_flags[cellId] = 1; 
    }
}

// 2. SMART PARALLEL COPY: 1 Thread Block per Cell.
// No binary search. Threads inside the block cooperate to copy the cell's primitives.
__global__ void blockParallelCompactPrimitivesKernel(
    int oldTotalActiveCells,
    const uint32_t* d_flags,
    const uint32_t* d_newIndices,
    const uint32_t* oldOffsets,
    const uint32_t* oldPrimIDs,
    const uint32_t* newOffsets,
    uint32_t* newPrimIDs
) {
    // blockIdx.x acts as our cell iterator. Grid-stride over blocks in case there are >65k cells
    for (int oldId = blockIdx.x; oldId < oldTotalActiveCells; oldId += gridDim.x) {
        if (d_flags[oldId] == 1) {
            uint32_t newId = d_newIndices[oldId];
            
            uint32_t oldStart = oldOffsets[oldId];
            uint32_t primCount = oldOffsets[oldId + 1] - oldStart;
            uint32_t newStart = newOffsets[newId];
            
            // All 256 threads in this block work together!
            // Thread 0 copies prim 0, 256, 512...
            // Thread 1 copies prim 1, 257, 513...
            for (uint32_t i = threadIdx.x; i < primCount; i += blockDim.x) {
                newPrimIDs[newStart + i] = oldPrimIDs[oldStart + i];
            }
        }
    }
}

__global__ void reindexPairsKernel_Opt(uint32_t* d_intersectPairs, int numIntersections, const uint32_t* d_newIndices) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < numIntersections) {
        uint32_t oldId = d_intersectPairs[idx];
        d_intersectPairs[idx] = d_newIndices[oldId];
    }
}

__global__ void computeCellSizesKernel_Opt(
    int oldTotalActiveCells,
    const uint32_t* d_flags,
    const uint32_t* d_newIndices,
    const uint32_t* oldOffsets,
    uint32_t* d_newCellPrimCounts
) {
    int oldId = blockIdx.x * blockDim.x + threadIdx.x;
    if (oldId < oldTotalActiveCells) {
        if (d_flags[oldId] == 1) {
            uint32_t newId = d_newIndices[oldId];
            d_newCellPrimCounts[newId] = oldOffsets[oldId + 1] - oldOffsets[oldId];
        }
    }
}

void parallelPruneAndReindexAll(
    uint32_t* d_intersectPairsA,
    int numIntersections,
    uint32_t*& d_outSortedPrimIDsA,
    uint32_t*& d_outNodeOffsetsA,
    uint32_t& outTotalActiveCellsA,
    int& outNumPrimsA,
    cudaStream_t s,
    cuBQL::DeviceMemoryResource& memResource
) {
    if (numIntersections == 0) {
        outTotalActiveCellsA = 0;
        outNumPrimsA = 0;
        return;
    }

    uint32_t oldCells = outTotalActiveCellsA;

    uint32_t* d_flagsA = nullptr;
    uint32_t* d_newIndicesA = nullptr;
    cudaMallocAsync(&d_flagsA, oldCells * sizeof(uint32_t), s);
    cudaMallocAsync(&d_newIndicesA, oldCells * sizeof(uint32_t), s);
    cudaMemsetAsync(d_flagsA, 0, oldCells * sizeof(uint32_t), s);

    int block = 256;
    int gridPairs = (numIntersections + block - 1) / block;
    int gridCells = (oldCells + block - 1) / block;

    markActiveCellsKernel_Opt<<<gridPairs, block, 0, s>>>(d_intersectPairsA, numIntersections, d_flagsA);

    // This scan creates the "unique indices" you mentioned, telling us exactly who goes where
    thrust::exclusive_scan(thrust::cuda::par.on(s), d_flagsA, d_flagsA + oldCells, d_newIndicesA);

    uint32_t* d_newCellPrimCounts = nullptr;
    uint32_t* d_packedOffsets = nullptr;
    cudaMallocAsync(&d_newCellPrimCounts, (oldCells + 1) * sizeof(uint32_t), s);
    cudaMallocAsync(&d_packedOffsets, (oldCells + 1) * sizeof(uint32_t), s);
    cudaMemsetAsync(d_newCellPrimCounts, 0, (oldCells + 1) * sizeof(uint32_t), s);

    computeCellSizesKernel_Opt<<<gridCells, block, 0, s>>>(
        oldCells, d_flagsA, d_newIndicesA, d_outNodeOffsetsA, d_newCellPrimCounts
    );

    // This scan calculates the exact destination start index for each cell's chunk of primitives
    thrust::exclusive_scan(thrust::cuda::par.on(s), d_newCellPrimCounts, d_newCellPrimCounts + oldCells + 1, d_packedOffsets);

    uint32_t* d_packedPrimIDs = nullptr;
    cudaMallocAsync(&d_packedPrimIDs, outNumPrimsA * sizeof(uint32_t), s);

    // Launch the smart block-parallel copy kernel
    // We launch up to min(oldCells, 32768) blocks to ensure we don't exceed max grid size, 
    // but give plenty of blocks to saturate the GPU.
    int copyGrid = (oldCells < 32768) ? oldCells : 32768;
    blockParallelCompactPrimitivesKernel<<<copyGrid, block, 0, s>>>(
        oldCells, d_flagsA, d_newIndicesA, d_outNodeOffsetsA, d_outSortedPrimIDsA, d_packedOffsets, d_packedPrimIDs
    );

    reindexPairsKernel_Opt<<<gridPairs, block, 0, s>>>(d_intersectPairsA, numIntersections, d_newIndicesA);

    cudaFreeAsync(d_outSortedPrimIDsA, s);
    cudaFreeAsync(d_outNodeOffsetsA, s);

    uint32_t finalCellsCount, finalPrimsCount;
    cudaMemcpyAsync(&finalCellsCount, &d_newIndicesA[oldCells - 1], sizeof(uint32_t), cudaMemcpyDeviceToHost, s);
    cudaMemcpyAsync(&finalPrimsCount, &d_packedOffsets[oldCells], sizeof(uint32_t), cudaMemcpyDeviceToHost, s);
    
    cudaStreamSynchronize(s); 

    d_outNodeOffsetsA = d_packedOffsets;
    d_outSortedPrimIDsA = d_packedPrimIDs;
    outTotalActiveCellsA = finalCellsCount + 1;
    outNumPrimsA = (int)finalPrimsCount;

    cudaFreeAsync(d_flagsA, s);
    cudaFreeAsync(d_newIndicesA, s);
    cudaFreeAsync(d_newCellPrimCounts, s);
}