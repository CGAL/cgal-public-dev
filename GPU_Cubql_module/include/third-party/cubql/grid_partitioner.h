// SPDX-FileCopyrightText: Copyright (c) 2026
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "builder_detail.h"

namespace cuBQL {
namespace ext {
namespace grid_forest {
namespace partitioner {

using namespace cuBQL::gpuBuilder_impl;
using cuBQL::builder_detail::atomicMaxFloat;
using cuBQL::builder_detail::atomicMinFloat;

struct CompactedCellState {
    uint32_t activeCellIdx;
    uint32_t primCount;
};

/**
 * @brief Flags non-empty grid cells for prefix sum stream compaction.
 */
__global__ inline void markOccupiedCellsKernel(uint32_t* flags, const uint32_t* cellCounts, uint32_t totalGridCells) {
    uint32_t idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < totalGridCells) {
        flags[idx] = (cellCounts[idx] > 0) ? 1 : 0;
    }
}

/**
 * @brief Initializes empty AABBs and transfers primitive counts for compacted active cells.
 */
template <typename T, int D>
__global__ void initializeCompactedCellBoxesKernel(cuBQL::box_t<T, D>* dstBoxes,
                                                   uint32_t* dstCounts,
                                                   const uint32_t* flags,
                                                   const uint32_t* scannedOffsets,
                                                   const uint32_t* srcCounts,
                                                   uint32_t totalGridCells) {
    uint32_t idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < totalGridCells && flags[idx]) {
        uint32_t compactedIdx   = scannedOffsets[idx];
        dstCounts[compactedIdx] = srcCounts[idx];

        for (int d = 0; d < D; ++d) {
            dstBoxes[compactedIdx].lower[d] = empty_box_lower<T>();
            dstBoxes[compactedIdx].upper[d] = empty_box_upper<T>();
        }
    }
}

/**
 * @brief Bins primitives into spatial cells based on their centroid position relative to the global bounding box.
 */
template <typename T, int D>
__global__ void assignPrimitivesToCellsKernel(uint32_t* primCellIDs,
                                             uint32_t* cellCounts,
                                             const box_t<T, D>* primBoxes,
                                             box_t<T, D> globalBox,
                                             uint32_t numPrims,
                                             int grid_x, int grid_y, int grid_z) {
    const int primID = threadIdx.x + blockIdx.x * blockDim.x;
    if (primID >= (int)numPrims) return;

    const box_t<T, D> box = primBoxes[primID];
    if (box.get_lower(0) > box.get_upper(0)) { // Invalid / degenerate box guard
        primCellIDs[primID] = (uint32_t)-1;
        return;
    }

    // Compute primitive centroid
    float cx = 0.5f * (box.get_lower(0) + box.get_upper(0));
    float cy = 0.5f * (box.get_lower(1) + box.get_upper(1));
    float cz = 0.5f * (box.get_lower(2) + box.get_upper(2));

    // Normalize centroid coordinates relative to global bounds [0.0, 1.0]
    float u = (cx - globalBox.get_lower(0)) / (globalBox.get_upper(0) - globalBox.get_lower(0));
    float v = (cy - globalBox.get_lower(1)) / (globalBox.get_upper(1) - globalBox.get_lower(1));
    float w = (cz - globalBox.get_lower(2)) / (globalBox.get_upper(2) - globalBox.get_lower(2));

    // Clamp cell coordinates to grid boundaries
    int x_idx = max(0, min((int)(u * grid_x), grid_x - 1));
    int y_idx = max(0, min((int)(v * grid_y), grid_y - 1));
    int z_idx = max(0, min((int)(w * grid_z), grid_z - 1));

    uint32_t cellID     = x_idx + (y_idx * grid_x) + (z_idx * grid_x * grid_y);
    primCellIDs[primID] = cellID;

    atomicAdd(&cellCounts[cellID], 1);
}

/**
 * @brief Scatters primitive indices into contiguous cell buffers and tightens cell AABBs.
 * 
 * Note: Cell bounding boxes are dynamically refitted around the full extent of member primitives.
 * Because primitives span non-zero volume across cell boundaries, the resulting cell AABBs can overlap.
 */
template <typename T, int D>
__global__ void scatterPrimitivesAndRefitCellBoxesKernel(uint32_t* dstSortedPrimIDs,
                                                        cuBQL::box_t<T, D>* dstBoxes,
                                                        uint32_t* dynamicCellCounters,
                                                        const uint32_t* nodeOffsets,
                                                        const uint32_t* flags,
                                                        const uint32_t* scannedOffsets,
                                                        const uint32_t* primCellIDs,
                                                        const box_t<T, D>* primBoxes,
                                                        uint32_t numPrims) {
    uint32_t primID = threadIdx.x + blockIdx.x * blockDim.x;
    if (primID >= numPrims) return;

    uint32_t cellID = primCellIDs[primID];
    if (cellID != (uint32_t)-1 && flags[cellID]) {
        uint32_t compactedCellIdx = scannedOffsets[cellID];

        // Scatter primitive index into destination array
        uint32_t localOffset    = atomicAdd(&dynamicCellCounters[compactedCellIdx], 1);
        uint32_t globalWritePos = nodeOffsets[compactedCellIdx] + localOffset;

        dstSortedPrimIDs[globalWritePos] = primID;

        // Atomically expand cell AABB to fit the primitive
        box_t<T, D> box  = primBoxes[primID];
        auto& targetBox = dstBoxes[compactedCellIdx];
        for (int d = 0; d < D; ++d) {
            atomicMinFloat((float*)&targetBox.lower[d], (float)box.get_lower(d));
            atomicMaxFloat((float*)&targetBox.upper[d], (float)box.get_upper(d));
        }
    }
}

/**
 * @brief Performs a direct level-cut spatial partition from root directly to level N.
 * 
 * Partitions primitives directly into a 3D grid layout by bisecting the longest axis across
 * `numIterations`. Occupied grid cells are stream-compacted, and primitive indices are reordered 
 * into contiguous memory. Cell bounding boxes are tightened around member primitives, creating an
 * overlapping spatial partition.
 */
template <typename T, int D>
void partitionPrimitivesToLevelCut(const box_t<T, D>* boxes,
                                   int numPrims,
                                   uint32_t numIterations,
                                   box_t<T, D> globalBox,
                                   cuBQL::box_t<T, D>*& outNodeBoxes,
                                   uint32_t*& outSortedPrimIDs,
                                   uint32_t*& outNodeOffsets,
                                   uint32_t& outTotalActiveCells,
                                   cudaStream_t s,
                                   cuBQL::DeviceMemoryResource& memResource) {
    // Determine split dimensions by repeatedly bisecting along the longest extent
    float sides[3] = {globalBox.get_upper(0) - globalBox.get_lower(0), 
                      globalBox.get_upper(1) - globalBox.get_lower(1),
                      globalBox.get_upper(2) - globalBox.get_lower(2)};
    int splits[3]  = {0, 0, 0};
    for (uint32_t i = 0; i < numIterations; ++i) {
        int longest = (sides[0] >= sides[1] && sides[0] >= sides[2]) ? 0 : (sides[1] >= sides[2]) ? 1 : 2;
        sides[longest] /= 2.0f;
        splits[longest]++;
    }

    int grid_x              = 1 << splits[0];
    int grid_y              = 1 << splits[1];
    int grid_z              = 1 << splits[2];
    uint32_t totalGridCells = grid_x * grid_y * grid_z;

    // Allocate temporary bins and compute centroid cell mappings
    uint32_t* d_initialCellCounts = nullptr;
    uint32_t* d_primCellIDs       = nullptr;
    _ALLOC(d_initialCellCounts, totalGridCells, s, memResource);
    _ALLOC(d_primCellIDs, numPrims, s, memResource);
    CUBQL_CUDA_CALL(MemsetAsync(d_initialCellCounts, 0, totalGridCells * sizeof(uint32_t), s));

    assignPrimitivesToCellsKernel<T, D><<<cuBQL::divRoundUp((uint32_t)numPrims, 1024u), 1024, 0, s>>>(
        d_primCellIDs, d_initialCellCounts, boxes, globalBox, numPrims, grid_x, grid_y, grid_z);

    // Compact occupied cells using CUB exclusive sum scan
    uint32_t* d_flags          = nullptr;
    uint32_t* d_scannedOffsets = nullptr;
    _ALLOC(d_flags, totalGridCells, s, memResource);
    _ALLOC(d_scannedOffsets, totalGridCells, s, memResource);

    markOccupiedCellsKernel<<<cuBQL::divRoundUp(totalGridCells, 1024u), 1024, 0, s>>>(d_flags, d_initialCellCounts, totalGridCells);

    size_t tempStorageBytes = 0;
    cub::DeviceScan::ExclusiveSum(nullptr, tempStorageBytes, d_flags, d_scannedOffsets, totalGridCells, s);
    uint8_t* d_tempStorage = nullptr;
    _ALLOC(d_tempStorage, tempStorageBytes, s, memResource);
    cub::DeviceScan::ExclusiveSum((void*)d_tempStorage, tempStorageBytes, d_flags, d_scannedOffsets, totalGridCells, s);

    uint32_t lastFlag = 0, lastOffset = 0;
    CUBQL_CUDA_CALL(MemcpyAsync(&lastFlag, d_flags + totalGridCells - 1, sizeof(uint32_t), cudaMemcpyDeviceToHost, s));
    CUBQL_CUDA_CALL(MemcpyAsync(&lastOffset, d_scannedOffsets + totalGridCells - 1, sizeof(uint32_t), cudaMemcpyDeviceToHost, s));
    CUBQL_CUDA_CALL(StreamSynchronize(s));

    outTotalActiveCells = lastFlag + lastOffset;

    // Allocate output buffers for compacted cell nodes and sorted primitive offsets
    _ALLOC(outNodeBoxes, outTotalActiveCells, s, memResource);
    _ALLOC(outNodeOffsets, outTotalActiveCells + 1, s, memResource);
    _ALLOC(outSortedPrimIDs, numPrims, s, memResource);

    uint32_t* d_compactedCounts = nullptr;
    uint32_t* d_dynamicCounters = nullptr;
    _ALLOC(d_compactedCounts, outTotalActiveCells, s, memResource);
    _ALLOC(d_dynamicCounters, outTotalActiveCells, s, memResource);
    CUBQL_CUDA_CALL(MemsetAsync(d_dynamicCounters, 0, outTotalActiveCells * sizeof(uint32_t), s));

    initializeCompactedCellBoxesKernel<T, D><<<cuBQL::divRoundUp(totalGridCells, 1024u), 1024, 0, s>>>(
        outNodeBoxes, d_compactedCounts, d_flags, d_scannedOffsets, d_initialCellCounts, totalGridCells);

    size_t scanOffsetBytes = 0;
    cub::DeviceScan::ExclusiveSum(nullptr, scanOffsetBytes, d_compactedCounts, outNodeOffsets, outTotalActiveCells + 1, s);
    uint8_t* d_scanOffsetTemp = nullptr;
    _ALLOC(d_scanOffsetTemp, scanOffsetBytes, s, memResource);
    cub::DeviceScan::ExclusiveSum((void*)d_scanOffsetTemp, scanOffsetBytes, d_compactedCounts, outNodeOffsets, outTotalActiveCells + 1, s);

    // Scatter primitive IDs and refit bounding boxes for active cells
    scatterPrimitivesAndRefitCellBoxesKernel<T, D><<<cuBQL::divRoundUp((uint32_t)numPrims, 1024u), 1024, 0, s>>>(
        outSortedPrimIDs, outNodeBoxes, d_dynamicCounters, outNodeOffsets, d_flags, d_scannedOffsets, d_primCellIDs, boxes, numPrims);

    CUBQL_CUDA_CALL(StreamSynchronize(s));

    // Cleanup workspace allocations
    _FREE(d_scanOffsetTemp, s, memResource);
    _FREE(d_compactedCounts, s, memResource);
    _FREE(d_dynamicCounters, s, memResource);
    _FREE(d_tempStorage, s, memResource);
    _FREE(d_flags, s, memResource);
    _FREE(d_scannedOffsets, s, memResource);
    _FREE(d_primCellIDs, s, memResource);
    _FREE(d_initialCellCounts, s, memResource);
}

} // namespace partitioner
} // namespace grid_forest
} // namespace ext
} // namespace cuBQL