// SPDX-FileCopyrightText: Copyright (c) 2026
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "cuBQL/builder/cuda/builder_common.h"
#include "cuBQL/builder/cuda/refit.h"
#include <iostream>

namespace cuBQL {
namespace gpuBuilder_v3 {

using namespace cuBQL::gpuBuilder_impl;
using cuBQL::gpuBuilder_impl::BuildState;
using cuBQL::gpuBuilder_impl::NodeState;

// Helper atomic functions for float types using standard CAS loops
__device__ inline float atomicMinFloat(float* address, float val) {
    int* address_as_int = (int*)address;
    int old = *address_as_int, assumed;
    while (val < __int_as_float(old)) {
        assumed = old;
        old = atomicCAS(address_as_int, assumed, __float_as_int(val));
        if (assumed == old) break;
    }
    return __int_as_float(old);
}

__device__ inline float atomicMaxFloat(float* address, float val) {
    int* address_as_int = (int*)address;
    int old = *address_as_int, assumed;
    while (val > __int_as_float(old)) {
        assumed = old;
        old = atomicCAS(address_as_int, assumed, __float_as_int(val));
        if (assumed == old) break;
    }
    return __int_as_float(old);
}

struct PrimState
{
  union {
    struct
    {
      uint64_t primID : 31; //!< prim we're talking about
      uint64_t done : 1;
      uint64_t nodeID : 32; //!< node the given prim is (currently) in.
    };
    uint64_t bits;
  };
};

template <typename T, int D> struct CUBQL_ALIGN(16) TempNode
{
  using box_t = cuBQL::box_t<T, D>;
  union {
    struct
    {
      AtomicBox<box_t> centBounds; // Used temporarily for cell discovery metrics if needed
      uint32_t count;
      uint32_t alreadyMarked;
    } openBranch;
    struct
    {
      uint32_t offset;
      int dim;
      uint32_t tieBreaker;
      float pos;
    } openNode;
    struct
    {
      uint32_t offset;
      uint32_t count;
      uint32_t unused[2];
    } doneNode;
  };
};

template <typename T, int D> struct CUBQL_ALIGN(16) FinalNormalNode
{
  using box_t = cuBQL::box_t<T, D>;

  struct OpenBranch
  {
    box_t bounds;
    uint32_t count;
    uint32_t alreadyMarked;
  };

  struct OpenNode
  {
    uint32_t offset;
    int dim;
    uint32_t tieBreaker;
    float pos;
  };

  struct DoneNode
  {
    uint32_t offset;
    uint32_t count;
    uint32_t unused[2];
  };

  __device__ __host__ FinalNormalNode() {}
  __device__ __host__ ~FinalNormalNode() {}

  union {
    OpenBranch openBranch;
    OpenNode   openNode;
    DoneNode   doneNode;
  };
};

__global__ void clearSkippedNodes(NodeState* nodeStates, uint32_t count) {
  uint32_t idx = threadIdx.x + blockIdx.x * blockDim.x;
  if(idx < count) {
    nodeStates[idx] = DONE_NODE;
  }
}

template <typename T, int D>
__global__ void
flagActiveNodes(uint32_t* flags, const TempNode<T, D>* tempNodes, uint32_t firstActiveNodeID, uint32_t totalGridCells) {
  uint32_t idx = threadIdx.x + blockIdx.x * blockDim.x;
  if(idx >= totalGridCells)
    return;

  uint32_t nodeID = firstActiveNodeID + idx;
  flags[idx] = (tempNodes[nodeID].openBranch.count > 0) ? 1 : 0;
}

// Step A & B: Remap indices and layouts
template <typename T, int D>
__global__ void compactNodesAndRemapPrims(FinalNormalNode<T, D>* dstNodes,
                                          const TempNode<T, D>* srcNodes,
                                          PrimState* primStates,
                                          const uint32_t* flags,
                                          const uint32_t* scannedOffsets,
                                          uint32_t numPrims,
                                          uint32_t firstActiveNodeID,
                                          uint32_t totalGridCells) {
  uint32_t idx = threadIdx.x + blockIdx.x * blockDim.x;

  if(idx < totalGridCells) {
    if(flags[idx]) {
      uint32_t oldNodeID = firstActiveNodeID + idx;
      uint32_t newNodeID = firstActiveNodeID + scannedOffsets[idx];

      auto& src = srcNodes[oldNodeID].openBranch;
      auto& dst = dstNodes[newNodeID].openBranch;

      // Using cuBQL's native type-specific bounding box min/max limits
      for(int d = 0; d < D; ++d) {
        dst.bounds.lower[d] = empty_box_lower<T>();
        dst.bounds.upper[d] = empty_box_upper<T>();
      }
      dst.count = src.count;
      dst.alreadyMarked = src.alreadyMarked;
    }
  }

  uint32_t primStride = blockDim.x * gridDim.x;
  for(uint32_t primID = idx; primID < numPrims; primID += primStride) {
    auto& prim = primStates[primID];
    if(!prim.done && prim.nodeID != (uint32_t)-1) {
      uint32_t cellIdx = prim.nodeID - firstActiveNodeID;
      prim.nodeID = firstActiveNodeID + scannedOffsets[cellIdx];
    }
  }
}

// NEW STEP C: Ultra-fast parallel bounding box calculation using floating point safe CAS loop
template <typename T, int D>
__global__ void tightenCompactedNodes(FinalNormalNode<T, D>* dstNodes,
                                      const PrimState* primStates,
                                      const box_t<T, D>* primBoxes,
                                      uint32_t numPrims) {
  uint32_t primID = threadIdx.x + blockIdx.x * blockDim.x;
  if(primID >= numPrims)
    return;

  const auto& prim = primStates[primID];
  if(!prim.done && prim.nodeID != (uint32_t)-1) {
    box_t<T, D> box = primBoxes[prim.primID];
    auto& node = dstNodes[prim.nodeID].openBranch;

    for(int d = 0; d < D; ++d) {
      // Swapped standard atomic calls with type-safe float CAS functions
      atomicMinFloat((float*)&node.bounds.lower[d], (float)box.get_lower(d));
      atomicMaxFloat((float*)&node.bounds.upper[d], (float)box.get_upper(d));
    }
  }
}

template <typename T, int D>
__global__ void mapPrimsToSubRoots(PrimState* primStates,
                                   TempNode<T, D>* tempNodes,
                                   const box_t<T, D>* primBoxes,
                                   box_t<T, D> globalBox,
                                   uint32_t numPrims,
                                   uint32_t firstActiveNodeID,
                                   int grid_x,
                                   int grid_y,
                                   int grid_z) {
  const int primID = threadIdx.x + blockIdx.x * blockDim.x;
  if(primID >= numPrims)
    return;

  auto& me = primStates[primID];
  me.primID = primID;

  const box_t<T, D> box = primBoxes[primID];
  if(box.get_lower(0) > box.get_upper(0)) {
    me.nodeID = (uint32_t)-1;
    me.done = true;
    return;
  }

  float cx = 0.5f * (box.get_lower(0) + box.get_upper(0));
  float cy = 0.5f * (box.get_lower(1) + box.get_upper(1));
  float cz = 0.5f * (box.get_lower(2) + box.get_upper(2));

  float u = (cx - globalBox.get_lower(0)) / (globalBox.get_upper(0) - globalBox.get_lower(0));
  float v = (cy - globalBox.get_lower(1)) / (globalBox.get_upper(1) - globalBox.get_lower(1));
  float w = (cz - globalBox.get_lower(2)) / (globalBox.get_upper(2) - globalBox.get_lower(2));

  int x_idx = max(0, min((int)(u * grid_x), grid_x - 1));
  int y_idx = max(0, min((int)(v * grid_y), grid_y - 1));
  int z_idx = max(0, min((int)(w * grid_z), grid_z - 1));

  uint32_t cellID = x_idx + (y_idx * grid_x) + (z_idx * grid_x * grid_y);
  uint32_t assignedNodeID = firstActiveNodeID + cellID;

  me.nodeID = assignedNodeID;
  me.done = false;

  atomicAdd(&tempNodes[assignedNodeID].openBranch.count, 1);
}

template <typename T, int D>
void test_speedrun_initialization(const box_t<T, D>* boxes,
                                  int numPrims,
                                  uint32_t numIterations,
                                  box_t<T, D> globalBox,
                                  uint32_t* d_nodeDescendantCounts,
                                  FinalNormalNode<T, D>*& outTempNodes,
                                  PrimState*& outPrimStates,
                                  uint32_t& outTotalAllocatedNodes,
                                  uint32_t& outFirstActiveNodeID,
                                  cudaStream_t s,
                                  GpuMemoryResource& memResource) {
  assert(sizeof(PrimState) == sizeof(uint64_t));

  // 1. Grid geometry generation
  float sides[3] = {globalBox.get_upper(0) - globalBox.get_lower(0), globalBox.get_upper(1) - globalBox.get_lower(1),
                    globalBox.get_upper(2) - globalBox.get_lower(2)};
  int splits[3] = {0, 0, 0};
  for(uint32_t i = 0; i < numIterations; ++i) {
    int longest = (sides[0] >= sides[1] && sides[0] >= sides[2]) ? 0 : (sides[1] >= sides[2]) ? 1 : 2;
    sides[longest] /= 2.0f;
    splits[longest]++;
  }

  int grid_x = 1 << splits[0];
  int grid_y = 1 << splits[1];
  int grid_z = 1 << splits[2];
  uint32_t totalGridCells = grid_x * grid_y * grid_z;

  outFirstActiveNodeID = 1 << numIterations;
  uint32_t initialTotalNodes = outFirstActiveNodeID + totalGridCells;

  // 2. Initial Setup Allocations
  TempNode<T, D>* d_uncompactedNodes = 0;
  NodeState* nodeStates = 0;
  BuildState* buildState = 0;

  size_t allocationSize = max((size_t)initialTotalNodes * 4, (size_t)2 * numPrims);

  _ALLOC(d_uncompactedNodes, allocationSize, s, memResource);
  _ALLOC(nodeStates, allocationSize, s, memResource);
  _ALLOC(outPrimStates, numPrims, s, memResource);
  _ALLOC(buildState, 1, s, memResource);

  CUBQL_CUDA_CALL(MemsetAsync(d_nodeDescendantCounts, 0, allocationSize * sizeof(uint32_t), s));
  CUBQL_CUDA_CALL(MemsetAsync(d_uncompactedNodes, 0, initialTotalNodes * sizeof(TempNode<T, D>), s));
  CUBQL_CUDA_CALL(MemsetAsync(nodeStates, 0, allocationSize * sizeof(NodeState), s));

  clearSkippedNodes<<<cuBQL::divRoundUp((int)outFirstActiveNodeID, 1024), 1024, 0, s>>>(nodeStates,
                                                                                        outFirstActiveNodeID);
  CUBQL_CUDA_CALL(MemsetAsync(nodeStates + outFirstActiveNodeID, OPEN_BRANCH, totalGridCells * sizeof(NodeState), s));

  mapPrimsToSubRoots<<<divRoundUp(numPrims, 1024), 1024, 0, s>>>(
      outPrimStates, d_uncompactedNodes, boxes, globalBox, numPrims, outFirstActiveNodeID, grid_x, grid_y, grid_z);

  // ==================================================================
  // HIGH-SPEED COMPACTION PHASE
  // ==================================================================
  uint32_t* d_flags = 0;
  uint32_t* d_scannedOffsets = 0;
  _ALLOC(d_flags, totalGridCells, s, memResource);
  _ALLOC(d_scannedOffsets, totalGridCells, s, memResource);

  flagActiveNodes<<<cuBQL::divRoundUp((int)totalGridCells, 1024), 1024, 0, s>>>(d_flags, d_uncompactedNodes,
                                                                                outFirstActiveNodeID, totalGridCells);

  size_t tempStorageBytes = 0;
  cub::DeviceScan::ExclusiveSum(nullptr, tempStorageBytes, d_flags, d_scannedOffsets, totalGridCells, s);

  uint8_t* d_tempStorage = nullptr;
  _ALLOC(d_tempStorage, tempStorageBytes, s, memResource);
  cub::DeviceScan::ExclusiveSum((void*)d_tempStorage, tempStorageBytes, d_flags, d_scannedOffsets, totalGridCells, s);

  CUBQL_CUDA_CALL(StreamSynchronize(s));
  uint32_t totalActiveCells = 0, lastFlag = 0, lastOffset = 0;
  CUBQL_CUDA_CALL(MemcpyAsync(&lastFlag, d_flags + totalGridCells - 1, sizeof(uint32_t), cudaMemcpyDeviceToHost, s));
  CUBQL_CUDA_CALL(
      MemcpyAsync(&lastOffset, d_scannedOffsets + totalGridCells - 1, sizeof(uint32_t), cudaMemcpyDeviceToHost, s));

  CUBQL_CUDA_CALL(StreamSynchronize(s));
  totalActiveCells = lastFlag + lastOffset;
  outTotalAllocatedNodes = outFirstActiveNodeID + totalActiveCells;

  _ALLOC(outTempNodes, outTotalAllocatedNodes, s, memResource);
  CUBQL_CUDA_CALL(MemsetAsync(outTempNodes, 0, outFirstActiveNodeID * sizeof(FinalNormalNode<T, D>), s));

  uint32_t executionThreads = max(totalGridCells, (uint32_t)numPrims);
  compactNodesAndRemapPrims<<<cuBQL::divRoundUp((int)executionThreads, 1024), 1024, 0, s>>>(
      outTempNodes, d_uncompactedNodes, outPrimStates, d_flags, d_scannedOffsets, numPrims, outFirstActiveNodeID,
      totalGridCells);

  // ==================================================================
  // RUN FRESH PARALLEL TIGHTENING PASS
  // ==================================================================
  tightenCompactedNodes<<<cuBQL::divRoundUp(numPrims, 1024), 1024, 0, s>>>(outTempNodes, outPrimStates, boxes,
                                                                           numPrims);

  CUBQL_CUDA_CALL(
      MemcpyAsync(&buildState->numNodes, &outTotalAllocatedNodes, sizeof(uint32_t), cudaMemcpyHostToDevice, s));
  CUBQL_CUDA_CALL(StreamSynchronize(s));

  _FREE(d_tempStorage, s, memResource);
  _FREE(d_flags, s, memResource);
  _FREE(d_scannedOffsets, s, memResource);
  _FREE(d_uncompactedNodes, s, memResource);
  _FREE(nodeStates, s, memResource);
  _FREE(buildState, s, memResource);
}

} // namespace gpuBuilder_v3
} // namespace cuBQL