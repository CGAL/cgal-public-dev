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

    struct PrimState {
      union {
        struct {
          uint64_t primID:31; //!< prim we're talking about
          uint64_t done  : 1;
          uint64_t nodeID:32; //!< node the given prim is (currently) in.
        };
        uint64_t bits;
      };
    };

    template<typename T, int D>
    struct CUBQL_ALIGN(16) TempNode {
      using box_t = cuBQL::box_t<T,D>;
      union {
        struct {
          AtomicBox<box_t> centBounds;
          uint32_t         count;
          uint32_t         alreadyMarked; 
        } openBranch;
        struct {
          uint32_t offset;
          int      dim;
          uint32_t tieBreaker;
          float    pos;
        } openNode;
        struct {
          uint32_t offset;
          uint32_t count;
          uint32_t unused[2];
        } doneNode;
      };
    };
    
    // Marks the skipped top levels of the virtual tree as DONE_NODE
    __global__ void clearSkippedNodes(
        NodeState* nodeStates,
        uint32_t   count) 
    {
        uint32_t idx = threadIdx.x + blockIdx.x * blockDim.x;
        if (idx < count) {
            nodeStates[idx] = DONE_NODE; 
        }
    }

    // Flag active vs empty grid sub-roots
    template<typename T, int D>
    __global__ void flagActiveNodes(
        uint32_t            *flags,
        const TempNode<T,D> *tempNodes,
        uint32_t             firstActiveNodeID,
        uint32_t             totalGridCells)
    {
        uint32_t idx = threadIdx.x + blockIdx.x * blockDim.x;
        if (idx >= totalGridCells) return;
        
        uint32_t nodeID = firstActiveNodeID + idx;
        flags[idx] = (tempNodes[nodeID].openBranch.count > 0) ? 1 : 0;
    }

    // High-speed parallel remapping kernel
    template<typename T, int D>
    __global__ void compactNodesAndRemapPrims(
        TempNode<T,D>       *dstNodes,
        const TempNode<T,D> *srcNodes,
        PrimState           *primStates,
        const uint32_t      *flags,
        const uint32_t      *scannedOffsets,
        uint32_t             numPrims,
        uint32_t             firstActiveNodeID,
        uint32_t             totalGridCells)
    {
        uint32_t idx = threadIdx.x + blockIdx.x * blockDim.x;
        
        // Step A: Compact the TempNodes (Grid Cells)
        if (idx < totalGridCells) {
            if (flags[idx]) {
                uint32_t oldNodeID = firstActiveNodeID + idx;
                uint32_t newNodeID = firstActiveNodeID + scannedOffsets[idx];
                dstNodes[newNodeID] = srcNodes[oldNodeID];
            }
        }

        // Step B: Update Primitive Map Links Concurrently
        uint32_t primStride = blockDim.x * gridDim.x;
        for (uint32_t primID = idx; primID < numPrims; primID += primStride) {
            auto &prim = primStates[primID];
            if (!prim.done && prim.nodeID != (uint32_t)-1) {
                uint32_t cellIdx = prim.nodeID - firstActiveNodeID;
                prim.nodeID = firstActiveNodeID + scannedOffsets[cellIdx];
            }
        }
    }

    // High-speed O(1) quantization pass mapping primitives into forest sub-roots
    template<typename T, int D>
    __global__ void mapPrimsToSubRoots(
        PrimState        *primStates,
        TempNode<T,D>    *tempNodes,
        const box_t<T,D> *primBoxes,
        box_t<T,D>        globalBox,
        uint32_t          numPrims,
        uint32_t          firstActiveNodeID,
        int               grid_x,
        int               grid_y,
        int               grid_z)
    {
        const int primID = threadIdx.x + blockIdx.x * blockDim.x;
        if (primID >= numPrims) return;

        auto &me = primStates[primID];
        me.primID = primID;

        const box_t<T,D> box = primBoxes[primID];
        if (box.get_lower(0) > box.get_upper(0)) {
            me.nodeID = (uint32_t)-1;
            me.done = true;
            return;
        }

        // Calculate primitive center
        float cx = 0.5f * (box.get_lower(0) + box.get_upper(0));
        float cy = 0.5f * (box.get_lower(1) + box.get_upper(1));
        float cz = 0.5f * (box.get_lower(2) + box.get_upper(2));

        // Normalize coordinates relative to the scene's global bounding box
        float u = (cx - globalBox.get_lower(0)) / (globalBox.get_upper(0) - globalBox.get_lower(0));
        float v = (cy - globalBox.get_lower(1)) / (globalBox.get_upper(1) - globalBox.get_lower(1));
        float w = (cz - globalBox.get_lower(2)) / (globalBox.get_upper(2) - globalBox.get_lower(2));

        // Clamp positions tightly into grid bounds
        int x_idx = max(0, min((int)(u * grid_x), grid_x - 1));
        int y_idx = max(0, min((int)(v * grid_y), grid_y - 1));
        int z_idx = max(0, min((int)(w * grid_z), grid_z - 1));

        // 3D Grid flattening to linear cell index
        uint32_t cellID = x_idx + (y_idx * grid_x) + (z_idx * grid_x * grid_y);
        uint32_t assignedNodeID = firstActiveNodeID + cellID;

        me.nodeID = assignedNodeID;
        me.done = false;

        // Atomically accumulate local bounding metrics directly on the sub-root slot
        atomicAdd(&tempNodes[assignedNodeID].openBranch.count, 1);
        gpuBuilder_impl::atomic_grow(tempNodes[assignedNodeID].openBranch.centBounds, box.center());
    }

    /**
     * Upgraded Initialization Pipeline with Extensive Debug Infrastructure.
     */
    template<typename T, int D>
    void test_speedrun_initialization(
               const box_t<T,D>  *boxes,
               int                numPrims,
               uint32_t           numIterations,            
               box_t<T,D>         globalBox,                
               uint32_t          *d_nodeDescendantCounts,   
               // --- Exposed Device References Out ---
               TempNode<T,D>    *&outTempNodes,         
               PrimState        *&outPrimStates,        
               uint32_t          &outTotalAllocatedNodes,
               uint32_t          &outFirstActiveNodeID,
               cudaStream_t       s,
               GpuMemoryResource &memResource)
    {
      //std::cout << "[DBUG] => Entering test_speedrun_initialization | Prims: " << numPrims << " | Iters: " << numIterations << std::endl;
      assert(sizeof(PrimState) == sizeof(uint64_t));

      // 1. Grid geometry generation
      float sides[3] = {
          globalBox.get_upper(0) - globalBox.get_lower(0),
          globalBox.get_upper(1) - globalBox.get_lower(1),
          globalBox.get_upper(2) - globalBox.get_lower(2)
      };
      int splits[3] = {0, 0, 0};
      for (uint32_t i = 0; i < numIterations; ++i) {
          int longest = (sides[0] >= sides[1] && sides[0] >= sides[2]) ? 0 : 
                        (sides[1] >= sides[2]) ? 1 : 2;
          sides[longest] /= 2.0f;
          splits[longest]++;
      }

      int grid_x = 1 << splits[0];
      int grid_y = 1 << splits[1];
      int grid_z = 1 << splits[2];
      uint32_t totalGridCells = grid_x * grid_y * grid_z;

      outFirstActiveNodeID = 1 << numIterations; 
      uint32_t initialTotalNodes = outFirstActiveNodeID + totalGridCells;
     // std::cout << "[DBUG] Grid Sizes: [" << grid_x << ", " << grid_y << ", " << grid_z << "] Total Cells: " << totalGridCells << std::endl;

      // 2. Initial Setup Allocations
      TempNode<T,D> *d_uncompactedNodes = 0;
      NodeState     *nodeStates = 0;
      BuildState    *buildState = 0;

      size_t allocationSize = max((size_t)initialTotalNodes * 4, (size_t)2 * numPrims);
    //  std::cout << "[DBUG] Requesting core base allocations (Size: " << allocationSize << ")..." << std::endl;
      
      _ALLOC(d_uncompactedNodes, allocationSize, s, memResource);
      _ALLOC(nodeStates, allocationSize, s, memResource);
      _ALLOC(outPrimStates, numPrims, s, memResource);
      _ALLOC(buildState, 1, s, memResource);
    //  std::cout << "[DBUG] Core allocations successfully complete." << std::endl;

      CUBQL_CUDA_CALL(MemsetAsync(d_nodeDescendantCounts, 0, allocationSize * sizeof(uint32_t), s));
      CUBQL_CUDA_CALL(MemsetAsync(d_uncompactedNodes, 0, initialTotalNodes * sizeof(TempNode<T,D>), s));
      CUBQL_CUDA_CALL(MemsetAsync(nodeStates, 0, allocationSize * sizeof(NodeState), s));

     // std::cout << "[DBUG] Launching clearSkippedNodes..." << std::endl;
      clearSkippedNodes<<<cuBQL::divRoundUp((int)outFirstActiveNodeID, 1024), 1024, 0, s>>>(nodeStates, outFirstActiveNodeID);
      CUBQL_CUDA_CALL(MemsetAsync(nodeStates + outFirstActiveNodeID, OPEN_BRANCH, totalGridCells * sizeof(NodeState), s));

   //   std::cout << "[DBUG] Launching mapPrimsToSubRoots..." << std::endl;
      mapPrimsToSubRoots<<<divRoundUp(numPrims, 1024), 1024, 0, s>>>(
          outPrimStates, d_uncompactedNodes, boxes, globalBox, numPrims, 
          outFirstActiveNodeID, grid_x, grid_y, grid_z
      );

      // ==================================================================
      // HIGH-SPEED COMPACTION & EMPTY NODE DELETION PHASE
      // ==================================================================
      uint32_t *d_flags = 0;
      uint32_t *d_scannedOffsets = 0;
      //std::cout << "[DBUG] Allocating structural flags/scans workspace vectors..." << std::endl;
      _ALLOC(d_flags, totalGridCells, s, memResource);
      _ALLOC(d_scannedOffsets, totalGridCells, s, memResource);

     // std::cout << "[DBUG] Launching flagActiveNodes..." << std::endl;
      flagActiveNodes<<<cuBQL::divRoundUp((int)totalGridCells, 1024), 1024, 0, s>>>(
          d_flags, d_uncompactedNodes, outFirstActiveNodeID, totalGridCells
      );

      //std::cout << "[DBUG] Computing size sizing metrics using cub::DeviceScan..." << std::endl;
      size_t tempStorageBytes = 0;
      cub::DeviceScan::ExclusiveSum(nullptr, tempStorageBytes, d_flags, d_scannedOffsets, totalGridCells, s);
      
      uint8_t *d_tempStorage = nullptr;
      _ALLOC(d_tempStorage, tempStorageBytes, s, memResource);
      cub::DeviceScan::ExclusiveSum((void*)d_tempStorage, tempStorageBytes, d_flags, d_scannedOffsets, totalGridCells, s);

      //std::cout << "[DBUG] Synchronizing explicitly prior to reading back device properties..." << std::endl;
      CUBQL_CUDA_CALL(StreamSynchronize(s));
      cudaError_t err = cudaDeviceSynchronize();
      if(err != cudaSuccess) { std::cout << "[FATAL] Pipeline broken before reading back metrics: " << cudaGetErrorString(err) << std::endl; }

      uint32_t totalActiveCells = 0;
      uint32_t lastFlag = 0, lastOffset = 0;
      //std::cout << "[DBUG] Asynchronously processing DeviceToHost metrics..." << std::endl;
      CUBQL_CUDA_CALL(MemcpyAsync(&lastFlag, d_flags + totalGridCells - 1, sizeof(uint32_t), cudaMemcpyDeviceToHost, s));
      CUBQL_CUDA_CALL(MemcpyAsync(&lastOffset, d_scannedOffsets + totalGridCells - 1, sizeof(uint32_t), cudaMemcpyDeviceToHost, s));
      
      //std::cout << "[DBUG] Synchronizing metrics readout streams..." << std::endl;
      CUBQL_CUDA_CALL(StreamSynchronize(s));
      totalActiveCells = lastFlag + lastOffset;
      //std::cout << "[DBUG] Metrics pulled => totalActiveCells: " << totalActiveCells << std::endl;

      outTotalAllocatedNodes = outFirstActiveNodeID + totalActiveCells;

      //std::cout << "[DBUG] Allocating final compacted outTempNodes array size: " << outTotalAllocatedNodes << std::endl;
      _ALLOC(outTempNodes, outTotalAllocatedNodes, s, memResource);
      CUBQL_CUDA_CALL(MemsetAsync(outTempNodes, 0, outFirstActiveNodeID * sizeof(TempNode<T,D>), s));

      uint32_t executionThreads = max(totalGridCells, (uint32_t)numPrims);
      //std::cout << "[DBUG] Executing compaction kernel across executionThreads: " << executionThreads << std::endl;
      compactNodesAndRemapPrims<<<cuBQL::divRoundUp((int)executionThreads, 1024), 1024, 0, s>>>(
          outTempNodes, d_uncompactedNodes, outPrimStates, d_flags, d_scannedOffsets,
          numPrims, outFirstActiveNodeID, totalGridCells
      );

      //std::cout << "[DBUG] Writing tracking states inside buildState payload..." << std::endl;
      CUBQL_CUDA_CALL(MemcpyAsync(&buildState->numNodes, &outTotalAllocatedNodes, sizeof(uint32_t), cudaMemcpyHostToDevice, s));
      
      //std::cout << "[DBUG] Syncing post compaction execution..." << std::endl;
      CUBQL_CUDA_CALL(StreamSynchronize(s));
      err = cudaDeviceSynchronize();
      if(err != cudaSuccess) { std::cout << "[FATAL] Compaction kernel triggered an internal fault: " << cudaGetErrorString(err) << std::endl; }

      //std::cout << "[DBUG] Freeing internal dynamic allocations..." << std::endl;
      _FREE(d_tempStorage, s, memResource);
      _FREE(d_flags, s, memResource);
      _FREE(d_scannedOffsets, s, memResource);
      _FREE(d_uncompactedNodes, s, memResource);
      _FREE(nodeStates, s, memResource);
      _FREE(buildState, s, memResource);

      //std::cout << "[DBUG] Exiting test_speedrun_initialization safely." << std::endl;
    }

  } // ::cuBQL::gpuBuilder_v3
} // ::cuBQL