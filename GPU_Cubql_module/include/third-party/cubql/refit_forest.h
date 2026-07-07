// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "cuBQL/builder/cuda/builder_common.h"
#include <iostream>

namespace cuBQL {
  namespace cuda_forest {

    using namespace cuBQL::gpuBuilder_impl;
    
    template<typename T, int D>
    __global__ void
    refit_init_forest(const typename BinaryBVH<T,D>::Node *nodes,
                      uint32_t                            *refitData,
                      int                                  numNodes,
                      uint32_t                             numInit)
    {
      const int nodeID = threadIdx.x + blockIdx.x * blockDim.x;
      if (nodeID >= numNodes) return;
      
      const auto &node = nodes[nodeID];
      
      // If it's a leaf, it has no children to map, so just exit.
      if (node.admin.count > 0 || node.admin.offset == (uint32_t)-1) {
        return;
      }

      // FIX: Skip uninitialized alignment padding slots (e.g., Node 259)
      if (node.admin.offset == 0) {
        return;
      }

      uint32_t c0 = node.admin.offset + 0;
      uint32_t c1 = node.admin.offset + 1;

      // DEVICE DIAGNOSTIC TRAP 1: Out of bounds child pointer initialization
      if (c0 >= (uint32_t)numNodes || c1 >= (uint32_t)numNodes) {
        printf("[DEVICE ERROR] refit_init_forest: Internal nodeID %d points to out-of-bounds children [%u, %u]. Max numNodes is %d\n", 
               nodeID, c0, c1, numNodes);
        return;
      }

      refitData[c0] = nodeID << 1;
      refitData[c1] = nodeID << 1;
    }
    
    template<typename T, int D>
    __global__ void 
    refit_run_forest(BinaryBVH<T,D>    bvh,
                     uint32_t         *refitData,
                     const box_t<T,D> *boxes,
                     uint32_t          numInit)
    {
      int nodeID = threadIdx.x + blockIdx.x * blockDim.x;
      if (nodeID >= bvh.numNodes) return;
      
      typename BinaryBVH<T,D>::Node *node = &bvh.nodes[nodeID];
      
      bool isLeaf = (node->admin.count > 0) || (node->admin.offset == (uint32_t)-1);
      if (!isLeaf) return;

      box_t<T,D> bounds; 
      bounds.set_empty();
      
      if (node->admin.count > 0) {
        uint32_t offset = node->admin.offset;
        
        // DEVICE DIAGNOSTIC TRAP 2: Item List indexing boundary check
        if (offset + node->admin.count > (uint32_t)bvh.numPrims) {
          printf("[DEVICE ERROR] refit_run_forest: Leaf nodeID %d has item offset %u + count %u exceeding bvh.numPrims %d\n",
                 nodeID, offset, node->admin.count, bvh.numPrims);
          return;
        }
        
        for (int i = 0; i < node->admin.count; i++) {
          uint32_t primID = bvh.primIDs[offset + i];
          const box_t<T,D> primBox = boxes[primID];
          bounds.lower = min(bounds.lower, primBox.lower);
          bounds.upper = max(bounds.upper, primBox.upper);
        }
      }

      int parentID = (refitData[nodeID] >> 1);
      while (true) {
        if (nodeID < 0 || nodeID >= bvh.numNodes) {
          printf("[DEVICE ERROR] refit_run_forest: Execution loop tracking nodeID %d went out of bounds (max %d)\n", nodeID, bvh.numNodes);
          break;
        }

        bvh.nodes[nodeID].bounds = bounds;
        __threadfence(); 
        
        if (nodeID < (int)numInit)
          break;

        // DEVICE DIAGNOSTIC TRAP 3: Vertical tree traversal parent index validation
        if (parentID < 0 || parentID >= bvh.numNodes) {
          printf("[DEVICE ERROR] refit_run_forest: NodeID %d attempted upward hop to invalid parentID %d (max %d). refitData base payload is %u\n",
                 nodeID, parentID, bvh.numNodes, refitData[nodeID]);
          break;
        }

        uint32_t refitBits = atomicAdd(&refitData[parentID], 1u);
        if ((refitBits & 1) == 0)
          break; 

        nodeID   = parentID;
        node     = &bvh.nodes[parentID];
        parentID = (refitBits >> 1);
        
        uint32_t c0 = node->admin.offset + 0;
        uint32_t c1 = node->admin.offset + 1;
        
        // DEVICE DIAGNOSTIC TRAP 4: Internal bounding box reduction index validation
        if (c0 >= (uint32_t)bvh.numNodes || c1 >= (uint32_t)bvh.numNodes) {
          printf("[DEVICE ERROR] refit_run_forest: Ascending parent nodeID %d features broken child tracking pointers [%u, %u] (max %d)\n",
                 nodeID, c0, c1, bvh.numNodes);
          break;
        }

        typename BinaryBVH<T,D>::Node l = bvh.nodes[c0];
        typename BinaryBVH<T,D>::Node r = bvh.nodes[c1];
        bounds.lower = min(l.bounds.lower, r.bounds.lower);
        bounds.upper = max(l.bounds.upper, r.bounds.upper);
      }
    }

    template<typename T, int D>
    void refit_forest(BinaryBVH<T,D>    &bvh,
                      const box_t<T,D>  *boxes,
                      uint32_t           numCells,
                      cudaStream_t       s,
                      GpuMemoryResource &memResource)
    {
      int numNodes = bvh.numNodes;

      if (numNodes <= 0 || bvh.nodes == nullptr || bvh.primIDs == nullptr) {
        std::cout << "   [HOST WARNING] refit_forest: Aborting due to invalid/empty BVH structural tracking values." << std::endl;
        return;
      }
      
      const uint32_t numInit = numCells + (numCells & 1u);
      
      uint32_t *refitData = 0;
      memResource.malloc((void**)&refitData, numNodes * sizeof(*refitData), s);
      
      cudaMemsetAsync(refitData, 0, numNodes * sizeof(*refitData), s);
      
      refit_init_forest<T,D><<<divRoundUp(numNodes, 1024), 1024, 0, s>>>(
        bvh.nodes, refitData, numNodes, numInit
      );

      refit_run_forest<T,D><<<divRoundUp(numNodes, 32), 32, 0, s>>>(
        bvh, refitData, boxes, numInit
      );

      memResource.free((void*)refitData, s);
    }
    
  } // namespace cuda_forest
} // namespace cuBQL