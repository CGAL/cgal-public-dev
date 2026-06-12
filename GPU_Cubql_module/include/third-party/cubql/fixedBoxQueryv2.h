// SPDX-FileCopyrightText: Copyright (c) 2025 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "cuBQL/traversal/fixedAnyShapeQuery.h"

namespace cuBQL {
  namespace fixedBoxQueryv2 {
    
    /*! This query finds all primitives within a given fixed (ie, never
      changing) axis-aligned cartesian box, starting from a user-specified 
      node index, and calls the provided callback-lambda for each such prim.
    */
    template<typename T, int D, typename Lambda>
    inline __cubql_both
    void forEachPrim(const Lambda &lambdaToCallOnEachPrim,
                     const BinaryBVH<T,D> bvh,
                     const box3f queryBox,
                     uint32_t startNodeIdx = 0,
                     bool dbg=false);
  
    template<typename T, int D, typename Lambda>
    inline __cubql_both
    void forEachLeaf(const Lambda &lambdaToCallOnEachLeaf,
                     const BinaryBVH<T,D> bvh,
                     const box3f queryBox,
                     uint32_t startNodeIdx = 0,
                     bool dbg=false);
  
    template<typename T, int D, int W, typename Lambda>
    inline __cubql_both
    void forEachPrim(const Lambda &lambdaToCallOnEachPrim,
                     const WideBVH<T,D,W> bvh,
                     const box3f queryBox,
                     uint64_t startNodeIdx = 0,
                     bool dbg=false);
  
    template<typename T, int D, int W, typename Lambda>
    inline __cubql_both
    void forEachLeaf(const Lambda &lambdaToCallOnEachLeaf,
                     const WideBVH<T,D,W> bvh,
                     const box3f queryBox,
                     uint64_t startNodeIdx = 0,
                     bool dbg=false);
  


    // ==================================================================
    // IMPLEMENTATION
    // ==================================================================

    template<typename T, int D, typename Lambda>
    inline __cubql_both
    void forEachLeaf(const Lambda &lambdaToCallOnEachLeaf,
                     const BinaryBVH<T,D> bvh,
                     const box3f queryBox,
                     uint32_t startNodeIdx,
                     bool dbg)
    {
      bvh3f::node_t::Admin traversalStack[64], *stackPtr = traversalStack;
      
      // Fixed: Modified to start traversal at user-provided startNodeIdx
      bvh3f::node_t::Admin node = bvh.nodes[startNodeIdx].admin;
      
      while (true) {
        while (true) {
          if (node.count != 0)
            break;

          uint32_t n0Idx = (uint32_t)node.offset+0;
          uint32_t n1Idx = (uint32_t)node.offset+1;
          bvh3f::node_t n0 = bvh.nodes[n0Idx];
          bvh3f::node_t n1 = bvh.nodes[n1Idx];
          bool o0 = queryBox.overlaps(n0.bounds);
          bool o1 = queryBox.overlaps(n1.bounds);
          
          if (o0) {
            if (o1) {
              *stackPtr++ = n1.admin;
            }
            node = n0.admin;
          } else {
            if (o1) {
              node = n1.admin;
            } else {
              node.count = 0;
              break;
            }
          }
        }

        if (node.count != 0) {
          int leafResult
            = lambdaToCallOnEachLeaf(bvh.primIDs+node.offset,(uint32_t)node.count);
          if (leafResult == CUBQL_TERMINATE_TRAVERSAL)
            return;
        }

        if (stackPtr == traversalStack)
          return;
        node = *--stackPtr;
      }
    }

    template<typename T, int D, typename Lambda>
    inline __cubql_both
    void forEachPrim(const Lambda &lambdaToCallOnEachPrim,
                     const BinaryBVH<T,D> bvh,
                     const box3f queryBox,
                     uint32_t startNodeIdx,
                     bool dbg)
    {
      auto leafCode = [&lambdaToCallOnEachPrim]
        (const uint32_t *primIDs, size_t numPrims) -> int
      {
        for (int i=0;i<(int)numPrims;i++)
          if (lambdaToCallOnEachPrim(primIDs[i]) == CUBQL_TERMINATE_TRAVERSAL)
            return CUBQL_TERMINATE_TRAVERSAL;
        return CUBQL_CONTINUE_TRAVERSAL;
      };
      forEachLeaf(leafCode,bvh,queryBox,startNodeIdx,dbg);
    }


    template<typename T, int D, int W, typename Lambda>
    inline __cubql_both
    void forEachLeaf(const Lambda &lambdaToCallOnEachLeaf,
                     const WideBVH<T,D,W> bvh,
                     const box3f queryBox,
                     uint64_t startNodeIdx,
                     bool dbg)
    {
      struct StackEntry {
        uint64_t nodeID:48;
        uint64_t childID:16;
      };
      StackEntry traversalStack[64], *stackPtr = traversalStack;
      StackEntry current;
      
      // Fixed: Initialize to user-provided node ID and look at child index 0
      current.nodeID = startNodeIdx;
      current.childID = 0;

#if 1
      typename WideBVH<T,D,W>::Node::Child child;
      while (true) {
        while (true) {
          child = bvh.nodes[current.nodeID].children[current.childID];
          if (!child.valid
              ||
              !queryBox.overlaps(child.bounds)) {
            if (current.childID+1 < W)
              current.childID = current.childID+1;
            else if (stackPtr > traversalStack)
              current = *--stackPtr;
            else
              return;
          } else if (child.count == 0) {
            if (current.childID+1 < W) {
              current.childID++;
              *stackPtr++ = current;
            }
            current.nodeID = child.offset;
            current.childID = 0;
          } else
            break;
        }

        int leafResult
          = lambdaToCallOnEachLeaf(bvh.primIDs+child.offset,(uint32_t)child.count);
        if (leafResult == CUBQL_TERMINATE_TRAVERSAL)
          return;
        if (current.childID+1 < W)
          current.childID = current.childID+1;
        else if (stackPtr > traversalStack)
          current = *--stackPtr;
        else
          return;
      }
#else
      while (true) {
        auto child = bvh.nodes[current.nodeID].children[current.childID];

        if (!child.valid)
          /* nothing to do here */;
        else if (!queryBox.overlaps(child.bounds))
          /* nothing to do here */;
        else if (child.count != 0) {
          int leafResult
            = lambdaToCallOnEachLeaf(bvh.primIDs+child.offset,(uint32_t)child.count);
          if (leafResult == CUBQL_TERMINATE_TRAVERSAL)
            return;
        } else {
          if (current.childID+1 < W) {
            current.childID++;
            *stackPtr++ = current;
          }
          current.nodeID = child.offset;
          current.childID = 0;
          continue;
        }
        if (current.childID+1 < W)
          current.childID = current.childID+1;
        else if (stackPtr > traversalStack)
          current = *--stackPtr;
        else
          return;
      }
#endif
    }

    template<typename T, int D, int W, typename Lambda>
    inline __cubql_both
    void forEachPrim(const Lambda &lambdaToCallOnEachPrim,
                     const WideBVH<T,D,W> bvh,
                     const box3f queryBox,
                     uint64_t startNodeIdx,
                     bool dbg)
    {
      auto leafCode = [&lambdaToCallOnEachPrim]
        (const uint32_t *primIDs, size_t numPrims) -> int
      {
        for (int i=0;i<(int)numPrims;i++)
          if (lambdaToCallOnEachPrim(primIDs[i]) == CUBQL_TERMINATE_TRAVERSAL)
            return CUBQL_TERMINATE_TRAVERSAL;
        return CUBQL_CONTINUE_TRAVERSAL;
      };
      forEachLeaf(leafCode,bvh,queryBox,startNodeIdx,dbg);
    }
    
  } // ::cubql::fixedBoxQueryv2
} // ::cubql