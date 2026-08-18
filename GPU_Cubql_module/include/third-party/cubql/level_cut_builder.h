// SPDX-FileCopyrightText: Copyright (c) 2026
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "builder_detail.h"

namespace cuBQL {
namespace ext {
namespace level_cut {

using namespace cuBQL::gpuBuilder_impl;
using cuBQL::builder_detail::PrimState;
using cuBQL::builder_detail::TempNode;

template <typename T, int D>
__global__ void initState(BuildState* buildState,
                           NodeState* nodeStates,
                           TempNode<T, D>* nodes) {
    buildState->numNodes = 2;

    nodeStates[0]             = OPEN_BRANCH;
    nodes[0].openBranch.count = 0;
    nodes[0].openBranch.level = 0;
    nodes[0].openBranch.centBounds.set_empty();

    nodeStates[1]            = DONE_NODE;
    nodes[1].doneNode.offset = 0;
    nodes[1].doneNode.count  = 0;
}

template <typename T, int D>
__global__ void initPrims(TempNode<T, D>* nodes,
                          PrimState* primState,
                          const box_t<T, D>* primBoxes,
                          uint32_t numPrims) {
    const int primID = threadIdx.x + blockIdx.x * blockDim.x;
    if (primID >= (int)numPrims) return;

    auto& me  = primState[primID];
    me.primID = primID;

    const box_t<T, D> box = primBoxes[primID];
    if (box.get_lower(0) <= box.get_upper(0)) {
        me.nodeID = 0;
        me.done   = false;
        atomicAdd(&nodes[0].openBranch.count, 1);
        gpuBuilder_impl::atomic_grow(nodes[0].openBranch.centBounds, box.center());
    } else {
        me.nodeID = (uint32_t)-1;
        me.done   = true;
    }
}

template <typename T, int D>
__global__ void selectSplits(BuildState* buildState,
                             NodeState* nodeStates,
                             TempNode<T, D>* nodes,
                             uint32_t numNodes,
                             BuildConfig buildConfig,
                             uint32_t targetLevel,
                             uint32_t* globalMarkedCounter,
                             uint32_t* tempMarkedIndices,
                             uint32_t* d_nodeDescendantCounts) {
    __shared__ int l_newNodeOfs;
    if (threadIdx.x == 0) l_newNodeOfs = 0;
    __syncthreads();

    int* t_nodeOffsetToWrite = 0;
    int  t_localOffsetToAdd  = 0;
    uint32_t currentLevel    = 0;

    while (true) {
        const int nodeID = threadIdx.x + blockIdx.x * blockDim.x;
        if (nodeID >= (int)numNodes) break;

        NodeState& nodeState = nodeStates[nodeID];
        if (nodeState == DONE_NODE) break;

        if (nodeState == OPEN_NODE) {
            nodeState   = DONE_NODE;
            int offset  = nodes[nodeID].openNode.offset;
            auto& done  = nodes[nodeID].doneNode;
            done.count  = 0;
            done.offset = offset;
            break;
        }

        auto in      = nodes[nodeID].openBranch;
        currentLevel = in.level;

        d_nodeDescendantCounts[nodeID] = in.count;

        bool shouldMark = (currentLevel == targetLevel) ||
                          (currentLevel < targetLevel && in.count <= buildConfig.makeLeafThreshold);

        if (shouldMark) {
            uint32_t slot           = atomicAdd(globalMarkedCounter, 1);
            tempMarkedIndices[slot] = nodeID;
        }

        if (in.count <= buildConfig.makeLeafThreshold) {
            auto& done  = nodes[nodeID].doneNode;
            done.count  = in.count;
            done.offset = (uint32_t)-1;
            nodeState   = DONE_NODE;
        } else {
            float widestWidth = 0.f;
            int   widestDim   = -1;
            float widestLo = 0.f, widestHi = 0.f, widestCtr = 0.f;
#pragma unroll
            for (int d = 0; d < D; d++) {
                float lo    = in.centBounds.get_lower(d);
                float hi    = in.centBounds.get_upper(d);
                float width = hi - lo;
                if (width <= widestWidth) continue;
                float ctr = 0.5f * (hi + lo);

                widestWidth = width;
                widestDim   = d;
                widestLo    = lo;
                widestHi    = hi;
                widestCtr   = ctr;
            }

            auto& open = nodes[nodeID].openNode;
            if (widestDim >= 0) { open.pos = widestCtr; }
            open.dim = (widestDim < 0 || widestCtr == widestLo || widestCtr == widestHi) ? -1 : widestDim;

            t_nodeOffsetToWrite = (int*)&open.offset;
            t_localOffsetToAdd  = atomicAdd(&l_newNodeOfs, 2);
            nodeState           = OPEN_NODE;
        }
        break;
    }
    __syncthreads();
    if (threadIdx.x == 0 && l_newNodeOfs > 0)
        l_newNodeOfs = atomicAdd(&buildState->numNodes, l_newNodeOfs);
    __syncthreads();
    if (t_nodeOffsetToWrite) {
        int openOffset = *t_nodeOffsetToWrite = l_newNodeOfs + t_localOffsetToAdd;
#pragma unroll
        for (int side = 0; side < 2; side++) {
            const int childID   = openOffset + side;
            auto& child         = nodes[childID].openBranch;
            child.centBounds.set_empty();
            child.count         = 0;
            child.level         = currentLevel + 1;
            nodeStates[childID] = OPEN_BRANCH;
        }
    }
}

template <typename T, int D>
__global__ void updatePrims_shm(NodeState* nodeStates,
                                TempNode<T, D>* nodes,
                                PrimState* primStates,
                                const box_t<T, D>* primBoxes,
                                int numPrims,
                                int nodeBegin,
                                int numPasses = 8) {
    enum { numShm = 512 };
    __shared__ AtomicBox<box_t<T, D>> l_boxes[numShm];
    __shared__ int l_count[numShm];
    for (int i = threadIdx.x; i < numShm; i += blockDim.x) {
        l_boxes[i].set_empty();
        l_count[i] = 0;
    }

    __syncthreads();
    for (int pass = 0; pass < numPasses; pass++) {
        while (true) {
            const int primID = threadIdx.x + pass * blockDim.x + numPasses * blockIdx.x * blockDim.x;
            if (primID >= numPrims) break;

            const auto me = primStates[primID];
            if (me.done) break;

            const auto ns = nodeStates[me.nodeID];
            if (ns == DONE_NODE) {
                primStates[primID].done = true;
                break;
            }

            const auto split  = nodes[me.nodeID].openNode;
            const box_t<T, D> primBox = primBoxes[me.primID];
            int side                  = 0;
            if (split.dim == -1) {
                side = (atomicAdd(&nodes[me.nodeID].openNode.tieBreaker, 1) & 1);
            } else {
                const float center = 0.5f * (primBox.get_lower(split.dim) + primBox.get_upper(split.dim));
                side               = (center >= split.pos);
            }
            int newNodeID  = split.offset + side;
            auto& myBranch = nodes[newNodeID].openBranch;
            if (newNodeID - nodeBegin < numShm) {
                gpuBuilder_impl::atomic_grow(l_boxes[newNodeID - nodeBegin], primBox.center());
                atomicAdd(&l_count[newNodeID - nodeBegin], 1);
            } else {
                gpuBuilder_impl::atomic_grow(myBranch.centBounds, primBox.center());
                atomicAdd(&myBranch.count, 1);
            }
            primStates[primID].nodeID = newNodeID;
            break;
        }
    }
    __syncthreads();
    for (int i = threadIdx.x; i < numShm; i += blockDim.x) {
        if (l_count[i] > 0) {
            atomicAdd(&nodes[nodeBegin + i].openBranch.count, l_count[i]);
            gpuBuilder_impl::atomic_grow(nodes[nodeBegin + i].openBranch.centBounds, l_boxes[i]);
        }
    }
}

template <typename T, int D>
__global__ void updatePrims(NodeState* nodeStates,
                            TempNode<T, D>* nodes,
                            PrimState* primStates,
                            const box_t<T, D>* primBoxes,
                            int numPrims) {
    const int primID = threadIdx.x + blockIdx.x * blockDim.x;
    if (primID >= numPrims) return;

    const auto me = primStates[primID];
    if (me.done) return;

    const auto ns = nodeStates[me.nodeID];
    if (ns == DONE_NODE) {
        primStates[primID].done = true;
        return;
    }

    auto& split               = nodes[me.nodeID].openNode;
    const box_t<T, D> primBox = primBoxes[me.primID];
    int side                  = 0;
    if (split.dim == -1) {
        side = (atomicAdd(&split.tieBreaker, 1) & 1);
    } else {
        const float center = 0.5f * (primBox.get_lower(split.dim) + primBox.get_upper(split.dim));
        side               = (center >= split.pos);
    }
    int newNodeID  = split.offset + side;
    auto& myBranch = nodes[newNodeID].openBranch;
    atomicAdd(&myBranch.count, 1);
    gpuBuilder_impl::atomic_grow(myBranch.centBounds, primBox.center());
    primStates[primID].nodeID = newNodeID;
}

template <typename T, int D>
__global__ void writePrimsAndLeafOffsets(TempNode<T, D>* nodes,
                                          uint32_t* bvhItemList,
                                          PrimState* primStates,
                                          int numPrims) {
    const int offset = threadIdx.x + blockIdx.x * blockDim.x;
    if (offset >= numPrims) return;

    auto& ps            = primStates[offset];
    bvhItemList[offset] = ps.primID;

    if ((int)ps.nodeID < 0) return;
    auto& node = nodes[ps.nodeID];
    atomicMin(&node.doneNode.offset, offset);
}

template <typename T, int D>
__global__ void writeNodes(typename BinaryBVH<T, D>::Node* finalNodes,
                            TempNode<T, D>* tempNodes,
                            int numNodes) {
    const int nodeID = threadIdx.x + blockIdx.x * blockDim.x;
    if (nodeID >= numNodes) return;

    finalNodes[nodeID].admin.offset = tempNodes[nodeID].doneNode.offset;
    finalNodes[nodeID].admin.count  = tempNodes[nodeID].doneNode.count;
}

template <typename T, int D>
void build_custom(BinaryBVH<T, D>& bvh,
                 const box_t<T, D>* boxes,
                 int numPrims,
                 BuildConfig buildConfig,
                 uint32_t targetLevel,
                 uint32_t* d_markedNodeIndices,
                 uint32_t* d_nodeDescendantCounts,
                 uint32_t* h_outMarkedCount,
                 cudaStream_t s,
                 GpuMemoryResource& memResource) {
    assert(sizeof(PrimState) == sizeof(uint64_t));

    TempNode<T, D>* tempNodes           = 0;
    NodeState* nodeStates              = 0;
    PrimState* primStates              = 0;
    BuildState* buildState             = 0;
    uint32_t* d_globalMarkedCounter    = 0;

    _ALLOC(tempNodes, 2 * numPrims, s, memResource);
    _ALLOC(nodeStates, 2 * numPrims, s, memResource);
    _ALLOC(primStates, numPrims, s, memResource);
    _ALLOC(buildState, 1, s, memResource);
    _ALLOC(d_globalMarkedCounter, 1, s, memResource);

    CUBQL_CUDA_CALL(MemsetAsync(d_nodeDescendantCounts, 0, 2 * numPrims * sizeof(uint32_t), s));
    CUBQL_CUDA_CALL(MemsetAsync(d_globalMarkedCounter, 0, sizeof(uint32_t), s));

    initState<<<1, 1, 0, s>>>(buildState, nodeStates, tempNodes);
    initPrims<<<divRoundUp(numPrims, 1024), 1024, 0, s>>>(tempNodes, primStates, boxes, numPrims);

    int numDone = 0;
    int numNodes;

    cudaEvent_t stateDownloadedEvent;
    CUBQL_CUDA_CALL(EventCreate(&stateDownloadedEvent));

    while (true) {
        CUBQL_CUDA_CALL(MemcpyAsync(&numNodes, &buildState->numNodes, sizeof(numNodes), cudaMemcpyDeviceToHost, s));
        if (numNodes == numDone) break;

        CUBQL_CUDA_CALL(EventRecord(stateDownloadedEvent, s));
        CUBQL_CUDA_CALL(EventSynchronize(stateDownloadedEvent));

        selectSplits<<<divRoundUp(numNodes, 1024), 1024, 0, s>>>(
            buildState, nodeStates, tempNodes, numNodes, buildConfig,
            targetLevel, d_globalMarkedCounter, d_markedNodeIndices, d_nodeDescendantCounts);

        numDone = numNodes;

        if (sizeof(T) * D <= sizeof(float3)) {
            updatePrims_shm<<<divRoundUp(numPrims, 512), 512, 0, s>>>(nodeStates, tempNodes, primStates, boxes, numPrims, numDone);
        } else {
            updatePrims<<<divRoundUp(numPrims, 1024), 1024, 0, s>>>(nodeStates, tempNodes, primStates, boxes, numPrims);
        }
    }
    CUBQL_CUDA_CALL(EventDestroy(stateDownloadedEvent));

    if (h_outMarkedCount) {
        CUBQL_CUDA_CALL(MemcpyAsync(h_outMarkedCount, d_globalMarkedCounter, sizeof(uint32_t), cudaMemcpyDeviceToHost, s));
        CUBQL_CUDA_CALL(StreamSynchronize(s));
    }

    uint8_t* d_temp_storage     = NULL;
    size_t temp_storage_bytes   = 0;
    PrimState* sortedPrimStates = 0;

    _ALLOC(sortedPrimStates, numPrims, s, memResource);
    cub::DeviceRadixSort::SortKeys((void*&)d_temp_storage, temp_storage_bytes, (uint64_t*)primStates, (uint64_t*)sortedPrimStates, numPrims, 32, 64, s);
    _ALLOC(d_temp_storage, temp_storage_bytes, s, memResource);
    cub::DeviceRadixSort::SortKeys((void*&)d_temp_storage, temp_storage_bytes, (uint64_t*)primStates, (uint64_t*)sortedPrimStates, numPrims, 32, 64, s);
    _FREE(d_temp_storage, s, memResource);

    bvh.numPrims = numPrims;
    _ALLOC(bvh.primIDs, numPrims, s, memResource);
    writePrimsAndLeafOffsets<<<divRoundUp(numPrims, 1024), 1024, 0, s>>>(tempNodes, bvh.primIDs, sortedPrimStates, numPrims);

    bvh.numNodes = numNodes;
    _ALLOC(bvh.nodes, numNodes, s, memResource);
    writeNodes<<<divRoundUp(numNodes, 1024), 1024, 0, s>>>(bvh.nodes, tempNodes, numNodes);

    _FREE(sortedPrimStates, s, memResource);
    _FREE(tempNodes, s, memResource);
    _FREE(nodeStates, s, memResource);
    _FREE(primStates, s, memResource);
    _FREE(buildState, s, memResource);
    _FREE(d_globalMarkedCounter, s, memResource);
}

} // namespace level_cut
} // namespace ext
} // namespace cuBQL