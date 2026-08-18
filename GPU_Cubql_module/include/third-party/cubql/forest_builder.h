// SPDX-FileCopyrightText: Copyright (c) 2026
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "builder_detail.h"

namespace cuBQL {
namespace ext {
namespace grid_forest {
namespace forest_builder {

using namespace cuBQL::gpuBuilder_impl;
using cuBQL::builder_detail::clearPrimStates;
using cuBQL::builder_detail::PrimState;
using cuBQL::builder_detail::TempNode;

template <typename T, int D>
__global__ void initForestState(BuildState* buildState,
                               NodeState* nodeStates,
                               TempNode<T, D>* nodes,
                               uint32_t numCells,
                               uint32_t numInit) {
    const int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= (int)numInit) return;

    if (tid == 0) { buildState->numNodes = numInit; }

    if (tid < (int)numCells) {
        nodeStates[tid]             = OPEN_BRANCH;
        nodes[tid].openBranch.count = 0;
        nodes[tid].openBranch.level = 0;
        nodes[tid].openBranch.centBounds.set_empty();
    } else {
        nodeStates[tid]            = DONE_NODE;
        nodes[tid].doneNode.offset = 0;
        nodes[tid].doneNode.count  = 0;
    }
}

template <typename T, int D>
__global__ void initForestPrims(TempNode<T, D>* nodes,
                                PrimState* primState,
                                const box_t<T, D>* primBoxes,
                                const uint32_t* d_outSortedPrimIDs,
                                const uint32_t* d_outNodeOffsets,
                                uint32_t numCells,
                                uint32_t maxNodes,
                                uint32_t* d_errorFlag) {
    const int cellID = blockIdx.x;
    if (cellID >= (int)numCells) return;

    if (cellID >= maxNodes) {
        if (atomicCAS(d_errorFlag, 0, 1) == 0) {
            printf("[DEVICE ERROR] initForestPrims: cellID %d >= maxNodes %u\n", cellID, maxNodes);
        }
        return;
    }

    uint32_t startPrimIdx = d_outNodeOffsets[cellID];
    uint32_t endPrimIdx   = d_outNodeOffsets[cellID + 1];

    for (uint32_t slot = startPrimIdx + threadIdx.x; slot < endPrimIdx; slot += blockDim.x) {
        uint32_t primID = d_outSortedPrimIDs[slot];
        auto& me        = primState[primID];
        me.primID       = primID;

        const box_t<T, D> box = primBoxes[primID];
        if (box.get_lower(0) <= box.get_upper(0)) {
            me.nodeID = cellID;
            me.done   = false;
            atomicAdd(&nodes[cellID].openBranch.count, 1);
            gpuBuilder_impl::atomic_grow(nodes[cellID].openBranch.centBounds, box.center());
        } else {
            me.nodeID = (uint32_t)-1;
            me.done   = true;
        }
    }
}

template <typename T, int D>
__global__ void selectSplits(BuildState* buildState,
                             NodeState* nodeStates,
                             TempNode<T, D>* nodes,
                             uint32_t numNodes,
                             BuildConfig buildConfig,
                             uint32_t* d_nodeDescendantCounts,
                             uint32_t maxNodes,
                             uint32_t* d_errorFlag) {
    __shared__ int l_newNodeOfs;
    if (threadIdx.x == 0) l_newNodeOfs = 0;
    __syncthreads();

    int* t_nodeOffsetToWrite = 0;
    int  t_localOffsetToAdd  = 0;
    uint32_t currentLevel    = 0;

    while (true) {
        const int nodeID = threadIdx.x + blockIdx.x * blockDim.x;
        if (nodeID >= (int)numNodes) break;

        if (nodeID >= (int)maxNodes) {
            if (atomicCAS(d_errorFlag, 0, 1) == 0) {
                printf("[DEVICE ERROR] selectSplits: Attempted read nodeID %d >= maxNodes %u\n", nodeID, maxNodes);
            }
            break;
        }

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

        if (in.count <= 1 || in.count <= buildConfig.makeLeafThreshold) {
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

            if (widestDim < 0 || widestCtr == widestLo || widestCtr == widestHi) {
                auto& done  = nodes[nodeID].doneNode;
                done.count  = in.count;
                done.offset = (uint32_t)-1;
                nodeState   = DONE_NODE;
            } else {
                auto& open = nodes[nodeID].openNode;
                open.pos   = widestCtr;
                open.dim   = widestDim;

                t_nodeOffsetToWrite = (int*)&open.offset;
                t_localOffsetToAdd  = atomicAdd(&l_newNodeOfs, 2);
                nodeState           = OPEN_NODE;
            }
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
            const int childID = openOffset + side;

            if (childID >= (int)maxNodes) {
                if (atomicCAS(d_errorFlag, 0, 1) == 0) {
                    printf("\n!!! [CRITICAL TRAP] !!!\n"
                           "[DEVICE ERROR] selectSplits: Exploded childID %d >= maxNodes %u!\n", childID, maxNodes);
                }
                break;
            }

            auto& child         = nodes[childID].openBranch;
            child.centBounds.set_empty();
            child.count         = 0;
            child.level         = currentLevel + 1;
            nodeStates[childID] = OPEN_BRANCH;
        }
    }
}

template <typename T, int D>
__global__ void propagateDescendantCounts(uint32_t numNodes,
                                          const NodeState* nodeStates,
                                          const TempNode<T, D>* nodes,
                                          uint32_t* d_nodeDescendantCounts) {
    int nodeID = numNodes - 1 - (threadIdx.x + blockIdx.x * blockDim.x);
    if (nodeID < 0) return;

    if (nodeStates[nodeID] == DONE_NODE) {
        auto& doneNode = nodes[nodeID].doneNode;
        if (doneNode.offset != (uint32_t)-1) {
            uint32_t leftChild  = doneNode.offset;
            uint32_t rightChild = leftChild + 1;

            if (leftChild < numNodes) {
                d_nodeDescendantCounts[nodeID] = d_nodeDescendantCounts[leftChild] + d_nodeDescendantCounts[rightChild];
            }
        }
    }
}

template <typename T, int D>
__global__ void updatePrims_shm(NodeState* nodeStates,
                                TempNode<T, D>* nodes,
                                PrimState* primStates,
                                const box_t<T, D>* primBoxes,
                                int totalAllocationPrims,
                                int nodeBegin,
                                int numPasses,
                                uint32_t maxNodes,
                                uint32_t* d_errorFlag) {
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
            if (primID >= totalAllocationPrims) break;

            const auto me = primStates[primID];
            if (me.done) break;

            if (me.nodeID >= maxNodes) {
                if (atomicCAS(d_errorFlag, 0, 1) == 0) {
                    printf("[DEVICE ERROR] updatePrims_shm: primID %d points to OOB nodeID %u\n", primID, (uint32_t)me.nodeID);
                }
                break;
            }

            const auto ns = nodeStates[me.nodeID];
            if (ns == DONE_NODE) {
                primStates[primID].done = true;
                break;
            }

            const auto split          = nodes[me.nodeID].openNode;
            const box_t<T, D> primBox = primBoxes[me.primID];
            int side                  = 0;
            if (split.dim == -1) {
                side = (atomicAdd(&nodes[me.nodeID].openNode.tieBreaker, 1) & 1);
            } else {
                box_t<T, D> centerBox;
                centerBox.set_empty();
                centerBox.grow(primBox.center());
                const float center = centerBox.get_lower(split.dim);
                side               = (center >= split.pos);
            }
            int newNodeID = split.offset + side;

            if (newNodeID >= (int)maxNodes) {
                if (atomicCAS(d_errorFlag, 0, 1) == 0) {
                    printf("[DEVICE ERROR] updatePrims_shm: primID %d generated OOB target newNodeID %d\n", primID, newNodeID);
                }
                break;
            }

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
            if ((nodeBegin + i) >= (int)maxNodes) {
                if (atomicCAS(d_errorFlag, 0, 1) == 0) {
                    printf("[DEVICE ERROR] updatePrims_shm Shared-Flush: Flush index %d >= maxNodes %u\n", nodeBegin + i, maxNodes);
                }
                return;
            }
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
                            int totalAllocationPrims,
                            uint32_t maxNodes,
                            uint32_t* d_errorFlag) {
    const int primID = threadIdx.x + blockIdx.x * blockDim.x;
    if (primID >= totalAllocationPrims) return;

    const auto me = primStates[primID];
    if (me.done) return;

    if (me.nodeID >= maxNodes) {
        if (atomicCAS(d_errorFlag, 0, 1) == 0) {
            printf("[DEVICE ERROR] updatePrims: primID %d points to OOB nodeID %u\n", primID, (uint32_t)me.nodeID);
        }
        return;
    }

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
        box_t<T, D> centerBox;
        centerBox.set_empty();
        centerBox.grow(primBox.center());
        const float center = centerBox.get_lower(split.dim);
        side               = (center >= split.pos);
    }
    int newNodeID = split.offset + side;

    if (newNodeID >= (int)maxNodes) {
        if (atomicCAS(d_errorFlag, 0, 1) == 0) {
            printf("[DEVICE ERROR] updatePrims: primID %d generated OOB target newNodeID %d\n", primID, newNodeID);
        }
        return;
    }

    auto& myBranch = nodes[newNodeID].openBranch;
    atomicAdd(&myBranch.count, 1);
    gpuBuilder_impl::atomic_grow(myBranch.centBounds, primBox.center());
    primStates[primID].nodeID = newNodeID;
}

template <typename T, int D>
__global__ void writePrimsAndLeafOffsets(TempNode<T, D>* nodes,
                                          uint32_t* bvhItemList,
                                          PrimState* primStates,
                                          int totalAllocationPrims,
                                          uint32_t maxNodes,
                                          uint32_t* d_errorFlag) {
    const int offset = threadIdx.x + blockIdx.x * blockDim.x;
    if (offset >= totalAllocationPrims) return;

    auto& ps            = primStates[offset];
    bvhItemList[offset] = ps.primID;

    if ((int)ps.nodeID < 0) return;

    if (ps.nodeID >= maxNodes) {
        if (atomicCAS(d_errorFlag, 0, 1) == 0) {
            printf("[DEVICE ERROR] writePrimsAndLeafOffsets: Target nodeID %u >= maxNodes %u\n", (uint32_t)ps.nodeID, maxNodes);
        }
        return;
    }

    auto& node = nodes[ps.nodeID];
    atomicMin(&node.doneNode.offset, offset);
}

template <typename T, int D>
__global__ void writeNodes(typename BinaryBVH<T, D>::Node* finalNodes,
                            TempNode<T, D>* tempNodes,
                            int numNodes,
                            uint32_t maxNodes,
                            uint32_t* d_errorFlag) {
    const int nodeID = threadIdx.x + blockIdx.x * blockDim.x;
    if (nodeID >= numNodes) return;

    if (nodeID >= (int)maxNodes) {
        if (atomicCAS(d_errorFlag, 0, 1) == 0) {
            printf("[DEVICE ERROR] writeNodes: nodeID %d >= maxNodes %u\n", nodeID, maxNodes);
        }
        return;
    }

    finalNodes[nodeID].admin.offset = tempNodes[nodeID].doneNode.offset;
    finalNodes[nodeID].admin.count  = tempNodes[nodeID].doneNode.count;
}

template <typename T, int D>
void build_forest(BinaryBVH<T, D>& bvh,
                  const box_t<T, D>* boxes,
                  int numPrims,
                  int originalMaxPrims,
                  uint32_t numCells,
                  const uint32_t* d_outSortedPrimIDs,
                  const uint32_t* d_outNodeOffsets,
                  BuildConfig buildConfig,
                  uint32_t* d_nodeDescendantCounts,
                  cudaStream_t s,
                  GpuMemoryResource& memResource) {
    assert(sizeof(PrimState) == sizeof(uint64_t));

    const uint32_t numInit  = numCells + (numCells & 1u);
    const uint32_t maxNodes = 2u * (uint32_t)originalMaxPrims + numInit + 2u;

    TempNode<T, D>* tempNodes = 0;
    NodeState* nodeStates     = 0;
    PrimState* primStates     = 0;
    BuildState* buildState   = 0;
    uint32_t* d_errorFlag     = 0;

    _ALLOC(tempNodes, maxNodes, s, memResource);
    _ALLOC(nodeStates, maxNodes, s, memResource);
    _ALLOC(primStates, originalMaxPrims, s, memResource);
    _ALLOC(buildState, 1, s, memResource);
    _ALLOC(d_errorFlag, 1, s, memResource);

    CUBQL_CUDA_CALL(MemsetAsync(d_nodeDescendantCounts, 0, maxNodes * sizeof(uint32_t), s));
    CUBQL_CUDA_CALL(MemsetAsync(d_errorFlag, 0, sizeof(uint32_t), s));

    clearPrimStates<<<divRoundUp((uint32_t)originalMaxPrims, 1024U), 1024, 0, s>>>(primStates, originalMaxPrims);

    initForestState<<<divRoundUp(numInit, 1024U), 1024, 0, s>>>(buildState, nodeStates, tempNodes, numCells, numInit);

    initForestPrims<<<(int)numCells, 256, 0, s>>>(
        tempNodes, primStates, boxes, d_outSortedPrimIDs, d_outNodeOffsets, numCells, maxNodes, d_errorFlag);

    int numDone            = 0;
    int numNodes           = 0;
    uint32_t hostErrorFlag = 0;

    cudaEvent_t stateDownloadedEvent;
    CUBQL_CUDA_CALL(EventCreate(&stateDownloadedEvent));

    while (true) {
        CUBQL_CUDA_CALL(MemcpyAsync(&numNodes, &buildState->numNodes, sizeof(numNodes), cudaMemcpyDeviceToHost, s));
        CUBQL_CUDA_CALL(MemcpyAsync(&hostErrorFlag, d_errorFlag, sizeof(hostErrorFlag), cudaMemcpyDeviceToHost, s));

        CUBQL_CUDA_CALL(EventRecord(stateDownloadedEvent, s));
        CUBQL_CUDA_CALL(EventSynchronize(stateDownloadedEvent));

        if (hostErrorFlag > 0 || numNodes >= (int)maxNodes || numNodes == numDone) { break; }

        selectSplits<<<divRoundUp((uint32_t)numNodes, 1024U), 1024, 0, s>>>(
            buildState, nodeStates, tempNodes, numNodes, buildConfig, d_nodeDescendantCounts, maxNodes, d_errorFlag);

        numDone = numNodes;

        if (sizeof(T) * D <= sizeof(float3)) {
            updatePrims_shm<<<divRoundUp((uint32_t)originalMaxPrims, 512U), 512, 0, s>>>(nodeStates, tempNodes, primStates, boxes, originalMaxPrims, numDone, 8, maxNodes, d_errorFlag);
        } else {
            updatePrims<<<divRoundUp((uint32_t)originalMaxPrims, 1024U), 1024, 0, s>>>(nodeStates, tempNodes, primStates, boxes, originalMaxPrims, maxNodes, d_errorFlag);
        }

        CUBQL_CUDA_CALL(StreamSynchronize(s));
    }
    CUBQL_CUDA_CALL(EventDestroy(stateDownloadedEvent));

    if (hostErrorFlag > 0 || numNodes >= (int)maxNodes) {
        _FREE(tempNodes, s, memResource);
        _FREE(nodeStates, s, memResource);
        _FREE(primStates, s, memResource);
        _FREE(buildState, s, memResource);
        _FREE(d_errorFlag, s, memResource);
        return;
    }

    if (numNodes > 0) {
        uint32_t tpb    = 256;
        uint32_t blocks = (numNodes + tpb - 1) / tpb;
        propagateDescendantCounts<<<blocks, tpb, 0, s>>>(numNodes, nodeStates, tempNodes, d_nodeDescendantCounts);
    }

    uint8_t* d_temp_storage     = NULL;
    size_t temp_storage_bytes   = 0;
    PrimState* sortedPrimStates = 0;

    _ALLOC(sortedPrimStates, originalMaxPrims, s, memResource);
    cub::DeviceRadixSort::SortKeys((void*&)d_temp_storage, temp_storage_bytes, (uint64_t*)primStates, (uint64_t*)sortedPrimStates, originalMaxPrims, 32, 64, s);
    _ALLOC(d_temp_storage, temp_storage_bytes, s, memResource);
    cub::DeviceRadixSort::SortKeys((void*&)d_temp_storage, temp_storage_bytes, (uint64_t*)primStates, (uint64_t*)sortedPrimStates, originalMaxPrims, 32, 64, s);
    _FREE(d_temp_storage, s, memResource);

    bvh.numPrims = originalMaxPrims;
    _ALLOC(bvh.primIDs, originalMaxPrims, s, memResource);
    writePrimsAndLeafOffsets<<<divRoundUp((uint32_t)originalMaxPrims, 1024U), 1024, 0, s>>>(tempNodes, bvh.primIDs, sortedPrimStates, originalMaxPrims, maxNodes, d_errorFlag);

    bvh.numNodes = numNodes;
    _ALLOC(bvh.nodes, numNodes, s, memResource);
    writeNodes<<<divRoundUp((uint32_t)numNodes, 1024U), 1024, 0, s>>>(bvh.nodes, tempNodes, numNodes, maxNodes, d_errorFlag);

    _FREE(sortedPrimStates, s, memResource);
    _FREE(tempNodes, s, memResource);
    _FREE(nodeStates, s, memResource);
    _FREE(primStates, s, memResource);
    _FREE(buildState, s, memResource);
    _FREE(d_errorFlag, s, memResource);
}

} // namespace forest_builder
} // namespace grid_forest
} // namespace ext
} // namespace cuBQL