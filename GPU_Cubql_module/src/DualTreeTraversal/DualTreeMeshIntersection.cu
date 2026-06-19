#include "DualTreeMeshIntersection.h"
#include <thrust/device_ptr.h>
#include <thrust/scan.h>
#include <thrust/execution_policy.h>
#include <thrust/fill.h>
#include <thrust/iterator/counting_iterator.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include "../custom_pipeline/GPUPredicatesCheck.h"

// Helper macro matching your project's naming convention
#ifndef CUBQL_CUDA_CALL
#define CUBQL_CUDA_CALL( call )                                                              \
{                                                                                            \
    cudaError_t err = cuda##call;                                                             \
    if( cudaSuccess != err ) {                                                               \
        fprintf(stderr, "CUDA error in file '%s' in line %i : %s.\n",                        \
                __FILE__, __LINE__, cudaGetErrorString( err ) );                             \
                exit(EXIT_FAILURE);                                                          \
    }                                                                                        \
}
#endif

// --------------------------------------------------------------------
// ARCHITECTURAL STRUCTURES
// --------------------------------------------------------------------

struct DualFrontierElement {
    uint32_t nodeID_A;
    uint32_t nodeID_B;
    uint32_t batchIdx; 
};

struct DualTreeScratchpad {
    uint32_t capacity;
    uint32_t maxLeaves;
    uint32_t candidateCapacity;
    
    DualFrontierElement* frontierA;
    DualFrontierElement* frontierB;
    
    int* childCounts;
    int* offsets;
    
    DualFrontierElement* leafPairs;
    int* d_leafCounter;

    int2* d_candidatePairs;
    int* d_pairStatuses;

    void init(uint32_t initialCapacity, uint32_t initialLeaves, uint32_t maxCandidates) {
        capacity = initialCapacity;
        maxLeaves = initialLeaves;
        candidateCapacity = maxCandidates;

        CUBQL_CUDA_CALL(Malloc(&frontierA,   capacity * sizeof(DualFrontierElement)));
        CUBQL_CUDA_CALL(Malloc(&frontierB,   capacity * sizeof(DualFrontierElement)));
        CUBQL_CUDA_CALL(Malloc(&childCounts, capacity * sizeof(int)));
        CUBQL_CUDA_CALL(Malloc(&offsets,     capacity * sizeof(int)));
        CUBQL_CUDA_CALL(Malloc(&leafPairs,   maxLeaves * sizeof(DualFrontierElement)));
        CUBQL_CUDA_CALL(Malloc(&d_leafCounter, sizeof(int)));

        CUBQL_CUDA_CALL(Malloc(&d_candidatePairs, candidateCapacity * sizeof(int2)));
        CUBQL_CUDA_CALL(Malloc(&d_pairStatuses,   candidateCapacity * sizeof(int)));
    }

    void reserve(uint32_t requiredCapacity) {
        if (requiredCapacity <= capacity) return;
        capacity = (uint32_t)(requiredCapacity * 1.5f);
        
        CUBQL_CUDA_CALL(Free(frontierB));
        CUBQL_CUDA_CALL(Free(childCounts));
        CUBQL_CUDA_CALL(Free(offsets));
        
        CUBQL_CUDA_CALL(Malloc(&frontierB,   capacity * sizeof(DualFrontierElement)));
        CUBQL_CUDA_CALL(Malloc(&childCounts, capacity * sizeof(int)));
        CUBQL_CUDA_CALL(Malloc(&offsets,     capacity * sizeof(int)));
    }

    void resizeLeafBuffer(uint32_t requiredLeaves) {
        maxLeaves = requiredLeaves;
        CUBQL_CUDA_CALL(Free(leafPairs));
        CUBQL_CUDA_CALL(Malloc(&leafPairs, maxLeaves * sizeof(DualFrontierElement)));
    }

    void resizeCandidateBuffer(uint32_t requiredCandidates) {
        candidateCapacity = requiredCandidates;
        CUBQL_CUDA_CALL(Free(d_candidatePairs));
        CUBQL_CUDA_CALL(Free(d_pairStatuses));
        CUBQL_CUDA_CALL(Malloc(&d_candidatePairs, candidateCapacity * sizeof(int2)));
        CUBQL_CUDA_CALL(Malloc(&d_pairStatuses,   candidateCapacity * sizeof(int)));
    }

    void freeAll() {
        CUBQL_CUDA_CALL(Free(frontierA));
        CUBQL_CUDA_CALL(Free(frontierB));
        CUBQL_CUDA_CALL(Free(childCounts));
        CUBQL_CUDA_CALL(Free(offsets));
        CUBQL_CUDA_CALL(Free(leafPairs));
        CUBQL_CUDA_CALL(Free(d_leafCounter));
        CUBQL_CUDA_CALL(Free(d_candidatePairs));
        CUBQL_CUDA_CALL(Free(d_pairStatuses));
    }
};

// --------------------------------------------------------------------
// KERNEL DEFINITIONS
// --------------------------------------------------------------------

__global__ void initializeFrontierSeedKernel(
    DualFrontierElement* d_frontierA,
    const uint32_t* outPairsA,
    const uint32_t* outPairsB,
    const uint32_t* markedNodeIndicesA,
    const uint32_t* markedNodeIndicesB,
    int chunkStartBatchIdx,
    int activeBatchesInChunk
) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= activeBatchesInChunk) return;

    int globalBatchIdx = chunkStartBatchIdx + idx;
    d_frontierA[idx].nodeID_A = markedNodeIndicesA[outPairsA[globalBatchIdx]];
    d_frontierA[idx].nodeID_B = markedNodeIndicesB[outPairsB[globalBatchIdx]];
    d_frontierA[idx].batchIdx = (uint32_t)globalBatchIdx;
}

__global__ void countDualTreeExpansionsKernel(
    const DualFrontierElement* d_inFrontier, 
    uint32_t size, 
    int* d_childCounts,
    cuBQL::bvh3f bvhA,
    cuBQL::bvh3f bvhB
) { 
    uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= size) return;

    DualFrontierElement curr = d_inFrontier[tid];
    auto nodeA = bvhA.nodes[curr.nodeID_A];
    auto nodeB = bvhB.nodes[curr.nodeID_B];

    if (nodeA.admin.count > 0 && nodeB.admin.count > 0) {
        d_childCounts[tid] = 0;
        return;
    }

    int validChildren = 0;
    if (nodeA.admin.count == 0 && nodeB.admin.count == 0) {
        uint32_t al = nodeA.admin.offset; uint32_t ar = al + 1;
        uint32_t bl = nodeB.admin.offset; uint32_t br = bl + 1;
        if (bvhA.nodes[al].bounds.overlaps(bvhB.nodes[bl].bounds)) validChildren++;
        if (bvhA.nodes[al].bounds.overlaps(bvhB.nodes[br].bounds)) validChildren++;
        if (bvhA.nodes[ar].bounds.overlaps(bvhB.nodes[bl].bounds)) validChildren++;
        if (bvhA.nodes[ar].bounds.overlaps(bvhB.nodes[br].bounds)) validChildren++;
    } 
    else if (nodeA.admin.count > 0) {
        uint32_t bl = nodeB.admin.offset; uint32_t br = bl + 1;
        if (nodeA.bounds.overlaps(bvhB.nodes[bl].bounds)) validChildren++;
        if (nodeA.bounds.overlaps(bvhB.nodes[br].bounds)) validChildren++;
    } 
    else {
        uint32_t al = nodeA.admin.offset; uint32_t ar = al + 1;
        if (bvhA.nodes[al].bounds.overlaps(nodeB.bounds)) validChildren++;
        if (bvhA.nodes[ar].bounds.overlaps(nodeB.bounds)) validChildren++;
    }
    d_childCounts[tid] = validChildren;
}

__global__ void populateDualTreeFrontierKernel(
    const DualFrontierElement* d_inFrontier, 
    uint32_t size, 
    const int* d_offsets, 
    DualFrontierElement* d_outFrontier, 
    DualFrontierElement* d_leafPairs, 
    int* d_leafCounter,
    uint32_t maxLeaves,
    cuBQL::bvh3f bvhA,
    cuBQL::bvh3f bvhB
) { 
    uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= size) return;

    DualFrontierElement curr = d_inFrontier[tid];
    auto nodeA = bvhA.nodes[curr.nodeID_A];
    auto nodeB = bvhB.nodes[curr.nodeID_B];

    if (nodeA.admin.count > 0 && nodeB.admin.count > 0) {
        int wPos = atomicAdd(d_leafCounter, 1);
        if (wPos < maxLeaves) {
            d_leafPairs[wPos] = curr;
        }
        return;
    }

    int writeIdx = d_offsets[tid];
    if (nodeA.admin.count == 0 && nodeB.admin.count == 0) {
        uint32_t al = nodeA.admin.offset; uint32_t ar = al + 1;
        uint32_t bl = nodeB.admin.offset; uint32_t br = bl + 1;
        if (bvhA.nodes[al].bounds.overlaps(bvhB.nodes[bl].bounds)) d_outFrontier[writeIdx++] = {al, bl, curr.batchIdx};
        if (bvhA.nodes[al].bounds.overlaps(bvhB.nodes[br].bounds)) d_outFrontier[writeIdx++] = {al, br, curr.batchIdx};
        if (bvhA.nodes[ar].bounds.overlaps(bvhB.nodes[bl].bounds)) d_outFrontier[writeIdx++] = {ar, bl, curr.batchIdx};
        if (bvhA.nodes[ar].bounds.overlaps(bvhB.nodes[br].bounds)) d_outFrontier[writeIdx++] = {ar, br, curr.batchIdx};
    } 
    else if (nodeA.admin.count > 0) {
        uint32_t bl = nodeB.admin.offset; uint32_t br = bl + 1;
        if (nodeA.bounds.overlaps(bvhB.nodes[bl].bounds)) d_outFrontier[writeIdx++] = {curr.nodeID_A, bl, curr.batchIdx};
        if (nodeA.bounds.overlaps(bvhB.nodes[br].bounds)) d_outFrontier[writeIdx++] = {curr.nodeID_A, br, curr.batchIdx};
    } 
    else {
        uint32_t al = nodeA.admin.offset; uint32_t ar = al + 1;
        if (bvhA.nodes[al].bounds.overlaps(nodeB.bounds)) d_outFrontier[writeIdx++] = {al, curr.nodeID_B, curr.batchIdx};
        if (bvhA.nodes[ar].bounds.overlaps(nodeB.bounds)) d_outFrontier[writeIdx++] = {ar, curr.nodeID_B, curr.batchIdx};
    }
}

__global__ void evaluateLeafPrimitivePairsKernel(
    const DualFrontierElement* d_leafPairs, 
    int leafCount, 
    int2* d_candidatePairs, 
    int* d_candCounter,
    uint32_t maxCandidates,
    cuBQL::bvh3f bvhA,
    cuBQL::bvh3f bvhB,
    const cuBQL::Triangle* dMeshA,
    const cuBQL::Triangle* dMeshB
) { 
    uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= leafCount) return;

    DualFrontierElement pair = d_leafPairs[tid];
    auto leafA = bvhA.nodes[pair.nodeID_A];
    auto leafB = bvhB.nodes[pair.nodeID_B];

    for (uint32_t i = 0; i < leafA.admin.count; ++i) {
        uint32_t primIdA = bvhA.primIDs[leafA.admin.offset + i];
        cuBQL::box3f boxA = dMeshA[primIdA].bounds();

        for (uint32_t j = 0; j < leafB.admin.count; ++j) {
            uint32_t primIdB = bvhB.primIDs[leafB.admin.offset + j];
            if (boxA.overlaps(dMeshB[primIdB].bounds())) {
                int writePos = atomicAdd(d_candCounter, 1);
                if (writePos < maxCandidates) {
                    d_candidatePairs[writePos] = make_int2((int)primIdA, (int)primIdB);
                }
            }
        }
    }
}

// --------------------------------------------------------------------
// CORE IMPLEMENTATION PIPELINE
// --------------------------------------------------------------------
uint64_t executeDualTreeTraversal(
    int batchMultiplier,
    int totalBatches,
    uint32_t maxDescendantsA,      // Max total primitives across subtrees in A
    uint32_t maxDescendantsB,      // Max total primitives across subtrees in B
    float expectedIntersectionDensity, // Value between (0.0f, 1.0f] matching expected overlap complexity
    const thrust::device_vector<uint32_t>& d_outPairsA,
    const thrust::device_vector<uint32_t>& d_outPairsB,
    const thrust::device_vector<uint32_t>& d_markedNodeIndicesA,
    const thrust::device_vector<uint32_t>& d_markedNodeIndicesB,
    const thrust::device_vector<uint32_t>& d_nodeDescendantCountsA, 
    const thrust::device_vector<uint32_t>& d_nodeDescendantCountsB, 
    uint32_t h_outMarkedCountA,
    uint32_t h_outMarkedCountB,
    const cuBQL::bvh3f& bvhA,
    const cuBQL::bvh3f& bvhB,
    const cuBQL::Triangle* dMeshA,
    const cuBQL::Triangle* dMeshB,
    std::vector<int2>& hGreenPairs,  
    std::vector<int2>& hYellowPairs, 
    IntersectionTimeTracker& outLoopTime
) {
    double tStart = cuBQL::getCurrentTime();

    const uint32_t* ptr_outPairsA          = thrust::raw_pointer_cast(d_outPairsA.data());
    const uint32_t* ptr_outPairsB          = thrust::raw_pointer_cast(d_outPairsB.data());
    const uint32_t* ptr_markedNodeIndicesA = thrust::raw_pointer_cast(d_markedNodeIndicesA.data());
    const uint32_t* ptr_markedNodeIndicesB = thrust::raw_pointer_cast(d_markedNodeIndicesB.data());

    // --- HEURISTIC MEMORY DIMENSIONING ---
    // Ensure density parameter sits in valid bounds
    float alpha = std::max(0.0001f, std::min(expectedIntersectionDensity, 1.0f));
    
    // Calculate adaptive bounds based on maximum subtree products in the chunk window
    uint64_t theoreticalMaxProduct = (uint64_t)maxDescendantsA * (uint64_t)maxDescendantsB;
    uint64_t scaleBoundedAllocation = (uint64_t)(batchMultiplier * theoreticalMaxProduct * alpha);

    // Apply baseline minimum floor bounds to avoid tiny or zero allocations
    uint32_t initialFrontierCapacity = (uint32_t)std::max<uint64_t>(4096, batchMultiplier * 4);
    uint32_t maxLeafCapacity         = (uint32_t)std::max<uint64_t>(4096, scaleBoundedAllocation);
    uint32_t maxCandidateCapacity    = (uint32_t)std::max<uint64_t>(8192, scaleBoundedAllocation * 2);
    
    DualTreeScratchpad allocator;
    allocator.init(initialFrontierCapacity, maxLeafCapacity, maxCandidateCapacity);
    
    outLoopTime.preallocateTimeMs = (cuBQL::getCurrentTime() - tStart) * 1000.0;
    uint64_t totalCandidatePairsProcessed = 0;

    for (int i = 0; i < totalBatches; i += batchMultiplier) {
        int activeBatchesInChunk = std::min(batchMultiplier, totalBatches - i);
        if (activeBatchesInChunk <= 0) continue;

        int zeroReset = 0;
        CUBQL_CUDA_CALL(Memcpy(allocator.d_leafCounter, &zeroReset, sizeof(int), cudaMemcpyHostToDevice));

        // --- SEED LAYER ZERO ---
        uint32_t currentFrontierSize = activeBatchesInChunk;
        int seedBlock = 256;
        int seedGrid  = (activeBatchesInChunk + seedBlock - 1) / seedBlock;
        initializeFrontierSeedKernel<<<seedGrid, seedBlock>>>(
            allocator.frontierA, ptr_outPairsA, ptr_outPairsB, 
            ptr_markedNodeIndicesA, ptr_markedNodeIndicesB, i, activeBatchesInChunk
        );
        
        DualFrontierElement* d_frontierCurrent = allocator.frontierA;
        DualFrontierElement* d_frontierNext    = allocator.frontierB;

        // --- LEVEL-BY-LEVEL TREE DESCENT ---
        while (currentFrontierSize > 0) {
            int block = 256;
            int grid  = (currentFrontierSize + block - 1) / block;

            countDualTreeExpansionsKernel<<<grid, block>>>(d_frontierCurrent, currentFrontierSize, allocator.childCounts, bvhA, bvhB);

            thrust::device_ptr<int> dev_counts(allocator.childCounts);
            thrust::device_ptr<int> dev_offsets(allocator.offsets);
            thrust::exclusive_scan(thrust::device, dev_counts, dev_counts + currentFrontierSize, dev_offsets);

            int lastCount = 0, lastOffset = 0;
            CUBQL_CUDA_CALL(MemcpyAsync(&lastCount, allocator.childCounts + currentFrontierSize - 1, sizeof(int), cudaMemcpyDeviceToHost));
            CUBQL_CUDA_CALL(Memcpy(&lastOffset, allocator.offsets + currentFrontierSize - 1, sizeof(int), cudaMemcpyDeviceToHost));
            uint32_t nextFrontierSize = lastCount + lastOffset;

            if (nextFrontierSize == 0) break;

            if (nextFrontierSize > allocator.capacity) {
                allocator.reserve(nextFrontierSize);
                d_frontierNext = allocator.frontierB; 
                thrust::exclusive_scan(thrust::device, dev_counts, dev_counts + currentFrontierSize, dev_offsets);
            }

            populateDualTreeFrontierKernel<<<grid, block>>>(
                d_frontierCurrent, currentFrontierSize, allocator.offsets, 
                d_frontierNext, allocator.leafPairs, allocator.d_leafCounter, allocator.maxLeaves, bvhA, bvhB
            );

            std::swap(d_frontierCurrent, d_frontierNext);
            currentFrontierSize = nextFrontierSize;
        }

        // --- HOST-SIDE LEAF BUFFER OVERFLOW CHECK ---
        int finalLeafPairsCount = 0;
        CUBQL_CUDA_CALL(Memcpy(&finalLeafPairsCount, allocator.d_leafCounter, sizeof(int), cudaMemcpyDeviceToHost));

        if (finalLeafPairsCount > allocator.maxLeaves) {
            allocator.resizeLeafBuffer((uint32_t)(finalLeafPairsCount * 1.5f));
            i -= batchMultiplier;
            continue;
        }

        // --- CHUNK PRIMITIVE RESOLUTION & COMPACTION ---
        if (finalLeafPairsCount > 0) {
            int leafBlock = 256;
            int leafGrid  = (finalLeafPairsCount + leafBlock - 1) / leafBlock;

            CUBQL_CUDA_CALL(Memcpy(allocator.d_leafCounter, &zeroReset, sizeof(int), cudaMemcpyHostToDevice));

            evaluateLeafPrimitivePairsKernel<<<leafGrid, leafBlock>>>(
                allocator.leafPairs, finalLeafPairsCount, allocator.d_candidatePairs, 
                allocator.d_leafCounter, allocator.candidateCapacity, bvhA, bvhB, dMeshA, dMeshB
            );

            int totalChunkPairs = 0;
            CUBQL_CUDA_CALL(Memcpy(&totalChunkPairs, allocator.d_leafCounter, sizeof(int), cudaMemcpyDeviceToHost));

            // --- HOST-SIDE TRIANGLE CANDIDATE OVERFLOW CHECK ---
            if (totalChunkPairs > allocator.candidateCapacity) {
                allocator.resizeCandidateBuffer((uint32_t)(totalChunkPairs * 1.5f));
                i -= batchMultiplier;
                continue;
            }

            if (totalChunkPairs > 0) {
                totalCandidatePairsProcessed += totalChunkPairs;

                double internalPredicateTimeDummy = 0.0;
                evaluateAndCompactPairs(
                    allocator.d_candidatePairs, allocator.d_pairStatuses, 
                    dMeshA, dMeshB, totalChunkPairs, internalPredicateTimeDummy
                );

                thrust::device_ptr<int2> dev_evaluated(allocator.d_candidatePairs);
                thrust::device_ptr<int> dev_statuses(allocator.d_pairStatuses);

                thrust::device_vector<int2> dev_green(totalChunkPairs);
                thrust::device_vector<int2> dev_yellow(totalChunkPairs);

                thrust::device_ptr<int2> dev_green_out(thrust::raw_pointer_cast(dev_green.data()));
                thrust::device_ptr<int2> dev_yellow_out(thrust::raw_pointer_cast(dev_yellow.data()));

                auto green_end  = thrust::copy_if(thrust::device, dev_evaluated, dev_evaluated + totalChunkPairs, dev_statuses, dev_green_out, IsTargetPairStatus{(int)PAIR_GREEN});
                auto yellow_end = thrust::copy_if(thrust::device, dev_evaluated, dev_evaluated + totalChunkPairs, dev_statuses, dev_yellow_out, IsTargetPairStatus{(int)PAIR_YELLOW});

                int totalGreen  = green_end - dev_green_out;
                int totalYellow = yellow_end - dev_yellow_out;

                if (totalGreen > 0) {
                    size_t oldSize = hGreenPairs.size();
                    hGreenPairs.resize(oldSize + totalGreen);
                    CUBQL_CUDA_CALL(Memcpy(hGreenPairs.data() + oldSize, thrust::raw_pointer_cast(dev_green.data()), totalGreen * sizeof(int2), cudaMemcpyDeviceToHost));
                }
                if (totalYellow > 0) {
                    size_t oldSize = hYellowPairs.size();
                    hYellowPairs.resize(oldSize + totalYellow);
                    CUBQL_CUDA_CALL(Memcpy(hYellowPairs.data() + oldSize, thrust::raw_pointer_cast(dev_yellow.data()), totalYellow * sizeof(int2), cudaMemcpyDeviceToHost));
                }
            }
        }
    }

    allocator.freeAll();
    outLoopTime.totalLoopTimeMs = (cuBQL::getCurrentTime() - tStart) * 1000.0;
    
    return totalCandidatePairsProcessed; 
}