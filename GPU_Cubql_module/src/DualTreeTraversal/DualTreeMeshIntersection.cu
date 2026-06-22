#include "DualTreeMeshIntersection.h"
#include <thrust/device_ptr.h>
#include <thrust/scan.h>
#include <thrust/execution_policy.h>
#include <thrust/fill.h>
#include <thrust/iterator/counting_iterator.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include "../custom_pipeline/GPUPredicatesCheck.h"

#ifndef CUBQL_CUDA_CALL
#define CUBQL_CUDA_CALL( call )                                                              \
{                                                                                            \
    cudaError_t err = cuda##call;                                                             \
    if( cudaSuccess != err ) {                                                               \
        fprintf(stderr, "\n[CUDA CRITICAL ERROR] in file '%s' at line %i : %s.\n",            \
                __FILE__, __LINE__, cudaGetErrorString( err ) );                             \
                exit(EXIT_FAILURE);                                                          \
    }                                                                                        \
}
#endif

#define CHECK_CUDA_KERNEL_STATE(kernelName)                                                   \
{                                                                                            \
    cudaDeviceSynchronize();                                                                 \
    cudaError_t error = cudaGetLastError();                                                   \
    if(error != cudaSuccess) {                                                               \
        fprintf(stderr, "\n[KERNEL LAUNCH FAILURE] %s dropped error at line %i: %s\n",       \
                kernelName, __LINE__, cudaGetErrorString(error));                            \
                exit(EXIT_FAILURE);                                                          \
    }                                                                                        \
}

struct DualFrontierElement {
    uint32_t nodeID_A;
    uint32_t nodeID_B;
    uint32_t batchIdx; 
};

struct IsTargetPairStatusTwo {
    int target;
    __host__ __device__ bool operator()(const int status) const {
        return status == target;
    }
};

// --------------------------------------------------------------------
// CORE traversal KERNELS IMPLEMENTATION
// --------------------------------------------------------------------

__global__ void initializeFrontierSeedKernel(
    DualFrontierElement* d_frontierA,
    const uint32_t* outPairsA, const uint32_t* outPairsB,
    const uint32_t* markedNodeIndicesA, const uint32_t* markedNodeIndicesB,
    int chunkStartBatchIdx, int activeBatchesInChunk
) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= activeBatchesInChunk) return;

    d_frontierA[idx].nodeID_A = markedNodeIndicesA[outPairsA[idx]];
    d_frontierA[idx].nodeID_B = markedNodeIndicesB[outPairsB[idx]];
    d_frontierA[idx].batchIdx = (uint32_t)(chunkStartBatchIdx + idx);
}

__global__ void countDualTreeExpansionsKernel(
    const DualFrontierElement* d_inFrontier, uint32_t size, int* d_childCounts,
    cuBQL::bvh3f bvhA, cuBQL::bvh3f bvhB
) { 
    uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= size) return;

    DualFrontierElement curr = d_inFrontier[tid];
    auto nodeA = bvhA.nodes[curr.nodeID_A];
    auto nodeB = bvhB.nodes[curr.nodeID_B];

    bool isLeafA = (nodeA.admin.count > 0);
    bool isLeafB = (nodeB.admin.count > 0);

    // If both are leaves, traversal stops here; no children to explore
    if (isLeafA && isLeafB) {
        d_childCounts[tid] = 0;
        return;
    }

    int validChildren = 0;
    if (!isLeafA && !isLeafB) {
        uint32_t al = nodeA.admin.offset; uint32_t ar = al + 1;
        uint32_t bl = nodeB.admin.offset; uint32_t br = bl + 1;
        if (bvhA.nodes[al].bounds.overlaps(bvhB.nodes[bl].bounds)) validChildren++;
        if (bvhA.nodes[al].bounds.overlaps(bvhB.nodes[br].bounds)) validChildren++;
        if (bvhA.nodes[ar].bounds.overlaps(bvhB.nodes[bl].bounds)) validChildren++;
        if (bvhA.nodes[ar].bounds.overlaps(bvhB.nodes[br].bounds)) validChildren++;
    } 
    else if (isLeafA) { 
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
    const DualFrontierElement* d_inFrontier, uint32_t size, const int* d_offsets, 
    DualFrontierElement* d_outFrontier, DualFrontierElement* d_leafPairs, 
    int* d_leafCounter, uint32_t maxLeaves, uint32_t frontierCapacity,
    cuBQL::bvh3f bvhA, cuBQL::bvh3f bvhB
) { 
    uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= size) return;

    DualFrontierElement curr = d_inFrontier[tid];
    auto nodeA = bvhA.nodes[curr.nodeID_A];
    auto nodeB = bvhB.nodes[curr.nodeID_B];

    bool isLeafA = (nodeA.admin.count > 0);
    bool isLeafB = (nodeB.admin.count > 0);

    // Hit the base intersection case: target collected for primitive processing
    if (isLeafA && isLeafB) {
        int wPos = atomicAdd(d_leafCounter, 1);
        if (wPos < maxLeaves) {
            d_leafPairs[wPos] = curr;
        }
        return;
    }

    int writeIdx = d_offsets[tid];
    
    auto push_frontier = [&](uint32_t a, uint32_t b) {
        if (writeIdx < frontierCapacity) {
            d_outFrontier[writeIdx++] = {a, b, curr.batchIdx};
        }
    };

    if (!isLeafA && !isLeafB) {
        uint32_t al = nodeA.admin.offset; uint32_t ar = al + 1;
        uint32_t bl = nodeB.admin.offset; uint32_t br = bl + 1;
        if (bvhA.nodes[al].bounds.overlaps(bvhB.nodes[bl].bounds)) push_frontier(al, bl);
        if (bvhA.nodes[al].bounds.overlaps(bvhB.nodes[br].bounds)) push_frontier(al, br);
        if (bvhA.nodes[ar].bounds.overlaps(bvhB.nodes[bl].bounds)) push_frontier(ar, bl);
        if (bvhA.nodes[ar].bounds.overlaps(bvhB.nodes[br].bounds)) push_frontier(ar, br);
    } 
    else if (isLeafA) { 
        uint32_t bl = nodeB.admin.offset; uint32_t br = bl + 1;
        if (nodeA.bounds.overlaps(bvhB.nodes[bl].bounds)) push_frontier(curr.nodeID_A, bl);
        if (nodeA.bounds.overlaps(bvhB.nodes[br].bounds)) push_frontier(curr.nodeID_A, br);
    } 
    else { 
        uint32_t al = nodeA.admin.offset; uint32_t ar = al + 1;
        if (bvhA.nodes[al].bounds.overlaps(nodeB.bounds)) push_frontier(al, curr.nodeID_B);
        if (bvhA.nodes[ar].bounds.overlaps(nodeB.bounds)) push_frontier(ar, curr.nodeID_B);
    }
}

__global__ void evaluateLeafPrimitivePairsKernel(
    const DualFrontierElement* d_leafPairs, int leafCount, 
    int2* d_candidatePairs, int* d_candCounter, uint32_t maxCandidates,
    cuBQL::bvh3f bvhA, cuBQL::bvh3f bvhB,
    const cuBQL::Triangle* dMeshA, const cuBQL::Triangle* dMeshB
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
// PIPELINE TRAVERSAL EXECUTION ENGINE
// --------------------------------------------------------------------
uint64_t executeDualTreeTraversal(
    int batchMultiplier, int totalBatches,
    uint32_t maxDescendantsA, uint32_t maxDescendantsB, float expectedIntersectionDensity,
    const thrust::device_vector<uint32_t>& d_outPairsA, const thrust::device_vector<uint32_t>& d_outPairsB,
    const thrust::device_vector<uint32_t>& d_markedNodeIndicesA, const thrust::device_vector<uint32_t>& d_markedNodeIndicesB,
    const thrust::device_vector<uint32_t>& d_nodeDescendantCountsA, const thrust::device_vector<uint32_t>& d_nodeDescendantCountsB, 
    uint32_t h_outMarkedCountA, uint32_t h_outMarkedCountB,
    const cuBQL::bvh3f& bvhA, const cuBQL::bvh3f& bvhB,
    const cuBQL::Triangle* dMeshA, const cuBQL::Triangle* dMeshB,
    std::vector<int2>& hGreenPairs, std::vector<int2>& hYellowPairs, 
    IntersectionTimeTracker& outLoopTime
) {
    double tStart = cuBQL::getCurrentTime();
    uint64_t totalCandidatePairsProcessed = 0;

    const uint32_t* base_outPairsA          = thrust::raw_pointer_cast(d_outPairsA.data());
    const uint32_t* base_outPairsB          = thrust::raw_pointer_cast(d_outPairsB.data());
    const uint32_t* base_markedNodeIndicesA = thrust::raw_pointer_cast(d_markedNodeIndicesA.data());
    const uint32_t* base_markedNodeIndicesB = thrust::raw_pointer_cast(d_markedNodeIndicesB.data());

    // Compute exact structured upper-bound sizing limits
    uint64_t safeMaxDescA = std::max<uint64_t>(1, maxDescendantsA);
    uint64_t safeMaxDescB = std::max<uint64_t>(1, maxDescendantsB);
    uint64_t maxLeafPairsPossible = (uint64_t)batchMultiplier * safeMaxDescA * safeMaxDescB;
    
    uint32_t leafCapacity = (uint32_t)std::max<uint64_t>(4000000, maxLeafPairsPossible);
    uint32_t candCapacity = (uint32_t)std::max<uint64_t>(16000000, 
        (uint64_t)(maxLeafPairsPossible * expectedIntersectionDensity));

    DualFrontierElement* d_leafPairs = nullptr;
    int* d_leafCounter               = nullptr;
    int2* d_candidatePairs           = nullptr;
    int* d_pairStatuses             = nullptr;

    CUBQL_CUDA_CALL(Malloc(&d_leafPairs,      leafCapacity * sizeof(DualFrontierElement)));
    CUBQL_CUDA_CALL(Malloc(&d_leafCounter,    sizeof(int)));
    CUBQL_CUDA_CALL(Malloc(&d_candidatePairs, candCapacity * sizeof(int2)));
    CUBQL_CUDA_CALL(Malloc(&d_pairStatuses,   candCapacity * sizeof(int)));

    outLoopTime.preallocateTimeMs = (cuBQL::getCurrentTime() - tStart) * 1000.0;

    for (int i = 0; i < totalBatches; i += batchMultiplier) {
        int activeBatchesInChunk = std::min(batchMultiplier, totalBatches - i);
        if (activeBatchesInChunk <= 0) continue;

        double tAssemblyStart = cuBQL::getCurrentTime();
        int zeroReset = 0;
        CUBQL_CUDA_CALL(Memcpy(d_leafCounter, &zeroReset, sizeof(int), cudaMemcpyHostToDevice));

        uint32_t currentFrontierSize = activeBatchesInChunk;
        DualFrontierElement* d_frontierCurrent = nullptr;
        CUBQL_CUDA_CALL(Malloc(&d_frontierCurrent, currentFrontierSize * sizeof(DualFrontierElement)));

        int seedBlock = 256;
        int seedGrid  = (activeBatchesInChunk + seedBlock - 1) / seedBlock;
        
        initializeFrontierSeedKernel<<<seedGrid, seedBlock>>>(
            d_frontierCurrent, 
            base_outPairsA + i, 
            base_outPairsB + i, 
            base_markedNodeIndicesA, 
            base_markedNodeIndicesB, 
            i, activeBatchesInChunk
        );
        CHECK_CUDA_KERNEL_STATE("initializeFrontierSeedKernel");
        
        outLoopTime.assemblyPhaseMs += (cuBQL::getCurrentTime() - tAssemblyStart) * 1000.0;
        double tExecutionStart = cuBQL::getCurrentTime();
        
        while (currentFrontierSize > 0) {
            int block = 256;
            int grid  = (currentFrontierSize + block - 1) / block;

            int* d_childCounts = nullptr;
            int* d_offsets     = nullptr;
            CUBQL_CUDA_CALL(Malloc(&d_childCounts, currentFrontierSize * sizeof(int)));
            CUBQL_CUDA_CALL(Malloc(&d_offsets,     currentFrontierSize * sizeof(int)));

            countDualTreeExpansionsKernel<<<grid, block>>>(d_frontierCurrent, currentFrontierSize, d_childCounts, bvhA, bvhB);
            CHECK_CUDA_KERNEL_STATE("countDualTreeExpansionsKernel");

            thrust::device_ptr<int> dev_counts(d_childCounts);
            thrust::device_ptr<int> dev_offsets(d_offsets);
            thrust::exclusive_scan(thrust::device, dev_counts, dev_counts + currentFrontierSize, dev_offsets);

            int lastCount = 0, lastOffset = 0;
            CUBQL_CUDA_CALL(Memcpy(&lastCount, d_childCounts + currentFrontierSize - 1, sizeof(int), cudaMemcpyDeviceToHost));
            CUBQL_CUDA_CALL(Memcpy(&lastOffset, d_offsets + currentFrontierSize - 1, sizeof(int), cudaMemcpyDeviceToHost));
            uint32_t nextFrontierSize = lastCount + lastOffset;

            if (nextFrontierSize == 0) {
                CUBQL_CUDA_CALL(Free(d_childCounts));
                CUBQL_CUDA_CALL(Free(d_offsets));
                CUBQL_CUDA_CALL(Free(d_frontierCurrent));
                break;
            }

            DualFrontierElement* d_frontierNext = nullptr;
            CUBQL_CUDA_CALL(Malloc(&d_frontierNext, nextFrontierSize * sizeof(DualFrontierElement)));

            populateDualTreeFrontierKernel<<<grid, block>>>(
                d_frontierCurrent, currentFrontierSize, d_offsets, 
                d_frontierNext, d_leafPairs, d_leafCounter, 
                leafCapacity, nextFrontierSize, bvhA, bvhB
            );
            CHECK_CUDA_KERNEL_STATE("populateDualTreeFrontierKernel");

            CUBQL_CUDA_CALL(Free(d_childCounts));
            CUBQL_CUDA_CALL(Free(d_offsets));
            CUBQL_CUDA_CALL(Free(d_frontierCurrent));

            d_frontierCurrent   = d_frontierNext;
            currentFrontierSize = nextFrontierSize;
        }

        outLoopTime.executionPhaseMs += (cuBQL::getCurrentTime() - tExecutionStart) * 1000.0;

        int finalLeafPairsCount = 0;
        CUBQL_CUDA_CALL(Memcpy(&finalLeafPairsCount, d_leafCounter, sizeof(int), cudaMemcpyDeviceToHost));

        if (finalLeafPairsCount > 0) {
            double tFineStart = cuBQL::getCurrentTime();
            int leafBlock = 256;
            int leafGrid  = (finalLeafPairsCount + leafBlock - 1) / leafBlock;

            CUBQL_CUDA_CALL(Memcpy(d_leafCounter, &zeroReset, sizeof(int), cudaMemcpyHostToDevice));

            evaluateLeafPrimitivePairsKernel<<<leafGrid, leafBlock>>>(
                d_leafPairs, finalLeafPairsCount, d_candidatePairs, 
                d_leafCounter, candCapacity, bvhA, bvhB, dMeshA, dMeshB
            );
            CHECK_CUDA_KERNEL_STATE("evaluateLeafPrimitivePairsKernel");

            int totalChunkPairs = 0;
            CUBQL_CUDA_CALL(Memcpy(&totalChunkPairs, d_leafCounter, sizeof(int), cudaMemcpyDeviceToHost));

            if (totalChunkPairs > 0) {
                totalCandidatePairsProcessed += totalChunkPairs;

                double internalPredicateTimeDummy = 0.0;
                evaluateAndCompactPairs(
                    d_candidatePairs, d_pairStatuses, 
                    dMeshA, dMeshB, totalChunkPairs, internalPredicateTimeDummy
                );
                CHECK_CUDA_KERNEL_STATE("evaluateAndCompactPairs");
                
                double tDownloadStart = cuBQL::getCurrentTime();
                thrust::device_ptr<int2> dev_evaluated(d_candidatePairs);
                thrust::device_ptr<int> dev_statuses(d_pairStatuses);

                thrust::device_vector<int2> dev_green(totalChunkPairs);
                thrust::device_vector<int2> dev_yellow(totalChunkPairs);

                thrust::device_ptr<int2> dev_green_out(thrust::raw_pointer_cast(dev_green.data()));
                thrust::device_ptr<int2> dev_yellow_out(thrust::raw_pointer_cast(dev_yellow.data()));

                auto green_end  = thrust::copy_if(thrust::device, dev_evaluated, dev_evaluated + totalChunkPairs, dev_statuses, dev_green_out, IsTargetPairStatusTwo{(int)PAIR_GREEN});
                auto yellow_end = thrust::copy_if(thrust::device, dev_evaluated, dev_evaluated + totalChunkPairs, dev_statuses, dev_yellow_out, IsTargetPairStatusTwo{(int)PAIR_YELLOW});
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
                
                outLoopTime.DownloadAndClean += (cuBQL::getCurrentTime() - tDownloadStart) * 1000.0;
            }
            outLoopTime.fineEvaluationPhaseMs += (cuBQL::getCurrentTime() - tFineStart) * 1000.0;
        }
    }

    CUBQL_CUDA_CALL(Free(d_leafPairs));
    CUBQL_CUDA_CALL(Free(d_leafCounter));
    CUBQL_CUDA_CALL(Free(d_candidatePairs));
    CUBQL_CUDA_CALL(Free(d_pairStatuses));

    outLoopTime.totalLoopTimeMs = (cuBQL::getCurrentTime() - tStart) * 1000.0;
    
    return totalCandidatePairsProcessed; 
}