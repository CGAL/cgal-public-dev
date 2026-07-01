#pragma once

#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/execution_policy.h>
#include "include/third-party/cubql/sm_builder_v3.h"

// 1. Reverted to standard, high-performance float box intersection helper
template<typename T, int D>
__device__ inline bool normalBoxesIntersect(
    const cuBQL::box_t<T,D>& a, 
    const cuBQL::box_t<T,D>& b) 
{
    for (int i = 0; i < D; ++i) {
        if (a.get_lower(i) > b.get_upper(i) || a.get_upper(i) < b.get_lower(i)) {
            return false;
        }
    }
    return true;
}

// 2. Pass 1: Count Intersections using FinalNormalNode
template<typename T, int D>
__global__ void countNormalNodeIntersections(
    const cuBQL::gpuBuilder_v3::FinalNormalNode<T,D>* d_nodesA, uint32_t startA, uint32_t endA,
    const cuBQL::gpuBuilder_v3::FinalNormalNode<T,D>* d_nodesB, uint32_t startB, uint32_t endB,
    uint32_t* d_perThreadCounts) 
{
    uint32_t threadIdxGlobal = threadIdx.x + blockIdx.x * blockDim.x;
    uint32_t totalActiveA = endA - startA;
    if (threadIdxGlobal >= totalActiveA) return;

    uint32_t nodeIDA = startA + threadIdxGlobal;
    auto& openBranchA = d_nodesA[nodeIDA].openBranch;
    
    if (openBranchA.count == 0) {
        d_perThreadCounts[threadIdxGlobal] = 0;
        return;
    }

    // Directly access the clean, native float box layout
    const auto& boxA = openBranchA.bounds; 
    uint32_t localIntersections = 0;

    for (uint32_t nodeIDB = startB; nodeIDB < endB; ++nodeIDB) {
        auto& openBranchB = d_nodesB[nodeIDB].openBranch;
        if (openBranchB.count == 0) continue;

        const auto& boxB = openBranchB.bounds;
        if (normalBoxesIntersect<T,D>(boxA, boxB)) { 
            localIntersections++; 
        }
    }
    d_perThreadCounts[threadIdxGlobal] = localIntersections;
}

// 3. Pass 2: Fill Intersections directly into pair arrays using FinalNormalNode
template<typename T, int D>
__global__ void fillNormalNodeIntersections(
    const cuBQL::gpuBuilder_v3::FinalNormalNode<T,D>* d_nodesA, uint32_t startA, uint32_t endA,
    const cuBQL::gpuBuilder_v3::FinalNormalNode<T,D>* d_nodesB, uint32_t startB, uint32_t endB,
    const uint32_t* d_offsets,
    uint32_t* d_outPairsA, uint32_t* d_outPairsB) 
{
    uint32_t threadIdxGlobal = threadIdx.x + blockIdx.x * blockDim.x;
    uint32_t totalActiveA = endA - startA;
    if (threadIdxGlobal >= totalActiveA) return;

    uint32_t nodeIDA = startA + threadIdxGlobal;
    auto& openBranchA = d_nodesA[nodeIDA].openBranch;
    if (openBranchA.count == 0) return;

    const auto& boxA = openBranchA.bounds;
    uint32_t writePos = d_offsets[threadIdxGlobal];

    for (uint32_t nodeIDB = startB; nodeIDB < endB; ++nodeIDB) {
        auto& openBranchB = d_nodesB[nodeIDB].openBranch;
        if (openBranchB.count == 0) continue;

        const auto& boxB = openBranchB.bounds;
        if (normalBoxesIntersect<T,D>(boxA, boxB)) { 
            d_outPairsA[writePos] = nodeIDA;
            d_outPairsB[writePos] = nodeIDB;
            writePos++;
        }
    }
}

// 4. Upgraded Orchestrator function mapping clean pointer streams
template<typename T, int D>
uint32_t executeNormalNodeCrossCheck(
    const cuBQL::gpuBuilder_v3::FinalNormalNode<T,D>* d_nodesA, uint32_t startA, uint32_t endA,
    const cuBQL::gpuBuilder_v3::FinalNormalNode<T,D>* d_nodesB, uint32_t startB, uint32_t endB,
    thrust::device_vector<uint32_t>& d_outPairsA,
    thrust::device_vector<uint32_t>& d_outPairsB,
    cudaStream_t s)
{
    uint32_t totalActiveA = endA - startA;
    uint32_t totalActiveB = endB - startB;
    if (totalActiveA == 0 || totalActiveB == 0) return 0;

    uint32_t threadsCross = 256;
    uint32_t blocksCross = (totalActiveA + threadsCross - 1) / threadsCross;

    auto exec_policy = thrust::cuda::par.on(s);

    // Pass 1: Count
    thrust::device_vector<uint32_t> d_perThreadCounts(totalActiveA, 0);
    countNormalNodeIntersections<T,D><<<blocksCross, threadsCross, 0, s>>>(
        d_nodesA, startA, endA,
        d_nodesB, startB, endB,
        thrust::raw_pointer_cast(d_perThreadCounts.data())
    );

    // Prefix Scan on stream
    thrust::device_vector<uint32_t> d_offsets(totalActiveA, 0);
    thrust::exclusive_scan(exec_policy, d_perThreadCounts.begin(), d_perThreadCounts.end(), d_offsets.begin());
    
    uint32_t lastCount = d_perThreadCounts.back();
    uint32_t lastOffset = d_offsets.back();
    uint32_t totalIntersections = lastOffset + lastCount;

    if (totalIntersections > 0) {
        d_outPairsA.resize(totalIntersections);
        d_outPairsB.resize(totalIntersections);

        // Pass 2: Fill
        fillNormalNodeIntersections<T,D><<<blocksCross, threadsCross, 0, s>>>(
            d_nodesA, startA, endA,
            d_nodesB, startB, endB,
            thrust::raw_pointer_cast(d_offsets.data()),
            thrust::raw_pointer_cast(d_outPairsA.data()),
            thrust::raw_pointer_cast(d_outPairsB.data())
        );
    }

    return totalIntersections;
}