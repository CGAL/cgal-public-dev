#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/execution_policy.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include "cuBQL/bvh.h"

struct DualFrontierPair {
    uint32_t nodeIdxA;
    uint32_t nodeIdxB;
};

// --------------------------------------------------------------------
// COMPILER-COMPATIBLE TH_FUNCTOR
// --------------------------------------------------------------------
struct ZipToPairFunctor {
    __device__ DualFrontierPair operator()(const thrust::tuple<uint32_t, uint32_t>& t) const {
        return DualFrontierPair{thrust::get<0>(t), thrust::get<1>(t)};
    }
};

// --------------------------------------------------------------------
// KERNELS
// --------------------------------------------------------------------

__global__ void countExpansionsKernel(
    const DualFrontierPair* d_inFrontier, 
    uint32_t frontierSize, 
    int* d_childCounts,
    cuBQL::bvh3f bvhA, 
    cuBQL::bvh3f bvhB) 
{
    uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= frontierSize) return;

    DualFrontierPair pair = d_inFrontier[tid];
    auto nodeA = bvhA.nodes[pair.nodeIdxA];
    auto nodeB = bvhB.nodes[pair.nodeIdxB];

    bool isLeafA = (nodeA.admin.count > 0);
    bool isLeafB = (nodeB.admin.count > 0);

    if (isLeafA && isLeafB) {
        d_childCounts[tid] = 1;
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

__global__ void populateFrontierKernel(
    const DualFrontierPair* d_inFrontier, 
    uint32_t frontierSize, 
    const int* d_offsets, 
    DualFrontierPair* d_outFrontier,
    cuBQL::bvh3f bvhA, 
    cuBQL::bvh3f bvhB) 
{
    uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= frontierSize) return;

    DualFrontierPair pair = d_inFrontier[tid];
    auto nodeA = bvhA.nodes[pair.nodeIdxA];
    auto nodeB = bvhB.nodes[pair.nodeIdxB];

    bool isLeafA = (nodeA.admin.count > 0);
    bool isLeafB = (nodeB.admin.count > 0);

    int writeIdx = d_offsets[tid];

    if (isLeafA && isLeafB) {
        d_outFrontier[writeIdx] = pair;
        return;
    }

    auto push_pair = [&](uint32_t a, uint32_t b) {
        d_outFrontier[writeIdx++] = {a, b};
    };

    if (!isLeafA && !isLeafB) {
        uint32_t al = nodeA.admin.offset; uint32_t ar = al + 1;
        uint32_t bl = nodeB.admin.offset; uint32_t br = bl + 1;

        if (bvhA.nodes[al].bounds.overlaps(bvhB.nodes[bl].bounds)) push_pair(al, bl);
        if (bvhA.nodes[al].bounds.overlaps(bvhB.nodes[br].bounds)) push_pair(al, br);
        if (bvhA.nodes[ar].bounds.overlaps(bvhB.nodes[bl].bounds)) push_pair(ar, bl);
        if (bvhA.nodes[ar].bounds.overlaps(bvhB.nodes[br].bounds)) push_pair(ar, br);
    } 
    else if (isLeafA) { 
        uint32_t bl = nodeB.admin.offset; uint32_t br = bl + 1;

        if (nodeA.bounds.overlaps(bvhB.nodes[bl].bounds)) push_pair(pair.nodeIdxA, bl);
        if (nodeA.bounds.overlaps(bvhB.nodes[br].bounds)) push_pair(pair.nodeIdxA, br);
    } 
    else { 
        uint32_t al = nodeA.admin.offset; uint32_t ar = al + 1;

        if (bvhA.nodes[al].bounds.overlaps(nodeB.bounds)) push_pair(al, pair.nodeIdxB);
        if (bvhA.nodes[ar].bounds.overlaps(nodeB.bounds)) push_pair(ar, pair.nodeIdxB);
    }
}

__global__ void finalUnpackKernel(
    const DualFrontierPair* d_frontier, 
    uint32_t totalPairs,
    uint32_t* d_outPairsA, 
    uint32_t* d_outPairsB) 
{
    uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= totalPairs) return;

    d_outPairsA[tid] = d_frontier[tid].nodeIdxA;
    d_outPairsB[tid] = d_frontier[tid].nodeIdxB;
}

// --------------------------------------------------------------------
// PIPELINE ORCHESTRATOR
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// PIPELINE ORCHESTRATOR (Optimized Memory Management)
// --------------------------------------------------------------------
uint64_t executeDualTreeStep(
    int numSteps,
    uint32_t maxDescendantsA,      
    uint32_t maxDescendantsB,      
    thrust::device_vector<uint32_t>& d_outPairsA,
    thrust::device_vector<uint32_t>& d_outPairsB,
    thrust::device_vector<uint32_t>& d_markedNodeIndicesA,      
    thrust::device_vector<uint32_t>& d_markedNodeIndicesB,      
    thrust::device_vector<uint32_t>& d_nodeDescendantCountsA,   
    thrust::device_vector<uint32_t>& d_nodeDescendantCountsB,   
    uint32_t& h_outMarkedCountA,                                
    uint32_t& h_outMarkedCountB,                                
    const cuBQL::bvh3f& bvhA,
    const cuBQL::bvh3f& bvhB,
    const void* dMeshA, 
    const void* dMeshB
)  {
    if (d_outPairsA.empty()) return 0;

    uint32_t currentFrontierSize = d_outPairsA.size();

    // Start with a reasonable initial capacity instead of maximum possible
    uint32_t currentCapacity = std::max((uint32_t)1000000, currentFrontierSize * 2);
    
    thrust::device_vector<DualFrontierPair> frontierBuffer1(currentCapacity);
    thrust::device_vector<DualFrontierPair> frontierBuffer2(currentCapacity);
    thrust::device_vector<int> d_childCounts(currentCapacity);
    thrust::device_vector<int> d_offsets(currentCapacity);

    // Track which vector is currently active
    thrust::device_vector<DualFrontierPair>* currentBuffer = &frontierBuffer1;
    thrust::device_vector<DualFrontierPair>* nextBuffer    = &frontierBuffer2;

    auto inputZipBegin = thrust::make_zip_iterator(thrust::make_tuple(d_outPairsA.begin(), d_outPairsB.begin()));
    auto inputZipEnd   = inputZipBegin + currentFrontierSize;
    
    thrust::transform(thrust::device, inputZipBegin, inputZipEnd, currentBuffer->begin(), ZipToPairFunctor());

    for (int step = 0; step < numSteps; ++step) {
        if (currentFrontierSize == 0) break;

        uint32_t threads = 256;
        uint32_t blocks  = (currentFrontierSize + threads - 1) / threads;

        // ALWAYS grab a fresh raw pointer before kernel launch
        DualFrontierPair* d_frontierCurrent = thrust::raw_pointer_cast(currentBuffer->data());

        countExpansionsKernel<<<blocks, threads>>>(
            d_frontierCurrent, currentFrontierSize, 
            thrust::raw_pointer_cast(d_childCounts.data()), 
            bvhA, bvhB
        );
        cudaDeviceSynchronize();

        thrust::exclusive_scan(
            thrust::device, 
            d_childCounts.begin(), 
            d_childCounts.begin() + currentFrontierSize, 
            d_offsets.begin()
        );

        int lastCount = d_childCounts[currentFrontierSize - 1];
        int lastOffset = d_offsets[currentFrontierSize - 1];
        uint32_t nextFrontierSize = lastCount + lastOffset;

        if (nextFrontierSize == 0) {
            currentFrontierSize = 0;
            break;
        }

        // --- DYNAMIC RESIZING LOGIC ---
        if (nextFrontierSize > currentCapacity) {
            // Grow by 1.5x to amortize allocation costs
            currentCapacity = std::max((uint32_t)(currentCapacity * 1.5), nextFrontierSize);
            
            currentBuffer->resize(currentCapacity);
            nextBuffer->resize(currentCapacity);
            d_childCounts.resize(currentCapacity);
            d_offsets.resize(currentCapacity);

            // Re-acquire pointer because resize() invalidated the old one
            d_frontierCurrent = thrust::raw_pointer_cast(currentBuffer->data());
        }

        DualFrontierPair* d_frontierNext = thrust::raw_pointer_cast(nextBuffer->data());

        populateFrontierKernel<<<blocks, threads>>>(
            d_frontierCurrent, currentFrontierSize, 
            thrust::raw_pointer_cast(d_offsets.data()), 
            d_frontierNext, 
            bvhA, bvhB
        );
        cudaDeviceSynchronize();

        // Swap the vector references for the next iteration
        std::swap(currentBuffer, nextBuffer);
        currentFrontierSize = nextFrontierSize;
    }

    if (currentFrontierSize > 0) {
        d_outPairsA.resize(currentFrontierSize);
        d_outPairsB.resize(currentFrontierSize);

        uint32_t threads = 256;
        uint32_t blocks  = (currentFrontierSize + threads - 1) / threads;

        DualFrontierPair* d_finalFrontier = thrust::raw_pointer_cast(currentBuffer->data());

        finalUnpackKernel<<<blocks, threads>>>(
            d_finalFrontier, currentFrontierSize,
            thrust::raw_pointer_cast(d_outPairsA.data()),
            thrust::raw_pointer_cast(d_outPairsB.data())
        );
        cudaDeviceSynchronize();

        thrust::device_vector<uint32_t> d_tempMarkedA = d_outPairsA;
        thrust::device_vector<uint32_t> d_tempMarkedB = d_outPairsB;

        thrust::sort(thrust::device, d_tempMarkedA.begin(), d_tempMarkedA.end());
        thrust::sort(thrust::device, d_tempMarkedB.begin(), d_tempMarkedB.end());

        auto endA = thrust::unique(thrust::device, d_tempMarkedA.begin(), d_tempMarkedA.end());
        auto endB = thrust::unique(thrust::device, d_tempMarkedB.begin(), d_tempMarkedB.end());

        h_outMarkedCountA = endA - d_tempMarkedA.begin();
        h_outMarkedCountB = endB - d_tempMarkedB.begin();

        d_markedNodeIndicesA.resize(h_outMarkedCountA);
        d_markedNodeIndicesB.resize(h_outMarkedCountB);
        thrust::copy_n(thrust::device, d_tempMarkedA.begin(), h_outMarkedCountA, d_markedNodeIndicesA.begin());
        thrust::copy_n(thrust::device, d_tempMarkedB.begin(), h_outMarkedCountB, d_markedNodeIndicesB.begin());
        cudaDeviceSynchronize();
    }

    return currentFrontierSize;
}