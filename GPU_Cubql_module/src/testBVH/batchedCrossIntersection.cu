#include "batchedCrossIntersection.h"
#include "include/third-party/cubql/fixedBoxQueryv2.h"
#include <thrust/device_ptr.h>
#include <thrust/scan.h>
#include <thrust/execution_policy.h>
#include <iostream>

// Helper macro matching your project's naming convention
#ifndef CUBQL_CUDA_CALL
#define CUBQL_CUDA_CALL( call )                                                              \
{                                                                                            \
    cudaError_t err = cuda##call;                                                             \
    if( cudaSuccess != err ) {                                                               \
        fprintf(stderr, "CUDA error in file '%s' in line %i : %s.\n",                        \
                __FILE__, __LINE__, cudaGetErrorString( err ) );                             \
                exit(EXIT_FAILURE);                                                                  \
    }                                                                                        \
}
#endif

// Functor to calculate size of batch bIdx completely on the GPU using offsets
struct GetBatchSizeFunctor {
    const uint32_t* outOffsetsB;
    uint32_t outMarkedCountB;
    uint32_t totalPrimsB;

    __host__ __device__
    GetBatchSizeFunctor(const uint32_t* _offsets, uint32_t _markedCount, uint32_t _totalPrims) 
        : outOffsetsB(_offsets), outMarkedCountB(_markedCount), totalPrimsB(_totalPrims) {}

    __host__ __device__
    int operator()(const uint32_t bIdx) const {
        uint32_t startOffset = outOffsetsB[bIdx];
        uint32_t endOffset = (bIdx + 1 < outMarkedCountB) ? outOffsetsB[bIdx + 1] : totalPrimsB;
        return (int)(endOffset - startOffset);
    }
};

__global__ void buildBobKernel(
    uint64_t* d_bob,
    const uint32_t* d_outPairsA,
    const uint32_t* d_outPairsB,
    const uint32_t* d_markedNodeIndicesA,
    const uint32_t* d_outOffsetsB,
    uint32_t totalBatches,
    uint32_t totalPrimsFlatBSize)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= totalBatches) return;

    // 1. Get the value to fill
    uint32_t pairAIdx = d_outPairsA[i];
    uint64_t valueToFill = (uint64_t)d_markedNodeIndicesA[pairAIdx];

    // 2. Get the range [start, end)
    uint32_t start = d_outOffsetsB[i];
    uint32_t end = (i + 1 < totalBatches) ? d_outOffsetsB[i + 1] : totalPrimsFlatBSize;

    // 3. Perform the fill for this batch's range
    for (uint32_t idx = start; idx < end; ++idx) {
        d_bob[idx] = valueToFill;
    }
}


// thrust::device_vector<uint64_t> d_bob(d_outPrimsFlatB.size());

// // 2. Launch once
// int blockSize = 256;
// int gridSize = (totalBatches + blockSize - 1) / blockSize;

// buildBobKernel<<<gridSize, blockSize>>>(
//     d_bob.data().get(),
//     d_outPairsA.data().get(),
//     d_outPairsB.data().get(),
//     d_markedNodeIndicesA.data().get(),
//     d_outOffsetsB.data().get(),
//     totalBatches,
//     d_outPrimsFlatB.size()
// );

// // 3. Sync and use
// cudaDeviceSynchronize();

// --------------------------------------------------------------------
// KERNEL DEFINITIONS
// --------------------------------------------------------------------
__global__ void countAABBOverlapsKernel_Indirected(
    int *pairCounts, 
    cuBQL::bvh3f bvhA, 
    const cuBQL::Triangle *triA, 
    const cuBQL::Triangle *triB, 
    const uint32_t* d_outPrimsFlatB, 
    uint32_t startOffsetB, 
    int numPrimsB, 
    const uint64_t* startNodeIdxAs) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= numPrimsB) return;

    uint32_t actualPrimIdB = d_outPrimsFlatB[startOffsetB + tid];
    uint64_t startNodeIdxA = startNodeIdxAs[startOffsetB + tid]; 
    cuBQL::Triangle b = triB[actualPrimIdB];
    cuBQL::box3f query = b.bounds();
    
    int count = 0;
    cuBQL::fixedBoxQueryv2::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        for (uint32_t i = 0; i < num; i++) {
            if (triA[ids[i]].bounds().overlaps(query)) { count++; }
        }
        return CUBQL_CONTINUE_TRAVERSAL;
    }, bvhA, query, startNodeIdxA);
    
    pairCounts[tid] = count;
}

__global__ void fillAABBOverlapsKernel_Indirected(
    int2 *candidatePairs, 
    const int *offsets, 
    cuBQL::bvh3f bvhA, 
    const cuBQL::Triangle *triA, 
    const cuBQL::Triangle *triB, 
    const uint32_t* d_outPrimsFlatB, 
    uint32_t startOffsetB, 
    int numPrimsB, 
    const uint64_t* startNodeIdxAs) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= numPrimsB) return;

    int wPos = offsets[tid];
    uint32_t actualPrimIdB = d_outPrimsFlatB[startOffsetB + tid];
    uint64_t startNodeIdxA = startNodeIdxAs[startOffsetB + tid]; 
    cuBQL::Triangle b = triB[actualPrimIdB];
    cuBQL::box3f query = b.bounds();
    
    cuBQL::fixedBoxQueryv2::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        for (uint32_t i = 0; i < num; i++) {
            if (triA[ids[i]].bounds().overlaps(query)) { 
                candidatePairs[wPos++] = make_int2((int)ids[i], (int)actualPrimIdB); 
            }
        }
        return CUBQL_CONTINUE_TRAVERSAL;
    }, bvhA, query, startNodeIdxA);
}

// --------------------------------------------------------------------
// HOST EXECUTABLE REGION
// --------------------------------------------------------------------
uint64_t executeBatchedCrossIntersectionLoop(
    int batchMultiplier,
    int totalBatches,
    const thrust::device_vector<uint32_t>& d_outPairsA,
    const thrust::device_vector<uint32_t>& d_outPairsB,
    const thrust::device_vector<uint32_t>& d_markedNodeIndicesA,
    const thrust::device_vector<uint32_t>& d_markedNodeIndicesB,
    const thrust::device_vector<uint32_t>& d_outOffsetsB,
    const thrust::device_vector<uint32_t>& d_outPrimsFlatB,
    const thrust::device_vector<uint32_t>& d_nodeDescendantCountsB, 
    uint32_t h_outMarkedCountB,
    const cuBQL::bvh3f& bvhA,
    const cuBQL::Triangle* dMeshA,
    const cuBQL::Triangle* dMeshB,
    double& outLoopTime
) {
    uint64_t finalCandidatePairs = 0;

    // 1. Get raw device pointers
    const uint32_t* ptr_outPairsA           = thrust::raw_pointer_cast(d_outPairsA.data());
    const uint32_t* ptr_outPairsB           = thrust::raw_pointer_cast(d_outPairsB.data());
    const uint32_t* ptr_markedNodeIndicesA  = thrust::raw_pointer_cast(d_markedNodeIndicesA.data());
    const uint32_t* ptr_outOffsetsB         = thrust::raw_pointer_cast(d_outOffsetsB.data());
    const uint32_t* ptr_outPrimsFlatB       = thrust::raw_pointer_cast(d_outPrimsFlatB.data());

    // Allocate the unified long flat root array matching d_outPrimsFlatB length 1-to-1
    thrust::device_vector<uint64_t> d_flatNodeRootsA(d_outPrimsFlatB.size());
    uint64_t* ptr_flatNodeRootsA = thrust::raw_pointer_cast(d_flatNodeRootsA.data());

    // 2. Find maximum primitive count using outOffsetsB directly
    int maxPrimsInBatch = 0;
    if (totalBatches > 0) {
        auto batchSizeTransformIterator = thrust::make_transform_iterator(
            d_outPairsB.begin(),
            GetBatchSizeFunctor(ptr_outOffsetsB, h_outMarkedCountB, d_outPrimsFlatB.size())
        );

        maxPrimsInBatch = thrust::reduce(
            thrust::device,
            batchSizeTransformIterator,
            batchSizeTransformIterator + totalBatches,
            0,
            thrust::maximum<int>()
        );
    }

    // 3. Preallocate our single reusable pair of intermediate count buffers
    int* d_pairCounts = nullptr;
    int* d_offsets = nullptr;
    if (maxPrimsInBatch > 0) {
        CUBQL_CUDA_CALL(Malloc(&d_pairCounts, maxPrimsInBatch * sizeof(int)));
        CUBQL_CUDA_CALL(Malloc(&d_offsets, maxPrimsInBatch * sizeof(int)));
    }

    double tBatchLoopStart = cuBQL::getCurrentTime();

    // 4. Pure Sequential Loop to Protect Memory 
    for (int i = 0; i < totalBatches; ++i) {
        uint32_t pairAVal, bIdx;
        CUBQL_CUDA_CALL(Memcpy(&pairAVal, ptr_outPairsA + i, sizeof(uint32_t), cudaMemcpyDeviceToHost));
        CUBQL_CUDA_CALL(Memcpy(&bIdx, ptr_outPairsB + i, sizeof(uint32_t), cudaMemcpyDeviceToHost));
        
        uint32_t startOffset;
        CUBQL_CUDA_CALL(Memcpy(&startOffset, ptr_outOffsetsB + bIdx, sizeof(uint32_t), cudaMemcpyDeviceToHost));
        
        uint32_t endOffset;
        if (bIdx + 1 < h_outMarkedCountB) {
            CUBQL_CUDA_CALL(Memcpy(&endOffset, ptr_outOffsetsB + bIdx + 1, sizeof(uint32_t), cudaMemcpyDeviceToHost));
        } else {
            endOffset = d_outPrimsFlatB.size();
        }
        
        int numPrims = endOffset - startOffset;
        if (numPrims == 0) continue;

        // Fetch the source root index value exactly like the original code
        uint32_t startNodeIdxA_32;
        CUBQL_CUDA_CALL(Memcpy(&startNodeIdxA_32, ptr_markedNodeIndicesA + pairAVal, sizeof(uint32_t), cudaMemcpyDeviceToHost));
        uint64_t startNodeIdxA = (uint64_t)startNodeIdxA_32;

        // FIXED: Dynamically write the scalar root index across this exact primitive window range 
        // inside the pre-allocated flat tracking array. This perfectly preserves the precise logic 
        // of your original version.
        thrust::fill_n(thrust::device, d_flatNodeRootsA.begin() + startOffset, numPrims, startNodeIdxA);

        // --- PHASE 1: Execution Count ---
        int blockSize = 128;
        int gridScale = (numPrims + blockSize - 1) / blockSize;

        countAABBOverlapsKernel_Indirected<<<gridScale, blockSize>>>(
            d_pairCounts, bvhA, const_cast<cuBQL::Triangle*>(dMeshA), const_cast<cuBQL::Triangle*>(dMeshB), 
            ptr_outPrimsFlatB, startOffset, numPrims, ptr_flatNodeRootsA
        );

        // --- PHASE 2: Reusable Device Inclusive/Exclusive Scan ---
        thrust::device_ptr<int> dev_counts(d_pairCounts);
        thrust::device_ptr<int> dev_offsets(d_offsets);
        thrust::exclusive_scan(thrust::device, dev_counts, dev_counts + numPrims, dev_offsets);

        // --- PHASE 3: Read Output Window Bounds ---
        int lastCount = 0;
        int lastOffset = 0;
        CUBQL_CUDA_CALL(Memcpy(&lastCount, d_pairCounts + numPrims - 1, sizeof(int), cudaMemcpyDeviceToHost));
        CUBQL_CUDA_CALL(Memcpy(&lastOffset, d_offsets + numPrims - 1, sizeof(int), cudaMemcpyDeviceToHost));
        int totalBatchPairs = lastCount + lastOffset;

        if (totalBatchPairs == 0) continue;
        
        finalCandidatePairs += totalBatchPairs;

        // --- PHASE 4: Local Write out and Sync ---
        int2* d_candidatePairs = nullptr;
        CUBQL_CUDA_CALL(Malloc(&d_candidatePairs, totalBatchPairs * sizeof(int2)));

        fillAABBOverlapsKernel_Indirected<<<gridScale, blockSize>>>(
            d_candidatePairs, d_offsets, bvhA, const_cast<cuBQL::Triangle*>(dMeshA), const_cast<cuBQL::Triangle*>(dMeshB), 
            ptr_outPrimsFlatB, startOffset, numPrims, ptr_flatNodeRootsA
        );

        cudaDeviceSynchronize();
        CUBQL_CUDA_CALL(Free(d_candidatePairs));
    }

    double tBatchLoopEnd = cuBQL::getCurrentTime();
    outLoopTime = (tBatchLoopEnd - tBatchLoopStart) * 1000.0;
    std::cout << "--------------------------------------------------\n\n";

    if (d_pairCounts) CUBQL_CUDA_CALL(Free(d_pairCounts));
    if (d_offsets) CUBQL_CUDA_CALL(Free(d_offsets));

    return finalCandidatePairs;
}