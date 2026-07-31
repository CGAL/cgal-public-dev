#include "batchedCrossIntersectionCommon.h"

__global__ void computeMaxChunkSizeKernel(
    const uint32_t* outPairsB,
    const uint32_t* outOffsetsB,
    const uint32_t* reverseMapB,
    uint32_t outMarkedCountB,
    uint32_t totalPrimsB,
    int totalBatches,
    int batchMultiplier,
    int numChunks,
    int* d_globalMax
) {
    int chunkIdx = blockIdx.x;
    if (chunkIdx >= numChunks) return;

    int startBatch = chunkIdx * batchMultiplier;
    int activeBatchesInChunk = (startBatch + batchMultiplier <= totalBatches) 
                               ? batchMultiplier 
                               : (totalBatches - startBatch);

    int threadSum = 0;

    for (int b = threadIdx.x; b < activeBatchesInChunk; b += blockDim.x) {
        uint32_t bIdx = outPairsB[startBatch + b];
        uint32_t markedIdxB = reverseMapB[bIdx];
        uint32_t startOffset = outOffsetsB[markedIdxB];
        uint32_t endOffset = (markedIdxB + 1 < outMarkedCountB) ? outOffsetsB[markedIdxB + 1] : totalPrimsB;
        threadSum += (int)(endOffset - startOffset);
    }

    extern __shared__ int sdata[];
    int tid = threadIdx.x;
    sdata[tid] = threadSum;
    __syncthreads();

    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sdata[tid] += sdata[tid + s];
        }
        __syncthreads();
    }

    if (tid == 0) {
        atomicMax(d_globalMax, sdata[0]);
    }
}

__global__ void assembleChunkBuffersByBatchKernel(
    uint32_t* d_BIter,
    uint64_t* d_AIter,
    const uint32_t* outPairsA,
    const uint32_t* outPairsB,
    const uint32_t* reverseMapB,
    const uint32_t* outOffsetsB,
    const uint32_t* outPrimsFlatB,
    uint32_t outMarkedCountB,
    uint32_t totalPrimsB,
    int chunkStartBatchIdx,
    int activeBatchesInChunk,
    const int* d_chunkBatchOffsets 
) {
    int b = threadIdx.x + blockIdx.x * blockDim.x;
    if (b >= activeBatchesInChunk) return;

    int currentBatchArrayIdx = chunkStartBatchIdx + b;
    
    uint32_t bIdx = outPairsB[currentBatchArrayIdx];
    uint32_t markedIdxB = reverseMapB[bIdx];
    
    uint32_t startOffsetB = outOffsetsB[markedIdxB];
    uint32_t endOffsetB = (markedIdxB + 1 < outMarkedCountB) ? outOffsetsB[markedIdxB + 1] : totalPrimsB;
    int numPrims = (int)(endOffsetB - startOffsetB);
    
    if (numPrims <= 0) return;

    int writeSandboxOffset = d_chunkBatchOffsets[b];
    uint64_t startNodeIdxA = (uint64_t)(outPairsA[currentBatchArrayIdx]);

    for (int p = 0; p < numPrims; ++p) {
        int targetWriteIdx = writeSandboxOffset + p;
        d_BIter[targetWriteIdx] = outPrimsFlatB[startOffsetB + p];
        d_AIter[targetWriteIdx] = startNodeIdxA;
    }
}

__global__ void countAABBOverlapsKernel_Indirected(
    int *pairCounts, 
    cuBQL::bvh3f bvhA, 
    const cuBQL::Triangle *triA, 
    const cuBQL::Triangle *triB, 
    const uint32_t* d_BIter,       
    uint32_t startOffsetB,         
    int numPrimsB, 
    const uint64_t* d_AIter)       
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= numPrimsB) return;

    uint32_t actualPrimIdB = d_BIter[startOffsetB + tid];
    uint64_t startNodeIdxA = d_AIter[startOffsetB + tid]; 
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

__global__ void countAABBOverlapsKernel_Indirected_boxes(
    int *pairCounts, 
    cuBQL::bvh3f bvhA, 
    const cuBQL::box3f *d_boxesA, 
    const cuBQL::box3f *d_boxesB, 
    const uint32_t* d_BIter,       
    uint32_t startOffsetB,         
    int numPrimsB, 
    const uint64_t* d_AIter)       
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= numPrimsB) return;

    uint32_t actualPrimIdB = d_BIter[startOffsetB + tid];
    uint64_t startNodeIdxA = d_AIter[startOffsetB + tid]; 

    // Direct 24-byte read of the query box
    cuBQL::box3f query = d_boxesB[actualPrimIdB];
    
    int count = 0;
    cuBQL::fixedBoxQueryv2::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        for (uint32_t i = 0; i < num; i++) {
            // Direct AABB overlap check on target box array
            if (d_boxesA[ids[i]].overlaps(query)) { 
                count++; 
            }
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
    const uint32_t* d_BIter,       
    uint32_t startOffsetB,         
    int numPrimsB, 
    const uint64_t* d_AIter)       
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= numPrimsB) return;

    int wPos = offsets[tid];
    uint32_t actualPrimIdB = d_BIter[startOffsetB + tid];
    uint64_t startNodeIdxA = d_AIter[startOffsetB + tid]; 
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

__global__ void fillAABBOverlapsKernel_Indirected_Boxes(
    int2 *candidatePairs, 
    const int *offsets, 
    cuBQL::bvh3f bvhA, 
    const cuBQL::box3f *d_boxesA, 
    const cuBQL::box3f *d_boxesB, 
    const uint32_t* d_BIter,       
    uint32_t startOffsetB,         
    int numPrimsB, 
    const uint64_t* d_AIter)       
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= numPrimsB) return;

    int wPos = offsets[tid];
    uint32_t actualPrimIdB = d_BIter[startOffsetB + tid];
    uint64_t startNodeIdxA = d_AIter[startOffsetB + tid]; 

    // Fetch bounding box directly for query
    cuBQL::box3f query = d_boxesB[actualPrimIdB];
    
    cuBQL::fixedBoxQueryv2::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        for (uint32_t i = 0; i < num; i++) {
            // Direct AABB-AABB overlap test against d_boxesA
            if (d_boxesA[ids[i]].overlaps(query)) { 
                candidatePairs[wPos++] = make_int2(static_cast<int>(ids[i]), static_cast<int>(actualPrimIdB)); 
            }
        }
        return CUBQL_CONTINUE_TRAVERSAL;
    }, bvhA, query, startNodeIdxA);
}