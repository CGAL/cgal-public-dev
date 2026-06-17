#define CUBQL_GPU_BUILDER_IMPLEMENTATION 1
#include "cuBQL/bvh.h"
#include "cuBQL/traversal/fixedBoxQuery.h"
#include <thrust/device_ptr.h>
#include <thrust/scan.h>
#include <thrust/copy.h>
#include <thrust/execution_policy.h>
#include <thrust/device_vector.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include "samples/common/loadOBJ.h"

#include "GPUPredicatesCheck.h"
#include "kernelsV2.h"

// struct GPUTimingBreakdown {
//     double uploadTime = 0.0;
//     double executionTime = 0.0; 
//     double downloadTime = 0.0;
//     double bvhBuildTime = 0.0;  
//     double queryTime = 0.0;     
//     double countAABBTime = 0.0;        
//     double fillAABBTime = 0.0;         
//     double evaluateGeometricTime = 0.0; 

//     long long totalCandidatePairs = 0;
// };


// --------------------------------------------------------------------
// KERNELS: Accelerated Vector Clean Structure Passes
// --------------------------------------------------------------------
__global__ void generateBoxes(cuBQL::box3f *boxes, const cuBQL::Triangle *tris, int N) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < N) boxes[i] = tris[i].bounds();
}

__global__ void countAABBOverlapsKernel(int *pairCounts, cuBQL::bvh3f bvhA, const cuBQL::Triangle *triA, const cuBQL::Triangle *triB, int currentBatchSize) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= currentBatchSize) return;

    cuBQL::Triangle b = triB[tid];
    cuBQL::box3f query = b.bounds();
    int count = 0;

    cuBQL::fixedBoxQuery::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        for (uint32_t i = 0; i < num; i++) {
            if (triA[ids[i]].bounds().overlaps(query)) {
                count++;
            }
        }
        return CUBQL_CONTINUE_TRAVERSAL;
    }, bvhA, query);

    pairCounts[tid] = count;
}

__global__ void fillAABBOverlapsKernel(int2 *candidatePairs, const int *offsets, cuBQL::bvh3f bvhA, const cuBQL::Triangle *triA, const cuBQL::Triangle *triB, int currentBatchSize) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= currentBatchSize) return;

    int wPos = offsets[tid];
    cuBQL::Triangle b = triB[tid];
    cuBQL::box3f query = b.bounds();

    cuBQL::fixedBoxQuery::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        for (uint32_t i = 0; i < num; i++) {
            if (triA[ids[i]].bounds().overlaps(query)) {
                candidatePairs[wPos++] = make_int2((int)ids[i], tid);
            }
        }
        return CUBQL_CONTINUE_TRAVERSAL;
    }, bvhA, query);
}

// --------------------------------------------------------------------
// THRUST COMPACTION PREDICATE STRUCTURE
// --------------------------------------------------------------------
struct IsTargetPairStatus {
    int targetStatus;
    __host__ __device__ bool operator()(const int &status) const {
        return status == targetStatus;
    }
};

// --------------------------------------------------------------------
// MAIN ENTRY PIPELINE
// --------------------------------------------------------------------
extern "C" void runPartitionedMeshIntersection(
    const cuBQL::Triangle* hA, int NA, 
    const cuBQL::Triangle* hB, int NB, 
    std::vector<IntersectionPair>& hGreenPairs,
    std::vector<IntersectionPair>& hYellowPairs,
    GPUTimingBreakdown& outTimings,
    int pipelineMode,
    int batchSize) 
{
    double t0, t_exec;

    int currentDevice = 0;
    if (cudaGetDevice(&currentDevice) == cudaSuccess) {
        cudaSetDevice(currentDevice);
    }

    outTimings.countAABBTime = 0.0;
    outTimings.fillAABBTime = 0.0;
    outTimings.evaluateGeometricTime = 0.0;

    outTimings.totalCandidatePairs = 0;

    cudaEvent_t startCount, stopCount;
    cudaEvent_t startFill, stopFill;
    cudaEvent_t startEval, stopEval;
    cudaEventCreate(&startCount); cudaEventCreate(&stopCount);
    cudaEventCreate(&startFill);  cudaEventCreate(&stopFill);
    cudaEventCreate(&startEval);  cudaEventCreate(&stopEval);

    // 1. UPGRADE: Allocation and immediate upload of BOTH complete mesh profiles upfront
    t0 = cuBQL::getCurrentTime();
    cuBQL::Triangle *dA;
    cuBQL::Triangle *dB_full;
    CUBQL_CUDA_CALL(Malloc(&dA, NA * sizeof(cuBQL::Triangle)));
    CUBQL_CUDA_CALL(Malloc(&dB_full, NB * sizeof(cuBQL::Triangle)));
    
    CUBQL_CUDA_CALL(Memcpy(dA, hA, NA * sizeof(cuBQL::Triangle), cudaMemcpyDefault));
    CUBQL_CUDA_CALL(Memcpy(dB_full, hB, NB * sizeof(cuBQL::Triangle), cudaMemcpyDefault));
    outTimings.uploadTime = cuBQL::getCurrentTime() - t0;

    // --- Start Total Execution Timing ---
    t_exec = cuBQL::getCurrentTime();
    
    // ==========================================
    // PHASE A: BVH Build
    // ==========================================
    double t_bvh_start = cuBQL::getCurrentTime();

    cuBQL::box3f *dBoxes;
    CUBQL_CUDA_CALL(Malloc(&dBoxes, NA * sizeof(cuBQL::box3f)));
    generateBoxes<<<cuBQL::divRoundUp(NA, 256), 256>>>(dBoxes, dA, NA);
    cudaDeviceSynchronize(); 
    
    cuBQL::bvh3f bvh;
    cuBQL::gpuBuilder(bvh, dBoxes, NA, cuBQL::BuildConfig());
    CUBQL_CUDA_CALL(Free(dBoxes));
    cudaDeviceSynchronize(); 

    outTimings.bvhBuildTime = cuBQL::getCurrentTime() - t_bvh_start;

    hGreenPairs.clear();
    hYellowPairs.clear();
    
    // Scratch execution space metrics scaled strictly to local batch sizing boundaries
    int *dPairCounts, *dOffsets;
    CUBQL_CUDA_CALL(Malloc(&dPairCounts, batchSize * sizeof(int)));
    CUBQL_CUDA_CALL(Malloc(&dOffsets, batchSize * sizeof(int)));

    // ==========================================
    // PHASE B: Query & Evaluation Loop
    // ==========================================
    double t_query_start = cuBQL::getCurrentTime();

    for (int bStart = 0; bStart < NB; bStart += batchSize) {
        int currentBatchSize = std::min(batchSize, NB - bStart);

        // Calculate the device memory offset address pointer dynamically
        cuBQL::Triangle* dB_batch = dB_full + bStart;

        // Time countAABBOverlapsKernel
        cudaEventRecord(startCount, 0);
        countAABBOverlapsKernel<<<cuBQL::divRoundUp(currentBatchSize, 128), 128>>>(
            dPairCounts, bvh, dA, dB_batch, currentBatchSize
        );
        cudaEventRecord(stopCount, 0);
        
        thrust::device_ptr<int> dev_counts(dPairCounts);
        thrust::device_ptr<int> dev_offsets(dOffsets);
        thrust::exclusive_scan(thrust::device, dev_counts, dev_counts + currentBatchSize, dev_offsets);
        
        int totalBatchPairs = 0, lastCount = 0;
        CUBQL_CUDA_CALL(Memcpy(&totalBatchPairs, dOffsets + currentBatchSize - 1, sizeof(int), cudaMemcpyDeviceToHost));
        CUBQL_CUDA_CALL(Memcpy(&lastCount, dPairCounts + currentBatchSize - 1, sizeof(int), cudaMemcpyDeviceToHost));
        totalBatchPairs += lastCount;

        outTimings.totalCandidatePairs += totalBatchPairs;

        cudaEventSynchronize(stopCount);
        float millisecondsCount = 0;
        cudaEventElapsedTime(&millisecondsCount, startCount, stopCount);
        outTimings.countAABBTime += (double)millisecondsCount / 1000.0;

        if (totalBatchPairs == 0) continue;

        int2 *dCandidatePairs;
       // IntersectionPair *dEvaluatedPairs;
        int *dPairStatuses;
        CUBQL_CUDA_CALL(Malloc(&dCandidatePairs, totalBatchPairs * sizeof(int2)));
       // CUBQL_CUDA_CALL(Malloc(&dEvaluatedPairs, totalBatchPairs * sizeof(IntersectionPair)));
        CUBQL_CUDA_CALL(Malloc(&dPairStatuses, totalBatchPairs * sizeof(int)));

        // Time fillAABBOverlapsKernel
        cudaEventRecord(startFill, 0);
        fillAABBOverlapsKernel<<<cuBQL::divRoundUp(currentBatchSize, 128), 128>>>(
            dCandidatePairs, dOffsets, bvh, dA, dB_batch, currentBatchSize
        );
        cudaEventRecord(stopFill, 0);

        // Sync and record Fill timing immediately before passing arrays to the evaluation stage
        cudaEventSynchronize(stopFill);
        float millisecondsFill = 0;
        cudaEventElapsedTime(&millisecondsFill, startFill, stopFill);
        outTimings.fillAABBTime += (double)millisecondsFill / 1000.0;

        // ====================================================================
        // MODULAR REFACTOR: Call your streamlined geometric evaluation step
        // ====================================================================
        evaluateAndCompactPairs(
            dCandidatePairs, 
            dPairStatuses, 
            dA, 
            dB_batch, 
            totalBatchPairs, 
            outTimings.evaluateGeometricTime
        );

        // ====================================================================
        // PHASE C: Thrust Compaction & Download (Streamlined)
        // ====================================================================
        // Reinterpret the int2* pointer as an IntersectionPair* pointer
        IntersectionPair* dCandidateAsEvaluated = reinterpret_cast<IntersectionPair*>(dCandidatePairs);

        thrust::device_ptr<IntersectionPair> dev_evaluated(dCandidateAsEvaluated);
        thrust::device_ptr<int> dev_statuses(dPairStatuses);

        thrust::device_vector<IntersectionPair> dev_green(totalBatchPairs);
        thrust::device_vector<IntersectionPair> dev_yellow(totalBatchPairs);

        thrust::device_ptr<IntersectionPair> dev_green_out(thrust::raw_pointer_cast(dev_green.data()));
        thrust::device_ptr<IntersectionPair> dev_yellow_out(thrust::raw_pointer_cast(dev_yellow.data()));

        auto green_end = thrust::copy_if(thrust::device, dev_evaluated, dev_evaluated + totalBatchPairs, dev_statuses, dev_green_out, IsTargetPairStatus{(int)PAIR_GREEN});
        auto yellow_end = thrust::copy_if(thrust::device, dev_evaluated, dev_evaluated + totalBatchPairs, dev_statuses, dev_yellow_out, IsTargetPairStatus{(int)PAIR_YELLOW});

        int totalGreen = green_end - dev_green_out;
        int totalYellow = yellow_end - dev_yellow_out;

        if (totalGreen > 0) {
            size_t oldSize = hGreenPairs.size();
            hGreenPairs.resize(oldSize + totalGreen);
            
            std::vector<IntersectionPair> hTmpGreen(totalGreen);
            CUBQL_CUDA_CALL(Memcpy(hTmpGreen.data(), thrust::raw_pointer_cast(dev_green.data()), totalGreen * sizeof(IntersectionPair), cudaMemcpyDefault));
            for (auto& pair : hTmpGreen) {
                pair.idB += bStart; 
            }
            std::copy(hTmpGreen.begin(), hTmpGreen.end(), hGreenPairs.begin() + oldSize);
        }
        if (totalYellow > 0) {
            size_t oldSize = hYellowPairs.size();
            hYellowPairs.resize(oldSize + totalYellow);
            
            std::vector<IntersectionPair> hTmpYellow(totalYellow);
            CUBQL_CUDA_CALL(Memcpy(hTmpYellow.data(), thrust::raw_pointer_cast(dev_yellow.data()), totalYellow * sizeof(IntersectionPair), cudaMemcpyDefault));
            for (auto& pair : hTmpYellow) {
                pair.idB += bStart; 
            }
            std::copy(hTmpYellow.begin(), hTmpYellow.end(), hYellowPairs.begin() + oldSize);
        }

        // Safe Cleanup (dEvaluatedPairs is completely gone!)
        CUBQL_CUDA_CALL(Free(dCandidatePairs));
        CUBQL_CUDA_CALL(Free(dPairStatuses));
    }
    cudaDeviceSynchronize(); 

    outTimings.queryTime = cuBQL::getCurrentTime() - t_query_start;

    // --- End Total Execution Timing ---
    outTimings.executionTime = cuBQL::getCurrentTime() - t_exec;
    outTimings.downloadTime = 0.0; 

    cudaEventDestroy(startCount); cudaEventDestroy(stopCount);
    cudaEventDestroy(startFill);  cudaEventDestroy(stopFill);
    cudaEventDestroy(startEval);  cudaEventDestroy(stopEval);

    CUBQL_CUDA_CALL(Free(dA));
    CUBQL_CUDA_CALL(Free(dB_full)); // Cleanup global array memory map
    CUBQL_CUDA_CALL(Free(dPairCounts));
    CUBQL_CUDA_CALL(Free(dOffsets));
    cuBQL::cuda::free(bvh);
}