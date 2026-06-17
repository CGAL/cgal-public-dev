#ifndef BATCHED_CROSS_INTERSECTION_H
#define BATCHED_CROSS_INTERSECTION_H

#include "cuBQL/bvh.h"
#include <thrust/device_vector.h>
#include <vector_types.h>
#include <iostream>
#include <vector>
#include "samples/common/loadOBJ.h"

struct IntersectionTimeTracker {
    double totalLoopTimeMs      = 0.0; 
    double preallocateTimeMs    = 0.0; 
    double assemblyPhaseMs      = 0.0; 
    double executionPhaseMs     = 0.0; 
    double fineEvaluationPhaseMs = 0.0; // Added tracker for the exact geometric predicates + compaction
    double cleanupTimeMs        = 0.0; 
    double DownloadAndClean = 0.0;

    void print() const {
        std::cout << "\n==================================================\n";
        std::cout << "          CROSS INTERSECTION TIME PROFILING       \n";
        std::cout << "==================================================\n";
        std::cout << " Preallocation Space Discovery : " << preallocateTimeMs << " ms\n";
        std::cout << " Sandbox Assembly & Host Reads : " << assemblyPhaseMs   << " ms\n";
        std::cout << " CUDA Kernels & Scan Compute   : " << executionPhaseMs  << " ms\n";
        std::cout << " Fine Geometric Evaluation     : " << fineEvaluationPhaseMs << " ms\n";
        std::cout << " Upload and clean up           : " << DownloadAndClean << " ms\n";
        std::cout << " Sandbox Cleanup & Free Cycles : " << cleanupTimeMs     << " ms\n";
        std::cout << "--------------------------------------------------\n";
        std::cout << " Total Tracked Pipeline Time   : " << totalLoopTimeMs   << " ms\n";
        std::cout << "==================================================\n\n";
    }
};

// Simple functor definition required for Thrust compaction filtering
struct IsTargetPairStatus {
    int target;
    __host__ __device__ bool operator()(const int status) const {
        return status == target;
    }
};

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
    std::vector<int2>& hGreenPairs,  // Output target for confirmed intersections
    std::vector<int2>& hYellowPairs, // Output target for coplanar / boundary elements
    IntersectionTimeTracker& outLoopTime
);

#endif // BATCHED_CROSS_INTERSECTION_H