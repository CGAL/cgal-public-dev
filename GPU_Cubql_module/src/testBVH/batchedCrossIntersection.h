#ifndef BATCHED_CROSS_INTERSECTION_H
#define BATCHED_CROSS_INTERSECTION_H

#include "cuBQL/bvh.h"
#include <thrust/device_vector.h>
#include <iostream> // CRITICAL: Required for std::cout and std::endl inside the print() method

#include "samples/common/loadOBJ.h"

/**
 * Executes the multi-stream batched bounding box cross-intersection loop.
 * Returns the final total number of AABB candidate pairs found.
 */
struct IntersectionTimeTracker {
    double totalLoopTimeMs   = 0.0; // Overall duration of the host execution function
    double preallocateTimeMs = 0.0; // Time spent on initial maximum bounds reduction and memory allocations
    double assemblyPhaseMs   = 0.0; // Total time spent stalling on host memory reads & filling sandbox buffers
    double executionPhaseMs  = 0.0; // Total time spent executing overlap counts, scans, and candidate writes on GPU
    double cleanupTimeMs     = 0.0; // Time taken to deallocate scratch allocations

    void print() const {
        std::cout << "\n==================================================\n";
        std::cout << "          CROSS INTERSECTION TIME PROFILING       \n";
        std::cout << "==================================================\n";
        std::cout << " Preallocation Space Discovery : " << preallocateTimeMs << " ms\n";
        std::cout << " Sandbox Assembly & Host Reads : " << assemblyPhaseMs   << " ms\n";
        std::cout << " CUDA Kernels & Scan Compute   : " << executionPhaseMs  << " ms\n";
        std::cout << " Sandbox Cleanup & Free Cycles : " << cleanupTimeMs     << " ms\n";
        std::cout << "--------------------------------------------------\n";
        std::cout << " Total Tracked Pipeline Time   : " << totalLoopTimeMs   << " ms\n";
        std::cout << "==================================================\n\n";
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
    IntersectionTimeTracker& outLoopTime
);

#endif // BATCHED_CROSS_INTERSECTION_H