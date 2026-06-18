#pragma once

#include <cstdint>
#include <vector>
#include <cuda_runtime.h>

#include "batchedCrossIntersection.h"

/**
 * @brief Aggregated pipeline metrics containing all printed structure data and timing records.
 */
struct ExecutionStats {
    // Structure Summaries
    uint32_t meshATotalNodes       = 0;
    uint32_t meshAExtractedTargets = 0;
    uint32_t meshBTotalNodes       = 0;
    uint32_t meshBExtractedTargets = 0;

    // Criss-Cross Bounding Box Cross-Check
    uint32_t totalIntersections     = 0;
    uint64_t totalPossiblePairs     = 0;
    double   intersectionPercentage = 0.0;

    // Timing Metrics Overview (in Milliseconds)
    double buildRefitMeshAMs   = 0.0;
    double buildRefitMeshBMs   = 0.0;
    double gpuCrossCheckEngineMs = 0.0;
    double parallelDfsDescentBMs = 0.0;

    // New 

    double GPUTotalTime = 0.0; // need to implement this

    // Dual-Tree Descent & Fine Evaluation Metrics
    int      totalCrissCrossBatches  = 0;
    uint64_t finalAabbCandidatePairs = 0;
    size_t   confirmedGreenPairs     = 0;
    size_t   confirmedYellowPairs    = 0;

    // Nested phase tracker
    IntersectionTimeTracker loopTracker;
};

// --------------------------------------------------------------------
// MAIN ENTRY TESTING PIPELINE INTERFACE
// --------------------------------------------------------------------
extern "C" {

    void kernelsTestBVH(
        const cuBQL::Triangle* hMeshA, int numTrianglesA, int maxCellSizeA,
        const cuBQL::Triangle* hMeshB, int numTrianglesB, int maxCellSizeB,
        int batchMultiplier,
        ExecutionStats& stats, std::vector<int2>& hGreenPairs,  // Output target for confirmed intersections
         std::vector<int2>& hYellowPairs // Output target for coplanar / boundary elements
    );

}