#ifndef ExecutionStats_H
#define ExecutionStats_H

#include <cstdint>
#include <vector>
#include <cuda_runtime.h>
#include "IntersectionTimeTracker.h"

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
    double buildRefitMeshAMs     = 0.0;
    double buildRefitMeshBMs     = 0.0;
    double gpuCrossCheckEngineMs = 0.0;
    double dualTreeStepMs        = 0.0; // <-- Added to track the down-level expansion step
    double parallelDfsDescentBMs = 0.0;

    // Detailed Allocation & Framework Overhead Tracks
    double initialAllocAndCopyMs = 0.0; // Captures initial raw cudaMalloc + cudaMemcpy
    double thrustInitOverheadMs  = 0.0; // Captures Thrust vector allocation + filling overhead
    double finalCleanupSyncMs    = 0.0; // Captures cudaFrees and destructor sync delays at the end
    
    double GPUTotalTime          = 0.0; // Comprehensive GPU Pipeline Time

    // Dual-Tree Descent & Fine Evaluation Metrics
    int      totalCrissCrossBatches  = 0;
    uint64_t finalAabbCandidatePairs = 0;
    size_t   confirmedGreenPairs     = 0;
    size_t   confirmedYellowPairs    = 0;

    // Nested phase tracker
    IntersectionTimeTracker loopTracker;
};

#endif