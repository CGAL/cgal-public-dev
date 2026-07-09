#pragma once

#include <cstdint>
#include <vector>
#include <cuda_runtime.h>

#include "batchedCrossIntersection.h"
#include "ExecutionStats.h"


// --------------------------------------------------------------------
// MAIN ENTRY TESTING PIPELINE INTERFACE
// --------------------------------------------------------------------
extern "C" void kernelsTestBVHV2(
    const float3* hVertsA, int numVertsA, const uint3* hIndicesA, int numTrianglesA, int maxCellSizeA,
    const float3* hVertsB, int numVertsB, const uint3* hIndicesB, int numTrianglesB, int maxCellSizeB,
    int batchMultiplier, int mode, int leafThreshold,
    ExecutionStats& stats,
    std::vector<int2>& hGreenPairs,
    std::vector<int2>& hYellowPairs
);