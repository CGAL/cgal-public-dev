#pragma once

#include <cstdint>
#include <vector>
#include <cuda_runtime.h>

#include "batchedCrossIntersection.h"
#include "ExecutionStats.h"


// --------------------------------------------------------------------
// MAIN ENTRY TESTING PIPELINE INTERFACE
// --------------------------------------------------------------------
extern "C" {

    void kernelsTestBVH(
        const cuBQL::Triangle* hMeshA, int numTrianglesA, int maxCellSizeA,
        const cuBQL::Triangle* hMeshB, int numTrianglesB, int maxCellSizeB,
        int batchMultiplier,int mode,
        ExecutionStats& stats, std::vector<int2>& hGreenPairs,  // Output target for confirmed intersections
         std::vector<int2>& hYellowPairs // Output target for coplanar / boundary elements
    );

}