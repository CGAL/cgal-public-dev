#pragma once

#include <cstdint>
#include <vector>
#include <cuda_runtime.h>

//#include "cuBQL/bvh.h"
//#include <thrust/device_vector.h>

//#include "../testBVH/batchedCrossIntersection.h"
#include "samples/common/loadOBJ.h"

#include "../testBVH/ExecutionStats.h"

// --------------------------------------------------------------------
// MAIN ENTRY TESTING PIPELINE INTERFACE
// --------------------------------------------------------------------
extern "C" {

    void dualTreeKernelLaunch(const cuBQL::Triangle* hMeshA, int numTrianglesA, int maxCellSizeA,
                               const cuBQL::Triangle* hMeshB, int numTrianglesB, int maxCellSizeB,
                               int batchMultiplier, float fillRatio,
                               ExecutionStats& stats, 
                               std::vector<int2>& hGreenPairs,  
                               std::vector<int2>& hYellowPairs);

}