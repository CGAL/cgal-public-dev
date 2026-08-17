#pragma once

#include <cstdint>
#include <vector>
#include <cuda_runtime.h>

#include "batchedCrossIntersection.h"
#include "ExecutionStats.h"
#include "../src/CPU/CgalDefinitions.h"
#include <tbb/concurrent_vector.h>



// --------------------------------------------------------------------
// MAIN ENTRY TESTING PIPELINE INTERFACE
// --------------------------------------------------------------------
// extern "C" void kernelsTestBVHV2(
//     const float3* hVertsA, int numVertsA, const uint3* hIndicesA, int numTrianglesA, int maxCellSizeA,
//     const float3* hVertsB, int numVertsB, const uint3* hIndicesB, int numTrianglesB, int maxCellSizeB,
//     int batchMultiplier, int mode, int leafThreshold,
//     ExecutionStats& stats,
//     std::vector<int2>& hGreenPairs,
//     std::vector<int2>& hYellowPairs
// );


extern "C" void kernelsTestBVHV2(Mesh & meshAcpu, Mesh & meshBcpu, const float3* hVertsA,
                                 int numVertsA,
                                 const uint3* hIndicesA,
                                 const float* hVertErrorsA,
                                 int numTrianglesA,
                                 int maxCellSizeA,
                                 const float3* hVertsB,
                                 int numVertsB,
                                 const uint3* hIndicesB,
                                 const float* hVertErrorsB,
                                 int numTrianglesB,
                                 int maxCellSizeB,
                                 int batchMultiplier,
                                 int mode,
                                 int leafThreshold,
                                 int activateAsyncDownload,
                                 ExecutionStats& stats,
                                 tbb::concurrent_vector<int2> & finalExactPairs) ;