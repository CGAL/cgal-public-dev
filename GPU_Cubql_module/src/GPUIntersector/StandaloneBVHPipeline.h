
#pragma once

#include <cstdint>
#include <vector>
#include <cuda_runtime.h>

// Project Utilities & CPU Definitions (relative to src/GPUIntersector/)
#include "common/ExecutionStats.h"
#include "CPU/CgalDefinitions.h"



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


extern "C" void kernelsTestBVHV3(Mesh & meshAcpu, Mesh & meshBcpu, const double3* hVertsA,
                                 int numVertsA,
                                 const uint3* hIndicesA,
                                 int numTrianglesA,
                                 int maxCellSizeA,
                                 const double3* hVertsB,
                                 int numVertsB,
                                 const uint3* hIndicesB,
                                 int numTrianglesB,
                                 int maxCellSizeB,
                                 int batchMultiplier,
                                 int mode,
                                 int gpuPredicateMode,
                                 int leafThreshold,
                                 ExecutionStats& stats,
                                 int2*& outFinalExactPairs,      
                                size_t& outFinalCount, Point3 m_centerA, Point3 m_centerB, double3 m_rotA, double3 m_transA, double3 m_rotB, double3 m_transB) ;