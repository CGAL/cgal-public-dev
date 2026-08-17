#ifndef PRUNE_PIPELINE_H
#define BPRUNE_PIPELINE_H

#include "cuBQL/bvh.h"
#include <thrust/device_vector.h>
#include <vector_types.h>
#include <iostream>
#include <vector>
#include "samples/common/loadOBJ.h"

#include "common/IntersectionTimeTracker.h"

// Simple functor definition required for Thrust compaction filtering

void parallelPruneAndReindexAll(
    uint32_t* d_intersectPairsA, // In/Out
    int numIntersections,
    
    // In/Out Structural Variables to update for Main Function
    uint32_t*& d_outSortedPrimIDsA,
    uint32_t*& d_outNodeOffsetsA,
    uint32_t& outTotalActiveCellsA,
    int& outNumPrimsA,
    
    cudaStream_t s,
    cuBQL::DeviceMemoryResource& memResource
);

#endif // PRUNE_PIPELINE