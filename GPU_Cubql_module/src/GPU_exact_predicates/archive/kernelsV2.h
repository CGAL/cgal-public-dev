#ifndef KERNELS_V2_H
#define KERNELS_V2_H

#include <vector>
#include <cuBQL/bvh.h> // Authoritative definition for cuBQL::Triangle

#include "samples/common/loadOBJ.h"

struct GPUTimingBreakdown {
    double uploadTime = 0.0;
    double executionTime = 0.0; 
    double downloadTime = 0.0;
    double bvhBuildTime = 0.0;  
    double queryTime = 0.0;     
    double countAABBTime = 0.0;        
    double fillAABBTime = 0.0;         
    double evaluateGeometricTime = 0.0; 

    long long totalCandidatePairs = 0;
};

struct IntersectionPair {
    int idA;
    int idB;
};

#ifdef __cplusplus
extern "C" {
#endif

void runPartitionedMeshIntersection(
    const cuBQL::Triangle* hA, int NA, 
    const cuBQL::Triangle* hB, int NB, 
    std::vector<IntersectionPair>& hGreenPairs,
    std::vector<IntersectionPair>& hYellowPairs,
    GPUTimingBreakdown& outTimings,
    int pipelineMode,
    int batchSize = 256000
);

#ifdef __cplusplus
}
#endif

#endif // KERNELS_V2_H