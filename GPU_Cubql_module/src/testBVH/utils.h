#ifndef UTILS_H
#define UTILS_H

#include <cuda_runtime.h>
#include "samples/common/loadOBJ.h"
#include "cuBQL/bvh.h"

// -----------------------------------------------------------------------------
// Kernel Declarations
// -----------------------------------------------------------------------------

__global__ void assembleTrianglesKernel(cuBQL::Triangle* dMesh,
                                        float2* dMetrics,
                                        const float3* dVerts,
                                        const float* dVertErrors,
                                        const uint3* dIndices,
                                        int numTriangles);

__global__ void generateBoxes(cuBQL::box3f* boxes, 
                              const cuBQL::Triangle* tris, 
                              int N);

__global__ void populateReverseMapBKernel(uint32_t* d_reverseMapB, 
                                           const uint32_t* d_markedNodeIndicesB, 
                                           uint32_t h_outMarkedCountB);

// -----------------------------------------------------------------------------
// Helper Launch Wrappers (Optional, for cleaner host code)
// -----------------------------------------------------------------------------

void launchAssembleTriangles(cuBQL::Triangle* dMesh,
                            float2* dMetrics,
                            const float3* dVerts,
                            const float* dVertErrors,
                            const uint3* dIndices,
                            int numTriangles,
                            cudaStream_t stream = 0);

void launchGenerateBoxes(cuBQL::box3f* dBoxes,
                         const cuBQL::Triangle* dTris,
                         int numTriangles,
                         cudaStream_t stream = 0);

#endif // UTILS_H