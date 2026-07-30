#ifndef UTILS_H
#define UTILS_H

#include <cuda_runtime.h>
#include "samples/common/loadOBJ.h"
#include "cuBQL/bvh.h"
#include "../custom_pipeline/TriangleDouble.h"

// =============================================================================
// Kernel Declarations (float3 - Legacy / Float Pipeline)
// =============================================================================

__global__ void assembleTrianglesKernel(cuBQL::Triangle* dMesh,
                                        float2* dMetrics,
                                        const float3* dVerts,
                                        const float* dVertErrors,
                                        const uint3* dIndices,
                                        int numTriangles);

__global__ void assembleTrianglesKernelTransformed(cuBQL::Triangle* dMesh,
                                                   float2* dMetrics,
                                                   const float3* dVerts,
                                                   const float* dVertErrors,
                                                   const uint3* dIndices,
                                                   int numTriangles, 
                                                   bool isTranslated,
                                                   bool isRotated);

// =============================================================================
// Kernel Declarations (double3 Pipeline)
// =============================================================================

// Primary Kernel (Handles both Float-Only and Dual Double+Float assembly via nullptr)
__global__ void assembleTrianglesKernel(TriangleDouble* dMeshDouble,
                                        cuBQL::Triangle* dMesh,
                                        float2* dMetrics,
                                        const double3* dVerts,
                                        const uint3* dIndices,
                                        int numTriangles);

// Legacy/Convenience Overload (for calls omitting dMeshDouble)
__global__ void assembleTrianglesKernel(cuBQL::Triangle* dMesh,
                                        float2* dMetrics,
                                        const double3* dVerts,
                                        const uint3* dIndices,
                                        int numTriangles);

// =============================================================================
// General Utility Kernels
// =============================================================================

__global__ void generateBoxes(cuBQL::box3f* boxes, 
                              const cuBQL::Triangle* tris, 
                              int N);

__global__ void populateReverseMapBKernel(uint32_t* d_reverseMapB, 
                                          const uint32_t* d_markedNodeIndicesB, 
                                          uint32_t h_outMarkedCountB);

// =============================================================================
// Host Launch Wrappers
// =============================================================================

void launchAssembleTriangles(cuBQL::Triangle* dMesh,
                             float2* dMetrics,
                             const float3* dVerts,
                             const float* dVertErrors,
                             const uint3* dIndices,
                             int numTriangles,
                             cudaStream_t stream = 0);

void launchAssembleTrianglesTransformed(cuBQL::Triangle* dMesh,
                                         float2* dMetrics,
                                         const float3* dVerts,
                                         const float* dVertErrors,
                                         const uint3* dIndices,
                                         int numTriangles,
                                         bool isTranslated,
                                         bool isRotated,
                                         cudaStream_t stream = 0);

// Host Launchers for double3
void launchAssembleTriangles(TriangleDouble* dMeshDouble,
                             cuBQL::Triangle* dMesh,
                             float2* dMetrics,
                             const double3* dVerts,
                             const uint3* dIndices,
                             int numTriangles,
                             cudaStream_t stream = 0);

void launchAssembleTriangles(cuBQL::Triangle* dMesh,
                             float2* dMetrics,
                             const double3* dVerts,
                             const uint3* dIndices,
                             int numTriangles,
                             cudaStream_t stream = 0);

void launchGenerateBoxes(cuBQL::box3f* dBoxes,
                         const cuBQL::Triangle* dTris,
                         int numTriangles,
                         cudaStream_t stream = 0);

#endif // UTILS_H