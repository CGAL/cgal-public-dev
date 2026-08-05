#ifndef UTILS_H
#define UTILS_H

#include <cuda_runtime.h>
#include "samples/common/loadOBJ.h"
#include "cuBQL/bvh.h"
#include "../custom_pipeline/TriangleDouble.h"

// =============================================================================
// Kernel Declarations (float3 - Legacy / Float Pipeline)
// =============================================================================

// -----------------------------------------------------------------------------
// Helper 3x3 Matrix Structure & Rigid Body Transformations
// -----------------------------------------------------------------------------
struct Mat3x3
{
  double m[3][3];
};

inline Mat3x3 makeRotationMatrixDeg(double rxDeg, double ryDeg, double rzDeg) {
  const double DEG_TO_RAD = 3.14159265358979323846 / 180.0;
  double radX = rxDeg * DEG_TO_RAD;
  double radY = ryDeg * DEG_TO_RAD;
  double radZ = rzDeg * DEG_TO_RAD;

  double cx = cos(radX), sx = sin(radX);
  double cy = cos(radY), sy = sin(radY);
  double cz = cos(radZ), sz = sin(radZ);

  Mat3x3 R;
  R.m[0][0] = cy * cz;
  R.m[0][1] = sx * sy * cz - cx * sz;
  R.m[0][2] = cx * sy * cz + sx * sz;

  R.m[1][0] = cy * sz;
  R.m[1][1] = sx * sy * sz + cx * cz;
  R.m[1][2] = cx * sy * sz - sx * cz;

  R.m[2][0] = -sy;
  R.m[2][1] = sx * cy;
  R.m[2][2] = cx * cy;

  return R;
}

// =============================================================================
// Vertex Transformation Kernels & Host Launchers
// =============================================================================

__global__ void transformVerticesKernel(double3* dVertsOut, 
                                        const double3* dVertsOrig, 
                                        int numVerts, 
                                        Mat3x3 R, 
                                        double3 center, 
                                        double3 trans);

void launchTransformVertices(double3* dVertsOut,
                             const double3* dVertsOrig,
                             int numVerts,
                             Mat3x3 R,
                             double3 center,
                             double3 trans,
                             cudaStream_t stream = 0);

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
// Box Generation Kernels
// =============================================================================

__global__ void generateBoxesTrisKernel(cuBQL::box3f* __restrict__ dBoxes,
                                        const double3* __restrict__ dVerts,
                                        const uint3* __restrict__ dIndices,
                                        int numTriangles);

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

void launchGenerateBoxesTris(cuBQL::box3f* dBoxes,
                             const double3* dVerts,
                             const uint3* dIndices, // Or const uint3* depending on your types
                             int numTriangles,
                             cudaStream_t stream = 0);

#endif // UTILS_H