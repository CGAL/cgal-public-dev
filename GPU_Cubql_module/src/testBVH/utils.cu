#include "utils.h"
#include "cuBQL/math/vec.h"
#include <cmath>
#include "../custom_pipeline/TriangleDouble.h"

// =============================================================================
// Internal Device Workhorse Functions
// =============================================================================

// -----------------------------------------------------------------------------
// 1. float3 Device Implementation (Legacy / Float Pipeline)
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// 2. double3 Device Implementation (Double -> Float Downcast Pipeline)
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// 2. double3 Device Implementation (Exact Downcast Drift & Metrics)
// -----------------------------------------------------------------------------
// __device__ __forceinline__ void assembleTrianglesDeviceImpl(
//     cuBQL::Triangle* dMesh, float2* dMetrics, const double3* dVerts, const uint3* dIndices, int numTriangles) {
//   int idx = threadIdx.x + blockIdx.x * blockDim.x;
//   if(idx >= numTriangles)
//     return;

//   uint3 triIdx = dIndices[idx];

//   double3 p0 = dVerts[triIdx.x];
//   double3 p1 = dVerts[triIdx.y];
//   double3 p2 = dVerts[triIdx.z];

//   // 1. Quantize double3 -> float components
//   float3 a = make_float3(static_cast<float>(p0.x), static_cast<float>(p0.y), static_cast<float>(p0.z));
//   float3 b = make_float3(static_cast<float>(p1.x), static_cast<float>(p1.y), static_cast<float>(p1.z));
//   float3 c = make_float3(static_cast<float>(p2.x), static_cast<float>(p2.y), static_cast<float>(p2.z));

//   // CORRECT
//   dMesh[idx].a.x = a.x;
//   dMesh[idx].a.y = a.y;
//   dMesh[idx].a.z = a.z;
//   dMesh[idx].b.x = b.x;
//   dMesh[idx].b.y = b.y;
//   dMesh[idx].b.z = b.z; // Fixed a.z -> b.z
//   dMesh[idx].c.x = c.x;
//   dMesh[idx].c.y = c.y;
//   dMesh[idx].c.z = c.z;

//   // 2. Compute L (Absolute Spatial Scale) and E (Max Edge Axis Extent) in double precision
//   double l0 = fmax(fmax(fabs(p0.x), fabs(p0.y)), fabs(p0.z));
//   double l1 = fmax(fmax(fabs(p1.x), fabs(p1.y)), fabs(p1.z));
//   double l2 = fmax(fmax(fabs(p2.x), fabs(p2.y)), fabs(p2.z));
//   float L = static_cast<float>(fmax(fmax(l0, l1), l2));

//   double ex = fmax(fmax(fabs(p0.x - p1.x), fabs(p1.x - p2.x)), fabs(p2.x - p0.x));
//   double ey = fmax(fmax(fabs(p0.y - p1.y), fabs(p1.y - p2.y)), fabs(p2.y - p0.y));
//   double ez = fmax(fmax(fabs(p0.z - p1.z), fabs(p1.z - p2.z)), fabs(p2.z - p0.z));
//   float E = static_cast<float>(fmax(fmax(ex, ey), ez));

//   // 3. Compute EXACT Quantization Drift per vertex: max(|double_pos - float_pos|)
//   double err0 = fmax(fmax(fabs(p0.x - static_cast<double>(a.x)), fabs(p0.y - static_cast<double>(a.y))),
//                      fabs(p0.z - static_cast<double>(a.z)));

//   double err1 = fmax(fmax(fabs(p1.x - static_cast<double>(b.x)), fabs(p1.y - static_cast<double>(b.y))),
//                      fabs(p1.z - static_cast<double>(b.z)));

//   double err2 = fmax(fmax(fabs(p2.x - static_cast<double>(c.x)), fabs(p2.y - static_cast<double>(c.y))),
//                      fabs(p2.z - static_cast<double>(c.z)));

//   float maxVertexError = static_cast<float>(fmax(fmax(err0, err1), err2));

//   // 5. Expand scale and edge metrics strictly by the exact error bound
//   L += maxVertexError;
//   E += (2.0f * maxVertexError);

//   // 6. Write out packed parameters
//   dMetrics[idx] = make_float2(L, E * E);
// }
__device__ __forceinline__ void assembleTrianglesDeviceImpl(cuBQL::Triangle* dMesh,
                                                            float2* dMetrics,
                                                            const float3* dVerts,
                                                            const float* dVertErrors,
                                                            const uint3* dIndices,
                                                            int numTriangles,
                                                            bool isTranslated,
                                                            bool isRotated) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if(idx >= numTriangles)
    return;

  uint3 triIdx = dIndices[idx];

  float3 p0 = dVerts[triIdx.x];
  float3 p1 = dVerts[triIdx.y];
  float3 p2 = dVerts[triIdx.z];

  dMesh[idx].a.x = p0.x;
  dMesh[idx].a.y = p0.y;
  dMesh[idx].a.z = p0.z;
  dMesh[idx].b.x = p1.x;
  dMesh[idx].b.y = p1.y;
  dMesh[idx].b.z = p1.z;
  dMesh[idx].c.x = p2.x;
  dMesh[idx].c.y = p2.y;
  dMesh[idx].c.z = p2.z;

  float l0 = fmaxf(fmaxf(fabsf(p0.x), fabsf(p0.y)), fabsf(p0.z));
  float l1 = fmaxf(fmaxf(fabsf(p1.x), fabsf(p1.y)), fabsf(p1.z));
  float l2 = fmaxf(fmaxf(fabsf(p2.x), fabsf(p2.y)), fabsf(p2.z));
  float L = fmaxf(fmaxf(l0, l1), l2);

  float ex = fmaxf(fmaxf(fabsf(p0.x - p1.x), fabsf(p1.x - p2.x)), fabsf(p2.x - p0.x));
  float ey = fmaxf(fmaxf(fabsf(p0.y - p1.y), fabsf(p1.y - p2.y)), fabsf(p2.y - p0.y));
  float ez = fmaxf(fmaxf(fabsf(p0.z - p1.z), fabsf(p1.z - p2.z)), fabsf(p2.z - p0.z));
  float E = fmaxf(fmaxf(ex, ey), ez);

  float e0, e1, e2;
  if(dVertErrors != nullptr) {
    e0 = dVertErrors[triIdx.x];
    e1 = dVertErrors[triIdx.y];
    e2 = dVertErrors[triIdx.z];
  } else {
    const float EPS_BOUND = 5.9604645e-8f;
    e0 = l0 * EPS_BOUND;
    e1 = l1 * EPS_BOUND;
    e2 = l2 * EPS_BOUND;
  }
  float maxVertexError = fmaxf(fmaxf(e0, e1), e2);

  const float eps_mach = 1.1920929e-7f;
  if(isTranslated) {
    maxVertexError += (L * eps_mach);
  }
  if(isRotated) {
    maxVertexError += (3.0f * eps_mach * L);
  }

  L += maxVertexError;
  E += (2.0f * maxVertexError);

  dMetrics[idx] = make_float2(L, E * E);
}

// -----------------------------------------------------------------------------
// 2. double3 Unified Device Implementation (Supports Float-only or Dual Double+Float)
// -----------------------------------------------------------------------------
__device__ __forceinline__ void assembleTrianglesDeviceImpl(TriangleDouble* dMeshDouble,
                                                            cuBQL::Triangle* dMesh,
                                                            float2* dMetrics,
                                                            const double3* dVerts,
                                                            const uint3* dIndices,
                                                            int numTriangles) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if(idx >= numTriangles)
    return;

  uint3 triIdx = dIndices[idx];

  double3 p0 = dVerts[triIdx.x];
  double3 p1 = dVerts[triIdx.y];
  double3 p2 = dVerts[triIdx.z];

  // 1. Store high-precision double triangle if buffer is requested
  if(dMeshDouble != nullptr) {
    dMeshDouble[idx].a = p0;
    dMeshDouble[idx].b = p1;
    dMeshDouble[idx].c = p2;
  }

  // 2. Downcast double3 -> float3 components for cuBQL mesh
  float3 a = make_float3(static_cast<float>(p0.x), static_cast<float>(p0.y), static_cast<float>(p0.z));
  float3 b = make_float3(static_cast<float>(p1.x), static_cast<float>(p1.y), static_cast<float>(p1.z));
  float3 c = make_float3(static_cast<float>(p2.x), static_cast<float>(p2.y), static_cast<float>(p2.z));

  dMesh[idx].a.x = a.x;
  dMesh[idx].a.y = a.y;
  dMesh[idx].a.z = a.z;
  dMesh[idx].b.x = b.x;
  dMesh[idx].b.y = b.y;
  dMesh[idx].b.z = b.z;
  dMesh[idx].c.x = c.x;
  dMesh[idx].c.y = c.y;
  dMesh[idx].c.z = c.z;

  // 3. Compute scale (L) and edge extent (E) in double precision
  // double l0 = fmax(fmax(fabs(p0.x), fabs(p0.y)), fabs(p0.z));
  // double l1 = fmax(fmax(fabs(p1.x), fabs(p1.y)), fabs(p1.z));
  // double l2 = fmax(fmax(fabs(p2.x), fabs(p2.y)), fabs(p2.z));
  // double L = fmax(fmax(l0, l1), l2);

  // double ex = fmax(fmax(fabs(p0.x - p1.x), fabs(p1.x - p2.x)), fabs(p2.x - p0.x));
  // double ey = fmax(fmax(fabs(p0.y - p1.y), fabs(p1.y - p2.y)), fabs(p2.y - p0.y));
  // double ez = fmax(fmax(fabs(p0.z - p1.z), fabs(p1.z - p2.z)), fabs(p2.z - p0.z));
  // double E = fmax(fmax(ex, ey), ez);

  // // 4. Compute exact quantization drift per vertex
  // double err0 = fmax(fmax(fabs(p0.x - static_cast<double>(a.x)), fabs(p0.y - static_cast<double>(a.y))),
  //                    fabs(p0.z - static_cast<double>(a.z)));
  // double err1 = fmax(fmax(fabs(p1.x - static_cast<double>(b.x)), fabs(p1.y - static_cast<double>(b.y))),
  //                    fabs(p1.z - static_cast<double>(b.z)));
  // double err2 = fmax(fmax(fabs(p2.x - static_cast<double>(c.x)), fabs(p2.y - static_cast<double>(c.y))),
  //                    fabs(p2.z - static_cast<double>(c.z)));

  // float maxVertexError = static_cast<float>(fmax(fmax(err0, err1), err2));

  // // 5. Convert maxVertexError into effective spatial scale L_eff
  // // FIX: Add maxVertexError directly to L to bound physical coordinate expansion
  // float L_eff = static_cast<float>(L) + maxVertexError;

  // // 6. Write out packed metrics
  // dMetrics[idx] = make_float2(L_eff, static_cast<float>(E * E));
}

// =============================================================================
// Global Entry-Point Kernels
// =============================================================================

__global__ void generateBoxesTrisKernel(cuBQL::box3f* __restrict__ dBoxes,
                                        const double3* __restrict__ dVerts,
                                        const uint3* __restrict__ dIndices,
                                        int numTriangles) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= numTriangles) return;

  uint3 triIdx = dIndices[idx];

  // Direct read: nvcc automatically uses read-only cache via const __restrict__
  double3 p0 = dVerts[triIdx.x];
  double3 p1 = dVerts[triIdx.y];
  double3 p2 = dVerts[triIdx.z];

  // Downcast to float for BVH bounding box
  float3 a = make_float3(static_cast<float>(p0.x), static_cast<float>(p0.y), static_cast<float>(p0.z));
  float3 b = make_float3(static_cast<float>(p1.x), static_cast<float>(p1.y), static_cast<float>(p1.z));
  float3 c = make_float3(static_cast<float>(p2.x), static_cast<float>(p2.y), static_cast<float>(p2.z));

  // Compute AABB
  cuBQL::box3f box;
  box.lower = cuBQL::vec3f(fminf(fminf(a.x, b.x), c.x),
                           fminf(fminf(a.y, b.y), c.y),
                           fminf(fminf(a.z, b.z), c.z));

  box.upper = cuBQL::vec3f(fmaxf(fmaxf(a.x, b.x), c.x),
                           fmaxf(fmaxf(a.y, b.y), c.y),
                           fmaxf(fmaxf(a.z, b.z), c.z));

  dBoxes[idx] = box;
}


// =============================================================================
// Global Entry-Point Kernels
// =============================================================================

// --- float3 Kernels ---
__global__ void assembleTrianglesKernel(cuBQL::Triangle* dMesh,
                                        float2* dMetrics,
                                        const float3* dVerts,
                                        const float* dVertErrors,
                                        const uint3* dIndices,
                                        int numTriangles) {
  assembleTrianglesDeviceImpl(dMesh, dMetrics, dVerts, dVertErrors, dIndices, numTriangles, false, false);
}

__global__ void assembleTrianglesKernelTransformed(cuBQL::Triangle* dMesh,
                                                   float2* dMetrics,
                                                   const float3* dVerts,
                                                   const float* dVertErrors,
                                                   const uint3* dIndices,
                                                   int numTriangles,
                                                   bool isTranslated,
                                                   bool isRotated) {
  assembleTrianglesDeviceImpl(dMesh, dMetrics, dVerts, dVertErrors, dIndices, numTriangles, isTranslated, isRotated);
}

// --- double3 Kernels ---
__global__ void assembleTrianglesKernel(TriangleDouble* dMeshDouble,
                                        cuBQL::Triangle* dMesh,
                                        float2* dMetrics,
                                        const double3* dVerts,
                                        const uint3* dIndices,
                                        int numTriangles) {
  assembleTrianglesDeviceImpl(dMeshDouble, dMesh, dMetrics, dVerts, dIndices, numTriangles);
}

__global__ void assembleTrianglesKernel(
    cuBQL::Triangle* dMesh, float2* dMetrics, const double3* dVerts, const uint3* dIndices, int numTriangles) {
  assembleTrianglesDeviceImpl(nullptr, dMesh, dMetrics, dVerts, dIndices, numTriangles);
}

// --- General Utility Kernels ---
__global__ void generateBoxes(cuBQL::box3f* boxes, const cuBQL::Triangle* tris, int N) {
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if(i < N) {
    boxes[i] = tris[i].bounds();
  }
}

__global__ void
populateReverseMapBKernel(uint32_t* d_reverseMapB, const uint32_t* d_markedNodeIndicesB, uint32_t h_outMarkedCountB) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if(idx < h_outMarkedCountB) {
    uint32_t directBvhNodeId = d_markedNodeIndicesB[idx];
    d_reverseMapB[directBvhNodeId] = idx;
  }
}

// =============================================================================
// Host Launch Wrappers
// =============================================================================

// --- float3 Launch Wrappers ---
void launchAssembleTriangles(cuBQL::Triangle* dMesh,
                             float2* dMetrics,
                             const float3* dVerts,
                             const float* dVertErrors,
                             const uint3* dIndices,
                             int numTriangles,
                             cudaStream_t stream) {
  int blockSize = 256;
  int gridSize = (numTriangles + blockSize - 1) / blockSize;
  assembleTrianglesKernel<<<gridSize, blockSize, 0, stream>>>(dMesh, dMetrics, dVerts, dVertErrors, dIndices,
                                                              numTriangles);
}

void launchAssembleTrianglesTransformed(cuBQL::Triangle* dMesh,
                                        float2* dMetrics,
                                        const float3* dVerts,
                                        const float* dVertErrors,
                                        const uint3* dIndices,
                                        int numTriangles,
                                        bool isTranslated,
                                        bool isRotated,
                                        cudaStream_t stream) {
  int blockSize = 256;
  int gridSize = (numTriangles + blockSize - 1) / blockSize;
  assembleTrianglesKernelTransformed<<<gridSize, blockSize, 0, stream>>>(dMesh, dMetrics, dVerts, dVertErrors, dIndices,
                                                                         numTriangles, isTranslated, isRotated);
}

// --- double3 Launch Wrappers ---
void launchAssembleTriangles(TriangleDouble* dMeshDouble,
                             cuBQL::Triangle* dMesh,
                             float2* dMetrics,
                             const double3* dVerts,
                             const uint3* dIndices,
                             int numTriangles,
                             cudaStream_t stream) {
  int blockSize = 256;
  int gridSize = (numTriangles + blockSize - 1) / blockSize;
  assembleTrianglesKernel<<<gridSize, blockSize, 0, stream>>>(dMeshDouble, dMesh, dMetrics, dVerts, dIndices,
                                                              numTriangles);
}

void launchAssembleTriangles(cuBQL::Triangle* dMesh,
                             float2* dMetrics,
                             const double3* dVerts,
                             const uint3* dIndices,
                             int numTriangles,
                             cudaStream_t stream) {
  launchAssembleTriangles(nullptr, dMesh, dMetrics, dVerts, dIndices, numTriangles, stream);
}

void launchGenerateBoxes(cuBQL::box3f* dBoxes, const cuBQL::Triangle* dTris, int numTriangles, cudaStream_t stream) {
  int blockSize = 256;
  int gridSize = (numTriangles + blockSize - 1) / blockSize;
  generateBoxes<<<gridSize, blockSize, 0, stream>>>(dBoxes, dTris, numTriangles);
}


// =============================================================================
// Host Launch Wrappers
// =============================================================================

void launchGenerateBoxesTris(cuBQL::box3f* dBoxes,
                             const double3* dVerts,
                             const uint3* dIndices,
                             int numTriangles,
                             cudaStream_t stream) {
  if (numTriangles <= 0) return;
  int blockSize = 256;
  int gridSize = (numTriangles + blockSize - 1) / blockSize;
  generateBoxesTrisKernel<<<gridSize, blockSize, 0, stream>>>(dBoxes, dVerts, dIndices, numTriangles);
}