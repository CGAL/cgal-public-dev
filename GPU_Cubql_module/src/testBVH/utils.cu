#include "utils.h"
#include "cuBQL/math/vec.h" // Ensures divRoundUp or equivalent CUDA math helpers are present

__global__ void assembleTrianglesKernel(cuBQL::Triangle* dMesh,
                                        float2* dMetrics,
                                        const float3* dVerts,
                                        const float* dVertErrors,
                                        const uint3* dIndices,
                                        int numTriangles) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= numTriangles)
    return;

  uint3 triIdx = dIndices[idx];

  float3 p0 = dVerts[triIdx.x];
  float3 p1 = dVerts[triIdx.y];
  float3 p2 = dVerts[triIdx.z];

  // 1. Map components directly into cuBQL fields
  dMesh[idx].a.x = p0.x; dMesh[idx].a.y = p0.y; dMesh[idx].a.z = p0.z;
  dMesh[idx].b.x = p1.x; dMesh[idx].b.y = p1.y; dMesh[idx].b.z = p1.z;
  dMesh[idx].c.x = p2.x; dMesh[idx].c.y = p2.y; dMesh[idx].c.z = p2.z;

  // 2. Compute L (Absolute Spatial Scale)
  float l0 = fmaxf(fmaxf(fabsf(p0.x), fabsf(p0.y)), fabsf(p0.z));
  float l1 = fmaxf(fmaxf(fabsf(p1.x), fabsf(p1.y)), fabsf(p1.z));
  float l2 = fmaxf(fmaxf(fabsf(p2.x), fabsf(p2.y)), fabsf(p2.z));
  float L  = fmaxf(fmaxf(l0, l1), l2);

  // 3. Compute E (Maximum Edge Axis Extent)
  float ex = fmaxf(fmaxf(fabsf(p0.x - p1.x), fabsf(p1.x - p2.x)), fabsf(p2.x - p0.x));
  float ey = fmaxf(fmaxf(fabsf(p0.y - p1.y), fabsf(p1.y - p2.y)), fabsf(p2.y - p0.y));
  float ez = fmaxf(fmaxf(fabsf(p0.z - p1.z), fabsf(p1.z - p2.z)), fabsf(p2.z - p0.z));
  float E  = fmaxf(fmaxf(ex, ey), ez);

  // 4. Absorb Host-Generated Rounding Error Bounds
  float e0, e1, e2;
  if (dVertErrors != nullptr) {
    e0 = dVertErrors[triIdx.x];
    e1 = dVertErrors[triIdx.y];
    e2 = dVertErrors[triIdx.z];
  } else {
    const float EPS_BOUND = 5.9604645e-8f; // 2^-24 upper bound
    e0 = l0 * EPS_BOUND;
    e1 = l1 * EPS_BOUND;
    e2 = l2 * EPS_BOUND;
  }
  float maxVertexError = fmaxf(fmaxf(e0, e1), e2);

  L += maxVertexError;
  E += (2.0f * maxVertexError);

  // 5. Write out packed parameters
  dMetrics[idx] = make_float2(L, E * E);
}

__global__ void generateBoxes(cuBQL::box3f* boxes, const cuBQL::Triangle* tris, int N) {
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < N) {
    boxes[i] = tris[i].bounds();
  }
}

__global__ void populateReverseMapBKernel(uint32_t* d_reverseMapB, 
                                           const uint32_t* d_markedNodeIndicesB, 
                                           uint32_t h_outMarkedCountB) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx < h_outMarkedCountB) {
    uint32_t directBvhNodeId = d_markedNodeIndicesB[idx];
    d_reverseMapB[directBvhNodeId] = idx;
  }
}

// Host Launch Wrappers
void launchAssembleTriangles(cuBQL::Triangle* dMesh,
                            float2* dMetrics,
                            const float3* dVerts,
                            const float* dVertErrors,
                            const uint3* dIndices,
                            int numTriangles,
                            cudaStream_t stream) {
  int blockSize = 256;
  int gridSize = (numTriangles + blockSize - 1) / blockSize;
  assembleTrianglesKernel<<<gridSize, blockSize, 0, stream>>>(
      dMesh, dMetrics, dVerts, dVertErrors, dIndices, numTriangles);
}

void launchGenerateBoxes(cuBQL::box3f* dBoxes,
                         const cuBQL::Triangle* dTris,
                         int numTriangles,
                         cudaStream_t stream) {
  int blockSize = 256;
  int gridSize = (numTriangles + blockSize - 1) / blockSize;
  generateBoxes<<<gridSize, blockSize, 0, stream>>>(dBoxes, dTris, numTriangles);
}