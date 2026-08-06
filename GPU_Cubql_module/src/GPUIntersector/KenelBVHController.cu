#define CUBQL_GPU_BUILDER_IMPLEMENTATION 1

#include <cuda_runtime.h>
#include <cmath>
#include <chrono>

#include "KernelBVHController.h"

#include <tbb/parallel_for.h>
#include <tbb/parallel_invoke.h>

// Core cuBQL headers
#include "cuBQL/builder/cuda.h"
#include "include/third-party/cubql/sm_builder_v2_2.h"
#include "include/third-party/cubql/refit_forest.h"
#include "cuBQL/bvh.h"

#include "include/third-party/cubql/fixedBoxQueryv2.h"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/scan.h>
#include <thrust/execution_policy.h>
#include <thrust/system/cuda/execution_policy.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/unique.h>

#include <vector>
#include <algorithm>
#include <iostream>
#include "samples/common/loadOBJ.h"

// Custom modules
#include "../testBVH/utils.h"
#include "../testBVH/crossCheck.h"
#include "../testBVH/DualTreeStep.h"
#include "../testBVH/rapidDescendKernel.h"
// #include "../testBVH/batchedCrossIntersection.h"
//#include "../testBVH/batchedCrossIntersectionV2.h"
#include "../testBVH/batchedCrossIntersectionDouble.h"
//#include "../testBVH/batchedCrossIntersectionDoubleOld.h"
//#include "../testBVH/batchedCrossIntersectionV3.h"
#include "../custom_pipeline/TriangleDouble.h"


#ifndef _ALLOC
#define _ALLOC(ptr, count, stream, memResource)                                                                        \
  CUBQL_CUDA_CALL(MallocAsync((void**)&(ptr), (size_t)(count) * sizeof(*(ptr)), stream))
#endif

#ifndef _FREE
#define _FREE(ptr, stream, memResource) CUBQL_CUDA_CALL(FreeAsync((ptr), stream))
#endif

// -----------------------------------------------------------------------------
// Helper 3x3 Matrix Structure & Rigid Body Device Transformations (Double Precision)
// -----------------------------------------------------------------------------
// struct Mat3x3
// {
//   double m[3][3];
// };

// inline Mat3x3 makeRotationMatrixDeg(double rxDeg, double ryDeg, double rzDeg) {
//   const double DEG_TO_RAD = 3.14159265358979323846 / 180.0;
//   double radX = rxDeg * DEG_TO_RAD;
//   double radY = ryDeg * DEG_TO_RAD;
//   double radZ = rzDeg * DEG_TO_RAD;

//   double cx = cos(radX), sx = sin(radX);
//   double cy = cos(radY), sy = sin(radY);
//   double cz = cos(radZ), sz = sin(radZ);

//   Mat3x3 R;
//   // Composite Rotation: R = Rz * Ry * Rx
//   R.m[0][0] = cy * cz;
//   R.m[0][1] = sx * sy * cz - cx * sz;
//   R.m[0][2] = cx * sy * cz + sx * sz;

//   R.m[1][0] = cy * sz;
//   R.m[1][1] = sx * sy * sz + cx * cz;
//   R.m[1][2] = cx * sy * sz - sx * cz;

//   R.m[2][0] = -sy;
//   R.m[2][1] = sx * cy;
//   R.m[2][2] = cx * cy;

//   return R;
// }

// __device__ inline double3 transformPoint(const double3& p, const Mat3x3& R, double3 center, double3 trans) {
//   double x = p.x - center.x;
//   double y = p.y - center.y;
//   double z = p.z - center.z;

//   double rx = R.m[0][0] * x + R.m[0][1] * y + R.m[0][2] * z;
//   double ry = R.m[1][0] * x + R.m[1][1] * y + R.m[1][2] * z;
//   double rz = R.m[2][0] * x + R.m[2][1] * y + R.m[2][2] * z;

//   return make_double3(rx + center.x + trans.x, ry + center.y + trans.y, rz + center.z + trans.z);
// }

// __global__ void transformVerticesKernel(
//     double3* dVertsOut, const double3* dVertsOrig, int numVerts, Mat3x3 R, double3 center, double3 trans) {
//   int idx = threadIdx.x + blockIdx.x * blockDim.x;
//   if(idx >= numVerts)
//     return;

//   dVertsOut[idx] = transformPoint(dVertsOrig[idx], R, center, trans);
// }

// __global__ void generateBoxesTrisKernel(cuBQL::box3f* __restrict__ dBoxes,
//                                      const double3* __restrict__ dVerts,
//                                      const uint3* __restrict__ dIndices,
//                                      int numTriangles) {
//   int idx = threadIdx.x + blockIdx.x * blockDim.x;
//   if (idx >= numTriangles) return;

//   uint3 triIdx = dIndices[idx];

//   // Direct read: nvcc automatically uses read-only cache via const __restrict__
//   double3 p0 = dVerts[triIdx.x];
//   double3 p1 = dVerts[triIdx.y];
//   double3 p2 = dVerts[triIdx.z];

//   // Downcast to float for BVH bounding box
//   float3 a = make_float3(static_cast<float>(p0.x), static_cast<float>(p0.y), static_cast<float>(p0.z));
//   float3 b = make_float3(static_cast<float>(p1.x), static_cast<float>(p1.y), static_cast<float>(p1.z));
//   float3 c = make_float3(static_cast<float>(p2.x), static_cast<float>(p2.y), static_cast<float>(p2.z));

//   // Compute AABB
//   cuBQL::box3f box;
//   box.lower = cuBQL::vec3f(fminf(fminf(a.x, b.x), c.x),
//                            fminf(fminf(a.y, b.y), c.y),
//                            fminf(fminf(a.z, b.z), c.z));

//   box.upper = cuBQL::vec3f(fmaxf(fmaxf(a.x, b.x), c.x),
//                            fmaxf(fmaxf(a.y, b.y), c.y),
//                            fmaxf(fmaxf(a.z, b.z), c.z));

//   dBoxes[idx] = box;
// }

// __global__ void generateBoxesAndTrisKernel(cuBQL::box3f* dBoxes,
//                                            TriangleDouble* dMeshDouble,
//                                            cuBQL::Triangle* dMesh,
//                                            const double3* dVerts,
//                                            const uint3* dIndices,
//                                            int numTriangles) {
//   int idx = threadIdx.x + blockIdx.x * blockDim.x;
//   if(idx >= numTriangles)
//     return;

//   uint3 triIdx = dIndices[idx];

//   double3 p0 = dVerts[triIdx.x];
//   double3 p1 = dVerts[triIdx.y];
//   double3 p2 = dVerts[triIdx.z];

//   // 1. Store high-precision double triangle (if requested)
//   if(dMeshDouble != nullptr) {
//     dMeshDouble[idx].a = p0;
//     dMeshDouble[idx].b = p1;
//     dMeshDouble[idx].c = p2;
//   }

//   // 2. Downcast to float3
//   float3 a = make_float3(static_cast<float>(p0.x), static_cast<float>(p0.y), static_cast<float>(p0.z));
//   float3 b = make_float3(static_cast<float>(p1.x), static_cast<float>(p1.y), static_cast<float>(p1.z));
//   float3 c = make_float3(static_cast<float>(p2.x), static_cast<float>(p2.y), static_cast<float>(p2.z));

//   // 3. Store cuBQL float triangle
//   dMesh[idx].a = cuBQL::vec3f(a.x, a.y, a.z);
//   dMesh[idx].b = cuBQL::vec3f(b.x, b.y, b.z);
//   dMesh[idx].c = cuBQL::vec3f(c.x, c.y, c.z);

//   // 4. Compute & write AABB directly from registers
//   cuBQL::box3f box;
//   box.lower = cuBQL::vec3f(fminf(fminf(a.x, b.x), c.x), fminf(fminf(a.y, b.y), c.y), fminf(fminf(a.z, b.z), c.z));

//   box.upper = cuBQL::vec3f(fmaxf(fmaxf(a.x, b.x), c.x), fmaxf(fmaxf(a.y, b.y), c.y), fmaxf(fmaxf(a.z, b.z), c.z));

//   dBoxes[idx] = box;
// }

// -----------------------------------------------------------------------------
// Parallel TBB CPU CGAL Transformations & Extraction Helpers
// -----------------------------------------------------------------------------
static void
transformCgalMesh(Mesh* mesh, const std::vector<Point3>& origPoints, Point3 center, double3 rotDeg, double3 trans) {
  if(!mesh || origPoints.empty())
    return;

  auto pmap = mesh->points();
  const size_t numPoints = mesh->num_vertices();

  const double DEG_TO_RAD = 3.14159265358979323846 / 180.0;
  double radX = rotDeg.x * DEG_TO_RAD;
  double radY = rotDeg.y * DEG_TO_RAD;
  double radZ = rotDeg.z * DEG_TO_RAD;

  double cx = cos(radX), sx = sin(radX);
  double cy = cos(radY), sy = sin(radY);
  double cz = cos(radZ), sz = sin(radZ);

  double r00 = cy * cz, r01 = sx * sy * cz - cx * sz, r02 = cx * sy * cz + sx * sz;
  double r10 = cy * sz, r11 = sx * sy * sz + cx * cz, r12 = cx * sy * sz - sx * cz;
  double r20 = -sy, r21 = sx * cy, r22 = cx * cy;

  double cx_p = center.x(), cy_p = center.y(), cz_p = center.z();
  double tx = trans.x;
  double ty = trans.y;
  double tz = trans.z;

  tbb::parallel_for(tbb::blocked_range<size_t>(0, numPoints),
                    [mesh, &pmap, &origPoints, r00, r01, r02, r10, r11, r12, r20, r21, r22, cx_p, cy_p, cz_p, tx, ty,
                     tz](const tbb::blocked_range<size_t>& range) {
                      for(size_t i = range.begin(); i != range.end(); ++i) {
                        Mesh::Vertex_index vd(static_cast<uint32_t>(i));
                        if(!mesh->is_removed(vd)) {
                          const Point3& pOrig = origPoints[i];
                          double x = pOrig.x() - cx_p;
                          double y = pOrig.y() - cy_p;
                          double z = pOrig.z() - cz_p;

                          double rx = r00 * x + r01 * y + r02 * z;
                          double ry = r10 * x + r11 * y + r12 * z;
                          double rz = r20 * x + r21 * y + r22 * z;

                          pmap[vd] = Point3(rx + cx_p + tx, ry + cy_p + ty, rz + cz_p + tz);
                        }
                      }
                    });
}

static void translateCgalAndExtractDouble3(
    Mesh* mesh, const std::vector<Point3>& origPoints, double3* hOutVerts, double xB, double yB, double zB) {
  if(!mesh || origPoints.empty() || !hOutVerts)
    return;

  auto pmap = mesh->points();
  const size_t numPoints = mesh->num_vertices();

  tbb::parallel_for(tbb::blocked_range<size_t>(0, numPoints),
                    [mesh, &pmap, &origPoints, hOutVerts, xB, yB, zB](const tbb::blocked_range<size_t>& range) {
                      for(size_t i = range.begin(); i != range.end(); ++i) {
                        Mesh::Vertex_index vd(static_cast<uint32_t>(i));
                        if(!mesh->is_removed(vd)) {
                          // 1. High-precision double translation from pristine baseline
                          const Point3& pOrig = origPoints[i];
                          Point3 pTranslated(pOrig.x() + xB, pOrig.y() + yB, pOrig.z() + zB);

                          // Update CGAL host mesh
                          pmap[vd] = pTranslated;

                          // 2. Extract double3 for GPU upload
                          hOutVerts[i] = make_double3(pTranslated.x(), pTranslated.y(), pTranslated.z());
                        }
                      }
                    });
}

// -----------------------------------------------------------------------------
// Controller Constructor & Destructor
// -----------------------------------------------------------------------------
KernelBVHController::KernelBVHController() {
  CUBQL_CUDA_CALL(StreamCreate(&m_stream));
}

KernelBVHController::~KernelBVHController() {
  cleanup();
  if(m_stream) {
    cudaStreamDestroy(m_stream);
    m_stream = nullptr;
  }
}

// -----------------------------------------------------------------------------
// Dynamic Rigid Body Transformation Engine
// -----------------------------------------------------------------------------
void KernelBVHController::setTransformBoth(double3 rotDegA,
                                           double3 transA,
                                           double3 rotDegB,
                                           double3 transB,
                                           float& timeGPU,
                                           float& timeCPU,
                                           float& timeTransformVerts,
                                           float& timeAssembleTris,
                                           float& timeGenBoxes) {
  timeGPU = 0.0f;
  timeCPU = 0.0f;
  timeTransformVerts = 0.0f;
  timeAssembleTris = 0.0f;
  timeGenBoxes = 0.0f;

  bool movedA = (rotDegA.x != m_rotA.x || rotDegA.y != m_rotA.y || rotDegA.z != m_rotA.z || transA.x != m_transA.x ||
                 transA.y != m_transA.y || transA.z != m_transA.z);

  bool movedB = (rotDegB.x != m_rotB.x || rotDegB.y != m_rotB.y || rotDegB.z != m_rotB.z || transB.x != m_transB.x ||
                 transB.y != m_transB.y || transB.z != m_transB.z);

  if(!movedA && !movedB)
    return;

  int block = 256;
  bool runA = (movedA && m_numVertsA > 0 && m_dVertsA && m_dVertsAOrig);
  bool runB = (movedB && m_numVertsB > 0 && m_dVertsB && m_dVertsBOrig);
  bool gpuWorkQueued = runA || runB;

  cudaEvent_t evStartGPU, evStopGPU;
  cudaEvent_t evStartTV, evStopTV;
  cudaEvent_t evStartAssemble, evStopAssemble;
  cudaEvent_t evStartRefit, evStopRefit;

  if(gpuWorkQueued) {
    cudaEventCreate(&evStartGPU);
    cudaEventCreate(&evStopGPU);
    cudaEventCreate(&evStartTV);
    cudaEventCreate(&evStopTV);
    cudaEventCreate(&evStartAssemble);
    cudaEventCreate(&evStopAssemble);
    cudaEventCreate(&evStartRefit);
    cudaEventCreate(&evStopRefit);

    cudaEventRecord(evStartGPU, m_stream);

    // --- STAGE 1: Transform Vertices (A + B) ---
    cudaEventRecord(evStartTV, m_stream);
    if(runA) {
      Mat3x3 RA = makeRotationMatrixDeg(rotDegA.x, rotDegA.y, rotDegA.z);
      double3 cA = make_double3(m_centerA.x(), m_centerA.y(), m_centerA.z());
      launchTransformVertices(m_dVertsA, m_dVertsAOrig, m_numVertsA, RA, cA, transA, m_stream);
    }
    if(runB) {
      Mat3x3 RB = makeRotationMatrixDeg(rotDegB.x, rotDegB.y, rotDegB.z);
      double3 cB = make_double3(m_centerB.x(), m_centerB.y(), m_centerB.z());
      launchTransformVertices(m_dVertsB, m_dVertsBOrig, m_numVertsB, RB, cB, transB, m_stream);
    }
    cudaEventRecord(evStopTV, m_stream);

    // --- STAGE 2: Fused Box & Triangle Assembly Kernel (A + B) ---

    cudaEventRecord(evStartAssemble, m_stream);
    // AFTER:
    if(runA) {
      launchGenerateBoxesTris(m_dBoxesA, m_dVertsA, m_dIndicesA, m_numTrianglesA, m_stream);
    }
    if(runB) {
      launchGenerateBoxesTris(m_dBoxesB, m_dVertsB, m_dIndicesB, m_numTrianglesB, m_stream);
    }
    cudaEventRecord(evStopAssemble, m_stream);

    // --- STAGE 3: BVH Tree Refit (A + B) ---
    cudaEventRecord(evStartRefit, m_stream);
    if(runA) {
      cuBQL::cuda::refit(m_bvhA, m_dBoxesA, m_stream);
    }
    if(runB) {
      cuBQL::cuda::refit(m_bvhB, m_dBoxesB, m_stream);
    }
    cudaEventRecord(evStopRefit, m_stream);

    cudaEventRecord(evStopGPU, m_stream);
  }

  // --- STEP 2: CONCURRENT CPU CGAL TRANSFORMS ---
  // auto cpuStart = std::chrono::high_resolution_clock::now();

  // tbb::parallel_invoke(
  //     [this, movedA, rotDegA, transA]() {
  //       if(movedA && m_meshAcpu) {
  //         transformCgalMesh(m_meshAcpu, m_origPointsA, m_centerA, rotDegA, transA);
  //       }
  //     },
  //     [this, movedB, rotDegB, transB]() {
  //       if(movedB && m_meshBcpu) {
  //         transformCgalMesh(m_meshBcpu, m_origPointsB, m_centerB, rotDegB, transB);
  //       }
  //     });

  // auto cpuEnd = std::chrono::high_resolution_clock::now();
  // timeCPU = std::chrono::duration<float, std::milli>(cpuEnd - cpuStart).count();

  timeCPU = 0;

  if(movedA) {
    m_rotA = rotDegA;
    m_transA = transA;
  }
  if(movedB) {
    m_rotB = rotDegB;
    m_transB = transB;
    shiftX = transB.x;
    shiftY = transB.y;
    shiftZ = transB.z;
  }

  // --- STEP 3: SYNCHRONIZE & TIMINGS ---
  if(gpuWorkQueued) {
    cudaEventSynchronize(evStopGPU);

    cudaEventElapsedTime(&timeGPU, evStartGPU, evStopGPU);
    cudaEventElapsedTime(&timeTransformVerts, evStartTV, evStopTV);
    cudaEventElapsedTime(&timeAssembleTris, evStartAssemble, evStopAssemble); // generateBoxesAndTrisKernel
    cudaEventElapsedTime(&timeGenBoxes, evStartRefit, evStopRefit);           // cuBQL::cuda::refit

    cudaEventDestroy(evStartGPU);
    cudaEventDestroy(evStopGPU);
    cudaEventDestroy(evStartTV);
    cudaEventDestroy(evStopTV);
    cudaEventDestroy(evStartAssemble);
    cudaEventDestroy(evStopAssemble);
    cudaEventDestroy(evStartRefit);
    cudaEventDestroy(evStopRefit);
  }
}
// Legacy Translation Handlers
void KernelBVHController::setTranslation(double xB, double yB, double zB) {
  float tGPU, tCPU, tTV, tAT, tGB;
  setTransformBoth(m_rotA, m_transA, m_rotB, make_double3(xB, yB, zB), tGPU, tCPU, tTV, tAT, tGB);
}


// -----------------------------------------------------------------------------
// Controller Construction Method
// -----------------------------------------------------------------------------
void KernelBVHController::construct(Mesh& meshAcpu,
                                    Mesh& meshBcpu,
                                    const double3* hVertsA,
                                    int numVertsA,
                                    const uint3* hIndicesA,
                                    int numTrianglesA,
                                    int levelA,
                                    const double3* hVertsB,
                                    int numVertsB,
                                    const uint3* hIndicesB,
                                    int numTrianglesB,
                                    int levelB,
                                    int leafThreshold,
                                    ExecutionStats& stats) {
  construct(meshAcpu, meshBcpu, Point3(0, 0, 0), Point3(0, 0, 0), hVertsA, numVertsA, hIndicesA, numTrianglesA, levelA,
            hVertsB, numVertsB, hIndicesB, numTrianglesB, levelB, leafThreshold, stats);
}

void KernelBVHController::construct(Mesh& meshAcpu,
                                    Mesh& meshBcpu,
                                    const Point3& centerA,
                                    const Point3& centerB,
                                    const double3* hVertsA,
                                    int numVertsA,
                                    const uint3* hIndicesA,
                                    int numTrianglesA,
                                    int levelA,
                                    const double3* hVertsB,
                                    int numVertsB,
                                    const uint3* hIndicesB,
                                    int numTrianglesB,
                                    int levelB,
                                    int leafThreshold,
                                    ExecutionStats& stats) {
  if(numTrianglesA <= 0 || numTrianglesB <= 0)
    return;

  cleanup();

  m_meshAcpu = &meshAcpu;
  m_meshBcpu = &meshBcpu;
  m_numTrianglesA = numTrianglesA;
  m_numVertsA = numVertsA;
  m_levelA = levelA;
  m_numTrianglesB = numTrianglesB;
  m_numVertsB = numVertsB;
  m_levelB = levelB;
  m_leafThreshold = leafThreshold;

  m_hVertsA = hVertsA;
  m_hIndicesA = hIndicesA;
  m_hVertsB = hVertsB;
  m_hIndicesB = hIndicesB;

  m_centerA = centerA;
  m_centerB = centerB;

  m_rotA = make_double3(0.0, 0.0, 0.0);
  m_transA = make_double3(0.0, 0.0, 0.0);
  m_rotB = make_double3(0.0, 0.0, 0.0);
  m_transB = make_double3(0.0, 0.0, 0.0);
  shiftX = 0.0;
  shiftY = 0.0;
  shiftZ = 0.0;

  cuBQL::BuildConfig buildConfig;
  buildConfig.makeLeafThreshold = m_leafThreshold;

  cudaEvent_t evAllocStart, evAllocStop;
  cudaEvent_t evBuildAStart, evBuildAStop;
  cudaEvent_t evBuildBStart, evBuildBStop;

  CUBQL_CUDA_CALL(EventCreate(&evAllocStart));
  CUBQL_CUDA_CALL(EventCreate(&evAllocStop));
  CUBQL_CUDA_CALL(EventCreate(&evBuildAStart));
  CUBQL_CUDA_CALL(EventCreate(&evBuildAStop));
  CUBQL_CUDA_CALL(EventCreate(&evBuildBStart));
  CUBQL_CUDA_CALL(EventCreate(&evBuildBStop));

  CUBQL_CUDA_CALL(EventRecord(evAllocStart, m_stream));

  // --- PERSISTENT DEVICE ALLOCATION ---
  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dVertsA, m_numVertsA * sizeof(double3), m_stream));
  CUBQL_CUDA_CALL(MemcpyAsync(m_dVertsA, hVertsA, m_numVertsA * sizeof(double3), cudaMemcpyHostToDevice, m_stream));

  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dVertsAOrig, m_numVertsA * sizeof(double3), m_stream));
  CUBQL_CUDA_CALL(MemcpyAsync(m_dVertsAOrig, hVertsA, m_numVertsA * sizeof(double3), cudaMemcpyHostToDevice, m_stream));

  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dIndicesA, m_numTrianglesA * sizeof(uint3), m_stream));
  CUBQL_CUDA_CALL(
      MemcpyAsync(m_dIndicesA, hIndicesA, m_numTrianglesA * sizeof(uint3), cudaMemcpyHostToDevice, m_stream));

  m_dVertErrorsA = nullptr;

  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dVertsB, m_numVertsB * sizeof(double3), m_stream));
  CUBQL_CUDA_CALL(MemcpyAsync(m_dVertsB, hVertsB, m_numVertsB * sizeof(double3), cudaMemcpyHostToDevice, m_stream));

  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dVertsBOrig, m_numVertsB * sizeof(double3), m_stream));
  CUBQL_CUDA_CALL(MemcpyAsync(m_dVertsBOrig, hVertsB, m_numVertsB * sizeof(double3), cudaMemcpyHostToDevice, m_stream));

  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dIndicesB, m_numTrianglesB * sizeof(uint3), m_stream));
  CUBQL_CUDA_CALL(
      MemcpyAsync(m_dIndicesB, hIndicesB, m_numTrianglesB * sizeof(uint3), cudaMemcpyHostToDevice, m_stream));

  m_dVertErrorsB = nullptr;

  // BVH Bounding Box Buffers
  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dBoxesA, m_numTrianglesA * sizeof(cuBQL::box3f), m_stream));
  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dBoxesB, m_numTrianglesB * sizeof(cuBQL::box3f), m_stream));

  // Generate initial bounding boxes directly from vertices and indices
  generateBoxesTrisKernel<<<cuBQL::divRoundUp(m_numTrianglesA, 256), 256, 0, m_stream>>>(m_dBoxesA, m_dVertsA,
                                                                                         m_dIndicesA, m_numTrianglesA);

  generateBoxesTrisKernel<<<cuBQL::divRoundUp(m_numTrianglesB, 256), 256, 0, m_stream>>>(m_dBoxesB, m_dVertsB,
                                                                                         m_dIndicesB, m_numTrianglesB);

  CUBQL_CUDA_CALL(EventRecord(evAllocStop, m_stream));

  // // --- CONCURRENT CPU CGAL BASELINE EXTRACTION ---
  // tbb::parallel_invoke(
  //     [this]() {
  //       if(m_meshAcpu) {
  //         size_t numVerts = m_meshAcpu->num_vertices();
  //         m_origPointsA.resize(numVerts);
  //         auto pmap = m_meshAcpu->points();
  //         tbb::parallel_for(tbb::blocked_range<size_t>(0, numVerts),
  //                           [&pmap, this](const tbb::blocked_range<size_t>& range) {
  //                             for(size_t i = range.begin(); i != range.end(); ++i) {
  //                               Mesh::Vertex_index vd(static_cast<uint32_t>(i));
  //                               m_origPointsA[i] = pmap[vd];
  //                             }
  //                           });
  //       }
  //     },
  //     [this]() {
  //       if(m_meshBcpu) {
  //         size_t numVerts = m_meshBcpu->num_vertices();
  //         m_origPointsB.resize(numVerts);
  //         auto pmap = m_meshBcpu->points();
  //         tbb::parallel_for(tbb::blocked_range<size_t>(0, numVerts),
  //                           [&pmap, this](const tbb::blocked_range<size_t>& range) {
  //                             for(size_t i = range.begin(); i != range.end(); ++i) {
  //                               Mesh::Vertex_index vd(static_cast<uint32_t>(i));
  //                               m_origPointsB[i] = pmap[vd];
  //                             }
  //                           });
  //       }
  //     });

  // --- THRUST DEVICE VECTOR INITIALIZATION ---
  double tThrustInitStart = cuBQL::getCurrentTime();

  uint32_t maxPossibleNodesA = 2 * m_numTrianglesA;
  uint32_t maxPossibleNodesB = 2 * m_numTrianglesB;

  m_dMarkedNodeIndicesA_Full.resize(maxPossibleNodesA, 0);
  m_dNodeDescendantCountsA.resize(maxPossibleNodesA, 0);
  m_dMarkedNodeIndicesB_Full.resize(maxPossibleNodesB, 0);
  m_dNodeDescendantCountsB.resize(maxPossibleNodesB, 0);
  m_dReverseMapB.resize(maxPossibleNodesB, 0);

  stats.thrustInitOverheadMs = (cuBQL::getCurrentTime() - tThrustInitStart) * 1000.0;

  // --- BUILD & REFIT BVHs ---
  CUBQL_CUDA_CALL(EventRecord(evBuildAStart, m_stream));
  cuBQL::gpuBuilder_v2_2::build_custom(m_bvhA, m_dBoxesA, m_numTrianglesA, buildConfig, (uint32_t)m_levelA,
                                       thrust::raw_pointer_cast(m_dMarkedNodeIndicesA_Full.data()),
                                       thrust::raw_pointer_cast(m_dNodeDescendantCountsA.data()),
                                       &m_hOutMarkedCountA_Full, m_stream, m_memResource);
  cuBQL::cuda::refit(m_bvhA, m_dBoxesA, m_stream);
  CUBQL_CUDA_CALL(EventRecord(evBuildAStop, m_stream));

  CUBQL_CUDA_CALL(EventRecord(evBuildBStart, m_stream));
  cuBQL::gpuBuilder_v2_2::build_custom(m_bvhB, m_dBoxesB, m_numTrianglesB, buildConfig, (uint32_t)m_levelB,
                                       thrust::raw_pointer_cast(m_dMarkedNodeIndicesB_Full.data()),
                                       thrust::raw_pointer_cast(m_dNodeDescendantCountsB.data()),
                                       &m_hOutMarkedCountB_Full, m_stream, m_memResource);
  cuBQL::cuda::refit(m_bvhB, m_dBoxesB, m_stream);
  CUBQL_CUDA_CALL(EventRecord(evBuildBStop, m_stream));

  CUBQL_CUDA_CALL(StreamSynchronize(m_stream));

  float msAlloc = 0.0f, msBuildA = 0.0f, msBuildB = 0.0f;
  CUBQL_CUDA_CALL(EventElapsedTime(&msAlloc, evAllocStart, evAllocStop));
  CUBQL_CUDA_CALL(EventElapsedTime(&msBuildA, evBuildAStart, evBuildAStop));
  CUBQL_CUDA_CALL(EventElapsedTime(&msBuildB, evBuildBStart, evBuildBStop));

  stats.initialAllocAndCopyMs = static_cast<double>(msAlloc);
  stats.buildRefitMeshAMs = static_cast<double>(msBuildA);
  stats.buildRefitMeshBMs = static_cast<double>(msBuildB);

  CUBQL_CUDA_CALL(EventDestroy(evAllocStart));
  CUBQL_CUDA_CALL(EventDestroy(evAllocStop));
  CUBQL_CUDA_CALL(EventDestroy(evBuildAStart));
  CUBQL_CUDA_CALL(EventDestroy(evBuildAStop));
  CUBQL_CUDA_CALL(EventDestroy(evBuildBStart));
  CUBQL_CUDA_CALL(EventDestroy(evBuildBStop));
}
// -----------------------------------------------------------------------------
// Execution Pipeline Function
// -----------------------------------------------------------------------------
void KernelBVHController::runIntersectionPipeline(int batchMultiplier,
                                                  int mode,
                                                  int activateAsyncDownload,
                                                  int2*& outFinalExactPairs, // Fast raw pointer return
                                                  size_t& outFinalCount,
                                                  ExecutionStats& stats, bool gpuDouble) {
  double tPipelineStart = cuBQL::getCurrentTime();

  m_dOutPairsA.clear();
  m_dOutPairsB.clear();

  uint32_t m_hOutMarkedCountA = m_hOutMarkedCountA_Full;
  uint32_t m_hOutMarkedCountB = m_hOutMarkedCountB_Full;
  thrust::device_vector<uint32_t> m_dMarkedNodeIndicesA = m_dMarkedNodeIndicesA_Full;
  thrust::device_vector<uint32_t> m_dMarkedNodeIndicesB = m_dMarkedNodeIndicesB_Full;

  // 1. CRISS-CROSS INTERSECTION
  double tCrossStart = cuBQL::getCurrentTime();

  uint32_t totalIntersections =
      executeCrissCrossIntersection(m_bvhA, m_dMarkedNodeIndicesA, m_hOutMarkedCountA, m_bvhB, m_dMarkedNodeIndicesB,
                                    m_hOutMarkedCountB, m_dOutPairsA, m_dOutPairsB);

  CUBQL_CUDA_CALL(StreamSynchronize(m_stream));
  double tCrossEnd = cuBQL::getCurrentTime();

  uint64_t totalPossiblePairs = (uint64_t)m_hOutMarkedCountA * m_hOutMarkedCountB;
  double intersectionPercentage =
      totalPossiblePairs > 0 ? ((double)totalIntersections / totalPossiblePairs) * 100.0 : 0.0;

  // 2. DUAL-TREE TRAVERSAL STEP
  double tDualStepStart = cuBQL::getCurrentTime();

  if(mode > 0) {
    executeDualTreeStep(mode, m_levelA, m_levelB, m_dOutPairsA, m_dOutPairsB, m_dMarkedNodeIndicesA,
                        m_dMarkedNodeIndicesB, m_dNodeDescendantCountsA, m_dNodeDescendantCountsB, m_hOutMarkedCountA,
                        m_hOutMarkedCountB, m_bvhA, m_bvhB);

    CUBQL_CUDA_CALL(StreamSynchronize(m_stream));
  }

  double tDualStepEnd = cuBQL::getCurrentTime();
  stats.dualTreeStepMs = (mode > 0) ? (tDualStepEnd - tDualStepStart) * 1000.0 : 0.0;

  int totalBatches = static_cast<int>(m_dOutPairsA.size());
  IntersectionTimeTracker tracker;
  uint64_t finalCandidatePairs = 0;

  // 3. RAPID DESCENT BFS & REVERSE MAP
  double tGpuBfsStart = cuBQL::getCurrentTime();

  executeRapidDescentBFS(m_bvhB, m_numTrianglesB, m_dMarkedNodeIndicesB, m_dNodeDescendantCountsB, m_hOutMarkedCountB,
                         m_dOutOffsetsB, m_dOutPrimsFlatB);

  if(m_hOutMarkedCountB > 0) {
    int blockSize = 256;
    int gridSize = (m_hOutMarkedCountB + blockSize - 1) / blockSize;
    populateReverseMapBKernel<<<gridSize, blockSize, 0, m_stream>>>(
        thrust::raw_pointer_cast(m_dReverseMapB.data()), thrust::raw_pointer_cast(m_dMarkedNodeIndicesB.data()),
        m_hOutMarkedCountB);
  }

  CUBQL_CUDA_CALL(StreamSynchronize(m_stream));
  double tGpuBfsEnd = cuBQL::getCurrentTime();

  // 4. BATCHED CROSS INTERSECTION LOOP (DOUBLE MODE ONLY)
  if(m_meshAcpu && m_meshBcpu) {
    finalCandidatePairs = executeBatchedCrossIntersectionLoopDouble(
        *m_meshAcpu, *m_meshBcpu, batchMultiplier, totalBatches, m_dOutPairsA, m_dOutPairsB, m_dReverseMapB,
        m_dMarkedNodeIndicesB, m_dOutOffsetsB, m_dOutPrimsFlatB, m_dNodeDescendantCountsB, m_hOutMarkedCountB, m_bvhA,
        m_dBoxesA, // <-- ADDED: Precomputed boxes for Mesh A
        m_dBoxesB, // <-- ADDED: Precomputed boxes for Mesh B
        m_dVertsA, m_dIndicesA, m_dVertsB, m_dIndicesB, outFinalExactPairs, outFinalCount, tracker, m_centerA,
        m_centerB, m_rotA, m_transA, m_rotB, m_transB, gpuDouble, m_stream);

    // if (gpuDouble)
    // {
    
    // }
    // else
    // {
    //   finalCandidatePairs = executeBatchedCrossIntersectionLoopDoubleOld(
    //     *m_meshAcpu, *m_meshBcpu, batchMultiplier, totalBatches, m_dOutPairsA, m_dOutPairsB, m_dReverseMapB,
    //     m_dMarkedNodeIndicesB, m_dOutOffsetsB, m_dOutPrimsFlatB, m_dNodeDescendantCountsB, m_hOutMarkedCountB, m_bvhA,
    //     m_dBoxesA, // <-- ADDED: Precomputed boxes for Mesh A
    //     m_dBoxesB, // <-- ADDED: Precomputed boxes for Mesh B
    //     m_dVertsA, m_dIndicesA, m_dVertsB, m_dIndicesB, outFinalExactPairs, outFinalCount, tracker, m_centerA,
    //     m_centerB, m_rotA, m_transA, m_rotB, m_transB, m_stream);


    // }
  }

  double tPipelineEnd = cuBQL::getCurrentTime();

  // 5. METRICS RECORDING
  stats.meshATotalNodes = m_bvhA.numNodes;
  stats.meshAExtractedTargets = m_hOutMarkedCountA;
  stats.meshBTotalNodes = m_bvhB.numNodes;
  stats.meshBExtractedTargets = m_hOutMarkedCountB;
  stats.totalIntersections = totalIntersections;
  stats.totalPossiblePairs = totalPossiblePairs;
  stats.intersectionPercentage = intersectionPercentage;
  stats.gpuCrossCheckEngineMs = (tCrossEnd - tCrossStart) * 1000.0;
  stats.parallelDfsDescentBMs = (tGpuBfsEnd - tGpuBfsStart) * 1000.0;
  stats.GPUTotalTime = (tPipelineEnd - tPipelineStart) * 1000.0;
  stats.totalCrissCrossBatches = totalBatches;
  stats.finalAabbCandidatePairs = finalCandidatePairs;
  stats.loopTracker = tracker;
}

// -----------------------------------------------------------------------------
// Cleanup Method
// -----------------------------------------------------------------------------
void KernelBVHController::cleanup() {
  clearGPU();

  m_hVertsA = nullptr;
  m_hIndicesA = nullptr;
  m_hVertsB = nullptr;
  m_hIndicesB = nullptr;

  m_meshAcpu = nullptr;
  m_meshBcpu = nullptr;

  m_origPointsA.clear();
  m_origPointsA.shrink_to_fit();
  m_origPointsB.clear();
  m_origPointsB.shrink_to_fit();

  m_centerA = Point3(0, 0, 0);
  m_centerB = Point3(0, 0, 0);
  m_rotA = make_double3(0.0, 0.0, 0.0);
  m_transA = make_double3(0.0, 0.0, 0.0);
  m_rotB = make_double3(0.0, 0.0, 0.0);
  m_transB = make_double3(0.0, 0.0, 0.0);
  shiftX = 0.0;
  shiftY = 0.0;
  shiftZ = 0.0;
}
void KernelBVHController::clearGPU() {
  if(m_stream) {
    CUBQL_CUDA_CALL(StreamSynchronize(m_stream)); // Removed 'cuda'
  }

  // Mesh A GPU Buffers
  if(m_dMeshA) {
    CUBQL_CUDA_CALL(FreeAsync(m_dMeshA, m_stream));
    m_dMeshA = nullptr;
  } // Removed 'cuda'
  if(m_dMeshMetricsA) {
    CUBQL_CUDA_CALL(FreeAsync(m_dMeshMetricsA, m_stream));
    m_dMeshMetricsA = nullptr;
  }
  if(m_dBoxesA) {
    CUBQL_CUDA_CALL(FreeAsync(m_dBoxesA, m_stream));
    m_dBoxesA = nullptr;
  }
  if(m_dVertsA) {
    CUBQL_CUDA_CALL(FreeAsync(m_dVertsA, m_stream));
    m_dVertsA = nullptr;
  }
  if(m_dVertsAOrig) {
    CUBQL_CUDA_CALL(FreeAsync(m_dVertsAOrig, m_stream));
    m_dVertsAOrig = nullptr;
  }
  if(m_dIndicesA) {
    CUBQL_CUDA_CALL(FreeAsync(m_dIndicesA, m_stream));
    m_dIndicesA = nullptr;
  }
  if(m_dVertErrorsA) {
    CUBQL_CUDA_CALL(FreeAsync(m_dVertErrorsA, m_stream));
    m_dVertErrorsA = nullptr;
  }

  // Mesh B GPU Buffers
  if(m_dMeshB) {
    CUBQL_CUDA_CALL(FreeAsync(m_dMeshB, m_stream));
    m_dMeshB = nullptr;
  }
  if(m_dMeshMetricsB) {
    CUBQL_CUDA_CALL(FreeAsync(m_dMeshMetricsB, m_stream));
    m_dMeshMetricsB = nullptr;
  }
  if(m_dBoxesB) {
    CUBQL_CUDA_CALL(FreeAsync(m_dBoxesB, m_stream));
    m_dBoxesB = nullptr;
  }
  if(m_dVertsB) {
    CUBQL_CUDA_CALL(FreeAsync(m_dVertsB, m_stream));
    m_dVertsB = nullptr;
  }
  if(m_dVertsBOrig) {
    CUBQL_CUDA_CALL(FreeAsync(m_dVertsBOrig, m_stream));
    m_dVertsBOrig = nullptr;
  }
  if(m_dIndicesB) {
    CUBQL_CUDA_CALL(FreeAsync(m_dIndicesB, m_stream));
    m_dIndicesB = nullptr;
  }
  if(m_dVertErrorsB) {
    CUBQL_CUDA_CALL(FreeAsync(m_dVertErrorsB, m_stream));
    m_dVertErrorsB = nullptr;
  }

  // cuBQL Acceleration Structures
  cuBQL::cuda::free(m_bvhA, m_stream, m_memResource);
  cuBQL::cuda::free(m_bvhB, m_stream, m_memResource);

  // Clear Thrust Device Vectors
  m_dMarkedNodeIndicesA_Full.clear();
  m_dMarkedNodeIndicesA_Full.shrink_to_fit();
  m_dNodeDescendantCountsA.clear();
  m_dNodeDescendantCountsA.shrink_to_fit();
  m_dMarkedNodeIndicesB_Full.clear();
  m_dMarkedNodeIndicesB_Full.shrink_to_fit();
  m_dNodeDescendantCountsB.clear();
  m_dNodeDescendantCountsB.shrink_to_fit();

  m_dReverseMapB.clear();
  m_dReverseMapB.shrink_to_fit();
  m_dOutPairsA.clear();
  m_dOutPairsA.shrink_to_fit();
  m_dOutPairsB.clear();
  m_dOutPairsB.shrink_to_fit();
  m_dOutOffsetsB.clear();
  m_dOutOffsetsB.shrink_to_fit();
  m_dOutPrimsFlatB.clear();
  m_dOutPrimsFlatB.shrink_to_fit();
  m_dOutOffsetsA.clear();
  m_dOutOffsetsA.shrink_to_fit();
  m_dOutPrimsFlatA.clear();
  m_dOutPrimsFlatA.shrink_to_fit();

  if(m_stream) {
    CUBQL_CUDA_CALL(StreamSynchronize(m_stream)); // Removed 'cuda'
  }
}

void KernelBVHController::reconstructGPU(ExecutionStats& stats) {
  if(!m_hVertsA || !m_hIndicesA || !m_hVertsB || !m_hIndicesB || m_numTrianglesA <= 0 || m_numTrianglesB <= 0) {
    std::cerr << "[KernelBVHController] Error: Cannot reconstruct GPU resources. Host pointers missing.\n";
    return;
  }

  // Clear existing GPU allocations safely
  clearGPU();

  cuBQL::BuildConfig buildConfig;
  buildConfig.makeLeafThreshold = m_leafThreshold;

  double tReconstructStart = cuBQL::getCurrentTime();

  // 1. Re-allocate & copy Mesh A
  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dVertsA, m_numVertsA * sizeof(double3), m_stream));
  CUBQL_CUDA_CALL(MemcpyAsync(m_dVertsA, m_hVertsA, m_numVertsA * sizeof(double3), cudaMemcpyHostToDevice, m_stream));

  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dVertsAOrig, m_numVertsA * sizeof(double3), m_stream));
  CUBQL_CUDA_CALL(
      MemcpyAsync(m_dVertsAOrig, m_hVertsA, m_numVertsA * sizeof(double3), cudaMemcpyHostToDevice, m_stream));

  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dIndicesA, m_numTrianglesA * sizeof(uint3), m_stream));
  CUBQL_CUDA_CALL(
      MemcpyAsync(m_dIndicesA, m_hIndicesA, m_numTrianglesA * sizeof(uint3), cudaMemcpyHostToDevice, m_stream));

  // 2. Re-allocate & copy Mesh B
  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dVertsB, m_numVertsB * sizeof(double3), m_stream));
  CUBQL_CUDA_CALL(MemcpyAsync(m_dVertsB, m_hVertsB, m_numVertsB * sizeof(double3), cudaMemcpyHostToDevice, m_stream));

  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dVertsBOrig, m_numVertsB * sizeof(double3), m_stream));
  CUBQL_CUDA_CALL(
      MemcpyAsync(m_dVertsBOrig, m_hVertsB, m_numVertsB * sizeof(double3), cudaMemcpyHostToDevice, m_stream));

  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dIndicesB, m_numTrianglesB * sizeof(uint3), m_stream));
  CUBQL_CUDA_CALL(
      MemcpyAsync(m_dIndicesB, m_hIndicesB, m_numTrianglesB * sizeof(uint3), cudaMemcpyHostToDevice, m_stream));

  // 3. Alloc & Compute Bounding Boxes
  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dBoxesA, m_numTrianglesA * sizeof(cuBQL::box3f), m_stream));
  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dBoxesB, m_numTrianglesB * sizeof(cuBQL::box3f), m_stream));

  generateBoxesTrisKernel<<<cuBQL::divRoundUp(m_numTrianglesA, 256), 256, 0, m_stream>>>(m_dBoxesA, m_dVertsA,
                                                                                         m_dIndicesA, m_numTrianglesA);
  generateBoxesTrisKernel<<<cuBQL::divRoundUp(m_numTrianglesB, 256), 256, 0, m_stream>>>(m_dBoxesB, m_dVertsB,
                                                                                         m_dIndicesB, m_numTrianglesB);

  // 4. Re-initialize Thrust Vectors
  uint32_t maxPossibleNodesA = 2 * m_numTrianglesA;
  uint32_t maxPossibleNodesB = 2 * m_numTrianglesB;

  m_dMarkedNodeIndicesA_Full.resize(maxPossibleNodesA, 0);
  m_dNodeDescendantCountsA.resize(maxPossibleNodesA, 0);
  m_dMarkedNodeIndicesB_Full.resize(maxPossibleNodesB, 0);
  m_dNodeDescendantCountsB.resize(maxPossibleNodesB, 0);
  m_dReverseMapB.resize(maxPossibleNodesB, 0);

  // 5. Build and Refit BVHs
  cuBQL::gpuBuilder_v2_2::build_custom(m_bvhA, m_dBoxesA, m_numTrianglesA, buildConfig, (uint32_t)m_levelA,
                                       thrust::raw_pointer_cast(m_dMarkedNodeIndicesA_Full.data()),
                                       thrust::raw_pointer_cast(m_dNodeDescendantCountsA.data()),
                                       &m_hOutMarkedCountA_Full, m_stream, m_memResource);
  cuBQL::cuda::refit(m_bvhA, m_dBoxesA, m_stream);

  cuBQL::gpuBuilder_v2_2::build_custom(m_bvhB, m_dBoxesB, m_numTrianglesB, buildConfig, (uint32_t)m_levelB,
                                       thrust::raw_pointer_cast(m_dMarkedNodeIndicesB_Full.data()),
                                       thrust::raw_pointer_cast(m_dNodeDescendantCountsB.data()),
                                       &m_hOutMarkedCountB_Full, m_stream, m_memResource);
  cuBQL::cuda::refit(m_bvhB, m_dBoxesB, m_stream);

  CUBQL_CUDA_CALL(StreamSynchronize(m_stream));

  // Re-apply transformations if any were previously active
  if(m_rotA.x != 0.0 || m_rotA.y != 0.0 || m_rotA.z != 0.0 || m_transA.x != 0.0 || m_transA.y != 0.0 ||
     m_transA.z != 0.0 || m_rotB.x != 0.0 || m_rotB.y != 0.0 || m_rotB.z != 0.0 || m_transB.x != 0.0 ||
     m_transB.y != 0.0 || m_transB.z != 0.0)
  {
    float tGPU, tCPU, tTV, tAT, tGB;
    double3 currentRotA = m_rotA, currentTransA = m_transA;
    double3 currentRotB = m_rotB, currentTransB = m_transB;

    // Reset internal tracker to force re-evaluation
    m_rotA = m_transA = m_rotB = m_transB = make_double3(0, 0, 0);
    setTransformBoth(currentRotA, currentTransA, currentRotB, currentTransB, tGPU, tCPU, tTV, tAT, tGB);
  }

  stats.initialAllocAndCopyMs = (cuBQL::getCurrentTime() - tReconstructStart) * 1000.0;
}