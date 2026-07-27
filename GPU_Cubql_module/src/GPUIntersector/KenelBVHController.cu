#define CUBQL_GPU_BUILDER_IMPLEMENTATION 1

#include <cuda_runtime.h>
#include <cmath>

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
#include "crossCheckFlexible.h"
#include "../testBVH/DualTreeStep.h"
#include "../testBVH/rapidDescendKernel.h"
#include "../testBVH/batchedCrossIntersection.h"
#include "../testBVH/batchedCrossIntersectionV2.h"
#include "../testBVH/batchedCrossIntersectionV3.h"

#ifndef _ALLOC
#define _ALLOC(ptr, count, stream, memResource)                                                                        \
  CUBQL_CUDA_CALL(MallocAsync((void**)&(ptr), (size_t)(count) * sizeof(*(ptr)), stream))
#endif

#ifndef _FREE
#define _FREE(ptr, stream, memResource) CUBQL_CUDA_CALL(FreeAsync((ptr), stream))
#endif

// -----------------------------------------------------------------------------
// Helper 3x3 Matrix Structure & Rigid Body Device Transformations
// -----------------------------------------------------------------------------
struct Mat3x3
{
  float m[3][3];
};

inline Mat3x3 makeRotationMatrixDeg(float rxDeg, float ryDeg, float rzDeg) {
  const float DEG_TO_RAD = 3.14159265358979323846f / 180.0f;
  float radX = rxDeg * DEG_TO_RAD;
  float radY = ryDeg * DEG_TO_RAD;
  float radZ = rzDeg * DEG_TO_RAD;

  float cx = cosf(radX), sx = sinf(radX);
  float cy = cosf(radY), sy = sinf(radY);
  float cz = cosf(radZ), sz = sinf(radZ);

  Mat3x3 R;
  // Composite Rotation: R = Rz * Ry * Rx
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

__device__ inline float3 transformPoint(const float3& p, const Mat3x3& R, float3 center, float3 trans) {
  float x = p.x - center.x;
  float y = p.y - center.y;
  float z = p.z - center.z;

  float rx = R.m[0][0] * x + R.m[0][1] * y + R.m[0][2] * z;
  float ry = R.m[1][0] * x + R.m[1][1] * y + R.m[1][2] * z;
  float rz = R.m[2][0] * x + R.m[2][1] * y + R.m[2][2] * z;

  return make_float3(rx + center.x + trans.x, ry + center.y + trans.y, rz + center.z + trans.z);
}

__global__ void transformVerticesKernel(
    float3* dVertsOut, const float3* dVertsOrig, int numVerts, Mat3x3 R, float3 center, float3 trans) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if(idx >= numVerts)
    return;

  dVertsOut[idx] = transformPoint(dVertsOrig[idx], R, center, trans);
}

// -----------------------------------------------------------------------------
// Parallel TBB CPU CGAL Transformations & Extraction Helpers
// -----------------------------------------------------------------------------
static void
transformCgalMesh(Mesh* mesh, const std::vector<Point3>& origPoints, Point3 center, float3 rotDeg, float3 trans) {
  if(!mesh || origPoints.empty())
    return;

  auto pmap = mesh->points();
  const size_t numPoints = mesh->num_vertices();

  const double DEG_TO_RAD = 3.14159265358979323846 / 180.0;
  double radX = static_cast<double>(rotDeg.x) * DEG_TO_RAD;
  double radY = static_cast<double>(rotDeg.y) * DEG_TO_RAD;
  double radZ = static_cast<double>(rotDeg.z) * DEG_TO_RAD;

  double cx = cos(radX), sx = sin(radX);
  double cy = cos(radY), sy = sin(radY);
  double cz = cos(radZ), sz = sin(radZ);

  double r00 = cy * cz, r01 = sx * sy * cz - cx * sz, r02 = cx * sy * cz + sx * sz;
  double r10 = cy * sz, r11 = sx * sy * sz + cx * cz, r12 = cx * sy * sz - sx * cz;
  double r20 = -sy, r21 = sx * cy, r22 = cx * cy;

  double cx_p = center.x(), cy_p = center.y(), cz_p = center.z();
  double tx = static_cast<double>(trans.x);
  double ty = static_cast<double>(trans.y);
  double tz = static_cast<double>(trans.z);

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

static void translateCgalAndExtractFloat3(
    Mesh* mesh, const std::vector<Point3>& origPoints, float3* hOutVerts, float xB, float yB, float zB) {
  if(!mesh || origPoints.empty() || !hOutVerts)
    return;

  auto pmap = mesh->points();
  const size_t numPoints = mesh->num_vertices();

  tbb::parallel_for(tbb::blocked_range<size_t>(0, numPoints), [mesh, &pmap, &origPoints, hOutVerts, xB, yB,
                                                               zB](const tbb::blocked_range<size_t>& range) {
    for(size_t i = range.begin(); i != range.end(); ++i) {
      Mesh::Vertex_index vd(static_cast<uint32_t>(i));
      if(!mesh->is_removed(vd)) {
        // 1. High-precision double translation from pristine baseline
        const Point3& pOrig = origPoints[i];
        Point3 pTranslated(pOrig.x() + static_cast<double>(xB), pOrig.y() + static_cast<double>(yB),
                           pOrig.z() + static_cast<double>(zB));

        // Update CGAL host mesh
        pmap[vd] = pTranslated;

        // 2. Quantize double -> float3 for GPU upload
        hOutVerts[i] = make_float3(static_cast<float>(pTranslated.x()), static_cast<float>(pTranslated.y()),
                                   static_cast<float>(pTranslated.z()));
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
void KernelBVHController::setTransformBoth(float3 rotDegA, float3 transA, float3 rotDegB, float3 transB) {
  // 1. Change detection per mesh
  bool movedA = (rotDegA.x != m_rotA.x || rotDegA.y != m_rotA.y || rotDegA.z != m_rotA.z || transA.x != m_transA.x ||
                 transA.y != m_transA.y || transA.z != m_transA.z);

  bool movedB = (rotDegB.x != m_rotB.x || rotDegB.y != m_rotB.y || rotDegB.z != m_rotB.z || transB.x != m_transB.x ||
                 transB.y != m_transB.y || transB.z != m_transB.z);

  if(!movedA && !movedB) {
    return;
  }

  bool gpuWorkQueued = false;
  int block = 256;

  // ---------------------------------------------------------------------------
  // STEP 1: QUEUE ALL GPU WORK FIRST (NON-BLOCKING COMMAND DISPATCH)
  // ---------------------------------------------------------------------------
  if(movedA && m_numVertsA > 0 && m_dVertsA && m_dVertsAOrig) {
    bool isTransA = (transA.x != 0.0f || transA.y != 0.0f || transA.z != 0.0f);
    bool isRotA = (rotDegA.x != 0.0f || rotDegA.y != 0.0f || rotDegA.z != 0.0f);

    Mat3x3 RA = makeRotationMatrixDeg(rotDegA.x, rotDegA.y, rotDegA.z);
    float3 cA = make_float3(static_cast<float>(m_centerA.x()), static_cast<float>(m_centerA.y()),
                            static_cast<float>(m_centerA.z()));

    int gridVertsA = cuBQL::divRoundUp(m_numVertsA, block);
    transformVerticesKernel<<<gridVertsA, block, 0, m_stream>>>(m_dVertsA, m_dVertsAOrig, m_numVertsA, RA, cA, transA);

    int gridTrisA = cuBQL::divRoundUp(m_numTrianglesA, block);
    assembleTrianglesKernelTransformed<<<gridTrisA, block, 0, m_stream>>>(
        m_dMeshA, m_dMeshMetricsA, m_dVertsA, m_dVertErrorsA, m_dIndicesA, m_numTrianglesA, isTransA, isRotA);

    generateBoxes<<<gridTrisA, block, 0, m_stream>>>(m_dBoxesA, m_dMeshA, m_numTrianglesA);
    cuBQL::cuda::refit(m_bvhA, m_dBoxesA, m_stream);

    gpuWorkQueued = true;
  }

  if(movedB && m_numVertsB > 0 && m_dVertsB && m_dVertsBOrig) {
    bool isTransB = (transB.x != 0.0f || transB.y != 0.0f || transB.z != 0.0f);
    bool isRotB = (rotDegB.x != 0.0f || rotDegB.y != 0.0f || rotDegB.z != 0.0f);

    Mat3x3 RB = makeRotationMatrixDeg(rotDegB.x, rotDegB.y, rotDegB.z);
    float3 cB = make_float3(static_cast<float>(m_centerB.x()), static_cast<float>(m_centerB.y()),
                            static_cast<float>(m_centerB.z()));

    int gridVertsB = cuBQL::divRoundUp(m_numVertsB, block);
    transformVerticesKernel<<<gridVertsB, block, 0, m_stream>>>(m_dVertsB, m_dVertsBOrig, m_numVertsB, RB, cB, transB);

    int gridTrisB = cuBQL::divRoundUp(m_numTrianglesB, block);
    assembleTrianglesKernelTransformed<<<gridTrisB, block, 0, m_stream>>>(
        m_dMeshB, m_dMeshMetricsB, m_dVertsB, m_dVertErrorsB, m_dIndicesB, m_numTrianglesB, isTransB, isRotB);

    generateBoxes<<<gridTrisB, block, 0, m_stream>>>(m_dBoxesB, m_dMeshB, m_numTrianglesB);
    cuBQL::cuda::refit(m_bvhB, m_dBoxesB, m_stream);

    gpuWorkQueued = true;
  }

  // ---------------------------------------------------------------------------
  // STEP 2: CONCURRENT CPU CGAL TRANSFORMS (TBB INVOKE OVERLAPPED WITH GPU)
  // ---------------------------------------------------------------------------
  tbb::parallel_invoke(
      [this, movedA, rotDegA, transA]() {
        if(movedA && m_meshAcpu) {
          transformCgalMesh(m_meshAcpu, m_origPointsA, m_centerA, rotDegA, transA);
        }
      },
      [this, movedB, rotDegB, transB]() {
        if(movedB && m_meshBcpu) {
          transformCgalMesh(m_meshBcpu, m_origPointsB, m_centerB, rotDegB, transB);
        }
      });

  // Update persistent state markers
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

  // ---------------------------------------------------------------------------
  // STEP 3: SYNCHRONIZE STREAM
  // ---------------------------------------------------------------------------
  if(gpuWorkQueued) {
    CUBQL_CUDA_CALL(StreamSynchronize(m_stream));
  }
}
// Legacy Translation Handlers
void KernelBVHController::setTranslation(float xB, float yB, float zB) {
  setTransformBoth(m_rotA, m_transA, m_rotB, make_float3(xB, yB, zB));
}

void KernelBVHController::setTranslationCPUHostUpload(float xB, float yB, float zB) {
  // 1. Exit if translation target hasn't changed
  if(xB == shiftX && yB == shiftY && zB == shiftZ)
    return;

  if(!m_meshBcpu || m_origPointsB.empty() || m_numVertsB <= 0)
    return;

  // 2. CPU: Update CGAL Mesh in double precision & prepare host float3 array
  std::vector<float3> h_vertsBUpdated(m_numVertsB);
  translateCgalAndExtractFloat3(m_meshBcpu, m_origPointsB, h_vertsBUpdated.data(), xB, yB, zB);

  // 3. GPU: Async upload updated float3 vertices to device memory (m_dVertsB)
  CUBQL_CUDA_CALL(
      MemcpyAsync(m_dVertsB, h_vertsBUpdated.data(), m_numVertsB * sizeof(float3), cudaMemcpyHostToDevice, m_stream));

  // 4. GPU: Triangle Assembly with Transformation Drift Flags
  int block = 256;
  int gridTris = cuBQL::divRoundUp(m_numTrianglesB, block);
  bool isTransB = (xB != 0.0f || yB != 0.0f || zB != 0.0f);
  bool isRotB = (m_rotB.x != 0.0f || m_rotB.y != 0.0f || m_rotB.z != 0.0f);

  assembleTrianglesKernelTransformed<<<gridTris, block, 0, m_stream>>>(
      m_dMeshB, m_dMeshMetricsB, m_dVertsB, m_dVertErrorsB, m_dIndicesB, m_numTrianglesB, isTransB, isRotB);

  // 5. GPU: Generate tight AABBs and refit BVH B
  generateBoxes<<<gridTris, block, 0, m_stream>>>(m_dBoxesB, m_dMeshB, m_numTrianglesB);
  cuBQL::cuda::refit(m_bvhB, m_dBoxesB, m_stream);

  // 6. Synchronize GPU stream
  CUBQL_CUDA_CALL(StreamSynchronize(m_stream));

  // Update persistent translation markers
  shiftX = xB;
  shiftY = yB;
  shiftZ = zB;
  m_transB = make_float3(xB, yB, zB);
}

// -----------------------------------------------------------------------------
// Controller Construction Method
// -----------------------------------------------------------------------------
void KernelBVHController::construct(Mesh& meshAcpu,
                                    Mesh& meshBcpu,
                                    const float3* hVertsA,
                                    int numVertsA,
                                    const uint3* hIndicesA,
                                    const float* hVertErrorsA,
                                    int numTrianglesA,
                                    int levelA,
                                    const float3* hVertsB,
                                    int numVertsB,
                                    const uint3* hIndicesB,
                                    const float* hVertErrorsB,
                                    int numTrianglesB,
                                    int levelB,
                                    int leafThreshold,
                                    ExecutionStats& stats) {
  construct(meshAcpu, meshBcpu, Point3(0, 0, 0), Point3(0, 0, 0), hVertsA, numVertsA, hIndicesA, hVertErrorsA,
            numTrianglesA, levelA, hVertsB, numVertsB, hIndicesB, hVertErrorsB, numTrianglesB, levelB, leafThreshold,
            stats);
}

void KernelBVHController::construct(Mesh& meshAcpu,
                                    Mesh& meshBcpu,
                                    const Point3& centerA,
                                    const Point3& centerB,
                                    const float3* hVertsA,
                                    int numVertsA,
                                    const uint3* hIndicesA,
                                    const float* hVertErrorsA,
                                    int numTrianglesA,
                                    int levelA,
                                    const float3* hVertsB,
                                    int numVertsB,
                                    const uint3* hIndicesB,
                                    const float* hVertErrorsB,
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

  m_centerA = centerA;
  m_centerB = centerB;

  m_rotA = make_float3(0.0f, 0.0f, 0.0f);
  m_transA = make_float3(0.0f, 0.0f, 0.0f);
  m_rotB = make_float3(0.0f, 0.0f, 0.0f);
  m_transB = make_float3(0.0f, 0.0f, 0.0f);
  shiftX = 0.0f;
  shiftY = 0.0f;
  shiftZ = 0.0f;

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

  // --------------------------------------------------------------------
  // 1. ASYNCHRONOUS PERSISTENT DEVICE ALLOCATION
  // --------------------------------------------------------------------
  // Mesh A Buffers
  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dVertsA, m_numVertsA * sizeof(float3), m_stream));
  CUBQL_CUDA_CALL(MemcpyAsync(m_dVertsA, hVertsA, m_numVertsA * sizeof(float3), cudaMemcpyHostToDevice, m_stream));

  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dVertsAOrig, m_numVertsA * sizeof(float3), m_stream));
  CUBQL_CUDA_CALL(MemcpyAsync(m_dVertsAOrig, hVertsA, m_numVertsA * sizeof(float3), cudaMemcpyHostToDevice, m_stream));

  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dIndicesA, m_numTrianglesA * sizeof(uint3), m_stream));
  CUBQL_CUDA_CALL(
      MemcpyAsync(m_dIndicesA, hIndicesA, m_numTrianglesA * sizeof(uint3), cudaMemcpyHostToDevice, m_stream));

  m_dVertErrorsA = nullptr;
  if(hVertErrorsA) {
    CUBQL_CUDA_CALL(MallocAsync((void**)&m_dVertErrorsA, m_numVertsA * sizeof(float), m_stream));
    CUBQL_CUDA_CALL(
        MemcpyAsync(m_dVertErrorsA, hVertErrorsA, m_numVertsA * sizeof(float), cudaMemcpyHostToDevice, m_stream));
  }

  // Mesh B Buffers
  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dVertsB, m_numVertsB * sizeof(float3), m_stream));
  CUBQL_CUDA_CALL(MemcpyAsync(m_dVertsB, hVertsB, m_numVertsB * sizeof(float3), cudaMemcpyHostToDevice, m_stream));

  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dVertsBOrig, m_numVertsB * sizeof(float3), m_stream));
  CUBQL_CUDA_CALL(MemcpyAsync(m_dVertsBOrig, hVertsB, m_numVertsB * sizeof(float3), cudaMemcpyHostToDevice, m_stream));

  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dIndicesB, m_numTrianglesB * sizeof(uint3), m_stream));
  CUBQL_CUDA_CALL(
      MemcpyAsync(m_dIndicesB, hIndicesB, m_numTrianglesB * sizeof(uint3), cudaMemcpyHostToDevice, m_stream));

  m_dVertErrorsB = nullptr;
  if(hVertErrorsB) {
    CUBQL_CUDA_CALL(MallocAsync((void**)&m_dVertErrorsB, m_numVertsB * sizeof(float), m_stream));
    CUBQL_CUDA_CALL(
        MemcpyAsync(m_dVertErrorsB, hVertErrorsB, m_numVertsB * sizeof(float), cudaMemcpyHostToDevice, m_stream));
  }

  // Primitive & BVH Bounding Box Allocation
  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dMeshA, m_numTrianglesA * sizeof(cuBQL::Triangle), m_stream));
  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dMeshMetricsA, m_numTrianglesA * sizeof(float2), m_stream));
  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dBoxesA, m_numTrianglesA * sizeof(cuBQL::box3f), m_stream));

  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dMeshB, m_numTrianglesB * sizeof(cuBQL::Triangle), m_stream));
  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dMeshMetricsB, m_numTrianglesB * sizeof(float2), m_stream));
  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dBoxesB, m_numTrianglesB * sizeof(cuBQL::box3f), m_stream));

  // Initial GPU Triangle Assembly & Box Generation
  assembleTrianglesKernel<<<cuBQL::divRoundUp(m_numTrianglesA, 256), 256, 0, m_stream>>>(
      m_dMeshA, m_dMeshMetricsA, m_dVertsA, m_dVertErrorsA, m_dIndicesA, m_numTrianglesA);
  generateBoxes<<<cuBQL::divRoundUp(m_numTrianglesA, 256), 256, 0, m_stream>>>(m_dBoxesA, m_dMeshA, m_numTrianglesA);

  assembleTrianglesKernel<<<cuBQL::divRoundUp(m_numTrianglesB, 256), 256, 0, m_stream>>>(
      m_dMeshB, m_dMeshMetricsB, m_dVertsB, m_dVertErrorsB, m_dIndicesB, m_numTrianglesB);
  generateBoxes<<<cuBQL::divRoundUp(m_numTrianglesB, 256), 256, 0, m_stream>>>(m_dBoxesB, m_dMeshB, m_numTrianglesB);

  CUBQL_CUDA_CALL(EventRecord(evAllocStop, m_stream));

  // --------------------------------------------------------------------
  // 2. CONCURRENT CPU CGAL BASELINE EXTRACTION
  // --------------------------------------------------------------------
  tbb::parallel_invoke(
      [this]() {
        if(m_meshAcpu) {
          size_t numVerts = m_meshAcpu->num_vertices();
          m_origPointsA.resize(numVerts);
          auto pmap = m_meshAcpu->points();
          tbb::parallel_for(tbb::blocked_range<size_t>(0, numVerts),
                            [&pmap, this](const tbb::blocked_range<size_t>& range) {
                              for(size_t i = range.begin(); i != range.end(); ++i) {
                                Mesh::Vertex_index vd(static_cast<uint32_t>(i));
                                m_origPointsA[i] = pmap[vd];
                              }
                            });
        }
      },
      [this]() {
        if(m_meshBcpu) {
          size_t numVerts = m_meshBcpu->num_vertices();
          m_origPointsB.resize(numVerts);
          auto pmap = m_meshBcpu->points();
          tbb::parallel_for(tbb::blocked_range<size_t>(0, numVerts),
                            [&pmap, this](const tbb::blocked_range<size_t>& range) {
                              for(size_t i = range.begin(); i != range.end(); ++i) {
                                Mesh::Vertex_index vd(static_cast<uint32_t>(i));
                                m_origPointsB[i] = pmap[vd];
                              }
                            });
        }
      });

  // --------------------------------------------------------------------
  // 3. THRUST DEVICE VECTOR INITIALIZATION
  // --------------------------------------------------------------------
  double tThrustInitStart = cuBQL::getCurrentTime();

  uint32_t maxPossibleNodesA = 2 * m_numTrianglesA;
  uint32_t maxPossibleNodesB = 2 * m_numTrianglesB;

  m_dMarkedNodeIndicesA_Full.resize(maxPossibleNodesA, 0);
  m_dNodeDescendantCountsA.resize(maxPossibleNodesA, 0);
  m_dMarkedNodeIndicesB_Full.resize(maxPossibleNodesB, 0);
  m_dNodeDescendantCountsB.resize(maxPossibleNodesB, 0);
  m_dReverseMapB.resize(maxPossibleNodesB, 0);

  stats.thrustInitOverheadMs = (cuBQL::getCurrentTime() - tThrustInitStart) * 1000.0;

  // --------------------------------------------------------------------
  // 4. BUILD & REFIT BVHs
  // --------------------------------------------------------------------
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

  // --------------------------------------------------------------------
  // 5. STREAM SYNCHRONIZATION & METRIC CALCULATION
  // --------------------------------------------------------------------
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
                                                  tbb::concurrent_vector<int2>& finalExactPairs,
                                                  ExecutionStats& stats) {
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
      executeCrossIntersectionFlexible(m_bvhA, m_dMarkedNodeIndicesA, m_hOutMarkedCountA, m_bvhB, m_dMarkedNodeIndicesB,
                                       m_hOutMarkedCountB, m_dOutPairsA, m_dOutPairsB, 0, 0, 0);

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
                        m_hOutMarkedCountB, m_bvhA, m_bvhB, m_dMeshA, m_dMeshB);

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

  // 4. BATCHED CROSS INTERSECTION LOOP
  if(m_meshAcpu && m_meshBcpu) {
    if(activateAsyncDownload == 0) {
      finalCandidatePairs = executeBatchedCrossIntersectionLoopV2(
          *m_meshAcpu, *m_meshBcpu, batchMultiplier, totalBatches, m_dOutPairsA, m_dOutPairsB, m_dReverseMapB,
          m_dMarkedNodeIndicesB, m_dOutOffsetsB, m_dOutPrimsFlatB, m_dNodeDescendantCountsB, m_hOutMarkedCountB, m_bvhA,
          m_dMeshA, m_dMeshB, m_dMeshMetricsA, m_dMeshMetricsB, finalExactPairs, tracker, m_stream);
    } else {
      finalCandidatePairs = executeBatchedCrossIntersectionLoopV3(
          *m_meshAcpu, *m_meshBcpu, batchMultiplier, totalBatches, m_dOutPairsA, m_dOutPairsB, m_dReverseMapB,
          m_dMarkedNodeIndicesB, m_dOutOffsetsB, m_dOutPrimsFlatB, m_dNodeDescendantCountsB, m_hOutMarkedCountB, m_bvhA,
          m_dMeshA, m_dMeshB, m_dMeshMetricsA, m_dMeshMetricsB, finalExactPairs, tracker, m_stream,
          activateAsyncDownload);
    }
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
  if(m_stream)
    CUBQL_CUDA_CALL(StreamSynchronize(m_stream));

  // Mesh A GPU Buffers
  if(m_dMeshA) {
    CUBQL_CUDA_CALL(FreeAsync(m_dMeshA, m_stream));
    m_dMeshA = nullptr;
  }
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

  // BVH Acceleration Structures
  cuBQL::cuda::free(m_bvhA, m_stream, m_memResource);
  cuBQL::cuda::free(m_bvhB, m_stream, m_memResource);

  // Thrust Vector Reset
  m_dMarkedNodeIndicesA_Full = thrust::device_vector<uint32_t>();
  m_dNodeDescendantCountsA = thrust::device_vector<uint32_t>();
  m_dMarkedNodeIndicesB_Full = thrust::device_vector<uint32_t>();
  m_dNodeDescendantCountsB = thrust::device_vector<uint32_t>();

  m_dReverseMapB = thrust::device_vector<uint32_t>();
  m_dOutPairsA = thrust::device_vector<uint32_t>();
  m_dOutPairsB = thrust::device_vector<uint32_t>();
  m_dOutOffsetsB = thrust::device_vector<uint32_t>();
  m_dOutPrimsFlatB = thrust::device_vector<uint32_t>();
  m_dOutOffsetsA = thrust::device_vector<uint32_t>();
  m_dOutPrimsFlatA = thrust::device_vector<uint32_t>();

  m_origPointsA.clear();
  m_origPointsA.shrink_to_fit();
  m_origPointsB.clear();
  m_origPointsB.shrink_to_fit();

  m_centerA = Point3(0, 0, 0);
  m_centerB = Point3(0, 0, 0);
  m_rotA = make_float3(0.0f, 0.0f, 0.0f);
  m_transA = make_float3(0.0f, 0.0f, 0.0f);
  m_rotB = make_float3(0.0f, 0.0f, 0.0f);
  m_transB = make_float3(0.0f, 0.0f, 0.0f);
  shiftX = 0.0f;
  shiftY = 0.0f;
  shiftZ = 0.0f;

  if(m_stream)
    CUBQL_CUDA_CALL(StreamSynchronize(m_stream));
}