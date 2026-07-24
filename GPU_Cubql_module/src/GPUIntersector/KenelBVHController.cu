#define CUBQL_GPU_BUILDER_IMPLEMENTATION 1

#include <cuda_runtime.h>

#include "KernelBVHController.h"

#include <tbb/parallel_for.h>

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
#include <iostream> // For deep debug printings
#include "samples/common/loadOBJ.h"

// Custom modules
#include "../testBVH/utils.h"
#include "crossCheckFlexible.h"
#include "../testBVH/DualTreeStep.h"
#include "../testBVH/rapidDescendKernel.h"
#include "../testBVH/batchedCrossIntersection.h"
#include "../testBVH/batchedCrossIntersectionV2.h"

#ifndef _ALLOC
#define _ALLOC(ptr, count, stream, memResource)                                                                        \
  CUBQL_CUDA_CALL(MallocAsync((void**)&(ptr), (size_t)(count) * sizeof(*(ptr)), stream))
#endif

#ifndef _FREE
#define _FREE(ptr, stream, memResource) CUBQL_CUDA_CALL(FreeAsync((ptr), stream))
#endif

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

// Kernel to translate both BVH AABBs and triangle mesh primitives in parallel by a delta offset
__global__ void translateMeshAndBoxesKernel(
    cuBQL::Triangle* dMesh, cuBQL::box3f* dBoxes, int numTriangles, float dx, float dy, float dz) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if(idx >= numTriangles)
    return;

  // 1. Shift Triangle Vertices by delta
  dMesh[idx].a.x += dx;
  dMesh[idx].a.y += dy;
  dMesh[idx].a.z += dz;
  dMesh[idx].b.x += dx;
  dMesh[idx].b.y += dy;
  dMesh[idx].b.z += dz;
  dMesh[idx].c.x += dx;
  dMesh[idx].c.y += dy;
  dMesh[idx].c.z += dz;

  // 2. Shift Bounding Box by delta
  dBoxes[idx].lower[0] += dx;
  dBoxes[idx].upper[0] += dx;
  dBoxes[idx].lower[1] += dy;
  dBoxes[idx].upper[1] += dy;
  dBoxes[idx].lower[2] += dz;
  dBoxes[idx].upper[2] += dz;
}

__global__ void translateVerticesKernel(
    float3* dVertsOut, const float3* dVertsOrig, int numVerts, float targetX, float targetY, float targetZ) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if(idx >= numVerts)
    return;

  float3 v = dVertsOrig[idx];
  dVertsOut[idx] = make_float3(v.x + targetX, v.y + targetY, v.z + targetZ);
}

void translateCgalMeshB(Mesh* mesh, const std::vector<Point3>& origPoints, float xB, float yB, float zB) {
  if(!mesh || origPoints.empty())
    return;

  auto pmap = mesh->points();
  const size_t numPoints = mesh->num_vertices();

  tbb::parallel_for(tbb::blocked_range<size_t>(0, numPoints),
                    [mesh, &pmap, &origPoints, xB, yB, zB](const tbb::blocked_range<size_t>& range) {
                      for(size_t i = range.begin(); i != range.end(); ++i) {
                        Mesh::Vertex_index vd(static_cast<uint32_t>(i));
                        if(!mesh->is_removed(vd)) {
                          // Compute absolute offset directly from pristine baseline point
                          const Point3& pOrig = origPoints[i];
                          pmap[vd] = Point3(pOrig.x() + xB, pOrig.y() + yB, pOrig.z() + zB);
                        }
                      }
                    });
}

void KernelBVHController::setTranslation(float xB, float yB, float zB) {
  // 1. Quick exit if target translation hasn't changed
  if(xB == shiftX && yB == shiftY && zB == shiftZ)
    return;

  bool gpuWorkQueued = false;

  // 2. GPU PIPELINE (Absolute shift from m_dVertsBOrig)
  if(m_numVertsB > 0 && m_dVertsB && m_dVertsBOrig) {
    int block = 256;

    int gridVerts = cuBQL::divRoundUp(m_numVertsB, block);
    translateVerticesKernel<<<gridVerts, block, 0, m_stream>>>(m_dVertsB, m_dVertsBOrig, m_numVertsB, xB, yB, zB);

    int gridTris = cuBQL::divRoundUp(m_numTrianglesB, block);
    assembleTrianglesKernelTranslated<<<gridTris, block, 0, m_stream>>>(
        m_dMeshB, m_dMeshMetricsB, m_dVertsB, m_dVertErrorsB, m_dIndicesB, m_numTrianglesB, true);

    generateBoxes<<<gridTris, block, 0, m_stream>>>(m_dBoxesB, m_dMeshB, m_numTrianglesB);

    cuBQL::cuda::refit(m_bvhB, m_dBoxesB, m_stream);

    gpuWorkQueued = true;
  }

  // 3. CPU CGAL SHIFT (Absolute shift from m_origPointsB)
  if(m_meshBcpu) {
    translateCgalMeshB(m_meshBcpu, m_origPointsB, xB, yB, zB);
  }

  // 4. SYNCHRONIZE GPU
  if(gpuWorkQueued) {
    CUBQL_CUDA_CALL(StreamSynchronize(m_stream));
  }

  shiftX = xB;
  shiftY = yB;
  shiftZ = zB;
}

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
  if(numTrianglesA <= 0 || numTrianglesB <= 0)
    return;

  // Cleanup old allocations if reconstruct is called
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

  shiftX = 0.0f;
  shiftY = 0.0f;
  shiftZ = 0.0f;

  cuBQL::BuildConfig buildConfig;
  buildConfig.makeLeafThreshold = m_leafThreshold;

  // --------------------------------------------------------------------
  // CUDA TIMING EVENTS (Non-blocking GPU timing)
  // --------------------------------------------------------------------
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
  // 1. ASYNCHRONOUS DEVICE ALLOCATION & STAGING (GPU QUEUE)
  // --------------------------------------------------------------------
  float3* dVertsA = nullptr;
  uint3* dIndicesA = nullptr;
  float* dVertErrorsA = nullptr;

  // Mesh A Staging
  CUBQL_CUDA_CALL(MallocAsync((void**)&dVertsA, m_numVertsA * sizeof(float3), m_stream));
  CUBQL_CUDA_CALL(MemcpyAsync(dVertsA, hVertsA, m_numVertsA * sizeof(float3), cudaMemcpyHostToDevice, m_stream));

  CUBQL_CUDA_CALL(MallocAsync((void**)&dIndicesA, m_numTrianglesA * sizeof(uint3), m_stream));
  CUBQL_CUDA_CALL(MemcpyAsync(dIndicesA, hIndicesA, m_numTrianglesA * sizeof(uint3), cudaMemcpyHostToDevice, m_stream));

  if(hVertErrorsA) {
    CUBQL_CUDA_CALL(MallocAsync((void**)&dVertErrorsA, m_numVertsA * sizeof(float), m_stream));
    CUBQL_CUDA_CALL(
        MemcpyAsync(dVertErrorsA, hVertErrorsA, m_numVertsA * sizeof(float), cudaMemcpyHostToDevice, m_stream));
  }

  // Mesh B Persistent Buffers
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

  // Mesh Primitive & Metric Buffers
  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dMeshA, m_numTrianglesA * sizeof(cuBQL::Triangle), m_stream));
  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dMeshMetricsA, m_numTrianglesA * sizeof(float2), m_stream));
  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dBoxesA, m_numTrianglesA * sizeof(cuBQL::box3f), m_stream));

  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dMeshB, m_numTrianglesB * sizeof(cuBQL::Triangle), m_stream));
  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dMeshMetricsB, m_numTrianglesB * sizeof(float2), m_stream));
  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dBoxesB, m_numTrianglesB * sizeof(cuBQL::box3f), m_stream));

  // Launch Assembly Kernels
  assembleTrianglesKernel<<<cuBQL::divRoundUp(m_numTrianglesA, 256), 256, 0, m_stream>>>(
      m_dMeshA, m_dMeshMetricsA, dVertsA, dVertErrorsA, dIndicesA, m_numTrianglesA);
  generateBoxes<<<cuBQL::divRoundUp(m_numTrianglesA, 256), 256, 0, m_stream>>>(m_dBoxesA, m_dMeshA, m_numTrianglesA);

  assembleTrianglesKernel<<<cuBQL::divRoundUp(m_numTrianglesB, 256), 256, 0, m_stream>>>(
      m_dMeshB, m_dMeshMetricsB, m_dVertsB, m_dVertErrorsB, m_dIndicesB, m_numTrianglesB);
  generateBoxes<<<cuBQL::divRoundUp(m_numTrianglesB, 256), 256, 0, m_stream>>>(m_dBoxesB, m_dMeshB, m_numTrianglesB);

  // Free staging memory asynchronously on stream
  CUBQL_CUDA_CALL(FreeAsync(dVertsA, m_stream));
  CUBQL_CUDA_CALL(FreeAsync(dIndicesA, m_stream));
  if(dVertErrorsA) {
    CUBQL_CUDA_CALL(FreeAsync(dVertErrorsA, m_stream));
  }

  CUBQL_CUDA_CALL(EventRecord(evAllocStop, m_stream));

  // --------------------------------------------------------------------
  // 2. PARALLEL CPU TASK (OVERLAPPED WITH GPU QUEUE EXECUTION)
  // --------------------------------------------------------------------
  // While the GPU runs allocations and triangle assembly in m_stream,
  // the CPU concurrently extracts pristine CGAL points using TBB threads.
  if(m_meshBcpu) {
    size_t numVerts = m_meshBcpu->num_vertices();
    m_origPointsB.resize(numVerts);
    auto pmap = m_meshBcpu->points();

    tbb::parallel_for(tbb::blocked_range<size_t>(0, numVerts), [&pmap, this](const tbb::blocked_range<size_t>& range) {
      for(size_t i = range.begin(); i != range.end(); ++i) {
        Mesh::Vertex_index vd(static_cast<uint32_t>(i));
        m_origPointsB[i] = pmap[vd];
      }
    });
  }

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
  // 4. BUILD & REFIT BVHs (ASYNCHRONOUS)
  // --------------------------------------------------------------------
  CUBQL_CUDA_CALL(EventRecord(evBuildAStart, m_stream));
  cuBQL::gpuBuilder_v2_2::build_custom(m_bvhA, m_dBoxesA, m_numTrianglesA, buildConfig, (uint32_t)m_levelA,
                                       thrust::raw_pointer_cast(m_dMarkedNodeIndicesA_Full.data()),
                                       thrust::raw_pointer_cast(m_dNodeDescendantCountsA.data()), &m_hOutMarkedCountA_Full,
                                       m_stream, m_memResource);
  cuBQL::cuda::refit(m_bvhA, m_dBoxesA, m_stream);
  CUBQL_CUDA_CALL(EventRecord(evBuildAStop, m_stream));

  CUBQL_CUDA_CALL(EventRecord(evBuildBStart, m_stream));
  cuBQL::gpuBuilder_v2_2::build_custom(m_bvhB, m_dBoxesB, m_numTrianglesB, buildConfig, (uint32_t)m_levelB,
                                       thrust::raw_pointer_cast(m_dMarkedNodeIndicesB_Full.data()),
                                       thrust::raw_pointer_cast(m_dNodeDescendantCountsB.data()), &m_hOutMarkedCountB_Full,
                                       m_stream, m_memResource);
  cuBQL::cuda::refit(m_bvhB, m_dBoxesB, m_stream);
  CUBQL_CUDA_CALL(EventRecord(evBuildBStop, m_stream));

  // --------------------------------------------------------------------
  // 5. SINGLE FINAL SYNCHRONIZATION & METRIC CALCULATION
  // --------------------------------------------------------------------
  CUBQL_CUDA_CALL(StreamSynchronize(m_stream));

  

  float msAlloc = 0.0f, msBuildA = 0.0f, msBuildB = 0.0f;
  CUBQL_CUDA_CALL(EventElapsedTime(&msAlloc, evAllocStart, evAllocStop));
  CUBQL_CUDA_CALL(EventElapsedTime(&msBuildA, evBuildAStart, evBuildAStop));
  CUBQL_CUDA_CALL(EventElapsedTime(&msBuildB, evBuildBStart, evBuildBStop));

  stats.initialAllocAndCopyMs = static_cast<double>(msAlloc);
  stats.buildRefitMeshAMs = static_cast<double>(msBuildA);
  stats.buildRefitMeshBMs = static_cast<double>(msBuildB);

  // Clean up timing events
  CUBQL_CUDA_CALL(EventDestroy(evAllocStart));
  CUBQL_CUDA_CALL(EventDestroy(evAllocStop));
  CUBQL_CUDA_CALL(EventDestroy(evBuildAStart));
  CUBQL_CUDA_CALL(EventDestroy(evBuildAStop));
  CUBQL_CUDA_CALL(EventDestroy(evBuildBStart));
  CUBQL_CUDA_CALL(EventDestroy(evBuildBStop));
}

void KernelBVHController::runIntersectionPipeline(int batchMultiplier,
                                                  int mode,
                                                  tbb::concurrent_vector<int2>& finalExactPairs,
                                                  ExecutionStats& stats) {
  double tPipelineStart = cuBQL::getCurrentTime();

  // RESET STALE RUN DATA FOR MULTIPLE RUN SAFETY
  m_dOutPairsA.clear();
  m_dOutPairsB.clear();

  uint32_t m_hOutMarkedCountA = m_hOutMarkedCountA_Full;
  uint32_t m_hOutMarkedCountB = m_hOutMarkedCountB_Full;
  thrust::device_vector<uint32_t> m_dMarkedNodeIndicesA = m_dMarkedNodeIndicesA_Full;
  thrust::device_vector<uint32_t> m_dMarkedNodeIndicesB = m_dMarkedNodeIndicesB_Full;


  // --------------------------------------------------------------------
  // 1. CRISS-CROSS INTERSECTION
  // --------------------------------------------------------------------
  double tCrossStart = cuBQL::getCurrentTime();

  uint32_t totalIntersections =
      executeCrossIntersectionFlexible(m_bvhA, m_dMarkedNodeIndicesA, m_hOutMarkedCountA, m_bvhB, m_dMarkedNodeIndicesB,
                                       m_hOutMarkedCountB, m_dOutPairsA, m_dOutPairsB, 0, 0, 0);

  CUBQL_CUDA_CALL(StreamSynchronize(m_stream));
  double tCrossEnd = cuBQL::getCurrentTime();

  uint64_t totalPossiblePairs = (uint64_t)m_hOutMarkedCountA * m_hOutMarkedCountB;
  double intersectionPercentage =
      totalPossiblePairs > 0 ? ((double)totalIntersections / totalPossiblePairs) * 100.0 : 0.0;

  // --------------------------------------------------------------------
  // 2. DUAL-TREE TRAVERSAL STEP
  // --------------------------------------------------------------------
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

  // --------------------------------------------------------------------
  // 3. RAPID DESCENT BFS & REVERSE MAP
  // --------------------------------------------------------------------
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

  // --------------------------------------------------------------------
  // 4. BATCHED CROSS INTERSECTION LOOP V2
  // --------------------------------------------------------------------
  if(m_meshAcpu && m_meshBcpu) {
    finalCandidatePairs = executeBatchedCrossIntersectionLoopV2(
        *m_meshAcpu, *m_meshBcpu, batchMultiplier, totalBatches, m_dOutPairsA, m_dOutPairsB, m_dReverseMapB,
        m_dMarkedNodeIndicesB, m_dOutOffsetsB, m_dOutPrimsFlatB, m_dNodeDescendantCountsB, m_hOutMarkedCountB, m_bvhA,
        m_dMeshA, m_dMeshB, m_dMeshMetricsA, m_dMeshMetricsB, finalExactPairs, tracker, m_stream);
  }

  double tPipelineEnd = cuBQL::getCurrentTime();

  // --------------------------------------------------------------------
  // 5. RECORD METRICS INTO EXECUTION STATS
  // --------------------------------------------------------------------
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

void KernelBVHController::cleanup() {
  // Explicit stream sync before deallocating memory
  if(m_stream)
    CUBQL_CUDA_CALL(StreamSynchronize(m_stream));

  // Free device memory allocated via CUDA
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

  // =========================================================================
  // FREE PERSISTENT RAW BUFFERS FOR MESH B RE-ASSEMBLY
  // =========================================================================
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

  // Release cuBQL BVH structures
  cuBQL::cuda::free(m_bvhA, m_stream, m_memResource);
  cuBQL::cuda::free(m_bvhB, m_stream, m_memResource);

  // Reclaim thrust memory correctly
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

  m_origPointsB.clear();
  m_origPointsB.shrink_to_fit();

  shiftX = 0.0f;
  shiftY = 0.0f;
  shiftZ = 0.0f;

  if(m_stream)
    CUBQL_CUDA_CALL(StreamSynchronize(m_stream));
}