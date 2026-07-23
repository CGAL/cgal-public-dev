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
    float3* dVertsOut, 
    const float3* dVertsOrig, 
    int numVerts, 
    float targetX, float targetY, float targetZ) 
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= numVerts) return;

    float3 v = dVertsOrig[idx];
    dVertsOut[idx] = make_float3(v.x + targetX, v.y + targetY, v.z + targetZ);
}

void translateCgalMeshB(Mesh* mesh, float dx, float dy, float dz) {
  if (!mesh) return;

  // Extract underlying point property map
  auto pmap = mesh->points();
  const size_t numPoints = mesh->num_vertices();

  tbb::parallel_for(tbb::blocked_range<size_t>(0, numPoints),
    [mesh, &pmap, dx, dy, dz](const tbb::blocked_range<size_t>& range) {
      for (size_t i = range.begin(); i != range.end(); ++i) {
        Mesh::Vertex_index vd(static_cast<uint32_t>(i));
        if (!mesh->is_removed(vd)) {
          const auto& p = pmap[vd];
          pmap[vd] = Point3(p.x() + dx, p.y() + dy, p.z() + dz);
        }
      }
    });
}

void KernelBVHController::setTranslation(float xB, float yB, float zB) {
  // 1. Quick exit if offset hasn't changed
  if (xB == shiftX && yB == shiftY && zB == shiftZ) return;

  bool gpuWorkQueued = false;

  // 2. LAUNCH GPU PIPELINE (Asynchronous on m_stream)
  if (m_numVertsB > 0 && m_dVertsB && m_dVertsBOrig) {
    int block = 256;

    // Step A: Shift raw vertices directly from pristine baseline
    int gridVerts = cuBQL::divRoundUp(m_numVertsB, block);
    translateVerticesKernel<<<gridVerts, block, 0, m_stream>>>(
        m_dVertsB, m_dVertsBOrig, m_numVertsB, xB, yB, zB);

    // Step B: Re-run assembly & metric recalculation on updated vertices
    int gridTris = cuBQL::divRoundUp(m_numTrianglesB, block);
    assembleTrianglesKernel<<<gridTris, block, 0, m_stream>>>(
        m_dMeshB, m_dMeshMetricsB, m_dVertsB, m_dVertErrorsB, m_dIndicesB, m_numTrianglesB);

    // Step C: Re-generate AABBs from assembled triangles
    generateBoxes<<<gridTris, block, 0, m_stream>>>(
        m_dBoxesB, m_dMeshB, m_numTrianglesB);

    // Step D: Refit BVH tree
    cuBQL::cuda::refit(m_bvhB, m_dBoxesB, m_stream);

    gpuWorkQueued = true;
  }

  // 3. CPU CGAL SHIFT (Runs concurrently on CPU using TBB)
  if (m_meshBcpu) {
    float dx = xB - shiftX;
    float dy = yB - shiftY;
    float dz = zB - shiftZ;
    translateCgalMeshB(m_meshBcpu, dx, dy, dz);
  }

  // 4. SYNCHRONIZE GPU WORK
  // FIX: Always sync m_stream if GPU kernels were queued, regardless of m_numVertsB checks,
  // guaranteeing the BVH refit is complete before any downstream CPU/GPU intersection call.
  if (gpuWorkQueued) {
    CUBQL_CUDA_CALL(StreamSynchronize(m_stream));
  }

  // Update tracked shift coordinates only after successful execution
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

  // IMPORTANT: Clean up any old allocations if construct() is called again
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

  // Reset shift trackers for newly constructed mesh pair
  shiftX = 0.0f;
  shiftY = 0.0f;
  shiftZ = 0.0f;

  cuBQL::BuildConfig buildConfig;
  buildConfig.makeLeafThreshold = m_leafThreshold;

  // --------------------------------------------------------------------
  // 1. ALLOCATION AND MESH ASSEMBLY (DEVICE STAGING)
  // --------------------------------------------------------------------
  double tAllocStart = cuBQL::getCurrentTime();

  // Mesh A: Temporary staging buffers
  float3* dVertsA = nullptr;
  uint3* dIndicesA = nullptr;
  float* dVertErrorsA = nullptr;

  CUBQL_CUDA_CALL(MallocAsync((void**)&dVertsA, m_numVertsA * sizeof(float3), m_stream));
  CUBQL_CUDA_CALL(MemcpyAsync(dVertsA, hVertsA, m_numVertsA * sizeof(float3), cudaMemcpyHostToDevice, m_stream));

  CUBQL_CUDA_CALL(MallocAsync((void**)&dIndicesA, m_numTrianglesA * sizeof(uint3), m_stream));
  CUBQL_CUDA_CALL(MemcpyAsync(dIndicesA, hIndicesA, m_numTrianglesA * sizeof(uint3), cudaMemcpyHostToDevice, m_stream));

  if(hVertErrorsA) {
    CUBQL_CUDA_CALL(MallocAsync((void**)&dVertErrorsA, m_numVertsA * sizeof(float), m_stream));
    CUBQL_CUDA_CALL(
        MemcpyAsync(dVertErrorsA, hVertErrorsA, m_numVertsA * sizeof(float), cudaMemcpyHostToDevice, m_stream));
  }

  // Mesh B: Persistent member buffers
  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dVertsB, m_numVertsB * sizeof(float3), m_stream));
  CUBQL_CUDA_CALL(MemcpyAsync(m_dVertsB, hVertsB, m_numVertsB * sizeof(float3), cudaMemcpyHostToDevice, m_stream));

  // Store pristine baseline vertices
  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dVertsBOrig, m_numVertsB * sizeof(float3), m_stream));
  CUBQL_CUDA_CALL(MemcpyAsync(m_dVertsBOrig, hVertsB, m_numVertsB * sizeof(float3), cudaMemcpyHostToDevice, m_stream));

  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dIndicesB, m_numTrianglesB * sizeof(uint3), m_stream));
  CUBQL_CUDA_CALL(MemcpyAsync(m_dIndicesB, hIndicesB, m_numTrianglesB * sizeof(uint3), cudaMemcpyHostToDevice, m_stream));

  // FIX 1: Nullify persistent pointer before conditional allocation!
  m_dVertErrorsB = nullptr;
  if(hVertErrorsB) {
    CUBQL_CUDA_CALL(MallocAsync((void**)&m_dVertErrorsB, m_numVertsB * sizeof(float), m_stream));
    CUBQL_CUDA_CALL(
        MemcpyAsync(m_dVertErrorsB, hVertErrorsB, m_numVertsB * sizeof(float), cudaMemcpyHostToDevice, m_stream));
  }

  // Allocate persistent triangle & box buffers
  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dMeshA, m_numTrianglesA * sizeof(cuBQL::Triangle), m_stream));
  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dMeshMetricsA, m_numTrianglesA * sizeof(float2), m_stream));
  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dBoxesA, m_numTrianglesA * sizeof(cuBQL::box3f), m_stream));

  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dMeshB, m_numTrianglesB * sizeof(cuBQL::Triangle), m_stream));
  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dMeshMetricsB, m_numTrianglesB * sizeof(float2), m_stream));
  CUBQL_CUDA_CALL(MallocAsync((void**)&m_dBoxesB, m_numTrianglesB * sizeof(cuBQL::box3f), m_stream));

  // Initial Assembly
  assembleTrianglesKernel<<<cuBQL::divRoundUp(m_numTrianglesA, 256), 256, 0, m_stream>>>(
      m_dMeshA, m_dMeshMetricsA, dVertsA, dVertErrorsA, dIndicesA, m_numTrianglesA);
  generateBoxes<<<cuBQL::divRoundUp(m_numTrianglesA, 256), 256, 0, m_stream>>>(m_dBoxesA, m_dMeshA, m_numTrianglesA);

  assembleTrianglesKernel<<<cuBQL::divRoundUp(m_numTrianglesB, 256), 256, 0, m_stream>>>(
      m_dMeshB, m_dMeshMetricsB, m_dVertsB, m_dVertErrorsB, m_dIndicesB, m_numTrianglesB);
  generateBoxes<<<cuBQL::divRoundUp(m_numTrianglesB, 256), 256, 0, m_stream>>>(m_dBoxesB, m_dMeshB, m_numTrianglesB);

  CUBQL_CUDA_CALL(StreamSynchronize(m_stream));
  stats.initialAllocAndCopyMs = (cuBQL::getCurrentTime() - tAllocStart) * 1000.0;

  // Free Mesh A's temporary staging buffers
  CUBQL_CUDA_CALL(FreeAsync(dVertsA, m_stream));
  CUBQL_CUDA_CALL(FreeAsync(dIndicesA, m_stream));
  if(dVertErrorsA) {
    CUBQL_CUDA_CALL(FreeAsync(dVertErrorsA, m_stream));
  }

  // --------------------------------------------------------------------
  // 2. THRUST DEVICE VECTOR INITIALIZATION
  // --------------------------------------------------------------------
  double tThrustInitStart = cuBQL::getCurrentTime();

  uint32_t maxPossibleNodesA = 2 * m_numTrianglesA;
  uint32_t maxPossibleNodesB = 2 * m_numTrianglesB;

  m_dMarkedNodeIndicesA.resize(maxPossibleNodesA, 0);
  m_dNodeDescendantCountsA.resize(maxPossibleNodesA, 0);
  m_dMarkedNodeIndicesB.resize(maxPossibleNodesB, 0);
  m_dNodeDescendantCountsB.resize(maxPossibleNodesB, 0);
  m_dReverseMapB.resize(maxPossibleNodesB, 0);

  // FIX 2: Ensure stream synchronization when Thrust resizes memory under the default stream/policy
  CUBQL_CUDA_CALL(StreamSynchronize(m_stream));
  stats.thrustInitOverheadMs = (cuBQL::getCurrentTime() - tThrustInitStart) * 1000.0;

  // --------------------------------------------------------------------
  // 3. BUILD & REFIT BVHs
  // --------------------------------------------------------------------
  double tInitAStart = cuBQL::getCurrentTime();

  cuBQL::gpuBuilder_v2_2::build_custom(m_bvhA, m_dBoxesA, m_numTrianglesA, buildConfig, (uint32_t)m_levelA,
                                       thrust::raw_pointer_cast(m_dMarkedNodeIndicesA.data()),
                                       thrust::raw_pointer_cast(m_dNodeDescendantCountsA.data()), &m_hOutMarkedCountA,
                                       m_stream, m_memResource);

  cuBQL::cuda::refit(m_bvhA, m_dBoxesA, m_stream);

  CUBQL_CUDA_CALL(StreamSynchronize(m_stream));
  stats.buildRefitMeshAMs = (cuBQL::getCurrentTime() - tInitAStart) * 1000.0;

  double tInitBStart = cuBQL::getCurrentTime();

  cuBQL::gpuBuilder_v2_2::build_custom(m_bvhB, m_dBoxesB, m_numTrianglesB, buildConfig, (uint32_t)m_levelB,
                                       thrust::raw_pointer_cast(m_dMarkedNodeIndicesB.data()),
                                       thrust::raw_pointer_cast(m_dNodeDescendantCountsB.data()), &m_hOutMarkedCountB,
                                       m_stream, m_memResource);

  cuBQL::cuda::refit(m_bvhB, m_dBoxesB, m_stream);

  CUBQL_CUDA_CALL(StreamSynchronize(m_stream));
  stats.buildRefitMeshBMs = (cuBQL::getCurrentTime() - tInitBStart) * 1000.0;
}

void KernelBVHController::runIntersectionPipeline(int batchMultiplier,
                                                  int mode,
                                                  tbb::concurrent_vector<int2>& finalExactPairs,
                                                  ExecutionStats& stats) {
  double tPipelineStart = cuBQL::getCurrentTime();

  // RESET STALE RUN DATA FOR MULTIPLE RUN SAFETY
  m_dOutPairsA.clear();
  m_dOutPairsB.clear();

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
  m_dMarkedNodeIndicesA = thrust::device_vector<uint32_t>();
  m_dNodeDescendantCountsA = thrust::device_vector<uint32_t>();
  m_dMarkedNodeIndicesB = thrust::device_vector<uint32_t>();
  m_dNodeDescendantCountsB = thrust::device_vector<uint32_t>();

  m_dReverseMapB = thrust::device_vector<uint32_t>();
  m_dOutPairsA = thrust::device_vector<uint32_t>();
  m_dOutPairsB = thrust::device_vector<uint32_t>();
  m_dOutOffsetsB = thrust::device_vector<uint32_t>();
  m_dOutPrimsFlatB = thrust::device_vector<uint32_t>();
  m_dOutOffsetsA = thrust::device_vector<uint32_t>();
  m_dOutPrimsFlatA = thrust::device_vector<uint32_t>();

  shiftX = 0.0f;
  shiftY = 0.0f;
  shiftZ = 0.0f;

  if(m_stream)
    CUBQL_CUDA_CALL(StreamSynchronize(m_stream));
}