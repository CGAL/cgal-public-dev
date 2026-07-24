#pragma once

#include "cuBQL/bvh.h"
#include "cuBQL/builder/cuda.h"
#include <thrust/device_vector.h>
#include <vector>
#include <tbb/concurrent_vector.h>

#include "samples/common/loadOBJ.h"
#include "../CPU/CgalDefinitions.h"
#include "../testBVH/ExecutionStats.h"


class KernelBVHController {
public:
    KernelBVHController();
    ~KernelBVHController();

    // Prevent accidental copying (CUDA resources / thrust vectors shouldn't be shallow-copied)
    KernelBVHController(const KernelBVHController&) = delete;
    KernelBVHController& operator=(const KernelBVHController&) = delete;

    // 1. Setup & Construction (Run Once)
    void construct(
        Mesh& meshAcpu, Mesh& meshBcpu,
        const float3* hVertsA, int numVertsA, const uint3* hIndicesA, const float* hVertErrorsA, int numTrianglesA, int levelA,
        const float3* hVertsB, int numVertsB, const uint3* hIndicesB, const float* hVertErrorsB, int numTrianglesB, int levelB,
        int leafThreshold, ExecutionStats& stats);

    // 2. Execution Pipeline (Safe to run multiple times)
    void runIntersectionPipeline(
        int batchMultiplier, int mode, 
        tbb::concurrent_vector<int2>& finalExactPairs, ExecutionStats& stats);

    // 3. Deallocates all GPU resources safely
    void cleanup();

    // Updates translation relative to original position recorded at construct()
    void setTranslation(float xB, float yB, float zB);

private:
    // Configuration Parameters
    int m_leafThreshold = 0;
    int m_levelA = 0;
    int m_levelB = 0;
    int m_numTrianglesA = 0;
    int m_numTrianglesB = 0;
    int m_numVertsA = 0;
    int m_numVertsB = 0;

    // Stores total applied translation offset relative to original mesh state
    float shiftX = 0.0f;
    float shiftY = 0.0f;
    float shiftZ = 0.0f;

    Mesh* m_meshAcpu = nullptr;
    Mesh* m_meshBcpu = nullptr;

    // CUDA Streams & Memory Resource
    cudaStream_t m_stream = nullptr;
    cuBQL::DeviceMemoryResource m_memResource;

    // Core Geometry & Topology Persistent Device Buffers
    cuBQL::Triangle* m_dMeshA = nullptr;
    float2*          m_dMeshMetricsA = nullptr;
    cuBQL::box3f*    m_dBoxesA = nullptr;

    cuBQL::Triangle* m_dMeshB = nullptr;
    float2*          m_dMeshMetricsB = nullptr;
    cuBQL::box3f*    m_dBoxesB = nullptr;

    // Persistent Raw GPU Buffers for Vertex-Level Shifting & Re-Assembly
    float3*          m_dVertsB = nullptr;       // Active transformed vertices
    float3*          m_dVertsBOrig = nullptr;   // Pristine baseline vertices (shift target)
    uint3*           m_dIndicesB = nullptr;     // Triangle vertex indices
    float*           m_dVertErrorsB = nullptr;  // Downcast/drift precision error bounds

    // Pristine baseline CGAL host points (for zero-accumulation CPU translations)
    std::vector<Point3> m_origPointsB;

    // BVH Structures (cuBQL v2_2 format)
    cuBQL::bvh3f m_bvhA;
    cuBQL::bvh3f m_bvhB;

    uint32_t m_hOutMarkedCountA_Full = 0;
    uint32_t m_hOutMarkedCountB_Full = 0;
    thrust::device_vector<uint32_t> m_dMarkedNodeIndicesA_Full;
    thrust::device_vector<uint32_t> m_dMarkedNodeIndicesB_Full;

    // Persistent Target & Extraction Counts
    //uint32_t m_hOutMarkedCountA = 0;
   // uint32_t m_hOutMarkedCountB = 0;


    // Persistent Thrust Vectors
   // thrust::device_vector<uint32_t> m_dMarkedNodeIndicesA;
    thrust::device_vector<uint32_t> m_dNodeDescendantCountsA;
   // thrust::device_vector<uint32_t> m_dMarkedNodeIndicesB;
    thrust::device_vector<uint32_t> m_dNodeDescendantCountsB;

    thrust::device_vector<uint32_t> m_dReverseMapB;
    thrust::device_vector<uint32_t> m_dOutPairsA;
    thrust::device_vector<uint32_t> m_dOutPairsB;

    thrust::device_vector<uint32_t> m_dOutOffsetsA;
    thrust::device_vector<uint32_t> m_dOutPrimsFlatA;
    thrust::device_vector<uint32_t> m_dOutOffsetsB;
    thrust::device_vector<uint32_t> m_dOutPrimsFlatB;
};