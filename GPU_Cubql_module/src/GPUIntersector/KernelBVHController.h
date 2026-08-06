#pragma once

#include "cuBQL/bvh.h"
#include "cuBQL/builder/cuda.h"
#include <thrust/device_vector.h>
#include <vector>
#include <tbb/concurrent_vector.h>

#include "samples/common/loadOBJ.h"
#include "../CPU/CgalDefinitions.h"
#include "../testBVH/ExecutionStats.h"

// REMOVED: #include "../custom_pipeline/TriangleDouble.h"

class KernelBVHController {
public:
    KernelBVHController();
    ~KernelBVHController();

    // Prevent accidental copying
    KernelBVHController(const KernelBVHController&) = delete;
    KernelBVHController& operator=(const KernelBVHController&) = delete;

    // 1. Setup & Construction - Double Precision Only
    void construct(
        Mesh& meshAcpu, Mesh& meshBcpu,
        const Point3& centerA, const Point3& centerB,
        const double3* hVertsA, int numVertsA, const uint3* hIndicesA, int numTrianglesA, int levelA,
        const double3* hVertsB, int numVertsB, const uint3* hIndicesB, int numTrianglesB, int levelB,
        int leafThreshold, ExecutionStats& stats);

    void construct(
        Mesh& meshAcpu, Mesh& meshBcpu,
        const double3* hVertsA, int numVertsA, const uint3* hIndicesA, int numTrianglesA, int levelA,
        const double3* hVertsB, int numVertsB, const uint3* hIndicesB, int numTrianglesB, int levelB,
        int leafThreshold, ExecutionStats& stats);

    // 2. Execution Pipeline - Double Precision Only
    void runIntersectionPipeline(
        int batchMultiplier, int mode, int activateAsyncDownload,
        int2*& outFinalExactPairs,       // Fast raw pointer return
    size_t& outFinalCount,  ExecutionStats& stats,bool gpuDouble);

    // 3. Deallocates all GPU resources safely
    void cleanup();
    void clearGPU();
    void reconstructGPU(ExecutionStats& stats);

    // Dynamic Dual Point Cloud Transformation
    void setTransformBoth(double3 rotDegA, double3 transA, 
                          double3 rotDegB, double3 transB,
                          float& timeGPU, float& timeCPU,
                          float& timeTransformVerts, 
                          float& timeAssembleTris, 
                          float& timeGenBoxes);

    // Legacy Translation Interfaces
    void setTranslation(double xB, double yB, double zB);
   // void setTranslationCPUHostUpload(double xB, double yB, double zB);

    // Centroid Getters
    Point3 getCenterA() const { return m_centerA; }
    Point3 getCenterB() const { return m_centerB; }

    bool isGPUAllocated() const { 
    return (m_dVertsA != nullptr && m_dVertsB != nullptr); 
}

private:
    // Configuration Parameters
    int m_leafThreshold = 0;
    int m_levelA = 0;
    int m_levelB = 0;
    int m_numTrianglesA = 0;
    int m_numTrianglesB = 0;
    int m_numVertsA = 0;
    int m_numVertsB = 0;

    // Pre-computed Mesh Centroids
    Point3 m_centerA{0, 0, 0};
    Point3 m_centerB{0, 0, 0};

    // Active transformations relative to pristine baseline states
    double3 m_rotA{0.0, 0.0, 0.0};
    double3 m_transA{0.0, 0.0, 0.0};
    double3 m_rotB{0.0, 0.0, 0.0};
    double3 m_transB{0.0, 0.0, 0.0};

    // Legacy shift markers
    double shiftX = 0.0f;
    double shiftY = 0.0f;
    double shiftZ = 0.0f;

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

    // REMOVED: TriangleDouble* m_dMeshDoubleA
    // REMOVED: TriangleDouble* m_dMeshDoubleB

        // ADD: Cached host data needed for GPU reconstruction
    const double3* m_hVertsA = nullptr;
    const uint3*   m_hIndicesA = nullptr;
    const double3* m_hVertsB = nullptr;
    const uint3*   m_hIndicesB = nullptr;

    // Persistent Raw GPU & CPU Buffers for Mesh A
    double3*             m_dVertsA = nullptr;       // Active transformed vertices
    double3*             m_dVertsAOrig = nullptr;   // Pristine baseline vertices
    uint3*              m_dIndicesA = nullptr;     // Triangle vertex indices
    float*              m_dVertErrorsA = nullptr;  // Precision error bounds
    std::vector<Point3> m_origPointsA;             // Baseline CGAL host points

    // Persistent Raw GPU & CPU Buffers for Mesh B
    double3*              m_dVertsB = nullptr;       // Active transformed vertices
    double3*              m_dVertsBOrig = nullptr;   // Pristine baseline vertices
    uint3*              m_dIndicesB = nullptr;     // Triangle vertex indices
    float*              m_dVertErrorsB = nullptr;  // Precision error bounds
    std::vector<Point3> m_origPointsB;             // Baseline CGAL host points

    // BVH Structures (cuBQL v2_2 format)
    cuBQL::bvh3f m_bvhA;
    cuBQL::bvh3f m_bvhB;

    uint32_t m_hOutMarkedCountA_Full = 0;
    uint32_t m_hOutMarkedCountB_Full = 0;
    thrust::device_vector<uint32_t> m_dMarkedNodeIndicesA_Full;
    thrust::device_vector<uint32_t> m_dMarkedNodeIndicesB_Full;

    // Persistent Thrust Vectors
    thrust::device_vector<uint32_t> m_dNodeDescendantCountsA;
    thrust::device_vector<uint32_t> m_dNodeDescendantCountsB;

    thrust::device_vector<uint32_t> m_dReverseMapB;
    thrust::device_vector<uint32_t> m_dOutPairsA;
    thrust::device_vector<uint32_t> m_dOutPairsB;

    thrust::device_vector<uint32_t> m_dOutOffsetsA;
    thrust::device_vector<uint32_t> m_dOutPrimsFlatA;
    thrust::device_vector<uint32_t> m_dOutOffsetsB;
    thrust::device_vector<uint32_t> m_dOutPrimsFlatB;



};