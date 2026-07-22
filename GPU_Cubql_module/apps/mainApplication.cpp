#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <chrono>

// CUDA & TBB Includes
#include <vector_types.h>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>

// CGAL & cuBQL Includes
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include "cuBQL/bvh.h"

// Custom Controller and Helpers
#include "../src/GPUIntersector/KernelBVHController.h"
#include "../src/CPU/CgalDefinitions.h"
#include "../src/Warmup/cuda_warmup.h"

// --- Minimal Double Box for Mesh Bounds ---
struct DoubleBox {
    double min_x = std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double min_z = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();
    double max_z = -std::numeric_limits<double>::infinity();

    void grow(double x, double y, double z) {
        min_x = std::min(min_x, x); max_x = std::max(max_x, x);
        min_y = std::min(min_y, y); max_y = std::max(max_y, y);
        min_z = std::min(min_z, z); max_z = std::max(max_z, z);
    }

    void grow(const DoubleBox& o) {
        min_x = std::min(min_x, o.min_x); max_x = std::max(max_x, o.max_x);
        min_y = std::min(min_y, o.min_y); max_y = std::max(max_y, o.max_y);
        min_z = std::min(min_z, o.min_z); max_z = std::max(max_z, o.max_z);
    }
};

DoubleBox computeMeshBoundsTBB(const Mesh& mesh) {
    size_t numVerts = num_vertices(mesh);
    auto coords = mesh.points();

    return tbb::parallel_reduce(
        tbb::blocked_range<size_t>(0, numVerts), DoubleBox(),
        [&](const tbb::blocked_range<size_t>& r, DoubleBox box) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                typename Mesh::Vertex_index v(static_cast<Mesh::size_type>(i));
                const auto& p = coords[v];
                box.grow(p.x(), p.y(), p.z());
            }
            return box;
        },
        [](DoubleBox a, const DoubleBox& b) { a.grow(b); return a; }
    );
}

void extractMeshTopologyNormalized(const Mesh& mesh, double cx, double cy, double cz,
                                   std::vector<float3>& outVerts,
                                   std::vector<float>& outVertexErrors,
                                   std::vector<uint3>& outIndices,
                                   bool usePreciseBounds) {
    size_t numVerts = num_vertices(mesh);
    size_t numFaces = num_faces(mesh);

    outVerts.resize(numVerts);
    outIndices.resize(numFaces);
    if (usePreciseBounds) outVertexErrors.resize(numVerts);
    else outVertexErrors.clear();

    auto coords = mesh.points();

    tbb::parallel_for(tbb::blocked_range<size_t>(0, numVerts), [&](const tbb::blocked_range<size_t>& r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
            typename Mesh::Vertex_index v(static_cast<Mesh::size_type>(i));
            const auto& p = coords[v];

            double sx = p.x() - cx, sy = p.y() - cy, sz = p.z() - cz;
            float fx = static_cast<float>(sx);
            float fy = static_cast<float>(sy);
            float fz = static_cast<float>(sz);

            outVerts[i] = {fx, fy, fz};

            if (usePreciseBounds) {
                double err_x = std::abs(sx - static_cast<double>(fx));
                double err_y = std::abs(sy - static_cast<double>(fy));
                double err_z = std::abs(sz - static_cast<double>(fz));
                double max_err = std::max({err_x, err_y, err_z});
                float max_mag = std::max({std::abs(fx), std::abs(fy), std::abs(fz)});
                outVertexErrors[i] = static_cast<float>(max_err) + max_mag * std::numeric_limits<float>::epsilon();
            }
        }
    });

    tbb::parallel_for(tbb::blocked_range<size_t>(0, numFaces), [&](const tbb::blocked_range<size_t>& r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
            typename Mesh::Face_index f(static_cast<Mesh::size_type>(i));
            auto h0 = mesh.halfedge(f);
            auto h1 = mesh.next(h0);
            auto h2 = mesh.next(h1);
            outIndices[i] = {static_cast<uint32_t>(mesh.target(h0)), 
                            static_cast<uint32_t>(mesh.target(h1)), 
                            static_cast<uint32_t>(mesh.target(h2))};
        }
    });
}

int main(int ac, char** av) {
    if (ac < 11) {
        std::cout << "Usage: " << av[0]
                  << " <meshA.off> <maxCellSizeA> <meshB.off> <maxCellSizeB> <batchmultiplier> <mode> "
                     "<leafThreshold> <tbbWorkers> <usePreciseBounds> <AsyncDownload>\n";
        return 1;
    }

    // 1. Warm up CUDA Context
    warmupCUDA();

    // Parse Arguments
    std::string meshPathA = av[1];
    int maxCellSizeA     = std::stoi(av[2]);
    std::string meshPathB = av[3];
    int maxCellSizeB     = std::stoi(av[4]);
    int batchMultiplier   = std::stoi(av[5]);
    int mode              = std::stoi(av[6]);
    int leafThreshold     = std::stoi(av[7]);
    int tbbWorkers        = std::stoi(av[8]);
    bool usePreciseBounds = (std::stoi(av[9]) != 0);

    tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, tbbWorkers);

    // 2. Load Meshes from Disk
    Mesh meshA, meshB;
    std::ifstream inFileA(meshPathA), inFileB(meshPathB);
    if (!(inFileA >> meshA) || !(inFileB >> meshB)) {
        std::cerr << "Error: Failed to load OFF input files.\n";
        return 1;
    }

    // 3. Centralize & Extract Topology
    DoubleBox box = computeMeshBoundsTBB(meshA);
    box.grow(computeMeshBoundsTBB(meshB));
    double cx = 0.5 * (box.min_x + box.max_x);
    double cy = 0.5 * (box.min_y + box.max_y);
    double cz = 0.5 * (box.min_z + box.max_z);

    std::vector<float3> hVertsA, hVertsB;
    std::vector<uint3> hIndicesA, hIndicesB;
    std::vector<float> hErrorsA, hErrorsB;

    extractMeshTopologyNormalized(meshA, cx, cy, cz, hVertsA, hErrorsA, hIndicesA, usePreciseBounds);
    extractMeshTopologyNormalized(meshB, cx, cy, cz, hVertsB, hErrorsB, hIndicesB, usePreciseBounds);

    KernelBVHController controller;
    ExecutionStats stats;

    // --------------------------------------------------------------------
    // MEASUREMENT 1: BVH CONSTRUCTION TIME
    // --------------------------------------------------------------------
    auto tConstructStart = std::chrono::high_resolution_clock::now();

    controller.construct(
        meshA, meshB,
        hVertsA.data(), static_cast<int>(hVertsA.size()), hIndicesA.data(),
        (usePreciseBounds ? hErrorsA.data() : nullptr), static_cast<int>(hIndicesA.size()), maxCellSizeA,
        hVertsB.data(), static_cast<int>(hVertsB.size()), hIndicesB.data(),
        (usePreciseBounds ? hErrorsB.data() : nullptr), static_cast<int>(hIndicesB.size()), maxCellSizeB,
        leafThreshold, stats
    );

    auto tConstructEnd = std::chrono::high_resolution_clock::now();
    double constructionTimeMs = std::chrono::duration<double, std::milli>(tConstructEnd - tConstructStart).count();

    // --------------------------------------------------------------------
    // MEASUREMENT 2: INTERSECTION PIPELINE TIME & PAIR COUNT
    // --------------------------------------------------------------------
    tbb::concurrent_vector<int2> finalExactPairs;

    auto tExecutionStart = std::chrono::high_resolution_clock::now();

    controller.runIntersectionPipeline(batchMultiplier, mode, finalExactPairs, stats);

    auto tExecutionEnd = std::chrono::high_resolution_clock::now();
    double executionTimeMs = std::chrono::duration<double, std::milli>(tExecutionEnd - tExecutionStart).count();

    // Free Controller Resources
    controller.cleanup();

    // --------------------------------------------------------------------
    // FINAL MINIMAL REPORTING
    // --------------------------------------------------------------------
    std::cout << "\n==================================================\n";
    std::cout << "1. BVH Construction Time: " << constructionTimeMs << " ms\n";
    std::cout << "2. Intersection Execution Time: " << executionTimeMs << " ms\n";
    std::cout << "3. Final Intersecting Pairs Count: " << finalExactPairs.size() << "\n";
    std::cout << "==================================================\n";

    return 0;
}