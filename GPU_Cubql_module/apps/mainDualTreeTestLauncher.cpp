#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <algorithm>

// --- CUDA VECTOR TYPES INCLUDE ---
#include <vector_types.h> // Pulls in int2 definition cleanly for the host compiler

// --- cuBQL INCLUDES ---
#include "include/loadOFFCGAL.h"
#include "include/loadOFF.h"

// --- CGAL INCLUDES ---
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/intersections.h>

// --- TBB INCLUDES ----
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/concurrent_vector.h>

// --- PIPELINE HEADER INCLUDE ---
// Relative path to resolve your source tree layout pointing to kernelsTestBVH.h
//#include "../src/testBVH/kernelsTestBVH.h"

#include "../src/DualTreeTraversal/DualTreeLauncher.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef Kernel::Point_3 Point3;
typedef Kernel::Triangle_3 Triangle3;

// Helper function to translate geometry from CGAL storage to raw floats
cuBQL::Triangle convertCgalFaceToCuBQL(const Mesh& mesh, face_descriptor face) {
    auto conn = mesh.vertices_around_face(mesh.halfedge(face));
    auto it = conn.begin();

    Point3 pA = mesh.point(*it++);
    Point3 pB = mesh.point(*it++);
    Point3 pC = mesh.point(*it);

    cuBQL::Triangle tri;
    tri.a = {(float)pA.x(), (float)pA.y(), (float)pA.z()};
    tri.b = {(float)pB.x(), (float)pB.y(), (float)pB.z()};
    tri.c = {(float)pC.x(), (float)pC.y(), (float)pC.z()};
    return tri;
}

// Reusable translation pipeline helper
std::vector<cuBQL::Triangle> processMeshLayout(const Mesh& mesh, std::vector<face_descriptor>& outFaces) {
    size_t numFaces = num_faces(mesh);
    outFaces.clear();
    outFaces.reserve(numFaces);
    for (auto f : mesh.faces()) {
        outFaces.push_back(f);
    }

    std::vector<cuBQL::Triangle> hMeshLayout(numFaces);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, numFaces), [&](const tbb::blocked_range<size_t>& r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
            hMeshLayout[i] = convertCgalFaceToCuBQL(mesh, outFaces[i]);
        }
    });
    return hMeshLayout;
}


int main(int ac, char** av) {
    if (ac < 7) {
        std::cout << "Usage: " << av[0] << " <meshA.off> <maxCellSizeA> <meshB.off> <maxCellSizeB> <batchmultiplier> <intersectionRatio>\n";
        return 1;
    }

    std::string meshPathA = av[1];
    int maxCellSizeA    = std::stoi(av[2]);
    std::string meshPathB = av[3];
    int maxCellSizeB    = std::stoi(av[4]);
    int batchmultipl    = std::stoi(av[5]); 
    float intersectionRatio = std::stof(av[6]);
    
    // Default to a solid baseline configuration for thread pooling
    tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, 8);

    double tStartTotal = cuBQL::getCurrentTime();

    // Arrays to store face descriptors for indexing down into standard narrow-phase intersections
    std::vector<face_descriptor> facesA;
    std::vector<face_descriptor> facesB;

    // --------------------------------------------------------------------
    // LOAD AND TRANSLATE MESH A
    // --------------------------------------------------------------------
    double tMeshAStart = cuBQL::getCurrentTime();
    Mesh meshA;
    std::ifstream inFileA(meshPathA);
    if (!(inFileA >> meshA)) {
        std::cerr << "Critical Error: Failed to open or parse mesh file A: " << meshPathA << "\n";
        return 1;
    }
    inFileA.close();
    
    std::vector<cuBQL::Triangle> hMeshLayoutA = processMeshLayout(meshA, facesA);
    std::cout << "[Step 1] Loaded and Processed Mesh A (" << hMeshLayoutA.size() << " tris): " 
              << cuBQL::prettyDouble(cuBQL::getCurrentTime() - tMeshAStart) << "s\n";

    // --------------------------------------------------------------------
    // LOAD AND TRANSLATE MESH B
    // --------------------------------------------------------------------
    double tMeshBStart = cuBQL::getCurrentTime();
    Mesh meshB;
    std::ifstream inFileB(meshPathB);
    if (!(inFileB >> meshB)) {
        std::cerr << "Critical Error: Failed to open or parse mesh file B: " << meshPathB << "\n";
        return 1;
    }
    inFileB.close();

    std::vector<cuBQL::Triangle> hMeshLayoutB = processMeshLayout(meshB, facesB);
    std::cout << "[Step 2] Loaded and Processed Mesh B (" << hMeshLayoutB.size() << " tris): " 
              << cuBQL::prettyDouble(cuBQL::getCurrentTime() - tMeshBStart) << "s\n";

    // --------------------------------------------------------------------
    // EXECUTE SPATIAL CROSS-INTERSECTION PIPELINE
    // --------------------------------------------------------------------
    std::cout << "[Step 3] Launching Dual-Mesh GPU Pipeline...\n";
    
    ExecutionStats stats;
    std::vector<int2> hGreenPairs;
    std::vector<int2> hYellowPairs;

  //  kernelsTestBVH(
        // hMeshLayoutA.data(), static_cast<int>(hMeshLayoutA.size()), maxCellSizeA,
        // hMeshLayoutB.data(), static_cast<int>(hMeshLayoutB.size()), maxCellSizeB, 
        // batchmultipl, stats, hGreenPairs, hYellowPairs
   // );

    dualTreeKernelLaunch( hMeshLayoutA.data(), static_cast<int>(hMeshLayoutA.size()), maxCellSizeA,
        hMeshLayoutB.data(), static_cast<int>(hMeshLayoutB.size()), maxCellSizeB, 
        batchmultipl, intersectionRatio , stats, hGreenPairs, hYellowPairs);
    

    // --------------------------------------------------------------------
    // HYBRID NARROW-PHASE FILTERING LOOP (Processes Only Yellow List via TBB)
    // --------------------------------------------------------------------
    std::cout << "[Step 4] Launching Parallel CGAL Narrow-Phase Filter over Yellow Candidates...\n";
    double t_inter_start = cuBQL::getCurrentTime();

    auto& coords1 = meshA.points();
    auto& coords2 = meshB.points();

    // Populate the thread-safe collection starting with your definite green hits
    tbb::concurrent_vector<int2> finalExactPairs(hGreenPairs.begin(), hGreenPairs.end());

    if (!hYellowPairs.empty()) {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, hYellowPairs.size()), [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                const auto& pair = hYellowPairs[i];

                // .x maps to Mesh A identifier, .y maps to Mesh B identifier
                face_descriptor fA = facesA[pair.x];
                face_descriptor fB = facesB[pair.y];

                auto halfedgeA = meshA.halfedge(fA);
                auto vA0 = meshA.target(halfedgeA);
                auto vA1 = meshA.target(meshA.next(halfedgeA));
                auto vA2 = meshA.target(meshA.next(meshA.next(halfedgeA)));
                Triangle3 triA(coords1[vA0], coords1[vA1], coords1[vA2]);

                auto halfedgeB = meshB.halfedge(fB);
                auto vB0 = meshB.target(halfedgeB);
                auto vB1 = meshB.target(meshB.next(halfedgeB));
                auto vB2 = meshB.target(meshB.next(meshB.next(halfedgeB)));
                Triangle3 triB(coords2[vB0], coords2[vB1], coords2[vB2]);

                // True strict mathematical intersection evaluation pass
                if (CGAL::do_intersect(triA, triB)) {
                    finalExactPairs.push_back(pair);
                }
            }
        });
    }
    double time_intersection = cuBQL::getCurrentTime() - t_inter_start;

    // --------------------------------------------------------------------
    // FORMATTED REPORTING OF METRICS FROM EXECUTIONSTATS
    // --------------------------------------------------------------------
    std::cout << "\n==================================================\n";
    std::cout << "             PIPELINE METRICS REPORT              \n";
    std::cout << "==================================================\n";
    std::cout << "STRUCTURE SUMMARY & PROPORTIONS:\n";
    std::cout << "  |- Mesh A Total Generated Nodes:   " << stats.meshATotalNodes << "\n";
    std::cout << "  |- Mesh A Extracted Targets (<" << maxCellSizeA << "): " << stats.meshAExtractedTargets << "\n";
    std::cout << "  |- Mesh B Total Generated Nodes:   " << stats.meshBTotalNodes << "\n";
    std::cout << "  |- Mesh B Extracted Targets (<" << maxCellSizeB << "): " << stats.meshBExtractedTargets << "\n\n";

    std::cout << "CRISS-CROSS BOUNDING BOX CROSS-CHECK:\n";
    std::cout << "  |- Intersections found:            " << stats.totalIntersections << " / " << stats.totalPossiblePairs << "\n";
    std::cout << "  |- Intersection Ratio:             " << std::fixed << std::setprecision(4) << stats.intersectionPercentage << "%\n\n";

    std::cout << "TIMING METRICS OVERVIEW:\n";
    std::cout << "  |- Build + Refit (Mesh A):         " << stats.buildRefitMeshAMs << " ms\n";
    std::cout << "  |- Build + Refit (Mesh B):         " << stats.buildRefitMeshBMs << " ms\n";
    std::cout << "  |- GPU Cross-Check Engine:         " << stats.gpuCrossCheckEngineMs << " ms\n";
    std::cout << "  |- Parallel DFS Descent (B):       " << stats.parallelDfsDescentBMs << " ms\n";
    std::cout << "  |- Comprehensive GPU Pipeline Time: " << stats.GPUTotalTime << " ms\n\n";

    std::cout << "DUAL-TREE DESCENT & FINE EVALUATION METRICS:\n";
    std::cout << "  |- Total Criss-Cross Batches:      " << stats.totalCrissCrossBatches << "\n";
    std::cout << "  |- Final AABB Candidate Pairs:     " << stats.finalAabbCandidatePairs << "\n";
    std::cout << "  |- Confirmed Green (Intersecting): " << stats.confirmedGreenPairs << " (Vector Size: " << hGreenPairs.size() << ")\n";
    std::cout << "  |- Confirmed Yellow (Coplanar):    " << stats.confirmedYellowPairs << " (Vector Size: " << hYellowPairs.size() << ")\n\n";

    std::cout << "FINE-GRAINED INTERSECTION LOOP BREAKDOWN:\n";
    std::cout << "  |- Total Loop Execution Time:      " << stats.loopTracker.totalLoopTimeMs << " ms\n";
    std::cout << "  |- Preallocation Phase Time:       " << stats.loopTracker.preallocateTimeMs << " ms\n";
    std::cout << "  |- Assembly Phase Time:            " << stats.loopTracker.assemblyPhaseMs << " ms\n";
    std::cout << "  |- Execution Phase Time:           " << stats.loopTracker.executionPhaseMs << " ms\n";
    std::cout << "  |- Fine Evaluation Phase Time:     " << stats.loopTracker.fineEvaluationPhaseMs << " ms\n";
    std::cout << "  |- Device Cleanup Time:            " << stats.loopTracker.cleanupTimeMs << " ms\n";
    std::cout << "  |- Host Download & Sync Time:      " << stats.loopTracker.DownloadAndClean << " ms\n\n";

    std::cout << "NARROW-PHASE HYBRID FILTER REPORT:\n";
    std::cout << "  |- CPU Narrow-Phase Compute Time:  " << cuBQL::prettyDouble(time_intersection) << " s\n";
    std::cout << "  |- Ambiguous Yellows Evaluated:    " << hYellowPairs.size() << "\n";
    std::cout << "  |- Total Final Intersecting Pairs: " << finalExactPairs.size() << " (Green Hit Count + Validated Yellows)\n";
    std::cout << "==================================================\n";

    std::cout << "Total Diagnostic App Execution Time: " 
              << cuBQL::prettyDouble(cuBQL::getCurrentTime() - tStartTotal) << "s\n";
    std::cout << "==================================================\n";

    return 0;
}