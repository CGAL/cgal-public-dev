#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <utility>

// --- CUDA VECTOR TYPES INCLUDE ---
#include <vector_types.h> 

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
#include "../src/testBVH/kernelsTestBVHAlternative.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef Kernel::Point_3 Point3;
typedef Kernel::Triangle_3 Triangle3;

// Helper function to translate geometry from CGAL storage to raw floats
// Direct, non-allocating conversion bypasses proxy iterators
inline cuBQL::Triangle convertCgalFaceToCuBQL(const Mesh& mesh, face_descriptor face) {
    // Manually unroll the halfedge triangle cycle to avoid iterator tracking overhead
    auto h0 = mesh.halfedge(face);
    auto h1 = mesh.next(h0);
    auto h2 = mesh.next(h1);

    const Point3& pA = mesh.point(mesh.target(h0));
    const Point3& pB = mesh.point(mesh.target(h1));
    const Point3& pC = mesh.point(mesh.target(h2));

    cuBQL::Triangle tri;
    tri.a = {static_cast<float>(pA.x()), static_cast<float>(pA.y()), static_cast<float>(pA.z())};
    tri.b = {static_cast<float>(pB.x()), static_cast<float>(pB.y()), static_cast<float>(pB.z())};
    tri.c = {static_cast<float>(pC.x()), static_cast<float>(pC.y()), static_cast<float>(pC.z())};
    return tri;
}

// Fully parallelized pipeline helper (Zero sequential bottlenecks)
void processMeshLayoutTimed(const Mesh& mesh, std::vector<face_descriptor>& outFaces, std::vector<cuBQL::Triangle>& outLayout) {
    size_t numFaces = num_faces(mesh);
    
    // Allocate up front to eliminate sequential push_back/reallocation locks
    outFaces.resize(numFaces);
    outLayout.resize(numFaces);

    // If your mesh contains deleted faces, use CGAL::Index_property_map instead.
    // For standard contiguous OFF loads, we can map indices directly in parallel.
    tbb::parallel_for(tbb::blocked_range<size_t>(0, numFaces), [&](const tbb::blocked_range<size_t>& r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
            // Reconstruct the face descriptor directly from the linear index
            face_descriptor f(static_cast<Mesh::size_type>(i));
            outFaces[i] = f;
            outLayout[i] = convertCgalFaceToCuBQL(mesh, f);
        }
    });
}
int main(int ac, char** av) {
    if (ac < 9) { 
        std::cout << "Usage: " << av[0] << " <meshA.off> <maxCellSizeA> <meshB.off> <maxCellSizeB> <batchmultiplier> <NumberOfDualTreeSteps> <Leaf Threshold> <tbb_workers>\n";
        return 1;
    }

    std::string meshPathA = av[1];
    int maxCellSizeA    = std::stoi(av[2]);
    std::string meshPathB = av[3];
    int maxCellSizeB    = std::stoi(av[4]);
    int batchmultipl    = std::stoi(av[5]); 
    int mode            = std::stoi(av[6]);
    int leafThreshhold  = std::stoi(av[7]); 
    int tbbWorkers      = std::stoi(av[8]); 
    
    tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, tbbWorkers);

    double tStartTotal = cuBQL::getCurrentTime();

    std::vector<face_descriptor> facesA;
    std::vector<face_descriptor> facesB;
    std::vector<cuBQL::Triangle> hMeshLayoutA;
    std::vector<cuBQL::Triangle> hMeshLayoutB;

    // --------------------------------------------------------------------
    // LOAD AND TRANSLATE MESH A
    // --------------------------------------------------------------------
    double tMeshAFileLoadStart = cuBQL::getCurrentTime();
    Mesh meshA;
    std::ifstream inFileA(meshPathA);
    if (!(inFileA >> meshA)) {
        std::cerr << "Critical Error: Failed to open or parse mesh file A: " << meshPathA << "\n";
        return 1;
    }
    inFileA.close();
    double tMeshAFileLoadTime = cuBQL::getCurrentTime() - tMeshAFileLoadStart;
    
    // Timer moved outside to precisely profile the full function duration
    double tMeshATranslationStart = cuBQL::getCurrentTime();
    processMeshLayoutTimed(meshA, facesA, hMeshLayoutA);
    double tMeshATranslationTime = cuBQL::getCurrentTime() - tMeshATranslationStart;

    // --------------------------------------------------------------------
    // LOAD AND TRANSLATE MESH B
    // --------------------------------------------------------------------
    double tMeshBFileLoadStart = cuBQL::getCurrentTime();
    Mesh meshB;
    std::ifstream inFileB(meshPathB);
    if (!(inFileB >> meshB)) {
        std::cerr << "Critical Error: Failed to open or parse mesh file B: " << meshPathB << "\n";
        return 1;
    }
    inFileB.close();
    double tMeshBFileLoadTime = cuBQL::getCurrentTime() - tMeshBFileLoadStart;

    // Timer moved outside to precisely profile the full function duration
    double tMeshBTranslationStart = cuBQL::getCurrentTime();
    processMeshLayoutTimed(meshB, facesB, hMeshLayoutB);
    double tMeshBTranslationTime = cuBQL::getCurrentTime() - tMeshBTranslationStart;

    // --------------------------------------------------------------------
    // EXECUTE SPATIAL CROSS-INTERSECTION PIPELINE
    // --------------------------------------------------------------------
    std::cout << "[Pipeline] Running Dual-Mesh GPU pipeline execution...\n";
    
    ExecutionStats stats;
    std::vector<int2> hGreenPairs;
    std::vector<int2> hYellowPairs;

    kernelsTestBVHV2(
        hMeshLayoutA.data(), static_cast<int>(hMeshLayoutA.size()), maxCellSizeA,
        hMeshLayoutB.data(), static_cast<int>(hMeshLayoutB.size()), maxCellSizeB, 
        batchmultipl, mode, leafThreshhold, stats, hGreenPairs, hYellowPairs 
    );

    // --------------------------------------------------------------------
    // HYBRID NARROW-PHASE FILTERING LOOP (Processes Only Yellow List via TBB)
    // --------------------------------------------------------------------
    double t_inter_start = cuBQL::getCurrentTime();

    auto& coords1 = meshA.points();
    auto& coords2 = meshB.points();

    tbb::concurrent_vector<int2> finalExactPairs(hGreenPairs.begin(), hGreenPairs.end());

    if (!hYellowPairs.empty()) {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, hYellowPairs.size()), [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                const auto& pair = hYellowPairs[i];

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
    std::cout << "  |- Pipeline Execution Mode:        " << mode << "\n";
    std::cout << "  |- TBB Max Parallelism Limit:      " << tbbWorkers << "\n";
    std::cout << "  |- Mesh A Total Generated Nodes:   " << stats.meshATotalNodes << "\n";
    std::cout << "  |- Mesh A Extracted Targets (<" << maxCellSizeA << "): " << stats.meshAExtractedTargets << "\n";
    std::cout << "  |- Mesh B Total Generated Nodes:   " << stats.meshBTotalNodes << "\n";
    std::cout << "  |- Mesh B Extracted Targets (<" << maxCellSizeB << "): " << stats.meshBExtractedTargets << "\n\n";

    std::cout << "HOST-SIDE INGESTION & TRANSLATION BREAKDOWN:\n";
    std::cout << "  |- Mesh A Disk Ingest (OFF Read):  " << std::fixed << std::setprecision(4) << tMeshAFileLoadTime << " s\n";
    std::cout << "  |- Mesh A Complete Host Translate: " << tMeshATranslationTime << " s\n";
    std::cout << "  |- Mesh B Disk Ingest (OFF Read):  " << tMeshBFileLoadTime << " s\n";
    std::cout << "  |- Mesh B Complete Host Translate: " << tMeshBTranslationTime << " s\n\n";

    std::cout << "CRISS-CROSS BOUNDING BOX CROSS-CHECK:\n";
    std::cout << "  |- Intersections found:            " << stats.totalIntersections << " / " << stats.totalPossiblePairs << "\n";
    std::cout << "  |- Intersection Ratio:             " << std::fixed << std::setprecision(4) << stats.intersectionPercentage << "%\n\n";

    std::cout << "TIMING METRICS OVERVIEW:\n";
    std::cout << "  |- Initial Alloc & Mesh Copy:      " << std::fixed << std::setprecision(4) << stats.initialAllocAndCopyMs << " ms\n";
    std::cout << "  |- Thrust Framework Init/Fill:     " << stats.thrustInitOverheadMs << " ms\n";
    std::cout << "  |- Build + Refit (Mesh A):         " << stats.buildRefitMeshAMs << " ms\n";
    std::cout << "  |- Build + Refit (Mesh B):         " << stats.buildRefitMeshBMs << " ms\n";
    std::cout << "  |- GPU Cross-Check Engine:         " << stats.gpuCrossCheckEngineMs << " ms\n";
    std::cout << "  |- Parallel DFS Descent (A & B):   " << stats.parallelDfsDescentBMs << " ms\n"; 
    std::cout << "  |- Dual Tree Expansion Step:       " << stats.dualTreeStepMs << " ms\n"; 
    std::cout << "  |- Explicit Device Cleanup Sync:   " << stats.finalCleanupSyncMs << " ms\n";
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