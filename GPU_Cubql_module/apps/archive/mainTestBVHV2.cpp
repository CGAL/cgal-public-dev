#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <utility>
#include <limits> // For numeric limits

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
#include <tbb/parallel_reduce.h> // Added for parallel reduction box computation
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>

// --- PIPELINE HEADER INCLUDE ---
#include "../src/testBVH/kernelsTestBVHAlternative.h"
#include "../src/CPU/YellowFilter.h"
#include "../src/CPU/CgalDefinitions.h"


#include "../src/Warmup/cuda_warmup.h"

// typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
// typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
// typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
// typedef Kernel::Point_3 Point3;
// typedef Kernel::Triangle_3 Triangle3;

// High-precision double bounding box structure for TBB reductions
struct DoubleBox
{
  double min_x = std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double min_z = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();
  double max_z = -std::numeric_limits<double>::infinity();

  void grow(double x, double y, double z) {
    if(x < min_x)
      min_x = x;
    if(x > max_x)
      max_x = x;
    if(y < min_y)
      min_y = y;
    if(y > max_y)
      max_y = y;
    if(z < min_z)
      min_z = z;
    if(z > max_z)
      max_z = z;
  }

  void grow(const DoubleBox& other) {
    if(other.min_x < min_x)
      min_x = other.min_x;
    if(other.max_x > max_x)
      max_x = other.max_x;
    if(other.min_y < min_y)
      min_y = other.min_y;
    if(other.max_y > max_y)
      max_y = other.max_y;
    if(other.min_z < min_z)
      min_z = other.min_z;
    if(other.max_z > max_z)
      max_z = other.max_z;
  }
};

// TBB Parallel Reduction worker to find double-precision bounds of a CGAL mesh
DoubleBox computeMeshBoundsTBB(const Mesh& mesh) {
  size_t numVerts = num_vertices(mesh);
  auto coords = mesh.points();

  return tbb::parallel_reduce(
      tbb::blocked_range<size_t>(0, numVerts), DoubleBox(),
      [&](const tbb::blocked_range<size_t>& r, DoubleBox localBox) -> DoubleBox {
        for(size_t i = r.begin(); i != r.end(); ++i) {
          typename Mesh::Vertex_index v(static_cast<Mesh::size_type>(i));
          const auto& p = coords[v];
          localBox.grow(p.x(), p.y(), p.z());
        }
        return localBox;
      },
      [](DoubleBox a, const DoubleBox& b) -> DoubleBox {
        a.grow(b);
        return a;
      });
}

// MODIFIED TRANSLATION: Applies double-precision local origin shift BEFORE casting to float
// UPGRADED TRANSLATION: Tracks precise float rounding error per vertex (dynamically optional)
void extractMeshTopologyNormalized(const Mesh& mesh,
                                   double cx,
                                   double cy,
                                   double cz,
                                   std::vector<float3>& outVerts,
                                   std::vector<float>& outVertexErrors, // Tracks L_inf error bounds
                                   std::vector<uint3>& outIndices,
                                   bool usePreciseBounds) {
  size_t numVerts = num_vertices(mesh);
  size_t numFaces = num_faces(mesh);

  outVerts.resize(numVerts);
  outIndices.resize(numFaces);

  if(usePreciseBounds) {
    outVertexErrors.resize(numVerts); // Allocate space for error terms
  } else {
    outVertexErrors.clear(); // Ensure it is empty
  }

  auto coords = mesh.points();

  // Pass 1: Shift geometry, downcast to float, and calculate exact rounding error (if requested)
  tbb::parallel_for(tbb::blocked_range<size_t>(0, numVerts), [&](const tbb::blocked_range<size_t>& r) {
    for(size_t i = r.begin(); i != r.end(); ++i) {
      typename Mesh::Vertex_index v(static_cast<Mesh::size_type>(i));
      const auto& p = coords[v];

      // 1. Calculate high-precision shifted coordinates
      double shifted_x = p.x() - cx;
      double shifted_y = p.y() - cy;
      double shifted_z = p.z() - cz;

      // 2. Downcast to single precision
      float fx = static_cast<float>(shifted_x);
      float fy = static_cast<float>(shifted_y);
      float fz = static_cast<float>(shifted_z);

      outVerts[i] = {fx, fy, fz};

      if(usePreciseBounds) {
        // 3. Compute the exact coordinate drift caused by the downcast
        double err_x = std::abs(shifted_x - static_cast<double>(fx));
        double err_y = std::abs(shifted_y - static_cast<double>(fy));
        double err_z = std::abs(shifted_z - static_cast<double>(fz));

        // 4. Evaluate tight L_infinity (Chebyshev) error bound
        // This is the absolute tightest bound for coordinate-aligned math in the GPU assembler
        double max_exact_err = std::max({err_x, err_y, err_z});

        // 5. Compute a relative safety guard for subsequent single-precision GPU reads/FMA operations
        // This scales dynamically: extremely small near the origin, larger at high coordinates
        float max_coord_mag = std::max({std::abs(fx), std::abs(fy), std::abs(fz)});
        float gpu_drift_guard = max_coord_mag * std::numeric_limits<float>::epsilon();

        // Store the final consolidated L_infinity error
        outVertexErrors[i] = static_cast<float>(max_exact_err) + gpu_drift_guard;
      }
    }
  });

  // Pass 2: Extract pure integer topology connections (Unchanged)
  tbb::parallel_for(tbb::blocked_range<size_t>(0, numFaces), [&](const tbb::blocked_range<size_t>& r) {
    for(size_t i = r.begin(); i != r.end(); ++i) {
      typename Mesh::Face_index f(static_cast<Mesh::size_type>(i));
      auto h0 = mesh.halfedge(f);
      auto h1 = mesh.next(h0);
      auto h2 = mesh.next(h1);

      outIndices[i] = {static_cast<unsigned int>(mesh.target(h0)), static_cast<unsigned int>(mesh.target(h1)),
                       static_cast<unsigned int>(mesh.target(h2))};
    }
  });
}
int main(int ac, char** av) {
  if(ac < 11) {
    std::cout << "Usage: " << av[0]
              << " <meshA.off> <maxCellSizeA> <meshB.off> <maxCellSizeB> <batchmultiplier> <NumberOfDualTreeSteps> "
                 "<Leaf Threshold> <tbb_workers> <usePreciseBounds><AsyncDownload>\n";
    return 1;
  }

  warmupCUDA();

  std::string meshPathA = av[1];
  int maxCellSizeA = std::stoi(av[2]);
  std::string meshPathB = av[3];
  int maxCellSizeB = std::stoi(av[4]);
  int batchmultipl = std::stoi(av[5]);
  int mode = std::stoi(av[6]);
  int leafThreshhold = std::stoi(av[7]);
  int tbbWorkers = std::stoi(av[8]);
  bool usePreciseBounds = (std::stoi(av[9]) != 0);
  int AsyncDownload = std::stoi(av[10]);

  tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, tbbWorkers);

  double tStartTotal = cuBQL::getCurrentTime();

  std::vector<float3> hVertsA;
  std::vector<uint3> hIndicesA;
  std::vector<float3> hVertsB;
  std::vector<uint3> hIndicesB;
  std::vector<float> outVertexErrorsA;
  std::vector<float> outVertexErrorsB;

  // --- LOAD MESH A ---
  double tMeshAFileLoadStart = cuBQL::getCurrentTime();
  Mesh meshA;
  std::ifstream inFileA(meshPathA);
  if(!(inFileA >> meshA)) {
    std::cerr << "Critical Error: Failed to open file A\n";
    return 1;
  }
  inFileA.close();
  double tMeshAFileLoadTime = cuBQL::getCurrentTime() - tMeshAFileLoadStart;

  // --- LOAD MESH B ---
  double tMeshBFileLoadStart = cuBQL::getCurrentTime();
  Mesh meshB;
  std::ifstream inFileB(meshPathB);
  if(!(inFileB >> meshB)) {
    std::cerr << "Critical Error: Failed to open file B\n";
    return 1;
  }
  inFileB.close();
  double tMeshBFileLoadTime = cuBQL::getCurrentTime() - tMeshBFileLoadStart;

  // --- TBB NORMALIZATION ---

  double tNormalizationStart0 = cuBQL::getCurrentTime();

  DoubleBox combinedBox = computeMeshBoundsTBB(meshA);
  combinedBox.grow(computeMeshBoundsTBB(meshB));

  double centerX = 0.5 * (combinedBox.min_x + combinedBox.max_x);
  double centerY = 0.5 * (combinedBox.min_y + combinedBox.max_y);
  double centerZ = 0.5 * (combinedBox.min_z + combinedBox.max_z);

    double tCentralizationTime = cuBQL::getCurrentTime() - tNormalizationStart0;

    double tNormalizationStart = cuBQL::getCurrentTime();
  extractMeshTopologyNormalized(meshA, centerX, centerY, centerZ, hVertsA, outVertexErrorsA, hIndicesA,
                                usePreciseBounds);
  extractMeshTopologyNormalized(meshB, centerX, centerY, centerZ, hVertsB, outVertexErrorsB, hIndicesB,
                                usePreciseBounds);
  double tNormalizationTime = cuBQL::getCurrentTime() - tNormalizationStart;

  // --------------------------------------------------------------------
  // EXECUTE SPATIAL CROSS-INTERSECTION PIPELINE
  // --------------------------------------------------------------------
  std::cout << "[Pipeline] Running Dual-Mesh GPU pipeline execution...\n";

  ExecutionStats stats;
 // std::vector<int2> hGreenPairs;
  //std::vector<int2> hYellowPairs;
tbb::concurrent_vector<int2>  finalExactPairs;
  // Phase 1: Call your GPU pipeline setup
  kernelsTestBVHV2(meshA, meshB, hVertsA.data(), static_cast<int>(hVertsA.size()), hIndicesA.data(),
                   (usePreciseBounds ? outVertexErrorsA.data() : nullptr), static_cast<int>(hIndicesA.size()),
                   maxCellSizeA, hVertsB.data(), static_cast<int>(hVertsB.size()), hIndicesB.data(),
                   (usePreciseBounds ? outVertexErrorsB.data() : nullptr), static_cast<int>(hIndicesB.size()),
                   maxCellSizeB, batchmultipl, mode, leafThreshhold, AsyncDownload, stats, finalExactPairs);

  // --------------------------------------------------------------------
  // HYBRID NARROW-PHASE FILTERING
  // --------------------------------------------------------------------
//   double t_inter_start = cuBQL::getCurrentTime();

//   // Directly initialize thread-safe storage with our definitive greens
//   tbb::concurrent_vector<int2> finalExactPairs(hGreenPairs.begin(), hGreenPairs.end());

//   // Phase 2: Filter out structural collisions using our clean helper
//   filterYellowPairsTBB(meshA, meshB, hYellowPairs.data(), hYellowPairs.size(),
//                                                          finalExactPairs);

//   double time_intersection = cuBQL::getCurrentTime() - t_inter_start;

  // --------------------------------------------------------------------
  // FORMATTED REPORTING OF METRICS FROM EXECUTIONSTATS
  // --------------------------------------------------------------------
  std::cout << "\n==================================================\n";
  std::cout << "             PIPELINE METRICS REPORT              \n";
  std::cout << "==================================================\n";
  std::cout << "STRUCTURE SUMMARY & PROPORTIONS:\n";
  std::cout << "  |- Pipeline Execution Mode:        " << mode << "\n";
  std::cout << "  |- TBB Max Parallelism Limit:      " << tbbWorkers << "\n";
  std::cout << "  |- Precise Error Bounds Active:    "
            << (usePreciseBounds ? "YES" : "NO (Using GPU Conservative Fallback)") << "\n";
  std::cout << "  |- Mesh A Total Generated Nodes:   " << stats.meshATotalNodes << "\n";
  std::cout << "  |- Mesh A Extracted Targets (<" << maxCellSizeA << "): " << stats.meshAExtractedTargets << "\n";
  std::cout << "  |- Mesh B Total Generated Nodes:   " << stats.meshBTotalNodes << "\n";
  std::cout << "  |- Mesh B Extracted Targets (<" << maxCellSizeB << "): " << stats.meshBExtractedTargets << "\n\n";

  std::cout << "HOST-SIDE INGESTION & TRANSLATION BREAKDOWN:\n";
  std::cout << "  |- Mesh A Disk Ingest (OFF Read):  " << std::fixed << std::setprecision(4) << tMeshAFileLoadTime
            << " s\n";
  std::cout << "  |- Mesh B Disk Ingest (OFF Read):  " << tMeshBFileLoadTime << " s\n";
  std::cout << "  |- TBB Double Bounds & Normalize:  " << tNormalizationTime << " s\n\n";
  std::cout << "  |- Centralization time:            " << tCentralizationTime << " s\n";

  std::cout << "CRISS-CROSS BOUNDING BOX CROSS-CHECK:\n";
  std::cout << "  |- Intersections found:            " << stats.totalIntersections << " / " << stats.totalPossiblePairs
            << "\n";
  std::cout << "  |- Intersection Ratio:             " << std::fixed << std::setprecision(4)
            << stats.intersectionPercentage << "%\n\n";

  std::cout << "TIMING METRICS OVERVIEW:\n";
  std::cout << "  |- Initial Alloc & Mesh Copy:      " << std::fixed << std::setprecision(4)
            << stats.initialAllocAndCopyMs << " ms\n";
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
  std::cout << "  |- Confirmed Green (Intersecting): " << stats.loopTracker.confirmedGreenPairs << "\n";
  std::cout << "  |- Confirmed Yellow (Coplanar):    " << stats.loopTracker.confirmedYellowPairs << "\n";

  std::cout << "FINE-GRAINED INTERSECTION LOOP BREAKDOWN:\n";
  std::cout << "  |- Total Loop Execution Time:      " << stats.loopTracker.totalLoopTimeMs << " ms\n";
  std::cout << "  |- Preallocation Phase Time:       " << stats.loopTracker.preallocateTimeMs << " ms\n";
  std::cout << "  |- Assembly Phase Time:            " << stats.loopTracker.assemblyPhaseMs << " ms\n";
  std::cout << "  |- Execution Phase Time:           " << stats.loopTracker.executionPhaseMs << " ms\n";
  std::cout << "  |- Fine Evaluation Phase Time:     " << stats.loopTracker.fineEvaluationPhaseMs << " ms\n";
  std::cout << "  |- Device Cleanup Time:            " << stats.loopTracker.cleanupTimeMs << " ms\n";
  std::cout << "  |- Host Download & Sync Time:      " << stats.loopTracker.DownloadAndClean << " ms\n\n";
  std::cout << "  |- Number of batches:              " << stats.loopTracker.numberOfBatchLoops << "\n";

  std::cout << "NARROW-PHASE HYBRID FILTER REPORT:\n";
  std::cout << "  |- CPU Narrow-Phase Compute Time:  " << stats.loopTracker.CPUPredicates << " s\n";
  std::cout << "  |- Ambiguous Yellows Evaluated:    " << stats.loopTracker.confirmedYellowPairs << "\n";
  std::cout << "  |- Total Final Intersecting Pairs: " << finalExactPairs.size()
            << " (Green Hit Count + Validated Yellows)\n";
  std::cout << "==================================================\n";

  std::cout << "Total Diagnostic App Execution Time: " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - tStartTotal)
            << "s\n";
  std::cout << "Total without file read :            " <<cuBQL::prettyDouble(cuBQL::getCurrentTime() - tNormalizationStart)  << "s\n";
  std::cout << "==================================================\n";

  return 0;
}