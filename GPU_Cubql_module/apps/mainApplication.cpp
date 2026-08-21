#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <thread>
#include <mutex>
#include <queue>
#include <atomic>
#include <condition_variable>
#include <algorithm>
#include <limits>

// GNU Readline headers for terminal arrow key & history support
#include <readline/readline.h>
#include <readline/history.h>

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
#include <cuBQL/bvh.h>

// Custom Controllers & CPU Utilities
#include "GPUIntersector/KernelBVHController.h"
#include "CPU/MeshIntersectionCurves.h"
#include "CPU/CgalDefinitions.h"
#include "CPU/PolyscopeBridge.h"
#include "CPU/ParallelCgalOffLoader.h"
#include "CPU/RotationTools.h"
#include "CPU/MeshTriangleDegeneracyVisualizer.h"
#include "CPU/CommandDispatcher.h"
#include "CPU/ExecutionTimingVisualizer.h"

#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/intersections.h>

#include "CPU/TestSuite.h"
#include "CPU/ApplicationState.h"

// Hardware Warmup & Pipeline Entry Points
#include "Warmup/cuda_warmup.h"
#include "GPUIntersector/StandaloneBVHPipeline.h"




static TestConfig g_testConfig;

/**
 * @brief Accumulator for computing mesh face and edge geometry metrics across parallel TBB threads.
 */
struct MeshStatsAccumulator
{
  double minArea = std::numeric_limits<double>::infinity();
  double maxArea = -1.0;
  double totalArea = 0.0;

  double minEdge = std::numeric_limits<double>::infinity();
  double maxEdge = -1.0;
  double totalEdge = 0.0;
  size_t totalEdgesCount = 0;
  size_t validFaces = 0;

  void merge(const MeshStatsAccumulator& rhs) {
    minArea = std::min(minArea, rhs.minArea);
    maxArea = std::max(maxArea, rhs.maxArea);
    totalArea += rhs.totalArea;

    minEdge = std::min(minEdge, rhs.minEdge);
    maxEdge = std::max(maxEdge, rhs.maxEdge);
    totalEdge += rhs.totalEdge;
    totalEdgesCount += rhs.totalEdgesCount;
    validFaces += rhs.validFaces;
  }
};

/**
 * @brief Double-precision Axis-Aligned Bounding Box (AABB) for spatial query boundaries.
 */
struct BoundingBox3D
{
  double min_x = std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double min_z = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();
  double max_z = -std::numeric_limits<double>::infinity();

  void expandByPoint(double x, double y, double z) {
    min_x = std::min(min_x, x);
    max_x = std::max(max_x, x);
    min_y = std::min(min_y, y);
    max_y = std::max(max_y, y);
    min_z = std::min(min_z, z);
    max_z = std::max(max_z, z);
  }

  void expandByBox(const BoundingBox3D& o) {
    min_x = std::min(min_x, o.min_x);
    max_x = std::max(max_x, o.max_x);
    min_y = std::min(min_y, o.min_y);
    max_y = std::max(max_y, o.max_y);
    min_z = std::min(min_z, o.min_z);
    max_z = std::max(max_z, o.max_z);
  }

  double3 getCenter() const {
    return make_double3(0.5 * (min_x + max_x), 0.5 * (min_y + max_y), 0.5 * (min_z + max_z));
  }
};

/**
 * @brief Calculates the global 3D AABB for a CGAL Surface_mesh in parallel using Intel TBB reduce.
 * @param mesh Input CGAL Surface_mesh.
 * @return BoundingBox3D Encapsulating boundary of all mesh vertices.
 */
BoundingBox3D computeMeshAABBParallel(const Mesh& mesh) {
  size_t numVerts = num_vertices(mesh);
  auto coords = mesh.points();

  return tbb::parallel_reduce(
      tbb::blocked_range<size_t>(0, numVerts), BoundingBox3D(),
      [&](const tbb::blocked_range<size_t>& r, BoundingBox3D box) {
        for(size_t i = r.begin(); i != r.end(); ++i) {
          typename Mesh::Vertex_index v(static_cast<Mesh::size_type>(i));
          const auto& p = coords[v];
          box.expandByPoint(p.x(), p.y(), p.z());
        }
        return box;
      },
      [](BoundingBox3D a, const BoundingBox3D& b) {
        a.expandByBox(b);
        return a;
      });
}

/**
 * @brief Translates mesh vertices to origin and scales coordinates in parallel using Intel TBB.
 * @param mesh Reference to CGAL mesh modified in-place.
 * @param cx Center X offset.
 * @param cy Center Y offset.
 * @param cz Center Z offset.
 * @param scaleFactor Normalization scale multiplier.
 */
void normalizeMeshCoordinatesParallel(Mesh& mesh, double cx, double cy, double cz, double scaleFactor) {
  size_t numVerts = num_vertices(mesh);
  auto coords = mesh.points();

  tbb::parallel_for(tbb::blocked_range<size_t>(0, numVerts), [&](const tbb::blocked_range<size_t>& r) {
    for(size_t i = r.begin(); i != r.end(); ++i) {
      typename Mesh::Vertex_index v(static_cast<Mesh::size_type>(i));
      const auto& p = coords[v];
      double nx = (p.x() - cx) * scaleFactor;
      double ny = (p.y() - cy) * scaleFactor;
      double nz = (p.z() - cz) * scaleFactor;
      coords[v] = Point3(nx, ny, nz);
    }
  });
}

void cmdTestConfig(ApplicationState& app, std::istringstream& iss) {
  int steps, prec, dualTreeSteps, queryLvl, refLvl, batch;
  double trans, rot;
  unsigned int seed;

  if(iss >> steps >> trans >> rot >> prec >> dualTreeSteps >> queryLvl >> refLvl >> batch >> seed) {
    g_testConfig.numSteps = steps;
    g_testConfig.maxTranslation = trans;
    g_testConfig.maxRotationDeg = rot;
    g_testConfig.enableGpuPrecision = prec;
    g_testConfig.dualTreeSteps = dualTreeSteps;
    g_testConfig.queryDescentLevel = queryLvl;
    g_testConfig.referenceDescentLevel = refLvl;
    g_testConfig.batchMultiplier = batch;
    g_testConfig.seed = seed;
    std::cout << "TestConfig updated successfully.\n";
  }

  std::cout << "\n=======================================================\n";
  std::cout << "            CURRENT TEST SUITE CONFIGURATION           \n";
  std::cout << "=======================================================\n";
  std::cout << "  Steps (numSteps)           : " << g_testConfig.numSteps << "\n";
  std::cout << "  Max Translation (maxTrans) : " << std::defaultfloat
            << std::setprecision(std::numeric_limits<double>::max_digits10) << g_testConfig.maxTranslation << "\n";
  std::cout << "  Max Rotation (maxRotDeg)   : " << g_testConfig.maxRotationDeg << "°\n";
  std::cout << "  GPU Precision Mode (prec)  : " << g_testConfig.enableGpuPrecision << "\n";
  std::cout << "  Dual tree steps performed  : " << g_testConfig.dualTreeSteps << "\n";
  std::cout << "  Query Descent Level        : " << g_testConfig.queryDescentLevel << "\n";
  std::cout << "  Reference Descent Level    : " << g_testConfig.referenceDescentLevel << "\n";
  std::cout << "  Batch Multiplier           : " << g_testConfig.batchMultiplier << "\n";
  std::cout << "  Random Seed                : " << g_testConfig.seed << "\n";
  std::cout << "=======================================================\n\n";
}


bool exportMeshToOff(const Mesh& mesh, const std::string& filepath) {
  std::ofstream outFile(filepath);
  if(!outFile.is_open())
    return false;
  outFile << std::setprecision(17);
  return static_cast<bool>(outFile << mesh);
}

/**
 * @brief Extracts normalized CGAL mesh points into a contiguous host vector of CUDA double3 elements.
 */
void extractVertexCoordinatesParallel(const Mesh& mesh, std::vector<double3>& outVerts) {
  size_t numVerts = num_vertices(mesh);
  outVerts.resize(numVerts);
  auto coords = mesh.points();

  tbb::parallel_for(tbb::blocked_range<size_t>(0, numVerts), [&](const tbb::blocked_range<size_t>& r) {
    for(size_t i = r.begin(); i != r.end(); ++i) {
      typename Mesh::Vertex_index v(static_cast<Mesh::size_type>(i));
      const auto& p = coords[v];
      outVerts[i] = make_double3(p.x(), p.y(), p.z());
    }
  });
}

/**
 * @brief Standard single-threaded stream parser for OFF mesh files using CGAL ifstream operator.
 */
bool loadMeshOffSequential(const std::string& filepath,
                           Mesh& mesh,
                           std::vector<double3>& outVerts,
                           std::vector<uint3>& outIndices) {
  std::ifstream inFile(filepath);
  if(!inFile.is_open() || !(inFile >> mesh))
    return false;

  size_t numVerts = num_vertices(mesh);
  outVerts.resize(numVerts);
  auto coords = mesh.points();
  for(size_t i = 0; i < numVerts; ++i) {
    typename Mesh::Vertex_index v(static_cast<Mesh::size_type>(i));
    const auto& p = coords[v];
    outVerts[i] = make_double3(p.x(), p.y(), p.z());
  }

  size_t numFaces = num_faces(mesh);
  outIndices.resize(numFaces);
  size_t fIdx = 0;
  for(auto f : mesh.faces()) {
    auto h = mesh.halfedge(f);
    outIndices[fIdx++] = {static_cast<uint32_t>(mesh.source(h)), static_cast<uint32_t>(mesh.target(h)),
                          static_cast<uint32_t>(mesh.target(mesh.next(h)))};
  }
  return true;
}


std::string parseArgument(std::istringstream& iss) {
  std::string arg;
  iss >> std::ws;
  if(iss.peek() == '"' || iss.peek() == '\'') {
    char quote = iss.get();
    std::getline(iss, arg, quote);
  } else {
    iss >> arg;
  }
  return arg;
}

// Thread-safe command queue
static std::mutex g_cmdMutex;
static std::condition_variable g_cmdCV;
static std::queue<std::string> g_cmdQueue;
static std::atomic<bool> g_running{true};
static bool g_readyForInput{true};

/**
 * @brief Asynchronous CLI reader thread maintaining interactive terminal input via GNU Readline.
 */
void inputThreadWorker() {
  char* rawInput = nullptr;
  while(g_running) {
    {
      std::unique_lock<std::mutex> lock(g_cmdMutex);
      g_cmdCV.wait(lock, [] { return g_readyForInput || !g_running; });
      if(!g_running)
        break;
    }

    rawInput = readline("> ");
    if(!rawInput) {
      g_running = false;
      break;
    }

    std::string line(rawInput);
    if(!line.empty())
      add_history(rawInput);
    free(rawInput);

    {
      std::lock_guard<std::mutex> lock(g_cmdMutex);
      g_cmdQueue.push(line);
      g_readyForInput = false;
    }

    std::istringstream iss(line);
    std::string cmd;
    iss >> cmd;
    if(cmd == "quit" || cmd == "exit") {
      g_running = false;
      break;
    }
  }
  g_running = false;
}

// --------------------------------------------------------------------
// APPLICATION STATE CONTEXT
// --------------------------------------------------------------------
// struct ApplicationState
// {
//   KernelBVHController controller;
//   ExecutionStats stats;
//   Mesh meshA, meshB;
//   std::vector<double3> hVertsA, hVertsB;
//   std::vector<uint3> hIndicesA, hIndicesB;
//   std::unique_ptr<tbb::global_control> tbbControl;
//   MeshTriangleDegeneracyVisualizer edgeVisualizer;
//   ExecutionTimingVisualizer timingVisualizer;

//   bool isLoaded = false;
//   double currentScaleFactor = 1.0;
//   double currentCenterX = 0.0, currentCenterY = 0.0, currentCenterZ = 0.0;
//   double currentMaxSpan = 0.0;
//   bool currentScaledToUnit = false;

//   double3 origCenterA{0.0, 0.0, 0.0}, origCenterB{0.0, 0.0, 0.0};
//   float3 normCenterA{0.0f, 0.0f, 0.0f}, normCenterB{0.0f, 0.0f, 0.0f};
// };

/**
 * @brief Computes scene bounds, centers geometry, scales coordinates, and updates CPU host buffers.
 * @param app Global application state context.
 * @param scaleToUnit If true, scales maximum scene bounding box dimension to fit within [-1, 1].
 */
void normalizeAndSynchronizeScene(ApplicationState& app, bool scaleToUnit) {
  BoundingBox3D boxA_orig = computeMeshAABBParallel(app.meshA);
  BoundingBox3D boxB_orig = computeMeshAABBParallel(app.meshB);
  app.origCenterA = boxA_orig.getCenter();
  app.origCenterB = boxB_orig.getCenter();

  BoundingBox3D sceneBox = boxA_orig;
  sceneBox.expandByBox(boxB_orig);

  double cx = 0.5 * (sceneBox.min_x + sceneBox.max_x);
  double cy = 0.5 * (sceneBox.min_y + sceneBox.max_y);
  double cz = 0.5 * (sceneBox.min_z + sceneBox.max_z);

  double spanX = sceneBox.max_x - sceneBox.min_x;
  double spanY = sceneBox.max_y - sceneBox.min_y;
  double spanZ = sceneBox.max_z - sceneBox.min_z;
  double maxSpan = std::max({spanX, spanY, spanZ});

  double scaleFactor = (scaleToUnit && maxSpan > 0.0) ? (2.0 / maxSpan) : 1.0;

  app.currentScaleFactor = scaleFactor;
  app.currentCenterX = cx;
  app.currentCenterY = cy;
  app.currentCenterZ = cz;
  app.currentMaxSpan = maxSpan;
  app.currentScaledToUnit = scaleToUnit;

  normalizeMeshCoordinatesParallel(app.meshA, cx, cy, cz, scaleFactor);
  normalizeMeshCoordinatesParallel(app.meshB, cx, cy, cz, scaleFactor);

  BoundingBox3D boxA_norm = computeMeshAABBParallel(app.meshA);
  BoundingBox3D boxB_norm = computeMeshAABBParallel(app.meshB);
  double3 cA_norm = boxA_norm.getCenter();
  double3 cB_norm = boxB_norm.getCenter();

  app.normCenterA =
      make_float3(static_cast<float>(cA_norm.x), static_cast<float>(cA_norm.y), static_cast<float>(cA_norm.z));
  app.normCenterB =
      make_float3(static_cast<float>(cB_norm.x), static_cast<float>(cB_norm.y), static_cast<float>(cB_norm.z));

  extractVertexCoordinatesParallel(app.meshA, app.hVertsA);
  extractVertexCoordinatesParallel(app.meshB, app.hVertsB);
}

// --------------------------------------------------------------------
// COMMAND HANDLERS
// --------------------------------------------------------------------

/**
 * @brief Fetches active Polyscope viewport transforms, executes GPU BVH intersection pipeline,
 *        and passes resulting intersecting polygon pairs back to visualizer.
 * @param app Application context state.
 * @param batchMultiplier Triangle candidate pair batch size modifier.
 * @param mode Pipeline execution mode bitmask.
 * @param enableGpuDoublePrecision Enables exact GPU double-precision floating point predicates when 1.
 */
void executeMeshIntersectionPipeline(ApplicationState& app,
                                     int batchMultiplier = std::numeric_limits<int>::max(),
                                     int mode = 0,
                                     int enableGpuDoublePrecision = 1,
                                     bool showCurves = false) {
  if(!app.isLoaded) {
    std::cout << "Error: You must 'load' meshes before computing.\n";
    return;
  }

  if(!app.controller.isGPUAllocated()) {
    std::cout << "[runCompute] GPU memory is cleared/empty. Automatically reconstructing...\n";
    app.controller.reconstructGPU(app.stats);
  }

  float3 rotA, transA, rotB, transB;
  if(!PolyscopeBridge::getCurrentTransforms(rotA, transA, rotB, transB)) {
    std::cout << "Error: Failed to fetch transform matrices from Polyscope.\n";
    return;
  }

  float tGPU, tCPU, tTV, tAT, tGB;
  app.controller.setTransformBoth(make_double3(rotA.x, rotA.y, rotA.z), make_double3(transA.x, transA.y, transA.z),
                                  make_double3(rotB.x, rotB.y, rotB.z), make_double3(transB.x, transB.y, transB.z),
                                  tGPU, tCPU, tTV, tAT, tGB);

  auto tStart = std::chrono::high_resolution_clock::now();
  int2* outIntersectionPairs = nullptr;
  size_t outFinalCount = 0;

  app.controller.runIntersectionPipeline(batchMultiplier, mode, 0, outIntersectionPairs, outFinalCount, app.stats,
                                         enableGpuDoublePrecision);

  app.timingVisualizer.updateAndShowCompute(app.stats);

  auto tEnd = std::chrono::high_resolution_clock::now();
  double elapsedMs = std::chrono::duration<double, std::milli>(tEnd - tStart).count();

  std::vector<int2> stdPairs;
  if(outIntersectionPairs && outFinalCount > 0) {
    stdPairs.assign(outIntersectionPairs, outIntersectionPairs + outFinalCount);
    PolyscopeBridge::highlightIntersections(stdPairs, num_faces(app.meshA), num_faces(app.meshB));
  } else {
    PolyscopeBridge::highlightIntersections({}, num_faces(app.meshA), num_faces(app.meshB));
  }

  // Handle curve generation via PolyscopeBridge abstraction
  double curveMs = 0.0;
  if (showCurves) {
    auto tCurveStart = std::chrono::high_resolution_clock::now();

    Point3 centerA(app.normCenterA.x, app.normCenterA.y, app.normCenterA.z);
    Point3 centerB(app.normCenterB.x, app.normCenterB.y, app.normCenterB.z);

    auto segments = computeExactIntersectionCurves(
        app.meshA, app.meshB, stdPairs,
        centerA, make_double3(rotA.x, rotA.y, rotA.z), make_double3(transA.x, transA.y, transA.z),
        centerB, make_double3(rotB.x, rotB.y, rotB.z), make_double3(transB.x, transB.y, transB.z));

    

    auto tCurveEnd = std::chrono::high_resolution_clock::now();

    PolyscopeBridge::showIntersectionCurves(segments);

    curveMs = std::chrono::duration<double, std::milli>(tCurveEnd - tCurveStart).count();
    
  } else {
    PolyscopeBridge::clearIntersectionCurves();
  }

  std::cout << "\n[Compute] Query completed in: " << std::fixed << std::setprecision(2) << elapsedMs << " ms. "
            << "Found " << outFinalCount << " intersections.\n";
  if (showCurves) {
    std::cout << "[Curves] Extracted exact curves in: " << std::fixed << std::setprecision(2) << curveMs << " ms.\n";
  }
  std::cout << "AABB Hits: " << app.stats.finalAabbCandidatePairs
            << " | Green: " << app.stats.loopTracker.confirmedGreenPairs
            << " | Yellow: " << app.stats.loopTracker.confirmedYellowPairs
            << " | Orange: " << app.stats.loopTracker.confirmedOrangePairs << std::endl;
  std::cout << "Time GPU Predicates (Yellow): " << app.stats.loopTracker.fineEvaluationPhaseMs
            << " ms | Time GPU Predicates (Orange): " << app.stats.loopTracker.gpuDoublePredicatesMs
            << " ms | Time CPU Predicates: " << app.stats.loopTracker.CPUPredicates << " ms\n";
  std::cout << "Time setup GPU: " << tGPU << " ms | Time setup CPU: " << tCPU << " ms" << std::endl;
  std::cout << "RotateTransformVerts: " << tTV << " timeAssembleTris: " << tAT << " timeGenBoxes: " << tGB << std::endl;

  if(outIntersectionPairs) {
    std::free(outIntersectionPairs);
  }
}

void cmdScale(ApplicationState& app, std::istringstream&) {
  if(!app.isLoaded) {
    std::cout << "Error: No meshes loaded. Run 'load' first.\n";
    return;
  }
  std::cout << "\n=======================================================\n";
  std::cout << "         MESH NORMALIZATION & SCALING INFO             \n";
  std::cout << "=======================================================\n";
  std::cout << "  Scale Factor Applied : " << std::setprecision(12) << app.currentScaleFactor << "\n";
  std::cout << "  Original Scene Center: (" << std::setprecision(8) << app.currentCenterX << ", " << app.currentCenterY
            << ", " << app.currentCenterZ << ")\n";
  std::cout << "  Max Original Extent  : " << app.currentMaxSpan << " units\n";
  std::cout << "  Coordinate Space     : "
            << (app.currentScaledToUnit ? "[-1.0, 1.0]^3 (Centered & Scaled)" : "Centered (Unscaled)") << "\n";
  std::cout << "  -----------------------------------------------------\n";
  std::cout << "  Original Centroid A  : (" << app.origCenterA.x << ", " << app.origCenterA.y << ", "
            << app.origCenterA.z << ")\n";
  std::cout << "  Original Centroid B  : (" << app.origCenterB.x << ", " << app.origCenterB.y << ", "
            << app.origCenterB.z << ")\n";
  std::cout << "  Normalized Centroid A: (" << app.normCenterA.x << ", " << app.normCenterA.y << ", "
            << app.normCenterA.z << ")\n";
  std::cout << "  Normalized Centroid B: (" << app.normCenterB.x << ", " << app.normCenterB.y << ", "
            << app.normCenterB.z << ")\n";
  std::cout << "=======================================================\n\n";
}

void cmdTbb(ApplicationState& app, std::istringstream& iss) {
  int numThreads = 0;
  if(iss >> numThreads && numThreads > 0) {
    app.tbbControl = std::make_unique<tbb::global_control>(tbb::global_control::max_allowed_parallelism, numThreads);
    std::cout << "TBB maximum worker limit updated to: " << numThreads << "\n";
  } else {
    std::cout << "Usage: tbb <num_threads> (must be > 0)\n";
  }
}

void cmdLoad(ApplicationState& app, std::istringstream& iss, bool useSequentialLoader) {
  std::string arg1 = parseArgument(iss);
  std::string arg2 = parseArgument(iss);

  int queryDescentLevel = 1, referenceDescentLevel = 1, leafThresh = 4, scaleToUnitInt = 0;

  if(arg1.empty()) {
    std::cout << "Usage: load" << (useSequentialLoader ? "Old" : "")
              << " <meshA> [meshB] [queryDescentLevel] [referenceDescentLevel] [leafThresh] [scaleToUnit(0/1)]\n";
    return;
  }

  std::string pathA, pathB;
  bool singleMeshMode = false;

  if(arg2.empty()) {
    singleMeshMode = true;
  } else {
    std::istringstream checkNum(arg2);
    int val;
    if(checkNum >> val) {
      singleMeshMode = true;
      queryDescentLevel = val;
    }
  }

  if(singleMeshMode) {
    pathA = arg1;
    pathB = arg1;
    iss >> referenceDescentLevel >> leafThresh >> scaleToUnitInt;
  } else {
    pathA = arg1;
    pathB = arg2;
    iss >> queryDescentLevel >> referenceDescentLevel >> leafThresh >> scaleToUnitInt;
  }

  bool scaleToUnit = (scaleToUnitInt != 0);

  app.controller.cleanup();
  app.meshA = Mesh();
  app.meshB = Mesh();
  app.hVertsA.clear();
  app.hIndicesA.clear();
  app.hVertsB.clear();
  app.hIndicesB.clear();
  app.isLoaded = false;

  std::cout << "Loading mesh(es) using " << (useSequentialLoader ? "Old CGAL Stream Loader" : "Fast Parallel IO Loader")
            << "...\n";

  auto tIoStart = std::chrono::high_resolution_clock::now();
  bool loadedSuccessfully = false;

  if(useSequentialLoader) {
    if(singleMeshMode) {
      if(loadMeshOffSequential(pathA, app.meshA, app.hVertsA, app.hIndicesA)) {
        app.meshB = app.meshA;
        app.hVertsB = app.hVertsA;
        app.hIndicesB = app.hIndicesA;
        loadedSuccessfully = true;
      }
    } else {
      if(loadMeshOffSequential(pathA, app.meshA, app.hVertsA, app.hIndicesA) &&
         loadMeshOffSequential(pathB, app.meshB, app.hVertsB, app.hIndicesB))
      {
        loadedSuccessfully = true;
      }
    }
  } else {
    if(singleMeshMode) {
      if(ParallelIO::loadOffToCgalMesh<Mesh, Point3>(pathA, app.meshA, app.hIndicesA)) {
        app.meshB = app.meshA;
        app.hIndicesB = app.hIndicesA;
        loadedSuccessfully = true;
      } else {
        std::cerr << "Error: Could not open/read Mesh file: " << pathA << "\n";
      }
    } else {
      if(ParallelIO::loadOffToCgalMesh<Mesh, Point3>(pathA, app.meshA, app.hIndicesA) &&
         ParallelIO::loadOffToCgalMesh<Mesh, Point3>(pathB, app.meshB, app.hIndicesB))
      {
        loadedSuccessfully = true;
      } else {
        std::cerr << "Error: Failed to load OFF input files via fast parallel IO loader.\n";
      }
    }
  }

  if(loadedSuccessfully) {
    auto tIoEnd = std::chrono::high_resolution_clock::now();
    double ioMs = std::chrono::duration<double, std::milli>(tIoEnd - tIoStart).count();

    auto tPrepStart = std::chrono::high_resolution_clock::now();

    normalizeAndSynchronizeScene(app, scaleToUnit);

    auto tPrepEnd = std::chrono::high_resolution_clock::now();
    double prepMs = std::chrono::duration<double, std::milli>(tPrepEnd - tPrepStart).count();

    auto tStart = std::chrono::high_resolution_clock::now();

    Point3 cA_cgal(app.normCenterA.x, app.normCenterA.y, app.normCenterA.z);
    Point3 cB_cgal(app.normCenterB.x, app.normCenterB.y, app.normCenterB.z);

    app.controller.construct(app.meshA, app.meshB, cA_cgal, cB_cgal, app.hVertsA.data(),
                             static_cast<int>(app.hVertsA.size()), app.hIndicesA.data(),
                             static_cast<int>(app.hIndicesA.size()), queryDescentLevel, app.hVertsB.data(),
                             static_cast<int>(app.hVertsB.size()), app.hIndicesB.data(),
                             static_cast<int>(app.hIndicesB.size()), referenceDescentLevel, leafThresh, app.stats);

    auto tEnd = std::chrono::high_resolution_clock::now();
    double elapsedMs = std::chrono::duration<double, std::milli>(tEnd - tStart).count();

    PolyscopeBridge::reset(app.meshA, app.meshB, app.normCenterA, app.normCenterB);

    std::cout << "\n--- Load & Construction Performance Breakdown ---\n";
    std::cout << "  Loader Used:               "
              << (useSequentialLoader ? "Standard CGAL Stream (loadOld)" : "Parallel IO (load)") << "\n";
    std::cout << "  1. Disk -> IO Load:        " << std::fixed << std::setprecision(2) << ioMs << " ms\n";
    std::cout << "  2. Topology/Layout Prep:   " << std::fixed << std::setprecision(2) << prepMs << " ms\n";
    std::cout << "  3. BVH Construction:       " << std::fixed << std::setprecision(2) << elapsedMs << " ms\n";
    std::cout << "  -----------------------------------------------\n";
    std::cout << "  Total Preparation Time:    " << std::fixed << std::setprecision(2) << (ioMs + prepMs + elapsedMs)
              << " ms\n";
    std::cout << "  Applied Scale Factor:      " << std::setprecision(10) << app.currentScaleFactor
              << (scaleToUnit ? " (Scaled to [-1, 1]^3)" : " (Unscaled)") << "\n";
    std::cout << "  Normalized Center A:       (" << app.normCenterA.x << ", " << app.normCenterA.y << ", "
              << app.normCenterA.z << ")\n";
    std::cout << "  Normalized Center B:       (" << app.normCenterB.x << ", " << app.normCenterB.y << ", "
              << app.normCenterB.z << ")\n\n";

    app.isLoaded = true;
  }
}

void cmdTransform(ApplicationState& app, std::istringstream& iss) {
  if(!app.isLoaded) {
    std::cout << "Error: You must 'load' meshes before transforming.\n";
    return;
  }
  double rotAx, rotAy, rotAz, transAx, transAy, transAz, rotBx, rotBy, rotBz, transBx, transBy, transBz;

  if(iss >> rotAx >> rotAy >> rotAz >> transAx >> transAy >> transAz >> rotBx >> rotBy >> rotBz >> transBx >> transBy >>
     transBz)
  {
    float3 fRotA = make_float3(static_cast<float>(rotAx), static_cast<float>(rotAy), static_cast<float>(rotAz));
    float3 fTransA = make_float3(static_cast<float>(transAx), static_cast<float>(transAy), static_cast<float>(transAz));
    float3 fRotB = make_float3(static_cast<float>(rotBx), static_cast<float>(rotBy), static_cast<float>(rotBz));
    float3 fTransB = make_float3(static_cast<float>(transBx), static_cast<float>(transBy), static_cast<float>(transBz));

    PolyscopeBridge::transformBoth(fRotA, fTransA, app.normCenterA, fRotB, fTransB, app.normCenterB);

    std::cout << "Transformations applied to viewport:\n";
    std::cout << "  Mesh A Rot: (" << rotAx << ", " << rotAy << ", " << rotAz << ") | Trans: (" << transAx << ", "
              << transAy << ", " << transAz << ")\n";
    std::cout << "  Mesh B Rot: (" << rotBx << ", " << rotBy << ", " << rotBz << ") | Trans: (" << transBx << ", "
              << transBy << ", " << transBz << ")\n";
  } else {
    std::cout << "Usage: transform <rotAx> <rotAy> <rotAz> <transAx> <transAy> <transAz> <rotBx> <rotBy> <rotBz> "
                 "<transBx> <transBy> <transBz>\n";
  }
}

void cmdTranslate(ApplicationState& app, std::istringstream& iss) {
  double x, y, z;
  if(iss >> x >> y >> z) {
    auto tStart = std::chrono::high_resolution_clock::now();
    app.controller.setTranslation(x, y, z);
    PolyscopeBridge::translateMeshB(static_cast<float>(x), static_cast<float>(y), static_cast<float>(z));
    auto tEnd = std::chrono::high_resolution_clock::now();
    double elapsedMs = std::chrono::duration<double, std::milli>(tEnd - tStart).count();
    std::cout << "Translation applied (GPU Kernel): (" << x << ", " << y << ", " << z << ") in " << elapsedMs
              << " ms\n";
  } else {
    std::cout << "Usage: translate <x> <y> <z>\n";
  }
}

void cmdClear(ApplicationState& app, std::istringstream&) {
  if(!app.isLoaded) {
    std::cout << "Warning: No meshes loaded, but clearing GPU memory allocations...\n";
  }
  app.controller.clearGPU();
  std::cout << "GPU memory buffers and cuBQL BVH acceleration structures cleared.\n";
}

void cmdReconstruct(ApplicationState& app, std::istringstream& iss) {
  if(!app.isLoaded) {
    std::cout << "Error: No mesh loaded on CPU. Use 'load' first before reconstructing.\n";
    return;
  }

  int levelA = 1;
  int levelB = 1;

  if(!(iss >> levelA >> levelB)) {
    std::cout << "Usage: reconstruct <queryDescentLevel> <referenceDescentLevel>\n";
    return;
  }

  std::cout << "Reconstructing GPU resources (LevelA: " << levelA << ", LevelB: " << levelB
            << ") and rebuilding cuBQL BVHs...\n";
  auto tStart = std::chrono::high_resolution_clock::now();
  app.controller.reconstructGPU(app.stats, levelA, levelB);
  auto tEnd = std::chrono::high_resolution_clock::now();
  double elapsedMs = std::chrono::duration<double, std::milli>(tEnd - tStart).count();
  std::cout << "GPU reconstruction and BVH rebuilding complete in " << std::fixed << std::setprecision(2) << elapsedMs
            << " ms.\n";
}

void cmdExport(ApplicationState& app, std::istringstream& iss) {
  if(!app.isLoaded) {
    std::cout << "Error: You must 'load' meshes before exporting.\n";
    return;
  }
  std::string meshTag = parseArgument(iss);
  std::string outPath = parseArgument(iss);

  if(meshTag.empty() || outPath.empty()) {
    std::cout << "Usage: export <A/B> <output.off>\n";
    return;
  }

  const Mesh* targetMesh = nullptr;
  if(meshTag == "A" || meshTag == "a")
    targetMesh = &app.meshA;
  else if(meshTag == "B" || meshTag == "b")
    targetMesh = &app.meshB;
  else {
    std::cout << "Error: Invalid mesh target '" << meshTag << "'. Choose 'A' or 'B'.\n";
    return;
  }

  auto tStart = std::chrono::high_resolution_clock::now();
  if(exportMeshToOff(*targetMesh, outPath)) {
    auto tEnd = std::chrono::high_resolution_clock::now();
    double elapsedMs = std::chrono::duration<double, std::milli>(tEnd - tStart).count();
    std::cout << "Successfully exported Mesh " << meshTag << " (" << num_vertices(*targetMesh) << " verts, "
              << num_faces(*targetMesh) << " faces) in normalized space to '" << outPath << "' in " << std::fixed
              << std::setprecision(2) << elapsedMs << " ms.\n";
  }
}

void cmdGizmo(ApplicationState& app, std::istringstream& iss) {
  int showA = 0, showB = -1;
  if(iss >> showA) {
    if(!(iss >> showB))
      showB = showA;
    PolyscopeBridge::setGizmosEnabled(showA != 0, showB != 0);
    std::cout << "Gizmo visibility updated: Mesh A = " << (showA != 0 ? "ON" : "OFF")
              << " | Mesh B = " << (showB != 0 ? "ON" : "OFF") << "\n";
  } else {
    std::cout << "Usage: gizmo <show(0/1)> [showB(0/1)]  (e.g., 'gizmo 0' hides both)\n";
  }
}

/**
 * @brief Executes classical CPU CGAL face-face intersection and updates the Polyscope viewport.
 */
void cmdComputeCGALClassic(ApplicationState& app, std::istringstream& iss) {
  if(!app.isLoaded) {
    std::cout << "Error: You must 'load' meshes before computing.\n";
    return;
  }

  // 1. Fetch current UI transformations from Polyscope
  float3 rotA, transA, rotB, transB;
  if(!PolyscopeBridge::getCurrentTransforms(rotA, transA, rotB, transB)) {
    std::cout << "Error: Failed to fetch transform matrices from Polyscope.\n";
    return;
  }

  // 2. Create transformed temporary mesh copies using RotationTools
  Mesh meshA_transformed, meshB_transformed;
  Point3 centerA(app.normCenterA.x, app.normCenterA.y, app.normCenterA.z);
  Point3 centerB(app.normCenterB.x, app.normCenterB.y, app.normCenterB.z);

  transformCgalMesh(app.meshA, meshA_transformed, centerA, make_double3(rotA.x, rotA.y, rotA.z),
                    make_double3(transA.x, transA.y, transA.z));

  transformCgalMesh(app.meshB, meshB_transformed, centerB, make_double3(rotB.x, rotB.y, rotB.z),
                    make_double3(transB.x, transB.y, transB.z));

  std::cout << "\nRunning Classical CGAL CPU Intersection...\n";

  using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;
  std::vector<std::pair<face_descriptor, face_descriptor>> cgal_intersected_tris;

  auto tStart = std::chrono::high_resolution_clock::now();

  // 3. Execute CGAL face-face intersection on transformed meshes
  CGAL::Polygon_mesh_processing::internal::compute_face_face_intersection(
      meshA_transformed, meshB_transformed, std::back_inserter(cgal_intersected_tris), CGAL::parameters::all_default(),
      CGAL::parameters::all_default());

  auto tEnd = std::chrono::high_resolution_clock::now();
  double elapsedMs = std::chrono::duration<double, std::milli>(tEnd - tStart).count();

  // Convert face descriptors to int2 pairs for Polyscope highlighting
  std::vector<int2> stdPairs;
  stdPairs.reserve(cgal_intersected_tris.size());
  for(const auto& pair : cgal_intersected_tris) {
    stdPairs.push_back(make_int2(static_cast<int>(pair.first), static_cast<int>(pair.second)));
  }

  PolyscopeBridge::highlightIntersections(stdPairs, num_faces(app.meshA), num_faces(app.meshB));

  std::cout << "=======================================================\n";
  std::cout << "         CLASSICAL CGAL CPU INTERSECTION               \n";
  std::cout << "=======================================================\n";
  std::cout << "  Execution Time      : " << std::fixed << std::setprecision(2) << elapsedMs << " ms\n";
  std::cout << "  Exact Intersections : " << cgal_intersected_tris.size() << " pairs\n";
  std::cout << "=======================================================\n\n";
}

void cmdTestBVH(ApplicationState& app, std::istringstream& iss) {
  if(!app.isLoaded) {
    std::cout << "Error: You must 'load' meshes first.\n";
    return;
  }
  int queryDescentLevel = 12, referenceDescentLevel = 12, mode = 0, leafThreshold = 4, predicateDoubleMode = 0;
  int batchMultiplier = std::numeric_limits<int>::max();
  std::string batchStr;

  iss >> queryDescentLevel >> referenceDescentLevel;

  if(iss >> batchStr) {
    if(batchStr == "INF" || batchStr == "inf") {
      batchMultiplier = std::numeric_limits<int>::max();
    } else {
      batchMultiplier = std::stoi(batchStr);
    }
    iss >> mode >> leafThreshold >> predicateDoubleMode;
  }

  float3 rotA, transA, rotB, transB;
  if(!PolyscopeBridge::getCurrentTransforms(rotA, transA, rotB, transB)) {
    std::cout << "!Error: Failed to fetch transform matrices from Polyscope.\n";
    return;
  }

  Point3 centerA(app.normCenterA.x, app.normCenterA.y, app.normCenterA.z);
  Point3 centerB(app.normCenterB.x, app.normCenterB.y, app.normCenterB.z);
  ExecutionStats testStats;
  int2* outIntersectionPairs = nullptr;
  size_t outFinalCount = 0;

  auto tStart = std::chrono::high_resolution_clock::now();
  runGridForestIntersectionPipeline(
      app.meshA, app.meshB, app.hVertsA.data(), static_cast<int>(app.hVertsA.size()), app.hIndicesA.data(),
      static_cast<int>(app.hIndicesA.size()), queryDescentLevel, app.hVertsB.data(),
      static_cast<int>(app.hVertsB.size()), app.hIndicesB.data(), static_cast<int>(app.hIndicesB.size()),
      referenceDescentLevel, batchMultiplier, mode, predicateDoubleMode, leafThreshold, testStats, outIntersectionPairs,
      outFinalCount, centerA, centerB, make_double3(rotA.x, rotA.y, rotA.z), make_double3(transA.x, transA.y, transA.z),
      make_double3(rotB.x, rotB.y, rotB.z), make_double3(transB.x, transB.y, transB.z));

  auto tEnd = std::chrono::high_resolution_clock::now();
  double elapsedMs = std::chrono::duration<double, std::milli>(tEnd - tStart).count();

  if(outIntersectionPairs && outFinalCount > 0) {
    std::vector<int2> stdPairs(outIntersectionPairs, outIntersectionPairs + outFinalCount);
    PolyscopeBridge::highlightIntersections(stdPairs, num_faces(app.meshA), num_faces(app.meshB));
    std::free(outIntersectionPairs);
  } else {
    PolyscopeBridge::highlightIntersections({}, num_faces(app.meshA), num_faces(app.meshB));
  }

  // Trigger dedicated testBVH visualizer window
  app.timingVisualizer.updateAndShowTestBVH(testStats);

  // Print dedicated terminal breakdown for testBVH
  std::cout << "\n=======================================================\n";
  std::cout << "        ALTERNATIVE BVH TEST PIPELINE STATS           \n";
  std::cout << "=======================================================\n";
  std::cout << "  Wall Clock Time        : " << std::fixed << std::setprecision(2) << elapsedMs << " ms\n";
  std::cout << "  GPU Pipeline Total     : " << std::fixed << std::setprecision(2) << testStats.GPUTotalTime << " ms\n";
  std::cout << "  -----------------------------------------------------\n";
  std::cout << "  [ALLOCATIONS & TRANSFERS]\n";
  std::cout << "  - Raw Alloc & Copy     : " << testStats.initialAllocAndCopyMs << " ms\n";
  std::cout << "  - Thrust Init Overhead : " << testStats.thrustInitOverheadMs << " ms\n";
  std::cout << "  [TREE BUILD & REFIT]\n";
  std::cout << "  - Mesh A Tree Build    : " << testStats.buildRefitMeshAMs << " ms\n";
  std::cout << "  - Mesh B Tree Build    : " << testStats.buildRefitMeshBMs << " ms\n";
  std::cout << "  [PIPELINE SEARCH & PREDICTS]\n";
  std::cout << "  - Cross-Check Engine   : " << testStats.gpuCrossCheckEngineMs << " ms\n";
  std::cout << "  - Dual-Tree Step      : " << testStats.dualTreeStepMs << " ms\n";
  std::cout << "  - Parallel DFS Descent : " << testStats.parallelDfsDescentBMs << " ms\n";
  std::cout << "  - Batched Loop Phases  : "
            << (testStats.loopTracker.assemblyPhaseMs + testStats.loopTracker.executionPhaseMs +
                testStats.loopTracker.fineEvaluationPhaseMs + testStats.loopTracker.gpuDoublePredicatesMs +
                testStats.loopTracker.DownloadAndClean + testStats.loopTracker.CPUPredicates)
            << " ms\n";
  std::cout << "  - Cleanup & Sync       : " << testStats.finalCleanupSyncMs << " ms\n";
  std::cout << "  -----------------------------------------------------\n";
  std::cout << "  Exact Intersections    : " << outFinalCount << "\n";
  std::cout << "  AABB Candidate Pairs   : " << testStats.finalAabbCandidatePairs << "\n";
  std::cout << "=======================================================\n\n";
}

void cmdStats(ApplicationState& app, std::istringstream&) {
  if(!app.isLoaded) {
    std::cout << "Error: No meshes loaded. Run 'load' first.\n";
    return;
  }
  app.edgeVisualizer.computeAndShow(app.meshA, app.meshB);
}

// --------------------------------------------------------------------
// MAIN
// --------------------------------------------------------------------
int main(int argc, char** argv) {
  std::cout << "Initializing CUDA environment...\n";
  warmupCUDA();
  std::cout << "Initialization complete.\n";

  std::cout << "Initializing Polyscope GUI...\n";
  PolyscopeBridge::init();

  ApplicationState app;
  app.edgeVisualizer.init();
  app.timingVisualizer.init();

  // Register Polyscope GUI button callback
  PolyscopeBridge::g_onFireCallback = [&app]() { executeMeshIntersectionPipeline(app); };

  // Register UI Commands
  CommandDispatcher ui;

  ui.registerCommand(
      "load", "<meshA.off> [meshB.off] [queryDescentLevel] [referenceDescentLevel] [leafThresh] [scaleToUnit(0/1)]",
      "Loads meshes via fast parallel IO, normalizes CGAL & GPU in-place, and constructs BVHs.",
      [&](std::istringstream& iss) { cmdLoad(app, iss, false); });

  ui.registerCommand(
      "loadOld", "<meshA.off> [meshB.off] [queryDescentLevel] [referenceDescentLevel] [leafThresh] [scaleToUnit(0/1)]",
      "Loads meshes via standard sequential CGAL stream loader.",
      [&](std::istringstream& iss) { cmdLoad(app, iss, true); });

 ui.registerCommand("compute", "[batchMultiplier] [DualTreeSteps] [GpuPredicates mode] [showCurves]",
                     "Syncs active viewport/gizmo transforms to GPU and executes intersection pipeline.",
                     [&](std::istringstream& iss) {
                       int batchMultiplier = std::numeric_limits<int>::max();
                       int mode = 0;
                       int enableGpuDoublePrecision = 1; // Default to 1 (true)
                       bool showCurves = false;          // Default to false

                       std::string batchStr;
                       if (iss >> batchStr) {
                         if (batchStr == "INF" || batchStr == "inf") {
                           batchMultiplier = std::numeric_limits<int>::max();
                         } else {
                           batchMultiplier = std::stoi(batchStr);
                         }
                         
                         if (iss >> mode >> enableGpuDoublePrecision) {
                           // Parse 0/1 or true/false for showCurves if provided
                           std::string curvesStr;
                           if (iss >> curvesStr) {
                             showCurves = (curvesStr == "1" || curvesStr == "true" || curvesStr == "TRUE");
                           }
                         }
                       }

                       executeMeshIntersectionPipeline(app, batchMultiplier, mode, enableGpuDoublePrecision, showCurves);
                     });

  ui.registerCommand(
      "transform",
      "<rotAx> <rotAy> <rotAz> <transAx> <transAy> <transAz> <rotBx> <rotBy> <rotBz> <transBx> <transBy> <transBz>",
      "Sets transformation matrices for Mesh A and Mesh B.", [&](std::istringstream& iss) { cmdTransform(app, iss); });

  ui.registerCommand("translate", "<x> <y> <z>", "Translates Mesh B directly via GPU kernel shifting.",
                     [&](std::istringstream& iss) { cmdTranslate(app, iss); });

  ui.registerCommand("gizmo", "<showA(0/1)> [showB(0/1)]", "Shows or hides 3D drag gizmos in the active viewport.",
                     [&](std::istringstream& iss) { cmdGizmo(app, iss); });

  ui.registerCommand(
      "standaloneCompute",
      "[queryDescentLevel] [referenceDescentLevel] [batchMultiplier] [NumDualTreeSteps] [leafThresh] [GpuPredicates "
      "mode]",
      "Computes intersection between two meshes using alternative routine by constructing new partial BVH",
      [&](std::istringstream& iss) { cmdTestBVH(app, iss); });

  ui.registerCommand("export", "<mesh_tag(A/B)> <out_filename.off>",
                     "Exports Mesh A or Mesh B in full 17-digit precision to OFF format.",
                     [&](std::istringstream& iss) { cmdExport(app, iss); });

  ui.registerCommand("clear", "", "Frees GPU memory buffers and cuBQL structures without wiping CPU mesh data.",
                     [&](std::istringstream& iss) { cmdClear(app, iss); });

  ui.registerCommand("reconstruct", "<queryDescentLevel> <referenceDescentLevel>",
                     "Re-allocates GPU buffers and rebuilds BVHs using cached CPU data.",
                     [&](std::istringstream& iss) { cmdReconstruct(app, iss); });

  ui.registerCommand("stats", "", "Displays high-precision geometric statistics and visualizes degeneracies.",
                     [&](std::istringstream& iss) { cmdStats(app, iss); });

  ui.registerCommand("scale", "", "Prints current scale factor, scene bounds center, and mesh centroids.",
                     [&](std::istringstream& iss) { cmdScale(app, iss); });

  ui.registerCommand("tbb", "<num_threads>", "Sets active CPU worker thread limit for TBB operations.",
                     [&](std::istringstream& iss) { cmdTbb(app, iss); });

  ui.registerCommand("help", "", "Displays this command usage reference.",
                     [&](std::istringstream&) { ui.printHelp(); });

  ui.registerCommand("profile", "", "Displays graphical stacked timing bar chart for execution phases.",
                     [&](std::istringstream&) {
                       if(app.isLoaded)
                         app.timingVisualizer.show();
                     });


  ui.registerCommand("ComputeCGALClassic", "",
                     "Computes intersections using classical CPU CGAL PMP and highlights them.",
                     [&](std::istringstream& iss) { cmdComputeCGALClassic(app, iss); });

  ui.registerCommand("testConfig",
                     "[steps] [trans] [rotDeg] [precision] [dualTreeSteps] [queryLvl] [refLvl] [batch] [seed]",
                     "Displays or updates test parameters for the automated test suite.",
                     [&](std::istringstream& iss) { cmdTestConfig(app, iss); });

  ui.registerCommand("runSuite", "", "Executes automated test suite sweep with active TestConfig parameters.",
                     [&](std::istringstream&) {
                       TestSuite suite(app);
                       suite.runSuite(g_testConfig);
                     });
  ui.registerCommand("HelpPredicateModes", "", 
                     "Displays the available GPU precision modes for exact predicates.",
                     [&](std::istringstream&) {
                       std::cout << "\n=======================================================\n";
                       std::cout << "               PREDICATE COMPUTE MODES                 \n";
                       std::cout << "=======================================================\n";
                       std::cout << "  Mode 0 : Float predicates -> CPU\n";
                       std::cout << "  Mode 1 : Float predicates -> Simulated double -> CPU\n";
                       std::cout << "  Mode 2 : Float predicates -> Inexact Big integer quantization (Experimental)\n";
                       std::cout << "  Mode 3 : Float predicates -> Exact Big integer quantization -> CPU (Experimental)\n";
                       std::cout << "  Mode 4 : Float predicates -> Real double -> CPU\n";
                       std::cout << "=======================================================\n\n";
                     });

  // Optional: Add a shorter alias so you can just type 'modes' in the terminal
  ui.registerAlias("modes", "HelpPredicateModes");

  ui.registerCommand("quit", "", "Exits the application.", [&](std::istringstream&) { g_running = false; });

  // Register shortcuts and aliases
  ui.registerAlias("cgalclassic", "ComputeCGALClassic");
  ui.registerAlias("cgal", "ComputeCGALClassic");

  // Register command aliases
  ui.registerAlias("loadold", "loadOld");
  ui.registerAlias("LoadOld", "loadOld");
  ui.registerAlias("rebuild", "reconstruct");
  ui.registerAlias("testbvh", "testBVH");
  ui.registerAlias("exit", "quit");

  std::cout << "Type 'help' for commands.\n";

  std::thread inputThread(inputThreadWorker);

  while(g_running) {
    std::string line;
    bool haveLine = false;
    {
      std::lock_guard<std::mutex> lock(g_cmdMutex);
      if(!g_cmdQueue.empty()) {
        line = g_cmdQueue.front();
        g_cmdQueue.pop();
        haveLine = true;
      }
    }

    if(haveLine) {
      ui.execute(line);

      {
        std::lock_guard<std::mutex> lock(g_cmdMutex);
        g_readyForInput = true;
      }
      g_cmdCV.notify_one();
    }

    PolyscopeBridge::drawFrame();
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
  }

  std::cout << "Cleaning up and exiting...\n";
  if(inputThread.joinable()) {
    if(g_running == false) {
      inputThread.detach();
    } else {
      inputThread.join();
    }
  }

  app.controller.cleanup();
  return 0;
}