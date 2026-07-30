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

// GNU Readline headers for Fedora terminal arrow key & history support
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
#include "cuBQL/bvh.h"

// Custom Controller and Helpers
#include "../src/GPUIntersector/KernelBVHController.h"
#include "../src/CPU/CgalDefinitions.h"
#include "../src/CPU/PolyscopeBridge.h"
#include "../src/CPU/ParallelCgalOffLoader.h"
#include "../src/Warmup/cuda_warmup.h"

// --- Minimal Double Box for Mesh Bounds ---
struct DoubleBox
{
  double min_x = std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double min_z = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();
  double max_z = -std::numeric_limits<double>::infinity();

  void grow(double x, double y, double z) {
    min_x = std::min(min_x, x);
    max_x = std::max(max_x, x);
    min_y = std::min(min_y, y);
    max_y = std::max(max_y, y);
    min_z = std::min(min_z, z);
    max_z = std::max(max_z, z);
  }

  void grow(const DoubleBox& o) {
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

DoubleBox computeMeshBoundsTBB(const Mesh& mesh) {
  size_t numVerts = num_vertices(mesh);
  auto coords = mesh.points();

  return tbb::parallel_reduce(
      tbb::blocked_range<size_t>(0, numVerts), DoubleBox(),
      [&](const tbb::blocked_range<size_t>& r, DoubleBox box) {
        for(size_t i = r.begin(); i != r.end(); ++i) {
          typename Mesh::Vertex_index v(static_cast<Mesh::size_type>(i));
          const auto& p = coords[v];
          box.grow(p.x(), p.y(), p.z());
        }
        return box;
      },
      [](DoubleBox a, const DoubleBox& b) {
        a.grow(b);
        return a;
      });
}

void normalizeCgalMeshTBB(Mesh& mesh, double cx, double cy, double cz, double scaleFactor) {
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

bool exportCgalMeshToOff(const Mesh& mesh, const std::string& filepath) {
  std::ofstream outFile(filepath);
  if(!outFile.is_open()) {
    std::cerr << "Error: Could not open file for writing: " << filepath << "\n";
    return false;
  }

  outFile << std::setprecision(17);

  if(!(outFile << mesh)) {
    std::cerr << "Error: Failed writing CGAL mesh data to: " << filepath << "\n";
    return false;
  }

  return true;
}

// Synchronizes normalized CGAL vertices directly to GPU double3 buffers without floating-point downcasting.
void syncNormalizedVertsTBB(const Mesh& mesh, std::vector<double3>& outVerts) {
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

// Old sequential CGAL loading method via standard std::ifstream and halfedge topology traversal
bool loadOffOldCGAL(const std::string& filepath,
                    Mesh& mesh,
                    std::vector<double3>& outVerts,
                    std::vector<uint3>& outIndices) {
  std::ifstream inFile(filepath);
  if(!inFile.is_open()) {
    std::cerr << "Error: Could not open file: " << filepath << "\n";
    return false;
  }

  if(!(inFile >> mesh)) {
    std::cerr << "Error: Failed reading CGAL mesh stream from: " << filepath << "\n";
    return false;
  }

  // Extract initial double3 vertex coordinates
  size_t numVerts = num_vertices(mesh);
  outVerts.resize(numVerts);
  auto coords = mesh.points();
  for(size_t i = 0; i < numVerts; ++i) {
    typename Mesh::Vertex_index v(static_cast<Mesh::size_type>(i));
    const auto& p = coords[v];
    outVerts[i] = make_double3(p.x(), p.y(), p.z());
  }

  // Extract face indices via CGAL halfedge navigation
  size_t numFaces = num_faces(mesh);
  outIndices.resize(numFaces);
  size_t fIdx = 0;
  for(auto f : mesh.faces()) {
    auto h = mesh.halfedge(f);
    uint3 tri;
    tri.x = static_cast<uint32_t>(mesh.source(h));
    tri.y = static_cast<uint32_t>(mesh.target(h));
    tri.z = static_cast<uint32_t>(mesh.target(mesh.next(h)));
    outIndices[fIdx++] = tri;
  }

  return true;
}

void printHelp() {
  std::cout << "\nAvailable Commands:\n";
  std::cout << "  load <meshA.off> [meshB.off] [maxCellA] [maxCellB] [leafThresh] [scaleToUnit(0/1)] [useDoublePreds(0/1)]\n";
  std::cout << "      - Loads meshes via fast parallel IO, normalizes CGAL & GPU in-place, and constructs BVHs.\n";
  std::cout << "  loadOld <meshA.off> [meshB.off] [maxCellA] [maxCellB] [leafThresh] [scaleToUnit(0/1)] [useDoublePreds(0/1)]\n";
  std::cout << "      - Loads meshes via standard CGAL stream loader (sequential std::ifstream).\n";
  std::cout << "  scale\n";
  std::cout << "      - Prints current scale factor, scene center, and individual mesh centroids.\n";
  std::cout << "  gizmo <show(0/1)> [showB(0/1)]\n";
  std::cout << "      - Shows or hides 3D drag gizmos in the viewport (e.g., 'gizmo 0' hides both).\n";
  std::cout << "  transform <rotAx> <rotAy> <rotAz> <rotBx> <rotBy> <rotBz> <transBx> <transBy> <transBz>\n";
  std::cout << "      - Transforms Mesh A (rotations) and Mesh B (rotations + translation) in GPU & Polyscope.\n";
  std::cout << "  translate <x> <y> <z>\n";
  std::cout << "      - Sets translation for Mesh B using GPU kernel shifting.\n";
  std::cout << "  translateCPU <x> <y> <z>\n";
  std::cout << "      - Sets translation for Mesh B via double-precision CPU recalculation & GPU upload.\n";
  std::cout << "  export <mesh_tag(A/B)> <out_filename.off>\n";
  std::cout << "      - Exports current Mesh A or B in full 17-digit precision to OFF format.\n";
  std::cout << "  compute [batchMultiplier] [DualTreeSteps] [async] [useDoublePreds(0/1)]\n";
  std::cout << "      - Runs intersection pipeline using explicit set transforms.\n";
  std::cout << "  fire [batchMultiplier] [DualTreeSteps] [async] [useDoublePreds(0/1)]\n";
  std::cout << "      - FireNow / Fire command: Syncs current viewport gizmo transforms directly to GPU & computes.\n";
  std::cout << "  tbb <num_threads>\n";
  std::cout << "      - Sets the active CPU thread limit for TBB worker operations.\n";
  std::cout << "  help\n";
  std::cout << "      - Displays this message.\n";
  std::cout << "  quit / exit\n";
  std::cout << "      - Exits the application.\n";
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

// Thread-safe command queue shared between input thread and main thread
static std::mutex g_cmdMutex;
static std::condition_variable g_cmdCV;
static std::queue<std::string> g_cmdQueue;
static std::atomic<bool> g_running{true};
static bool g_readyForInput{true};

void inputThreadFunc() {
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
    if(!line.empty()) {
      add_history(rawInput);
    }
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

int main(int argc, char** argv) {
  std::cout << "Initializing CUDA environment...\n";
  warmupCUDA();
  std::cout << "Initialization complete.\n";

  std::cout << "Initializing Polyscope GUI...\n";
  PolyscopeBridge::init();
  std::cout << "Type 'help' for commands.\n";

  std::unique_ptr<tbb::global_control> tbbControl;

  KernelBVHController controller;
  ExecutionStats stats;

  Mesh meshA, meshB;
  std::vector<double3> hVertsA, hVertsB;
  std::vector<uint3> hIndicesA, hIndicesB;

  bool isLoaded = false;

  double currentScaleFactor = 1.0;
  double currentCenterX = 0.0, currentCenterY = 0.0, currentCenterZ = 0.0;
  double currentMaxSpan = 0.0;
  bool currentScaledToUnit = false;

  double3 origCenterA{0.0, 0.0, 0.0}, origCenterB{0.0, 0.0, 0.0};
  float3 normCenterA{0.0f, 0.0f, 0.0f}, normCenterB{0.0f, 0.0f, 0.0f};

  // Helper lambda to run the pipeline with current drag-and-drop / gizmo positions
  // Helper lambda to run the pipeline with current drag-and-drop / gizmo positions
  auto executeFireNow = [&](int batchMultiplier = std::numeric_limits<int>::max(), int mode = 0,
                            int activateAsyncDownload = 0,
                            bool useDoublePreds = false) { // <--- Added parameter
    if(!isLoaded) {
      std::cout << "Error: You must 'load' meshes before computing.\n";
      return;
    }

    float3 rotA, transA, rotB, transB;
    if(!PolyscopeBridge::getCurrentTransforms(rotA, transA, rotB, transB)) {
      std::cout << "Error: Failed to fetch transform matrices from Polyscope.\n";
      return;
    }

    std::cout << "\n[FireNow] Syncing Viewport Transforms to GPU...\n";
    std::cout << "  Mesh A Rot: (" << rotA.x << ", " << rotA.y << ", " << rotA.z << ") | Trans: (" << transA.x << ", "
              << transA.y << ", " << transA.z << ")\n";
    std::cout << "  Mesh B Rot: (" << rotB.x << ", " << rotB.y << ", " << rotB.z << ") | Trans: (" << transB.x << ", "
              << transB.y << ", " << transB.z << ")\n";

    // 1. Dispatch rotation and translation to double3 GPU controller
    auto tStart = std::chrono::high_resolution_clock::now();
    controller.setTransformBoth(make_double3(rotA.x, rotA.y, rotA.z), make_double3(transA.x, transA.y, transA.z),
                                make_double3(rotB.x, rotB.y, rotB.z), make_double3(transB.x, transB.y, transB.z));

    // 2. Execute GPU pipeline with double predicate flag
    tbb::concurrent_vector<int2> finalExactPairs;

    controller.runIntersectionPipeline(batchMultiplier, mode, activateAsyncDownload, finalExactPairs, stats,
                                       useDoublePreds); // <--- Passed flag

    auto tEnd = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(tEnd - tStart).count();

    // 3. Highlight intersections in viewport
    std::vector<int2> stdPairs(finalExactPairs.begin(), finalExactPairs.end());
    PolyscopeBridge::highlightIntersections(stdPairs, num_faces(meshA), num_faces(meshB));

    std::cout << "Query completed in time: " << std::fixed << std::setprecision(2) << ms << " ms. "
              << "We got " << finalExactPairs.size() << " intersections.\n";
    std::cout << "AABB Hits: " << stats.finalAabbCandidatePairs
              << " | Greenlist : " << stats.loopTracker.confirmedGreenPairs
              << " | Yellowlist: " << stats.loopTracker.confirmedYellowPairs
              << " | Orangelist: " << stats.loopTracker.confirmedOrangePairs << std::endl;
  };

  // Connect GUI Button callback in Polyscope window to FireNow logic
  PolyscopeBridge::g_onFireCallback = [&]() { executeFireNow(); };

  std::thread inputThread(inputThreadFunc);

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
      std::istringstream iss(line);
      std::string cmd;
      if(!(iss >> cmd)) {
        // empty line
      } else if(cmd == "quit" || cmd == "exit") {
        g_running = false;
      } else if(cmd == "help") {
        printHelp();
      } else if(cmd == "scale") {
        if(!isLoaded) {
          std::cout << "Error: No meshes loaded. Run 'load' first.\n";
        } else {
          std::cout << "\n=======================================================\n";
          std::cout << "         MESH NORMALIZATION & SCALING INFO             \n";
          std::cout << "=======================================================\n";
          std::cout << "  Scale Factor Applied : " << std::setprecision(12) << currentScaleFactor << "\n";
          std::cout << "  Original Scene Center: (" << std::setprecision(8) << currentCenterX << ", " << currentCenterY
                    << ", " << currentCenterZ << ")\n";
          std::cout << "  Max Original Extent  : " << currentMaxSpan << " units\n";
          std::cout << "  Coordinate Space     : "
                    << (currentScaledToUnit ? "[-1.0, 1.0]^3 (Centered & Scaled)" : "Centered (Unscaled)") << "\n";
          std::cout << "  -----------------------------------------------------\n";
          std::cout << "  Original Centroid A  : (" << origCenterA.x << ", " << origCenterA.y << ", " << origCenterA.z
                    << ")\n";
          std::cout << "  Original Centroid B  : (" << origCenterB.x << ", " << origCenterB.y << ", " << origCenterB.z
                    << ")\n";
          std::cout << "  Normalized Centroid A: (" << normCenterA.x << ", " << normCenterA.y << ", " << normCenterA.z
                    << ")\n";
          std::cout << "  Normalized Centroid B: (" << normCenterB.x << ", " << normCenterB.y << ", " << normCenterB.z
                    << ")\n";
          std::cout << "=======================================================\n\n";
        }
      } else if(cmd == "tbb") {
        int numThreads = 0;
        if(iss >> numThreads && numThreads > 0) {
          tbbControl = std::make_unique<tbb::global_control>(tbb::global_control::max_allowed_parallelism, numThreads);
          std::cout << "TBB maximum worker limit updated to: " << numThreads << "\n";
        } else {
          std::cout << "Usage: tbb <num_threads> (must be > 0)\n";
        }
      } else if(cmd == "load" || cmd == "loadOld" || cmd == "loadold" || cmd == "LoadOld") {
        bool useOldLoader = (cmd == "loadOld" || cmd == "loadold" || cmd == "LoadOld");

        std::string arg1 = parseArgument(iss);
        std::string arg2 = parseArgument(iss);

        int maxCellA = 1, maxCellB = 1, leafThresh = 4, scaleToUnitInt = 0,
            useDoublePredsInt = 0; // <--- ADDED useDoublePredsInt

        if(arg1.empty()) {
          std::cout << "Usage: " << cmd
                    << " <meshA> [meshB] [maxCellA] [maxCellB] [leafThresh] [scaleToUnit(0/1)] [useDoublePreds(0/1)]\n";
        } else {
          std::string pathA, pathB;

          bool singleMeshMode = false;
          if(arg2.empty()) {
            singleMeshMode = true;
          } else {
            std::istringstream checkNum(arg2);
            int val;
            if(checkNum >> val) {
              singleMeshMode = true;
              maxCellA = val;
            }
          }

          if(singleMeshMode) {
            pathA = arg1;
            pathB = arg1;
            iss >> maxCellB >> leafThresh >> scaleToUnitInt >> useDoublePredsInt; // <--- Read flag
          } else {
            pathA = arg1;
            pathB = arg2;
            iss >> maxCellA >> maxCellB >> leafThresh >> scaleToUnitInt >> useDoublePredsInt; // <--- Read flag
          }

          bool scaleToUnit = (scaleToUnitInt != 0);
          bool useDoublePreds = (useDoublePredsInt != 0); // <--- Bool toggle

          controller.cleanup();
          meshA = Mesh();
          meshB = Mesh();
          hVertsA.clear();
          hIndicesA.clear();
          hVertsB.clear();
          hIndicesB.clear();
          isLoaded = false;

          std::cout << "Loading mesh(es) using "
                    << (useOldLoader ? "Old CGAL Stream Loader" : "Fast Parallel IO Loader") << "...\n";

          auto tIoStart = std::chrono::high_resolution_clock::now();
          bool loadedSuccessfully = false;

          if(useOldLoader) {
            // Old sequential std::ifstream >> CGAL mesh loader
            if(singleMeshMode) {
              if(loadOffOldCGAL(pathA, meshA, hVertsA, hIndicesA)) {
                meshB = meshA;
                hVertsB = hVertsA;
                hIndicesB = hIndicesA;
                loadedSuccessfully = true;
              }
            } else {
              if(loadOffOldCGAL(pathA, meshA, hVertsA, hIndicesA) && loadOffOldCGAL(pathB, meshB, hVertsB, hIndicesB)) {
                loadedSuccessfully = true;
              }
            }
          } else {
            // Fast POSIX mmap + TBB parsing directly into CGAL Mesh
            if(singleMeshMode) {
              if(ParallelIO::loadOffToCgalMesh<Mesh, Point3>(pathA, meshA, hIndicesA)) {
                meshB = meshA;
                hIndicesB = hIndicesA;
                loadedSuccessfully = true;
              } else {
                std::cerr << "Error: Could not open/read Mesh file: " << pathA << "\n";
              }
            } else {
              if(ParallelIO::loadOffToCgalMesh<Mesh, Point3>(pathA, meshA, hIndicesA) &&
                 ParallelIO::loadOffToCgalMesh<Mesh, Point3>(pathB, meshB, hIndicesB))
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

            DoubleBox boxA_orig = computeMeshBoundsTBB(meshA);
            DoubleBox boxB_orig = computeMeshBoundsTBB(meshB);
            origCenterA = boxA_orig.getCenter();
            origCenterB = boxB_orig.getCenter();

            DoubleBox sceneBox = boxA_orig;
            sceneBox.grow(boxB_orig);

            double cx = 0.5 * (sceneBox.min_x + sceneBox.max_x);
            double cy = 0.5 * (sceneBox.min_y + sceneBox.max_y);
            double cz = 0.5 * (sceneBox.min_z + sceneBox.max_z);

            double spanX = sceneBox.max_x - sceneBox.min_x;
            double spanY = sceneBox.max_y - sceneBox.min_y;
            double spanZ = sceneBox.max_z - sceneBox.min_z;
            double maxSpan = std::max({spanX, spanY, spanZ});

            double scaleFactor = 1.0;
            if(scaleToUnit && maxSpan > 0.0) {
              scaleFactor = 2.0 / maxSpan;
            }

            currentScaleFactor = scaleFactor;
            currentCenterX = cx;
            currentCenterY = cy;
            currentCenterZ = cz;
            currentMaxSpan = maxSpan;
            currentScaledToUnit = scaleToUnit;

            // Normalize CGAL mesh points
            normalizeCgalMeshTBB(meshA, cx, cy, cz, scaleFactor);
            normalizeCgalMeshTBB(meshB, cx, cy, cz, scaleFactor);

            DoubleBox boxA_norm = computeMeshBoundsTBB(meshA);
            DoubleBox boxB_norm = computeMeshBoundsTBB(meshB);
            double3 cA_norm = boxA_norm.getCenter();
            double3 cB_norm = boxB_norm.getCenter();

            normCenterA = make_float3(static_cast<float>(cA_norm.x), static_cast<float>(cA_norm.y),
                                      static_cast<float>(cA_norm.z));
            normCenterB = make_float3(static_cast<float>(cB_norm.x), static_cast<float>(cB_norm.y),
                                      static_cast<float>(cB_norm.z));

            // Sync normalized double3 vertex coordinates directly
            syncNormalizedVertsTBB(meshA, hVertsA);
            syncNormalizedVertsTBB(meshB, hVertsB);

            auto tPrepEnd = std::chrono::high_resolution_clock::now();
            double prepMs = std::chrono::duration<double, std::milli>(tPrepEnd - tPrepStart).count();

            auto tStart = std::chrono::high_resolution_clock::now();

            Point3 cA_cgal(cA_norm.x, cA_norm.y, cA_norm.z);
            Point3 cB_cgal(cB_norm.x, cB_norm.y, cB_norm.z);

            // Execute double-precision BVH construction overload
            controller.construct(meshA, meshB, cA_cgal, cB_cgal, hVertsA.data(), static_cast<int>(hVertsA.size()),
                                 hIndicesA.data(), static_cast<int>(hIndicesA.size()), maxCellA, hVertsB.data(),
                                 static_cast<int>(hVertsB.size()), hIndicesB.data(), static_cast<int>(hIndicesB.size()),
                                 maxCellB, leafThresh, stats, useDoublePreds);

            auto tEnd = std::chrono::high_resolution_clock::now();
            double ms = std::chrono::duration<double, std::milli>(tEnd - tStart).count();

            PolyscopeBridge::reset(meshA, meshB, normCenterA, normCenterB);

            std::cout << "\n--- Load & Construction Performance Breakdown ---\n";
            std::cout << "  Loader Used:               "
                      << (useOldLoader ? "Standard CGAL Stream (loadOld)" : "Parallel IO (load)") << "\n";
            std::cout << "  1. Disk -> IO Load:        " << std::fixed << std::setprecision(2) << ioMs << " ms\n";
            std::cout << "  2. Topology/Layout Prep:   " << std::fixed << std::setprecision(2) << prepMs << " ms\n";
            std::cout << "  3. BVH Construction:       " << std::fixed << std::setprecision(2) << ms << " ms\n";
            std::cout << "  -----------------------------------------------\n";
            std::cout << "  Total Preparation Time:    " << std::fixed << std::setprecision(2) << (ioMs + prepMs + ms)
                      << " ms\n";
            std::cout << "  Applied Scale Factor:      " << std::setprecision(10) << scaleFactor
                      << (scaleToUnit ? " (Scaled to [-1, 1]^3)" : " (Unscaled)") << "\n";
            std::cout << "  Normalized Center A:       (" << normCenterA.x << ", " << normCenterA.y << ", "
                      << normCenterA.z << ")\n";
            std::cout << "  Normalized Center B:       (" << normCenterB.x << ", " << normCenterB.y << ", "
                      << normCenterB.z << ")\n\n";

            isLoaded = true;
          }
        }
      } else if(cmd == "transform") {
        if(!isLoaded) {
          std::cout << "Error: You must 'load' meshes before transforming.\n";
        } else {
          double rotAx, rotAy, rotAz;
          double rotBx, rotBy, rotBz;
          double transBx, transBy, transBz;

          if(iss >> rotAx >> rotAy >> rotAz >> rotBx >> rotBy >> rotBz >> transBx >> transBy >> transBz) {
            double3 rotA = make_double3(rotAx, rotAy, rotAz);
            double3 transA = make_double3(0.0, 0.0, 0.0);
            double3 rotB = make_double3(rotBx, rotBy, rotBz);
            double3 transB = make_double3(transBx, transBy, transBz);

            auto tStart = std::chrono::high_resolution_clock::now();

            controller.setTransformBoth(rotA, transA, rotB, transB);

            float3 fRotA = make_float3(static_cast<float>(rotAx), static_cast<float>(rotAy), static_cast<float>(rotAz));
            float3 fTransA = make_float3(0.0f, 0.0f, 0.0f);
            float3 fRotB = make_float3(static_cast<float>(rotBx), static_cast<float>(rotBy), static_cast<float>(rotBz));
            float3 fTransB =
                make_float3(static_cast<float>(transBx), static_cast<float>(transBy), static_cast<float>(transBz));

            PolyscopeBridge::transformBoth(fRotA, fTransA, normCenterA, fRotB, fTransB, normCenterB);

            auto tEnd = std::chrono::high_resolution_clock::now();
            double ms = std::chrono::duration<double, std::milli>(tEnd - tStart).count();

            std::cout << "Transformations applied in " << std::fixed << std::setprecision(2) << ms << " ms:\n";
            std::cout << "  Mesh A Rot: (" << rotAx << ", " << rotAy << ", " << rotAz << ")\n";
            std::cout << "  Mesh B Rot: (" << rotBx << ", " << rotBy << ", " << rotBz << ") | Trans: (" << transBx
                      << ", " << transBy << ", " << transBz << ")\n";
          } else {
            std::cout
                << "Usage: transform <rotAx> <rotAy> <rotAz> <rotBx> <rotBy> <rotBz> <transBx> <transBy> <transBz>\n";
          }
        }
      } else if(cmd == "translate") {
        double x, y, z;
        if(iss >> x >> y >> z) {
          auto tStart = std::chrono::high_resolution_clock::now();

          controller.setTranslation(x, y, z);
          PolyscopeBridge::translateMeshB(static_cast<float>(x), static_cast<float>(y), static_cast<float>(z));

          auto tEnd = std::chrono::high_resolution_clock::now();
          double elapsedMs = std::chrono::duration<double, std::milli>(tEnd - tStart).count();

          std::cout << "Translation applied (GPU Kernel): (" << x << ", " << y << ", " << z << ") "
                    << "in " << elapsedMs << " ms\n";
        } else {
          std::cout << "Usage: translate <x> <y> <z>\n";
        }
      } else if(cmd == "translateCPU") {
        double x, y, z;
        if(iss >> x >> y >> z) {
          auto tStart = std::chrono::high_resolution_clock::now();

          controller.setTranslationCPUHostUpload(x, y, z);
          PolyscopeBridge::translateMeshB(static_cast<float>(x), static_cast<float>(y), static_cast<float>(z));

          auto tEnd = std::chrono::high_resolution_clock::now();
          double elapsedMs = std::chrono::duration<double, std::milli>(tEnd - tStart).count();

          std::cout << "Translation applied (CPU Host Upload): (" << x << ", " << y << ", " << z << ") "
                    << "in " << elapsedMs << " ms\n";
        } else {
          std::cout << "Usage: translateCPU <x> <y> <z>\n";
        }
      } else if(cmd == "compute") {
        if(!isLoaded) {
          std::cout << "Error: You must 'load' meshes before computing.\n";
        } else {
          int batchMultiplier = std::numeric_limits<int>::max();
          int mode = 0;
          int activateAsyncDownload = 0;
          int useDoublePredsInt = 0; // <--- ADDED

          iss >> batchMultiplier >> mode >> activateAsyncDownload >> useDoublePredsInt; // <--- READ FLAG

          tbb::concurrent_vector<int2> finalExactPairs;

          auto tStart = std::chrono::high_resolution_clock::now();

          controller.runIntersectionPipeline(batchMultiplier, mode, activateAsyncDownload, finalExactPairs, stats,
                                             (useDoublePredsInt != 0)); // <--- PASSED FLAG

          auto tEnd = std::chrono::high_resolution_clock::now();
          double ms = std::chrono::duration<double, std::milli>(tEnd - tStart).count();

          std::vector<int2> stdPairs(finalExactPairs.begin(), finalExactPairs.end());
          PolyscopeBridge::highlightIntersections(stdPairs, num_faces(meshA), num_faces(meshB));

          std::cout << "Query completed in time: " << std::fixed << std::setprecision(2) << ms << " ms. "
                    << "We got " << finalExactPairs.size() << " intersections.\n";
          std::cout << "AABB Hits: " << stats.finalAabbCandidatePairs
                    << " | Greenlist : " << stats.loopTracker.confirmedGreenPairs
                    << " | Yellowlist: " << stats.loopTracker.confirmedYellowPairs 
                    << " | Orangelist: " << stats.loopTracker.confirmedOrangePairs 
                    << std::endl;
        }
      } else if(cmd == "fire" || cmd == "FireNow") {
        int batchMultiplier = std::numeric_limits<int>::max();
        int mode = 0;
        int activateAsyncDownload = 0;
        int useDoublePredsInt = 0; // <--- ADDED

        iss >> batchMultiplier >> mode >> activateAsyncDownload >> useDoublePredsInt;           // <--- READ FLAG
        executeFireNow(batchMultiplier, mode, activateAsyncDownload, (useDoublePredsInt != 0)); // <--- PASSED FLAG
      } else if(cmd == "export") {
        if(!isLoaded) {
          std::cout << "Error: You must 'load' meshes before exporting.\n";
        } else {
          std::string meshTag = parseArgument(iss);
          std::string outPath = parseArgument(iss);

          if(meshTag.empty() || outPath.empty()) {
            std::cout << "Usage: export <A/B> <output.off>\n";
          } else {
            const Mesh* targetMesh = nullptr;
            if(meshTag == "A" || meshTag == "a") {
              targetMesh = &meshA;
            } else if(meshTag == "B" || meshTag == "b") {
              targetMesh = &meshB;
            } else {
              std::cout << "Error: Invalid mesh target '" << meshTag << "'. Choose 'A' or 'B'.\n";
            }

            if(targetMesh) {
              auto tStart = std::chrono::high_resolution_clock::now();

              if(exportCgalMeshToOff(*targetMesh, outPath)) {
                auto tEnd = std::chrono::high_resolution_clock::now();
                double ms = std::chrono::duration<double, std::milli>(tEnd - tStart).count();

                std::cout << "Successfully exported Mesh " << meshTag << " (" << num_vertices(*targetMesh) << " verts, "
                          << num_faces(*targetMesh) << " faces) in normalized space to '" << outPath << "' in "
                          << std::fixed << std::setprecision(2) << ms << " ms.\n";
              }
            }
          }
        }
      } else if(cmd == "gizmo") {
        int showA = 0, showB = -1;
        if(iss >> showA) {
          if(!(iss >> showB)) {
            showB = showA; // If only one argument given (e.g. 'gizmo 0'), apply to both
          }
          PolyscopeBridge::setGizmosEnabled(showA != 0, showB != 0);
          std::cout << "Gizmo visibility updated: Mesh A = " << (showA != 0 ? "ON" : "OFF")
                    << " | Mesh B = " << (showB != 0 ? "ON" : "OFF") << "\n";
        } else {
          std::cout << "Usage: gizmo <show(0/1)> [showB(0/1)]  (e.g., 'gizmo 0' hides both)\n";
        }
      } else {
        std::cout << "Unknown command: '" << cmd << "'. Type 'help' for a list of commands.\n";
      }

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

  controller.cleanup();
  return 0;
}