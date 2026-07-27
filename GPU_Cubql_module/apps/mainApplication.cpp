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

// Mutates CGAL mesh points in place to center and optionally scale them
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

// Extracts host float3 arrays from an already-normalized CGAL mesh
void extractMeshTopology(const Mesh& mesh,
                         std::vector<float3>& outVerts,
                         std::vector<float>& outVertexErrors,
                         std::vector<uint3>& outIndices,
                         bool usePreciseBounds) {
  size_t numVerts = num_vertices(mesh);
  size_t numFaces = num_faces(mesh);

  outVerts.resize(numVerts);
  outIndices.resize(numFaces);
  if(usePreciseBounds)
    outVertexErrors.resize(numVerts);
  else
    outVertexErrors.clear();

  auto coords = mesh.points();

  tbb::parallel_for(tbb::blocked_range<size_t>(0, numVerts), [&](const tbb::blocked_range<size_t>& r) {
    for(size_t i = r.begin(); i != r.end(); ++i) {
      typename Mesh::Vertex_index v(static_cast<Mesh::size_type>(i));
      const auto& p = coords[v];

      float fx = static_cast<float>(p.x());
      float fy = static_cast<float>(p.y());
      float fz = static_cast<float>(p.z());

      outVerts[i] = {fx, fy, fz};

      if(usePreciseBounds) {
        // Quantization error bound resulting from double -> float conversion
        double err_x = std::abs(p.x() - static_cast<double>(fx));
        double err_y = std::abs(p.y() - static_cast<double>(fy));
        double err_z = std::abs(p.z() - static_cast<double>(fz));
        double max_err = std::max({err_x, err_y, err_z});
        float max_mag = std::max({std::abs(fx), std::abs(fy), std::abs(fz)});
        outVertexErrors[i] = static_cast<float>(max_err) + max_mag * std::numeric_limits<float>::epsilon();
      }
    }
  });

  tbb::parallel_for(tbb::blocked_range<size_t>(0, numFaces), [&](const tbb::blocked_range<size_t>& r) {
    for(size_t i = r.begin(); i != r.end(); ++i) {
      typename Mesh::Face_index f(static_cast<Mesh::size_type>(i));
      auto h0 = mesh.halfedge(f);
      auto h1 = mesh.next(h0);
      auto h2 = mesh.next(h1);
      outIndices[i] = {static_cast<uint32_t>(mesh.target(h0)), static_cast<uint32_t>(mesh.target(h1)),
                       static_cast<uint32_t>(mesh.target(h2))};
    }
  });
}

void printHelp() {
  std::cout << "\nAvailable Commands:\n";
  std::cout
      << "  load <meshA.off> [meshB.off] [maxCellA] [maxCellB] [leafThresh] [usePrecise(0/1)] [scaleToUnit(0/1)]\n";
  std::cout << "      - Loads meshes, normalizes CGAL & GPU in-place, and constructs BVHs.\n";
  std::cout << "  scale\n";
  std::cout << "      - Prints current scale factor, scene center, and individual mesh centroids.\n";
  std::cout << "  transform <rotAx> <rotAy> <rotAz> <rotBx> <rotBy> <rotBz> <transBx> <transBy> <transBz>\n";
  std::cout << "      - Transforms Mesh A (rotations) and Mesh B (rotations + translation) in GPU & Polyscope.\n";
  std::cout << "  translate <x> <y> <z>\n";
  std::cout << "      - Sets translation for Mesh B using GPU kernel shifting.\n";
  std::cout << "  translateCPU <x> <y> <z>\n";
  std::cout << "      - Sets translation for Mesh B via double-precision CPU recalculation & GPU upload.\n";
  std::cout << "  export <mesh_tag(A/B)> <out_filename.off>\n";
  std::cout << "      - Exports current Mesh A or B in full 17-digit precision to OFF format.\n";
  std::cout << "  compute [batchMultiplier] [DualTreeSteps] [0 or size of batch for async computation]\n";
  std::cout << "      - Runs the intersection pipeline and returns exact pairs.\n";
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
  // 1. Initial Warmup
  std::cout << "Initializing CUDA environment...\n";
  warmupCUDA();
  std::cout << "Initialization complete.\n";

  // 2. Initialize Polyscope GUI Window
  std::cout << "Initializing Polyscope GUI...\n";
  PolyscopeBridge::init();
  std::cout << "Type 'help' for commands.\n";

  // Global TBB Worker Control
  std::unique_ptr<tbb::global_control> tbbControl;

  // Application State
  KernelBVHController controller;
  ExecutionStats stats;

  Mesh meshA, meshB;
  std::vector<float3> hVertsA, hVertsB;
  std::vector<uint3> hIndicesA, hIndicesB;
  std::vector<float> hErrorsA, hErrorsB;

  bool isLoaded = false;

  // Scaling / Normalization & Centroid State Tracking
  double currentScaleFactor = 1.0;
  double currentCenterX = 0.0, currentCenterY = 0.0, currentCenterZ = 0.0;
  double currentMaxSpan = 0.0;
  bool currentScaledToUnit = false;

  // Centroids in both Original (Unscaled) and Normalized Coordinate Spaces
  double3 origCenterA{0.0, 0.0, 0.0}, origCenterB{0.0, 0.0, 0.0};
  float3 normCenterA{0.0f, 0.0f, 0.0f}, normCenterB{0.0f, 0.0f, 0.0f};

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
          std::cout << "  Original Centroid A  : (" << origCenterA.x << ", " << origCenterA.y << ", " << origCenterA.z << ")\n";
          std::cout << "  Original Centroid B  : (" << origCenterB.x << ", " << origCenterB.y << ", " << origCenterB.z << ")\n";
          std::cout << "  Normalized Centroid A: (" << normCenterA.x << ", " << normCenterA.y << ", " << normCenterA.z << ")\n";
          std::cout << "  Normalized Centroid B: (" << normCenterB.x << ", " << normCenterB.y << ", " << normCenterB.z << ")\n";
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
      } else if(cmd == "load") {
        std::string arg1 = parseArgument(iss);
        std::string arg2 = parseArgument(iss);

        int maxCellA = 1, maxCellB = 1, leafThresh = 4, usePreciseInt = 1, scaleToUnitInt = 0;

        if(arg1.empty()) {
          std::cout << "Usage: load <meshA> [meshB] [maxCellA] [maxCellB] [leafThresh] [usePrecise(0/1)] "
                       "[scaleToUnit(0/1)]\n";
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
            iss >> maxCellB >> leafThresh >> usePreciseInt >> scaleToUnitInt;
          } else {
            pathA = arg1;
            pathB = arg2;
            iss >> maxCellA >> maxCellB >> leafThresh >> usePreciseInt >> scaleToUnitInt;
          }

          bool usePreciseBounds = (usePreciseInt != 0);
          bool scaleToUnit = (scaleToUnitInt != 0);

          controller.cleanup();
          meshA = Mesh();
          meshB = Mesh();
          isLoaded = false;

          std::cout << "Loading mesh(es)...\n";

          auto tIoStart = std::chrono::high_resolution_clock::now();
          bool loadedSuccessfully = false;

          if(singleMeshMode) {
            std::ifstream inFileA(pathA);
            if(inFileA.is_open() && (inFileA >> meshA)) {
              meshB = meshA; // In-memory topology copy
              loadedSuccessfully = true;
            } else {
              std::cerr << "Error: Could not open/read Mesh file: " << pathA << "\n";
            }
          } else {
            std::ifstream inFileA(pathA), inFileB(pathB);
            if(!inFileA.is_open() || !inFileB.is_open() || !(inFileA >> meshA) || !(inFileB >> meshB)) {
              std::cerr << "Error: Failed to load OFF input files from provided paths.\n";
              if(!inFileA.is_open())
                std::cerr << "  -> Could not open Mesh A: " << pathA << "\n";
              if(!inFileB.is_open())
                std::cerr << "  -> Could not open Mesh B: " << pathB << "\n";
            } else {
              loadedSuccessfully = true;
            }
          }

          if(loadedSuccessfully) {
            auto tIoEnd = std::chrono::high_resolution_clock::now();
            double ioMs = std::chrono::duration<double, std::milli>(tIoEnd - tIoStart).count();

            auto tPrepStart = std::chrono::high_resolution_clock::now();

            // 1. Compute ORIGINAL individual bounding boxes & centroids before normalization
            DoubleBox boxA_orig = computeMeshBoundsTBB(meshA);
            DoubleBox boxB_orig = computeMeshBoundsTBB(meshB);
            origCenterA = boxA_orig.getCenter();
            origCenterB = boxB_orig.getCenter();

            // Compute unified scene bounding box
            DoubleBox sceneBox = boxA_orig;
            sceneBox.grow(boxB_orig);

            double cx = 0.5 * (sceneBox.min_x + sceneBox.max_x);
            double cy = 0.5 * (sceneBox.min_y + sceneBox.max_y);
            double cz = 0.5 * (sceneBox.min_z + sceneBox.max_z);

            double spanX = sceneBox.max_x - sceneBox.min_x;
            double spanY = sceneBox.max_y - sceneBox.min_y;
            double spanZ = sceneBox.max_z - sceneBox.min_z;
            double maxSpan = std::max({spanX, spanY, spanZ});

            // Evaluate uniform scale factor
            double scaleFactor = 1.0;
            if(scaleToUnit && maxSpan > 0.0) {
              scaleFactor = 2.0 / maxSpan; // Fits exactly into [-1.0, 1.0]^3
            }

            // Save normalization info to state
            currentScaleFactor = scaleFactor;
            currentCenterX = cx;
            currentCenterY = cy;
            currentCenterZ = cz;
            currentMaxSpan = maxSpan;
            currentScaledToUnit = scaleToUnit;

            // 2. Mutate CGAL Meshes directly so CPU space matches GPU space perfectly
            normalizeCgalMeshTBB(meshA, cx, cy, cz, scaleFactor);
            normalizeCgalMeshTBB(meshB, cx, cy, cz, scaleFactor);

            // 3. Compute NORMALIZED centroids (used as rotation pivot centers)
            DoubleBox boxA_norm = computeMeshBoundsTBB(meshA);
            DoubleBox boxB_norm = computeMeshBoundsTBB(meshB);
            double3 cA_norm = boxA_norm.getCenter();
            double3 cB_norm = boxB_norm.getCenter();

            normCenterA = make_float3(static_cast<float>(cA_norm.x), static_cast<float>(cA_norm.y), static_cast<float>(cA_norm.z));
            normCenterB = make_float3(static_cast<float>(cB_norm.x), static_cast<float>(cB_norm.y), static_cast<float>(cB_norm.z));

            // 4. Extract GPU topologies directly from transformed CGAL meshes
            extractMeshTopology(meshA, hVertsA, hErrorsA, hIndicesA, usePreciseBounds);
            extractMeshTopology(meshB, hVertsB, hErrorsB, hIndicesB, usePreciseBounds);

            auto tPrepEnd = std::chrono::high_resolution_clock::now();
            double prepMs = std::chrono::duration<double, std::milli>(tPrepEnd - tPrepStart).count();

            auto tStart = std::chrono::high_resolution_clock::now();

            controller.construct(meshA, meshB, hVertsA.data(), static_cast<int>(hVertsA.size()), hIndicesA.data(),
                                 (usePreciseBounds ? hErrorsA.data() : nullptr), static_cast<int>(hIndicesA.size()),
                                 maxCellA, hVertsB.data(), static_cast<int>(hVertsB.size()), hIndicesB.data(),
                                 (usePreciseBounds ? hErrorsB.data() : nullptr), static_cast<int>(hIndicesB.size()),
                                 maxCellB, leafThresh, stats);

            auto tEnd = std::chrono::high_resolution_clock::now();
            double ms = std::chrono::duration<double, std::milli>(tEnd - tStart).count();

            // 5. Register identical normalized CGAL meshes with Polyscope
            PolyscopeBridge::reset(meshA, meshB);

            std::cout << "\n--- Load & Construction Performance Breakdown ---\n";
            std::cout << "  1. Disk -> CPU Load:       " << std::fixed << std::setprecision(2) << ioMs << " ms\n";
            std::cout << "  2. Topology/Layout Prep:   " << std::fixed << std::setprecision(2) << prepMs << " ms\n";
            std::cout << "  3. BVH Construction:       " << std::fixed << std::setprecision(2) << ms << " ms\n";
            std::cout << "  -----------------------------------------------\n";
            std::cout << "  Total Preparation Time:    " << std::fixed << std::setprecision(2) << (ioMs + prepMs + ms) << " ms\n";
            std::cout << "  Applied Scale Factor:      " << std::setprecision(10) << scaleFactor
                      << (scaleToUnit ? " (Scaled to [-1, 1]^3)" : " (Unscaled)") << "\n";
            std::cout << "  Normalized Center A:       (" << normCenterA.x << ", " << normCenterA.y << ", " << normCenterA.z << ")\n";
            std::cout << "  Normalized Center B:       (" << normCenterB.x << ", " << normCenterB.y << ", " << normCenterB.z << ")\n\n";

            isLoaded = true;
          }
        }
      } else if(cmd == "transform") {
        if(!isLoaded) {
          std::cout << "Error: You must 'load' meshes before transforming.\n";
        } else {
          float rotAx, rotAy, rotAz;
          float rotBx, rotBy, rotBz;
          float transBx, transBy, transBz;

          if(iss >> rotAx >> rotAy >> rotAz >> rotBx >> rotBy >> rotBz >> transBx >> transBy >> transBz) {
            float3 rotA = make_float3(rotAx, rotAy, rotAz);
            float3 transA = make_float3(0.0f, 0.0f, 0.0f); // Maintained at zero
            float3 rotB = make_float3(rotBx, rotBy, rotBz);
            float3 transB = make_float3(transBx, transBy, transBz);

            auto tStart = std::chrono::high_resolution_clock::now();

            // 1. Dispatch GPU & CPU transformations in CUDA Controller
            controller.setTransformBoth(rotA, transA, rotB, transB);

            // 2. Render transformations visually in PolyscopeBridge
            PolyscopeBridge::transformBoth(rotA, transA, normCenterA,
                                           rotB, transB, normCenterB);

            auto tEnd = std::chrono::high_resolution_clock::now();
            double ms = std::chrono::duration<double, std::milli>(tEnd - tStart).count();

            std::cout << "Transformations applied in " << std::fixed << std::setprecision(2) << ms << " ms:\n";
            std::cout << "  Mesh A Rot: (" << rotAx << ", " << rotAy << ", " << rotAz << ")\n";
            std::cout << "  Mesh B Rot: (" << rotBx << ", " << rotBy << ", " << rotBz << ") | Trans: (" 
                      << transBx << ", " << transBy << ", " << transBz << ")\n";
          } else {
            std::cout << "Usage: transform <rotAx> <rotAy> <rotAz> <rotBx> <rotBy> <rotBz> <transBx> <transBy> <transBz>\n";
          }
        }
      } else if(cmd == "translate") {
        float x, y, z;
        if(iss >> x >> y >> z) {
          auto tStart = std::chrono::high_resolution_clock::now();

          // GPU kernel shift
          controller.setTranslation(x, y, z);
          PolyscopeBridge::translateMeshB(x, y, z);

          auto tEnd = std::chrono::high_resolution_clock::now();
          double elapsedMs = std::chrono::duration<double, std::milli>(tEnd - tStart).count();

          std::cout << "Translation applied (GPU Kernel): (" << x << ", " << y << ", " << z << ") "
                    << "in " << elapsedMs << " ms\n";
        } else {
          std::cout << "Usage: translate <x> <y> <z>\n";
        }
      } else if(cmd == "translateCPU") {
        float x, y, z;
        if(iss >> x >> y >> z) {
          auto tStart = std::chrono::high_resolution_clock::now();

          // CPU double-precision recalculation & host upload
          controller.setTranslationCPUHostUpload(x, y, z);
          PolyscopeBridge::translateMeshB(x, y, z);

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
          int batchMultiplier = 1;
          int mode = 0;
          int activateAsyncDownload = 0;

          iss >> batchMultiplier >> mode >> activateAsyncDownload;

          tbb::concurrent_vector<int2> finalExactPairs;

          auto tStart = std::chrono::high_resolution_clock::now();

          controller.runIntersectionPipeline(batchMultiplier, mode, activateAsyncDownload, finalExactPairs, stats);

          auto tEnd = std::chrono::high_resolution_clock::now();
          double ms = std::chrono::duration<double, std::milli>(tEnd - tStart).count();

          std::vector<int2> stdPairs(finalExactPairs.begin(), finalExactPairs.end());
          PolyscopeBridge::highlightIntersections(stdPairs, num_faces(meshA), num_faces(meshB));

          std::cout << "Query completed in time: " << std::fixed << std::setprecision(2) << ms << " ms. "
                    << "We got " << finalExactPairs.size() << " intersections.\n";
          std::cout << "AABB Hits: " << stats.finalAabbCandidatePairs
                    << " | Greenlist : " << stats.loopTracker.confirmedGreenPairs
                    << " | Yellowlist: " << stats.loopTracker.confirmedYellowPairs << std::endl;
        }
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