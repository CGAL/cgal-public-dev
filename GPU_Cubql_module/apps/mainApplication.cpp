#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>
#include <iomanip>

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

bool exportCgalMeshToOff(const Mesh& mesh, const std::string& filepath) {
  std::ofstream outFile(filepath);
  if (!outFile.is_open()) {
    std::cerr << "Error: Could not open file for writing: " << filepath << "\n";
    return false;
  }

  // Set max double precision (17 digits) so floating point values export losslessly
  outFile << std::setprecision(17);

  // Direct CGAL Surface_mesh OFF stream operator export
  if (!(outFile << mesh)) {
    std::cerr << "Error: Failed writing CGAL mesh data to: " << filepath << "\n";
    return false;
  }

  return true;
}

void extractMeshTopologyNormalized(const Mesh& mesh,
                                   double cx,
                                   double cy,
                                   double cz,
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

      double sx = p.x() - cx, sy = p.y() - cy, sz = p.z() - cz;
      float fx = static_cast<float>(sx);
      float fy = static_cast<float>(sy);
      float fz = static_cast<float>(sz);

      outVerts[i] = {fx, fy, fz};

      if(usePreciseBounds) {
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
  std::cout << "  load <meshA.off> <meshB.off> [maxCellA] [maxCellB] [leafThresh] [usePrecise(0/1)]\n";
  std::cout << "      - Loads meshes, normalizes them, and constructs the BVHs.\n";
  std::cout << "  translate <x> <y> <z>\n";
  std::cout << "      - Sets the translation vector for Mesh B.\n";
  std::cout << "  export <mesh_tag(A/B)> <out_filename.off>\n";
  std::cout << "      - Exports Mesh A or B in full 17-digit precision to OFF format.\n";
  std::cout << "  compute [batchMultiplier] [mode]\n";
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
  iss >> std::ws; // Consume leading whitespace
  if(iss.peek() == '"' || iss.peek() == '\'') {
    char quote = iss.get();        // Read the opening quote
    std::getline(iss, arg, quote); // Read until the matching closing quote
  } else {
    iss >> arg; // Normal word parsing if no quotes
  }
  return arg;
}

int main(int argc, char** argv) {
  // 1. Initial Warmup
  std::cout << "Initializing CUDA environment...\n";
  warmupCUDA();
  std::cout << "Initialization complete. Type 'help' for commands.\n";

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

  // Readline CLI Loop
  char* rawInput = nullptr;
  while((rawInput = readline("> ")) != nullptr) {
    std::string line(rawInput);

    // Add non-empty lines to history so UP/DOWN arrow navigation works
    if(!line.empty()) {
      add_history(rawInput);
    }
    free(rawInput); // Free memory allocated by readline

    std::istringstream iss(line);
    std::string cmd;
    if(!(iss >> cmd))
      continue;

    if(cmd == "quit" || cmd == "exit") {
      break;
    } else if(cmd == "help") {
      printHelp();
    } else if(cmd == "tbb") {
      int numThreads = 0;
      if(iss >> numThreads && numThreads > 0) {
        tbbControl = std::make_unique<tbb::global_control>(tbb::global_control::max_allowed_parallelism, numThreads);
        std::cout << "TBB maximum worker limit updated to: " << numThreads << "\n";
      } else {
        std::cout << "Usage: tbb <num_threads> (must be > 0)\n";
      }
    } else if(cmd == "load") {
      // Use our helper function to safely strip single or double quotes
      std::string pathA = parseArgument(iss);
      std::string pathB = parseArgument(iss);

      int maxCellA = 1, maxCellB = 1, leafThresh = 4, usePreciseInt = 1;

      if(pathA.empty() || pathB.empty()) {
        std::cout << "Usage: load <meshA> <meshB> [maxCellA] [maxCellB] [leafThresh] [usePrecise(0/1)]\n";
        continue;
      }

      // Read optional remaining integer arguments
      iss >> maxCellA >> maxCellB >> leafThresh >> usePreciseInt;
      bool usePreciseBounds = (usePreciseInt != 0);

      // Clean up previous allocations & state
      controller.cleanup();
      meshA = Mesh();
      meshB = Mesh();
      isLoaded = false;

      std::cout << "Loading meshes...\n";

      // ---------------------------------------------------------------------
      // PHASE 1: Disk -> CPU Load Time
      // ---------------------------------------------------------------------
      auto tIoStart = std::chrono::high_resolution_clock::now();

      std::ifstream inFileA(pathA), inFileB(pathB);
      if(!inFileA.is_open() || !inFileB.is_open() || !(inFileA >> meshA) || !(inFileB >> meshB)) {
        std::cerr << "Error: Failed to load OFF input files from provided paths.\n";
        if(!inFileA.is_open())
          std::cerr << "  -> Could not open Mesh A: " << pathA << "\n";
        if(!inFileB.is_open())
          std::cerr << "  -> Could not open Mesh B: " << pathB << "\n";
        continue;
      }

      auto tIoEnd = std::chrono::high_resolution_clock::now();
      double ioMs = std::chrono::duration<double, std::milli>(tIoEnd - tIoStart).count();

      // ---------------------------------------------------------------------
      // PHASE 2: Layout Conversion & Normalization Time
      // ---------------------------------------------------------------------
      auto tPrepStart = std::chrono::high_resolution_clock::now();

      DoubleBox box = computeMeshBoundsTBB(meshA);
      box.grow(computeMeshBoundsTBB(meshB));
      double cx = 0.5 * (box.min_x + box.max_x);
      double cy = 0.5 * (box.min_y + box.max_y);
      double cz = 0.5 * (box.min_z + box.max_z);

      extractMeshTopologyNormalized(meshA, cx, cy, cz, hVertsA, hErrorsA, hIndicesA, usePreciseBounds);
      extractMeshTopologyNormalized(meshB, cx, cy, cz, hVertsB, hErrorsB, hIndicesB, usePreciseBounds);

      auto tPrepEnd = std::chrono::high_resolution_clock::now();
      double prepMs = std::chrono::duration<double, std::milli>(tPrepEnd - tPrepStart).count();

      // ---------------------------------------------------------------------
      // PHASE 3: BVH Construction Time (Original Block Unchanged)
      // ---------------------------------------------------------------------
      auto tStart = std::chrono::high_resolution_clock::now();

      controller.construct(meshA, meshB, hVertsA.data(), static_cast<int>(hVertsA.size()), hIndicesA.data(),
                           (usePreciseBounds ? hErrorsA.data() : nullptr), static_cast<int>(hIndicesA.size()), maxCellA,
                           hVertsB.data(), static_cast<int>(hVertsB.size()), hIndicesB.data(),
                           (usePreciseBounds ? hErrorsB.data() : nullptr), static_cast<int>(hIndicesB.size()), maxCellB,
                           leafThresh, stats);

      auto tEnd = std::chrono::high_resolution_clock::now();
      double ms = std::chrono::duration<double, std::milli>(tEnd - tStart).count();

      // ---------------------------------------------------------------------
      // Print Timings Breakdown
      // ---------------------------------------------------------------------
      std::cout << "\n--- Load & Construction Performance Breakdown ---\n";
      std::cout << "  1. Disk -> CPU Load:       " << std::fixed << std::setprecision(2) << ioMs << " ms\n";
      std::cout << "  2. Topology/Layout Prep:   " << std::fixed << std::setprecision(2) << prepMs << " ms\n";
      std::cout << "  3. BVH Construction:       " << std::fixed << std::setprecision(2) << ms << " ms\n";
      std::cout << "  -----------------------------------------------\n";
      std::cout << "  Total Preparation Time:    " << std::fixed << std::setprecision(2) << (ioMs + prepMs + ms)
                << " ms\n\n";

      isLoaded = true;
    
  }
  else if(cmd == "translate") {
    float x, y, z;
    if(iss >> x >> y >> z) {
      auto tStart = std::chrono::high_resolution_clock::now();
      
      controller.setTranslation(x, y, z);
      
      auto tEnd = std::chrono::high_resolution_clock::now();
      double elapsedMs = std::chrono::duration<double, std::milli>(tEnd - tStart).count();

      std::cout << "Translation applied: (" << x << ", " << y << ", " << z << ") "
                << "in " << elapsedMs << " ms\n";
    } else {
      std::cout << "Usage: translate <x> <y> <z>\n";
    }
  }
  else if(cmd == "compute") {
    if(!isLoaded) {
      std::cout << "Error: You must 'load' meshes before computing.\n";
      continue;
    }

    int batchMultiplier = 1;
    int mode = 0; // Default execution mode

    // Read optional configuration
    iss >> batchMultiplier >> mode;

    tbb::concurrent_vector<int2> finalExactPairs;

    // Execute Pipeline and measure time
    auto tStart = std::chrono::high_resolution_clock::now();

    controller.runIntersectionPipeline(batchMultiplier, mode, finalExactPairs, stats);

    auto tEnd = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(tEnd - tStart).count();

    std::cout << "Query completed in time: " << std::fixed << std::setprecision(2) << ms << " ms. "
              << "We got " << finalExactPairs.size() << " intersections.\n";
  }
  else if(cmd == "export") {
    if(!isLoaded) {
      std::cout << "Error: You must 'load' meshes before exporting.\n";
      continue;
    }

    std::string meshTag = parseArgument(iss);
    std::string outPath = parseArgument(iss);

    if(meshTag.empty() || outPath.empty()) {
      std::cout << "Usage: export <A/B> <output.off>\n";
      continue;
    }

    const Mesh* targetMesh = nullptr;
    if(meshTag == "A" || meshTag == "a") {
      targetMesh = &meshA;
    } else if(meshTag == "B" || meshTag == "b") {
      targetMesh = &meshB;
    } else {
      std::cout << "Error: Invalid mesh target '" << meshTag << "'. Choose 'A' or 'B'.\n";
      continue;
    }

    auto tStart = std::chrono::high_resolution_clock::now();

    if(exportCgalMeshToOff(*targetMesh, outPath)) {
      auto tEnd = std::chrono::high_resolution_clock::now();
      double ms = std::chrono::duration<double, std::milli>(tEnd - tStart).count();

      std::cout << "Successfully exported Mesh " << meshTag << " (" 
                << num_vertices(*targetMesh) << " verts, " 
                << num_faces(*targetMesh) << " faces) with 17-digit precision to '" 
                << outPath << "' in " << std::fixed << std::setprecision(2) << ms << " ms.\n";
    }
  }

  else {
    std::cout << "Unknown command: '" << cmd << "'. Type 'help' for a list of commands.\n";
  }
}

std::cout << "Cleaning up and exiting...\n";
controller.cleanup();
return 0;
}