#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm> // For std::transform

// --- cuBQL INCLUDES ---
#include "cuBQL/bvh.h"
#include "include/loadOFFCGAL.h"
#include "include/loadOFF.h"

// --- CGAL INCLUDES ---
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef Kernel::Point_3 Point3;

// Structured metadata mapping data per voxel to host-side profiling
struct VoxelContentMetrics {
  int voxelIndex;
  int countA;
  int countB;
  int firstPairOffset; // Start position in the flat sequential candidates array
};

struct VoxelTimingBreakdown {
  double uploadTime = 0.0;
  double bvhBuildTimeA = 0.0;
  double bvhBuildTimeB = 0.0;
  double countingPassTime = 0.0;
  double prefixScanTime = 0.0;
  double fillingPassTime = 0.0;
  double downloadTime = 0.0;
  double totalExecutionTime = 0.0;
};

// Extern declaration referencing our newly created dual-pass voxelizer routine
extern "C" void runDualPassVoxelizer(const cuBQL::Triangle* hA,
                                     int NA,
                                     const cuBQL::Triangle* hB,
                                     int NB,
                                     cuBQL::vec3i dims,
                                     std::vector<VoxelContentMetrics>& hVoxelMetrics,
                                     std::vector<int>& hFlatTrianglesA,
                                     std::vector<int>& hFlatTrianglesB,
                                     VoxelTimingBreakdown& outTimings);

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

bool should_save_file(const std::string& path) {
  if(path.empty())
    return false;
  std::string lower = path;
  std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c) { return std::tolower(c); });
  return (lower != "none" && lower != "skip" && lower != "\"\"");
}

int main(int ac, char** av) {
  if(ac < 5) {
    std::cout << "Usage: " << av[0]
              << " <A.off> <B.off> <grid_resolution_n> <voxel_metrics_output.csv|none>\n";
    std::cout << "Hint: Pass 'none' or 'skip' to bypass writing spatial voxel records to disk.\n";
    return 1;
  }

  int n = std::stoi(av[3]);
  std::string csvOutputPath = av[4];
  bool writeMetricsCSV = should_save_file(csvOutputPath);

  double t_start = cuBQL::getCurrentTime();

  // 1. UNIFIED MESH LOADING
  double t0 = cuBQL::getCurrentTime();
  Mesh m1, m2;
  std::ifstream in1(av[1]);
  std::ifstream in2(av[2]);
  if(!(in1 >> m1) || !(in2 >> m2)) {
    std::cerr << "Critical Error: Failed to open or parse input mesh files.\n";
    return 1;
  }
  in1.close();
  in2.close();
  std::cout << "[Step 1] Load Meshes from Disk via CGAL: " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - t0) << "\n";

  // 2. DATA TRANSLATION & UNIFIED WORLD BOUNDS COMPUTATION
  double t_conv_start = cuBQL::getCurrentTime();
  cuBQL::box3f worldBounds;

  std::vector<cuBQL::Triangle> hA;
  hA.reserve(num_faces(m1));
  for(auto f : m1.faces()) {
    cuBQL::Triangle t = convertCgalFaceToCuBQL(m1, f);
    hA.push_back(t);
    worldBounds.extend(t.bounds());
  }

  std::vector<cuBQL::Triangle> hB;
  hB.reserve(num_faces(m2));
  for(auto f : m2.faces()) {
    cuBQL::Triangle t = convertCgalFaceToCuBQL(m2, f);
    hB.push_back(t);
    worldBounds.extend(t.bounds());
  }

  int NA = hA.size();
  int NB = hB.size();
  double time_conversion = cuBQL::getCurrentTime() - t_conv_start;
  std::cout << "[Step 2] Convert Mesh Elements to cuBQL Formats: " << cuBQL::prettyDouble(time_conversion) << "\n";

  // STRICT UNIFORM PARAMETER MAPPING: Force N x N x N Matrix Space
  cuBQL::vec3i dims(n, n, n);
  int numCells = dims.x * dims.y * dims.z;

  std::cout << "--> Computed Voxel Matrix Grid Dims: [" << dims.x << " x " << dims.y << " x " << dims.z << "]\n";
  std::cout << "--> Global Simulation Bounding Box Limits: Min(" 
            << worldBounds.lower.x << "," << worldBounds.lower.y << "," << worldBounds.lower.z << ") Max("
            << worldBounds.upper.x << "," << worldBounds.upper.y << "," << worldBounds.upper.z << ")\n";

  // 3. BROAD-PHASE GPU VOXELIZATION EXECUTION (Dual Pass Track)
  std::vector<VoxelContentMetrics> voxelMetrics;
  std::vector<int> flatTrianglesA;
  std::vector<int> flatTrianglesB;
  VoxelTimingBreakdown voxelTimings;

  runDualPassVoxelizer(hA.data(), NA, hB.data(), NB, dims, voxelMetrics, flatTrianglesA, flatTrianglesB, voxelTimings);

  // 4. WRITE SPATIAL DECOMPOSITION VOXEL OUTPUT WITH EXPLICIT TRACKING TIMER
  if(writeMetricsCSV) {
    double t_csv_start = cuBQL::getCurrentTime();
    std::ofstream csvFile(csvOutputPath);
    csvFile << "VoxelIdx,GridX,GridY,GridZ,XMin,YMin,ZMin,XMax,YMax,ZMax,CountA,CountB,FlatListOffset\n";

    for(int i = 0; i < numCells; ++i) {
      const auto& cell = voxelMetrics[i];

      // Reverse map flat tracking index to structured 3D grid matrix locations
      int ix = i % dims.x;
      int iy = (i / dims.x) % dims.y;
      int iz = i / (dims.x * dims.y);

      // Interpolate local box bounds 
      cuBQL::vec3f f0 = cuBQL::vec3f(ix, iy, iz) / cuBQL::vec3f(dims);
      cuBQL::vec3f f1 = cuBQL::vec3f(ix + 1, iy + 1, iz + 1) / cuBQL::vec3f(dims);
      cuBQL::box3f cellBounds{worldBounds.lerp(f0), worldBounds.lerp(f1)};

      csvFile << cell.voxelIndex << "," 
              << ix << "," << iy << "," << iz << ","
              << cellBounds.lower.x << "," << cellBounds.lower.y << "," << cellBounds.lower.z << ","
              << cellBounds.upper.x << "," << cellBounds.upper.y << "," << cellBounds.upper.z << ","
              << cell.countA << "," << cell.countB << "," 
              << cell.firstPairOffset << "\n";
    }
    csvFile.close();
    double time_csv = cuBQL::getCurrentTime() - t_csv_start;
    std::cout << "[Step 5] Write Spatial Voxel Metrics CSV to Disk: " << cuBQL::prettyDouble(time_csv) << "\n";
  } else {
    std::cout << "[Step 5] Write Spatial Voxel Metrics CSV to Disk: SKIPPED\n";
  }

  // OUTPUT DETAILED BREAKDOWNS (CLEANED FROM TRAILING LITERALS)
  std::cout << "----------------------------------------------------\n";
  std::cout << "[Step 6s] GPU Voxelizer Pipeline Profiler Breakdown:\n";
  std::cout << "  |---> Requested Grid Base Parameter (n): " << n << "\n";
  std::cout << "  |---> Host->Device Upload Time:          " << cuBQL::prettyDouble(voxelTimings.uploadTime) << "\n";
  std::cout << "  |---> Total Combined Execution Time:     " << cuBQL::prettyDouble(voxelTimings.totalExecutionTime) << "\n";
  std::cout << "  |        |---> Mesh A BVH Build Time:    " << cuBQL::prettyDouble(voxelTimings.bvhBuildTimeA) << "\n"; 
  std::cout << "  |        |---> Mesh B BVH Build Time:    " << cuBQL::prettyDouble(voxelTimings.bvhBuildTimeB) << "\n"; 
  std::cout << "  |        |---> Query & Evaluation Loop:  " << cuBQL::prettyDouble(voxelTimings.countingPassTime + voxelTimings.prefixScanTime + voxelTimings.fillingPassTime) << "\n";    
  std::cout << "  |        |      |---> Pass 1 (Counting): " << cuBQL::prettyDouble(voxelTimings.countingPassTime) << "\n";
  std::cout << "  |        |      |---> Thrust Scan Stage: " << cuBQL::prettyDouble(voxelTimings.prefixScanTime) << "\n";
  std::cout << "  |        |      |---> Pass 2 (Filling):  " << cuBQL::prettyDouble(voxelTimings.fillingPassTime) << "\n";
  std::cout << "  |---> Device->Host Download Time:        " << cuBQL::prettyDouble(voxelTimings.downloadTime) << "\n";
  std::cout << "----------------------------------------------------\n";

  std::cout << "====================================================\n";
  std::cout << "Macro Statistics Summary:\n";
  std::cout << "  |-- Total Process Runtime:              " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - t_start) << "\n";
  std::cout << "  |-- Evaluated Grid Elements (Voxels):   " << voxelMetrics.size() << "\n";
  std::cout << "  |-- Total Mesh A Contained References:  " << flatTrianglesA.size() << "\n";
  std::cout << "  |-- Total Mesh B Contained References:  " << flatTrianglesB.size() << "\n";
  std::cout << "====================================================\n";

  return 0;
}