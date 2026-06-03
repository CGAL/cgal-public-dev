#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm> // For std::transform

// --- cuBQL INCLUDES ---
#include "include/loadOFFCGAL.h" 
#include "include/loadOFF.h"

// --- CGAL INCLUDES ---
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/intersections.h>

// --- TBB INCLUDES ----
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/concurrent_vector.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;

typedef Kernel::Point_3                             Point3;
typedef Kernel::Triangle_3                          Triangle3;

struct IntersectionPair {
    int idA;
    int idB;
};

struct GPUTimingBreakdown {
    double uploadTime = 0.0;
    double executionTime = 0.0; 
    double downloadTime = 0.0;
};

// Updated signature pointing to our newly compiled partition routine inside kernelsV2
extern "C" void runPartitionedMeshIntersection(
    const cuBQL::Triangle* hA, int NA, 
    const cuBQL::Triangle* hB, int NB, 
    std::vector<IntersectionPair>& hGreenPairs,
    std::vector<IntersectionPair>& hYellowPairs,
    GPUTimingBreakdown& outTimings,
    int pipelineMode);

cuBQL::Triangle convertCgalFaceToCuBQL(const Mesh& mesh, face_descriptor face) {
    auto conn = mesh.vertices_around_face(mesh.halfedge(face));
    auto it = conn.begin();
    
    Point3 pA = mesh.point(*it++);
    Point3 pB = mesh.point(*it++);
    Point3 pC = mesh.point(*it);

    cuBQL::Triangle tri;
    tri.a = { (float)pA.x(), (float)pA.y(), (float)pA.z() };
    tri.b = { (float)pB.x(), (float)pB.y(), (float)pB.z() };
    tri.c = { (float)pC.x(), (float)pC.y(), (float)pC.z() };
    return tri;
}

bool should_save_file(const std::string& path) {
    if (path.empty()) return false;
    std::string lower = path;
    std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c){ return std::tolower(c); });
    return (lower != "none" && lower != "skip" && lower != "\"\"");
}

int main(int ac, char** av) {
    if (ac < 7) { 
        std::cout << "Usage: " << av[0] << " <A.off> <B.off> <gpu_fine_output.csv|none> <final_exact_output.csv|none> <pipelineMode (0=Default/1/2)> <numThreads>\n";
        std::cout << "Hint: Pass 'none' or 'skip' to bypass writing intermediate or final CSV records to disk.\n";
        return 1;
    }

    std::string gpuOutputPath   = av[3]; 
    std::string exactOutputPath = av[4]; 
    int pipelineMode            = std::stoi(av[5]); 
    int numThreads              = std::stoi(av[6]); 

    bool writeGpuCSV   = should_save_file(gpuOutputPath);
    bool writeExactCSV = should_save_file(exactOutputPath);

    tbb::global_control global_limit(
        tbb::global_control::max_allowed_parallelism, numThreads
    );

    double t_start = cuBQL::getCurrentTime();

    // 1. UNIFIED MESH LOADING
    double t0 = cuBQL::getCurrentTime();
    Mesh m1, m2;
    std::ifstream in1(av[1]);
    std::ifstream in2(av[2]);
    if (!(in1 >> m1) || !(in2 >> m2)) {
        std::cerr << "Critical Error: Failed to open or parse input mesh files.\n";
        return 1;
    }
    in1.close();
    in2.close();
    std::cout << "[Step 1] Load Meshes from Disk via CGAL: " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - t0) << "s\n";

    // 2. DATA TRANSLATION
    double t_conv_start = cuBQL::getCurrentTime();
    
    std::vector<cuBQL::Triangle> hA; 
    hA.reserve(num_faces(m1));
    std::vector<face_descriptor> facesA; 
    facesA.reserve(num_faces(m1));
    for (auto f : m1.faces()) {
        hA.push_back(convertCgalFaceToCuBQL(m1, f));
        facesA.push_back(f);
    }

    std::vector<cuBQL::Triangle> hB; 
    hB.reserve(num_faces(m2));
    std::vector<face_descriptor> facesB; 
    facesB.reserve(num_faces(m2));
    for (auto f : m2.faces()) {
        hB.push_back(convertCgalFaceToCuBQL(m2, f));
        facesB.push_back(f);
    }
    
    int NA = hA.size();
    int NB = hB.size();
    double time_conversion = cuBQL::getCurrentTime() - t_conv_start;
    std::cout << "[Step 2] Convert Mesh Elements to cuBQL Formats: " << cuBQL::prettyDouble(time_conversion) << "s\n";

    // 3. BROAD-PHASE GPU INTERSECTION EXECUTION (Dual Lists Allocation Track)
    std::vector<IntersectionPair> confirmedIntersections; // Green List
    std::vector<IntersectionPair> ambiguousPairs;         // Yellow List
    GPUTimingBreakdown gpuTimings; 
    
    // Call our newly partitioned pipeline logic mapping dynamic configuration profile tokens
    runPartitionedMeshIntersection(
        hA.data(), NA, 
        hB.data(), NB, 
        confirmedIntersections, 
        ambiguousPairs, 
        gpuTimings, 
        pipelineMode
    );

    // 4. WRITE INTERMEDIATE GPU BVH OUTPUT (Saves everything the GPU found interesting: Green + Yellow)
    t0 = cuBQL::getCurrentTime();
    if (writeGpuCSV) {
        std::ofstream gpuFile(gpuOutputPath);
        gpuFile << "ID_A,ID_B,Type\n";
        for (const auto& p : confirmedIntersections) {
            gpuFile << p.idA << "," << p.idB << ",Green\n";
        }
        for (const auto& p : ambiguousPairs) {
            gpuFile << p.idA << "," << p.idB << ",Yellow\n";
        }
        std::cout << "[Step 5] Write GPU Fine CSV to Disk: " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - t0) << "s\n";
    } else {
        std::cout << "[Step 5] Write GPU Fine CSV to Disk: SKIPPED\n";
    }

    // 5. HYBRID NARROW-PHASE FILTERING LOOP (Processes Only Yellow List)
    double t_inter_start = cuBQL::getCurrentTime();
    
    auto& coords1 = m1.points();
    auto& coords2 = m2.points();

    // Initialize our results thread vector with the known safe elements from the green list
    tbb::concurrent_vector<IntersectionPair> finalExactPairs(confirmedIntersections.begin(), confirmedIntersections.end());

    // TBB runs parallel execution over ONLY the micro-ambiguous yellow elements!
    if (!ambiguousPairs.empty()) {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, ambiguousPairs.size()), [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                const auto& pair = ambiguousPairs[i];
                
                face_descriptor fA = facesA[pair.idA];
                face_descriptor fB = facesB[pair.idB];

                auto halfedgeA = m1.halfedge(fA);
                auto vA0 = m1.target(halfedgeA);
                auto vA1 = m1.target(m1.next(halfedgeA));
                auto vA2 = m1.target(m1.next(m1.next(halfedgeA)));
                Triangle3 triA(coords1[vA0], coords1[vA1], coords1[vA2]);

                auto halfedgeB = m2.halfedge(fB);
                auto vB0 = m2.target(halfedgeB);
                auto vB1 = m2.target(m2.next(halfedgeB));
                auto vB2 = m2.target(m2.next(m2.next(halfedgeB)));
                Triangle3 triB(coords2[vB0], coords2[vB1], coords2[vB2]);

                // The true mathematical filter
                if (CGAL::do_intersect(triA, triB)) {
                    finalExactPairs.push_back(pair);
                }
            }
        });
    }

    double time_intersection = cuBQL::getCurrentTime() - t_inter_start;
    
    // OUTPUT DETAILED BREAKDOWNS
    std::cout << "----------------------------------------------------\n";
    std::cout << "[Step 6] GPU Pipeline Profiler Breakdown:\n";
    std::cout << "  |---> Control Variant Mode Parameter: " << pipelineMode << "\n";
    std::cout << "  |---> Active TBB Worker Limit:        " << numThreads << "\n";
    std::cout << "  |---> Host->Device Upload Time:       " << cuBQL::prettyDouble(gpuTimings.uploadTime) << "s\n";
    std::cout << "  |---> BVH Build & Execution Time:     " << cuBQL::prettyDouble(gpuTimings.executionTime) << "s\n";
    std::cout << "  |---> Device->Host Download Time:     " << cuBQL::prettyDouble(gpuTimings.downloadTime) << "s\n";
    std::cout << "  |---> Definite Intersections (Green): " << confirmedIntersections.size() << "\n";
    std::cout << "  |---> Uncertain Candidates (Yellow):  " << ambiguousPairs.size() << "\n";
    
    std::cout << "[Step 7] Streamlined Parallel CGAL Narrow-Phase Processing:\n";
    std::cout << "  |---> Total CPU Intersect Routine:     " << cuBQL::prettyDouble(time_intersection) << "s\n";
    std::cout << "  |---> (CPU processed " << ambiguousPairs.size() << " elements instead of " << (confirmedIntersections.size() + ambiguousPairs.size()) << ")\n";
    std::cout << "----------------------------------------------------\n";

    // 6. FINAL EXACT OUTPUT GENERATION (Conditional)
    t0 = cuBQL::getCurrentTime();
    if (writeExactCSV && !finalExactPairs.empty()) {
        std::ofstream outFile(exactOutputPath);
        outFile << "ID_A,ID_B\n";
        for (const auto& p : finalExactPairs) {
            outFile << p.idA << "," << p.idB << "\n";
        }
        std::cout << "[Step 8] Write Final CSV to Disk:     " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - t0) << "s\n";
    } else {
        std::cout << "[Step 8] Write Final CSV to Disk:     SKIPPED\n";
    }

    std::cout << "====================================================\n";
    std::cout << "Macro Statistics Summary:\n";
    std::cout << "  |-- Total Process Runtime:      " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - t_start) << " s\n";
    std::cout << "  |-- GPU Definite/Uncertain Track: " << (confirmedIntersections.size() + ambiguousPairs.size()) << "\n";
    std::cout << "  |-- Total Confirmed Intersects:  " << finalExactPairs.size() << (writeExactCSV ? " (Saved)" : " (Not Saved)") << "\n";
    std::cout << "====================================================\n";

    return 0;
}