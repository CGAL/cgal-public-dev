#include <iostream>
#include <vector>
#include <fstream>
#include <string>

// --- cuBQL INCLUDES ---
#include "include/loadOFFCGAL.h" 
#include "include/loadOFF.h"

// --- CGAL INCLUDES ---
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/intersections.h>

// Use the robust Exact Kernel configuration requested
//typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
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
    double executionTime = 0.0; // BVH Build + Query Passes
    double downloadTime = 0.0;
};

// Linked external symbol explicitly matching your updated backend signature
extern "C" long long runMeshIntersection(
    const cuBQL::Triangle* hA, int NA, 
    const cuBQL::Triangle* hB, int NB, 
    std::vector<IntersectionPair>& hPairs,
    GPUTimingBreakdown& outTimings,
    bool doMoreWork);

// Helper helper to convert a topological CGAL face to raw cuBQL layout safely
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


int main(int ac, char** av) {
    if (ac < 6) {
        std::cout << "Usage: " << av[0] << " <A.off> <B.off> <gpu_fine_output.csv> <final_exact_output.csv> <doMoreWork (1/0)>\n";
        return 1;
    }

    std::string gpuOutputPath   = av[3]; 
    std::string exactOutputPath = av[4]; 
    bool doMoreWork             = (std::stoi(av[5]) != 0); // Convert 5th parameter into execution strategy
    
    double t_start = cuBQL::getCurrentTime();

    // 1. UNIFIED MESH LOADING (Using CGAL Engine)
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

    // 2. DATA TRANSLATION (CGAL mesh structures translated into cuBQL Linear Buffers)
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

    // 3. BROAD-PHASE GPU INTERSECTION EXECUTION
    std::vector<IntersectionPair> hPairs;
    GPUTimingBreakdown gpuTimings; 
    
    long long totalAABBHits = runMeshIntersection(hA.data(), NA, hB.data(), NB, hPairs, gpuTimings, doMoreWork);

    // 4. WRITE INTERMEDIATE GPU BVH OUTPUT
    t0 = cuBQL::getCurrentTime();
    if (!hPairs.empty()) {
        std::ofstream gpuFile(gpuOutputPath);
        gpuFile << "ID_A,ID_B\n";
        for (const auto& p : hPairs) {
            gpuFile << p.idA << "," << p.idB << "\n";
        }
    }
    std::cout << "[Step 5] Write GPU Fine CSV to Disk: " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - t0) << "s\n";

    // 5. COMBINED NARROW-PHASE FILTERING LOOP (On-The-Fly Construction)
    double t_inter_start = cuBQL::getCurrentTime();
    
    auto& coords1 = m1.points();
    auto& coords2 = m2.points();

    std::vector<IntersectionPair> exactPairs;
    exactPairs.reserve(hPairs.size());

    for (const auto& pair : hPairs) {
        face_descriptor fA = facesA[pair.idA];
        face_descriptor fB = facesB[pair.idB];

        // Zero-allocation topological pointer evaluation for Triangle A
        auto halfedgeA = m1.halfedge(fA);
        auto vA0 = m1.target(halfedgeA);
        auto vA1 = m1.target(m1.next(halfedgeA));
        auto vA2 = m1.target(m1.next(m1.next(halfedgeA)));
        Triangle3 triA(coords1[vA0], coords1[vA1], coords1[vA2]);

        // Zero-allocation topological pointer evaluation for Triangle B
        auto halfedgeB = m2.halfedge(fB);
        auto vB0 = m2.target(halfedgeB);
        auto vB1 = m2.target(m2.next(halfedgeB));
        auto vB2 = m2.target(m2.next(m2.next(halfedgeB)));
        Triangle3 triB(coords2[vB0], coords2[vB1], coords2[vB2]);

        // Evaluate intersection predicate immediately on hot stack data
        if (CGAL::do_intersect(triA, triB)) {
            exactPairs.push_back(pair);
        }
    }
    double time_intersection = cuBQL::getCurrentTime() - t_inter_start;
    
    // OUTPUT DETAILED BREAKDOWNS
    std::cout << "----------------------------------------------------\n";
    std::cout << "[Step 6] GPU Pipeline Profiler Breakdown:\n";
    std::cout << "  |---> Control Parameter (doMoreWork): " << (doMoreWork ? "TRUE" : "FALSE") << "\n";
    std::cout << "  |---> Host->Device Upload Time:       " << cuBQL::prettyDouble(gpuTimings.uploadTime) << "s\n";
    std::cout << "  |---> BVH Build & Execution Time:     " << cuBQL::prettyDouble(gpuTimings.executionTime) << "s\n";
    std::cout << "  |---> Device->Host Download Time:     " << cuBQL::prettyDouble(gpuTimings.downloadTime) << "s\n";
    
    std::cout << "[Step 7] Streamlined CGAL Narrow-Phase Processing:\n";
    std::cout << "  |---> Total Filtering Routine Wallclock: " << cuBQL::prettyDouble(time_intersection) << "s\n";
    std::cout << "  |---> Data Conversion Loop Overhead:    0.000000s (Fused into Filtering)\n";
    std::cout << "  |---> Exact do_intersect Verification:  " << cuBQL::prettyDouble(time_intersection) << "s\n";
    std::cout << "----------------------------------------------------\n";

    // 6. FINAL EXACT OUTPUT GENERATION
    t0 = cuBQL::getCurrentTime();
    if (!exactPairs.empty()) {
        std::ofstream outFile(exactOutputPath);
        outFile << "ID_A,ID_B\n";
        for (const auto& p : exactPairs) {
            outFile << p.idA << "," << p.idB << "\n";
        }
    }
    std::cout << "[Step 8] Write Final CSV to Disk:     " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - t0) << "s\n";

    std::cout << "====================================================\n";
    std::cout << "Macro Statistics Summary Summary:\n";
    std::cout << "  |-- Total Process Runtime: " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - t_start) << " s\n";
    std::cout << "  |-- GPU Broad-Phase Hits:  " << totalAABBHits << " (Saved: " << gpuOutputPath << ")\n";
    std::cout << "  |-- CGAL Narrow-Phase Hits: " << exactPairs.size() << " (Pruned " << (totalAABBHits - exactPairs.size()) << " duplicates)\n";
    std::cout << "====================================================\n";

    return 0;
}