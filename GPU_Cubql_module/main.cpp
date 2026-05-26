#include <iostream>
#include <vector>
#include <fstream>
#include "include/loadOFFCGAL.h" // Clean header-only include!
#include "include/loadOFF.h"

// --- CGAL INCLUDES ---
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/intersections.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3                             Point3;
typedef Kernel::Triangle_3                          Triangle3;

struct IntersectionPair {
    int idA;
    int idB;
};

// Simple struct to hold our pre-converted CGAL triangle pairs
struct CGALTrianglePair {
    Triangle3 triA;
    Triangle3 triB;
    IntersectionPair originalPair; // Keep track of original IDs for the output
};

// Mirror structural declaration from kernels.cu to prevent ABI mismatch
struct GPUTimingBreakdown {
    double uploadTime = 0.0;
    double executionTime = 0.0; // BVH Build + Query Passes
    double downloadTime = 0.0;
};

// Declare the external function compiled by nvcc (Updated with 6th argument)
extern "C" long long runMeshIntersection(
    const cuBQL::Triangle* hA, int NA, 
    const cuBQL::Triangle* hB, int NB, 
    std::vector<IntersectionPair>& hPairs,
    GPUTimingBreakdown& outTimings); // <--- Synchronized with kernels.cu

// Helper helper to convert your cuBQL::Triangle layout to a CGAL Triangle3 object
Triangle3 convertToCGAL(const cuBQL::Triangle& t) {
    return Triangle3(
        Point3(t.a.x, t.a.y, t.a.z),
        Point3(t.b.x, t.b.y, t.b.z),
        Point3(t.c.x, t.c.y, t.c.z)
    );
}

int main(int ac, char** av) {
    if (ac < 4) {
        std::cout << "Usage: " << av[0] << " <A.off> <B.off> <output.csv>\n";
        return 1;
    }

    std::string outputPath = av[3];
    double t_start = cuBQL::getCurrentTime();

    // 1. MESH LOADING
    double t0 = cuBQL::getCurrentTime();

    std::vector<cuBQL::Triangle> hA = cuBQL::samples::loadOFF(av[1]);
    std::vector<cuBQL::Triangle> hB = cuBQL::samples::loadOFF(av[2]);

    int NA = hA.size();
    int NB = hB.size();
    std::cout << "[Step 1] Load Meshes from Disk: " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - t0) << "s\n";

    // Call out into our CUDA module translation unit (Now collecting metrics)
    std::vector<IntersectionPair> hPairs;
    GPUTimingBreakdown gpuTimings; // <--- Target container instantiation
    long long totalAABBHits = runMeshIntersection(hA.data(), NA, hB.data(), NB, hPairs, gpuTimings);

    // --- NEW: TWO-LOOP CGAL FILTERING & CLEAN TIMING ---
    
    // LOOP 1: Convert all candidates to CGAL formats
    double t_conv_start = cuBQL::getCurrentTime();
    std::vector<CGALTrianglePair> cgalPairs;
    cgalPairs.reserve(hPairs.size());

    for (const auto& pair : hPairs) {
        cgalPairs.push_back({
            convertToCGAL(hA[pair.idA]),
            convertToCGAL(hB[pair.idB]),
            pair
        });
    }
    double time_conversion = cuBQL::getCurrentTime() - t_conv_start;

    // LOOP 2: Perform the robust geometric predicate checks
    double t_inter_start = cuBQL::getCurrentTime();
    std::vector<IntersectionPair> exactPairs;
    exactPairs.reserve(hPairs.size());

    for (const auto& cgalPair : cgalPairs) {
        if (CGAL::do_intersect(cgalPair.triA, cgalPair.triB)) {
            exactPairs.push_back(cgalPair.originalPair);
        }
    }
    double time_intersection = cuBQL::getCurrentTime() - t_inter_start;
    
    std::cout << "[Step 7.0] GPU Pipeline Breakdown:\n";
    std::cout << "  |---> Host->Device Upload:     " << cuBQL::prettyDouble(gpuTimings.uploadTime) << "s\n";
    std::cout << "  |---> BVH Build & Search:      " << cuBQL::prettyDouble(gpuTimings.executionTime) << "s\n";
    std::cout << "  |---> Device->Host Download:   " << cuBQL::prettyDouble(gpuTimings.downloadTime) << "s\n";
    
    std::cout << "[Step 7.5] CGAL Exact Filtering Total: " << cuBQL::prettyDouble(time_conversion + time_intersection) << "s\n";
    std::cout << "  |---> Data Conversion Loop:    " << cuBQL::prettyDouble(time_conversion) << "s\n";
    std::cout << "  |---> CGAL do_intersect Loop:  " << cuBQL::prettyDouble(time_intersection) << "s\n";
    // ----------------------------------------------------

    // 8. DISK I/O (Modified to write the verified exact pairs)
    t0 = cuBQL::getCurrentTime();
    if (!exactPairs.empty()) {
        std::ofstream outFile(outputPath);
        outFile << "ID_A,ID_B\n";
        for (const auto& p : exactPairs) outFile << p.idA << "," << p.idB << "\n";
    }
    std::cout << "[Step 8] Write CSV to Disk:     " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - t0) << "s\n";

    std::cout << "-------------------------------------------\n";
    std::cout << "Total runtime:      " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - t_start) << " s\n";
    std::cout << "GPU AABB Overlaps:  " << totalAABBHits << "\n";
    std::cout << "CGAL Exact Overlaps: " << exactPairs.size() << " (Dropped " << (totalAABBHits - exactPairs.size()) << " false positives)\n";

    return 0;
}