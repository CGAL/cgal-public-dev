#include <iostream>
#include <vector>
#include <fstream>
#include "include/loadOFFCGAL.h" // Clean header-only include!

struct IntersectionPair {
    int idA;
    int idB;
};

// Declare the external function compiled by nvcc
extern "C" long long runMeshIntersection(
    const cuBQL::Triangle* hA, int NA, 
    const cuBQL::Triangle* hB, int NB, 
    std::vector<IntersectionPair>& hPairs);

int main(int ac, char** av) {
    if (ac < 4) {
        std::cout << "Usage: " << av[0] << " <A.off> <B.off> <output.csv>\n";
        return 1;
    }

    std::string outputPath = av[3];
    double t_start = cuBQL::getCurrentTime();

    // 1. MESH LOADING
    double t0 = cuBQL::getCurrentTime();
    std::vector<cuBQL::Triangle> hA = cuBQL::samples::loadOFFCGAL(av[1]);
    std::vector<cuBQL::Triangle> hB = cuBQL::samples::loadOFFCGAL(av[2]);
    int NA = hA.size();
    int NB = hB.size();
    std::cout << "[Step 1] Load Meshes from Disk: " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - t0) << "s\n";

    // Call out into our CUDA module translation unit
    std::vector<IntersectionPair> hPairs;
    long long totalIntersections = runMeshIntersection(hA.data(), NA, hB.data(), NB, hPairs);

    // 8. DISK I/O
    t0 = cuBQL::getCurrentTime();
    if (totalIntersections > 0) {
        std::ofstream outFile(outputPath);
        outFile << "ID_A,ID_B\n";
        for (const auto& p : hPairs) outFile << p.idA << "," << p.idB << "\n";
    }
    std::cout << "[Step 8] Write CSV to Disk:     " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - t0) << "s\n";

    std::cout << "-------------------------------------------\n";
    std::cout << "Total runtime:  " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - t_start) << " s\n";
    std::cout << "Total overlaps: " << totalIntersections << "\n";

    return 0;
}