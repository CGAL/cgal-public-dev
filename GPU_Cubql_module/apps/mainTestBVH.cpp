#include <iostream>
#include <vector>
#include <fstream>
#include <string>

// --- cuBQL INCLUDES ---
#include "include/loadOFFCGAL.h"
#include "include/loadOFF.h"

// --- CGAL INCLUDES ---
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

// --- TBB INCLUDES ----
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef Kernel::Point_3 Point3;

// --- EXTERNAL C LINKAGE FOR YOUR UPGRADED TEST KERNEL ---
//extern "C" void kernelsTestBVH(const cuBQL::Triangle* hMeshA, int numTrianglesA, int maxCellSizeA,
 //                              const cuBQL::Triangle* hMeshB, int numTrianglesB, int maxCellSizeB);
 extern "C" void kernelsTestBVH(const cuBQL::Triangle* hMeshA, int numTrianglesA, int maxCellSizeA,
                               const cuBQL::Triangle* hMeshB, int numTrianglesB, int maxCellSizeB, int batchMultiplier);

// Helper function to translate geometry from CGAL storage to raw floats
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

// Reusable translation pipeline helper
std::vector<cuBQL::Triangle> processMeshLayout(const Mesh& mesh) {
    size_t numFaces = num_faces(mesh);
    std::vector<face_descriptor> faces;
    faces.reserve(numFaces);
    for (auto f : mesh.faces()) {
        faces.push_back(f);
    }

    std::vector<cuBQL::Triangle> hMeshLayout(numFaces);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, numFaces), [&](const tbb::blocked_range<size_t>& r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
            hMeshLayout[i] = convertCgalFaceToCuBQL(mesh, faces[i]);
        }
    });
    return hMeshLayout;
}

int main(int ac, char** av) {
    if (ac < 6) {
        std::cout << "Usage: " << av[0] << " <meshA.off> <maxCellSizeA> <meshB.off> <maxCellSizeB><batchmultiplier)\n";
        return 1;
    }

    std::string meshPathA = av[1];
    int maxCellSizeA    = std::stoi(av[2]);
    std::string meshPathB = av[3];
    int maxCellSizeB    = std::stoi(av[4]);
    int batchmultipl = std::stoi(av[5]); 
    
    // Default to a sane baseline for thread pooling during parsing
    tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, 8);

    double tStartTotal = cuBQL::getCurrentTime();

    // --------------------------------------------------------------------
    // LOAD AND TRANSLATE MESH A
    // --------------------------------------------------------------------
    double tMeshAStart = cuBQL::getCurrentTime();
    Mesh meshA;
    std::ifstream inFileA(meshPathA);
    if (!(inFileA >> meshA)) {
        std::cerr << "Critical Error: Failed to open or parse mesh file A: " << meshPathA << "\n";
        return 1;
    }
    inFileA.close();
    
    std::vector<cuBQL::Triangle> hMeshLayoutA = processMeshLayout(meshA);
    std::cout << "[Step 1] Loaded and Processed Mesh A (" << hMeshLayoutA.size() << " tris): " 
              << cuBQL::prettyDouble(cuBQL::getCurrentTime() - tMeshAStart) << "s\n";

    // --------------------------------------------------------------------
    // LOAD AND TRANSLATE MESH B
    // --------------------------------------------------------------------
    double tMeshBStart = cuBQL::getCurrentTime();
    Mesh meshB;
    std::ifstream inFileB(meshPathB);
    if (!(inFileB >> meshB)) {
        std::cerr << "Critical Error: Failed to open or parse mesh file B: " << meshPathB << "\n";
        return 1;
    }
    inFileB.close();

    std::vector<cuBQL::Triangle> hMeshLayoutB = processMeshLayout(meshB);
    std::cout << "[Step 2] Loaded and Processed Mesh B (" << hMeshLayoutB.size() << " tris): " 
              << cuBQL::prettyDouble(cuBQL::getCurrentTime() - tMeshBStart) << "s\n";

    // --------------------------------------------------------------------
    // EXECUTE SPATIAL CROSS-INTERSECTION PIPELINE
    // --------------------------------------------------------------------
    std::cout << "[Step 3] Launching Dual-Mesh GPU Pipeline...\n";
    kernelsTestBVH(
        hMeshLayoutA.data(), static_cast<int>(hMeshLayoutA.size()), maxCellSizeA,
        hMeshLayoutB.data(), static_cast<int>(hMeshLayoutB.size()), maxCellSizeB, batchmultipl
    );

    std::cout << "Total Diagnostic App Execution Time: " 
              << cuBQL::prettyDouble(cuBQL::getCurrentTime() - tStartTotal) << "s\n";
    std::cout << "==================================================\n";

    return 0;
}