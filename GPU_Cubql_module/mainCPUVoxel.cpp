#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <mutex>

// TBB Parallel Core Infrastructure
#include <tbb/task_group.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_vector.h>  // <-- Added for thread-safe value packing

// Minimal vector and bounding box math infrastructure to mirror cuBQL layout
struct Vec3f {
    float x, y, z;
    Vec3f() : x(0), y(0), z(0) {}
    Vec3f(float x, float y, float z) : x(x), y(y), z(z) {}
    
    Vec3f operator+(const Vec3f& o) const { return {x+o.x, y+o.y, z+o.z}; }
    Vec3f operator-(const Vec3f& o) const { return {x-o.x, y-o.y, z-o.z}; }
    Vec3f operator*(float s) const { return {x*s, y*s, z*s}; }
    Vec3f operator/(const Vec3f& o) const { return {x/o.x, y/o.y, z/o.z}; }
};

struct Box3f {
    Vec3f lower;
    Vec3f upper;
    Box3f() {
        lower = { 1e20f,  1e20f,  1e20f};
        upper = {-1e20f, -1e20f, -1e20f};
    }
    Box3f(Vec3f l, Vec3f u) : lower(l), upper(u) {}

    void extend(const Vec3f& p) {
        lower.x = std::min(lower.x, p.x); lower.y = std::min(lower.y, p.y); lower.z = std::min(lower.z, p.z);
        upper.x = std::max(upper.x, p.x); upper.y = std::max(upper.y, p.y); upper.z = std::max(upper.z, p.z);
    }
    void extend(const Box3f& b) {
        extend(b.lower);
        extend(b.upper);
    }
    Vec3f size() const { return upper - lower; }
};

struct CPUMeshTriangle {
    Vec3f a, b, c;
    Box3f bounds() const {
        Box3f box;
        box.extend(a); box.extend(b); box.extend(c);
        return box;
    }
};

// Simple Hash Combiner function for 3D Discrete Grid Locations
struct GridKey {
    int x, y, z;
    bool operator==(const GridKey& o) const { return x == o.x && y == o.y && z == o.z; }
};

struct GridKeyHash {
    std::size_t operator()(const GridKey& k) const {
        // High-performance spatial hashing murmur-style shift mixes
        std::size_t h1 = std::hash<int>()(k.x);
        std::size_t h2 = std::hash<int>()(k.y);
        std::size_t h3 = std::hash<int>()(k.z);
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

// FIX: Value storage must be concurrent_vector to prevent data-races during bucket execution
using SpatialHashMap = tbb::concurrent_unordered_map<GridKey, tbb::concurrent_vector<int>, GridKeyHash>;

// Mock OFF Parser (Reads basic geometry headers and triangle lists)
bool loadDummyOrOffMesh(const std::string& path, std::vector<CPUMeshTriangle>& outTris) {
    std::ifstream f(path);
    if (!f.is_open()) return false;
    
    std::string head;
    f >> head; // Expect "OFF"
    int nVerts = 0, nFaces = 0, nEdges = 0;
    f >> nVerts >> nFaces >> nEdges;
    
    std::vector<Vec3f> verts(nVerts);
    for (int i = 0; i < nVerts; ++i) {
        f >> verts[i].x >> verts[i].y >> verts[i].z;
    }
    for (int i = 0; i < nFaces; ++i) {
        int sideCount = 0;
        f >> sideCount;
        if (sideCount == 3) {
            int ia, ib, ic;
            f >> ia >> ib >> ic;
            outTris.push_back({verts[ia], verts[ib], verts[ic]});
        } else {
            int dummy;
            for(int j=0; j<sideCount; ++j) f >> dummy;
        }
    }
    return true;
}

int main(int ac, char** av) {
    if (ac < 4) {
        std::cout << "Usage: " << av[0] << " <A.off> <B.off> <grid_resolution_n>\n";
        return 1;
    }

    std::string pathA = av[1];
    std::string pathB = av[2];
    int n = std::stoi(av[3]);

    auto t_start = std::chrono::high_resolution_clock::now();

    // 1. Load Meshes Parallelly from storage disks
    std::vector<CPUMeshTriangle> meshA, meshB;
    bool loadedA = false, loadedB = false;
    
    tbb::task_group tg;
    tg.run([&]() { loadedA = loadDummyOrOffMesh(pathA, meshA); });
    tg.run([&]() { loadedB = loadDummyOrOffMesh(pathB, meshB); });
    tg.wait(); 

    if (!loadedA || !loadedB) {
        std::cerr << "Fatal Error parsing input files.\n";
        return 1;
    }

    auto t_loaded = std::chrono::high_resolution_clock::now();
    std::cout << "[Step 1] CPU Mesh Loading Phase Completed: " 
              << std::chrono::duration<double>(t_loaded - t_start).count() << "s\n";

    // 2. Uniform Bounds Calculation via reductions
    Box3f worldBounds;
    for (const auto& tri : meshA) worldBounds.extend(tri.bounds());
    for (const auto& tri : meshB) worldBounds.extend(tri.bounds());

    Vec3f worldSize = worldBounds.size();
    Vec3f invWorldSize(1.0f / worldSize.x, 1.0f / worldSize.y, 1.0f / worldSize.z);

    int dimX = n, dimY = n, dimZ = n;
    std::cout << "--> Target Core Matrix Grid Allocation Size: [" << dimX << " x " << dimY << " x " << dimZ << "]\n";

    SpatialHashMap mapA;
    SpatialHashMap mapB;

    auto t_voxel_start = std::chrono::high_resolution_clock::now();

    // 3. Parallel TBB Map Insertion Loop: Mesh A Voxel Classification
    tbb::parallel_for(tbb::blocked_range<size_t>(0, meshA.size()), [&](const tbb::blocked_range<size_t>& range) {
        for (size_t i = range.begin(); i != range.end(); ++i) {
            Box3f triBox = meshA[i].bounds();
            
            int minX = std::max(0, std::min(dimX - 1, static_cast<int>((triBox.lower.x - worldBounds.lower.x) * invWorldSize.x * dimX)));
            int maxX = std::max(0, std::min(dimX - 1, static_cast<int>((triBox.upper.x - worldBounds.lower.x) * invWorldSize.x * dimX)));
            int minY = std::max(0, std::min(dimY - 1, static_cast<int>((triBox.lower.y - worldBounds.lower.y) * invWorldSize.y * dimY)));
            int maxY = std::max(0, std::min(dimY - 1, static_cast<int>((triBox.upper.y - worldBounds.lower.y) * invWorldSize.y * dimY)));
            int minZ = std::max(0, std::min(dimZ - 1, static_cast<int>((triBox.lower.z - worldBounds.lower.z) * invWorldSize.z * dimZ)));
            int maxZ = std::max(0, std::min(dimZ - 1, static_cast<int>((triBox.upper.z - worldBounds.lower.z) * invWorldSize.z * dimZ)));

            for (int x = minX; x <= maxX; ++x) {
                for (int y = minY; y <= maxY; ++y) {
                    for (int z = minZ; z <= maxZ; ++z) {
                        mapA[{x, y, z}].push_back(static_cast<int>(i));
                    }
                }
            }
        }
    });

    // 4. Parallel TBB Map Insertion Loop: Mesh B Voxel Classification
    tbb::parallel_for(tbb::blocked_range<size_t>(0, meshB.size()), [&](const tbb::blocked_range<size_t>& range) {
        for (size_t i = range.begin(); i != range.end(); ++i) {
            Box3f triBox = meshB[i].bounds();
            
            int minX = std::max(0, std::min(dimX - 1, static_cast<int>((triBox.lower.x - worldBounds.lower.x) * invWorldSize.x * dimX)));
            int maxX = std::max(0, std::min(dimX - 1, static_cast<int>((triBox.upper.x - worldBounds.lower.x) * invWorldSize.x * dimX)));
            int minY = std::max(0, std::min(dimY - 1, static_cast<int>((triBox.lower.y - worldBounds.lower.y) * invWorldSize.y * dimY)));
            int maxY = std::max(0, std::min(dimY - 1, static_cast<int>((triBox.upper.y - worldBounds.lower.y) * invWorldSize.y * dimY)));
            int minZ = std::max(0, std::min(dimZ - 1, static_cast<int>((triBox.lower.z - worldBounds.lower.z) * invWorldSize.z * dimZ)));
            int maxZ = std::max(0, std::min(dimZ - 1, static_cast<int>((triBox.upper.z - worldBounds.lower.z) * invWorldSize.z * dimZ)));

            for (int x = minX; x <= maxX; ++x) {
                for (int y = minY; y <= maxY; ++y) {
                    for (int z = minZ; z <= maxZ; ++z) {
                        mapB[{x, y, z}].push_back(static_cast<int>(i));
                    }
                }
            }
        }
    });

    auto t_voxel_end = std::chrono::high_resolution_clock::now();
    double voxelizationTime = std::chrono::duration<double>(t_voxel_end - t_voxel_start).count();

    // 5. Gather statistics
    size_t totalRefA = 0;
    for (auto const& [key, list] : mapA) totalRefA += list.size();
    size_t totalRefB = 0;
    for (auto const& [key, list] : mapB) totalRefB += list.size();

    std::cout << "----------------------------------------------------\n";
    std::cout << "[Step 3 Execution] TBB CPU Multi-Thread Voxelizer Profiler Breakdown:\n";
    std::cout << "  |---> Parallel Hashing Voxelization Loop: " << voxelizationTime * 1000.0 << " ms\n";
    std::cout << "----------------------------------------------------\n";
    std::cout << "====================================================\n";
    std::cout << "Macro Execution Statistics Summary:\n";
    std::cout << "  |-- Total Process Runtime:              " << std::chrono::duration<double>(t_voxel_end - t_start).count() << " s\n";
    std::cout << "  |-- Active Buckets Contained (Mesh A):   " << mapA.size() << "\n";
    std::cout << "  |-- Active Buckets Contained (Mesh B):   " << mapB.size() << "\n";
    std::cout << "  |-- Total Mesh A Contained References:  " << totalRefA << "\n";
    std::cout << "  |-- Total Mesh B Contained References:  " << totalRefB << "\n";
    std::cout << "====================================================\n";

    return 0;
}