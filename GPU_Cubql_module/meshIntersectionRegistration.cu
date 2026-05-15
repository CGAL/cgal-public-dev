#define CUBQL_GPU_BUILDER_IMPLEMENTATION 1
#include "cuBQL/bvh.h"
#include "cuBQL/traversal/fixedBoxQuery.h"
#include "samples/common/loadOBJ.h"
#include "loadOFF.h"
#include <thrust/device_ptr.h>
#include <thrust/scan.h>
#include <thrust/reduce.h>
#include <fstream>
#include <vector>
#include <iomanip>

struct IntersectionPair {
    int idA;
    int idB;
};

// ------------------------------------------------------------
__global__ void generateBoxes(cuBQL::box3f *boxes, const cuBQL::Triangle *tris, int N) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < N) boxes[i] = tris[i].bounds();
}

__device__ bool intersects(const cuBQL::Triangle &a, const cuBQL::Triangle &b) {
    return a.bounds().overlaps(b.bounds());
}

__global__ void countKernel(int *counts, const cuBQL::Triangle *triA, cuBQL::bvh3f bvhA, const cuBQL::Triangle *triB, int N) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= N) return;
    cuBQL::Triangle b = triB[tid];
    cuBQL::box3f query = b.bounds();
    int count = 0;
    cuBQL::fixedBoxQuery::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        for (int i = 0; i < num; i++)
            if (intersects(triA[ids[i]], b)) count++;
        return CUBQL_CONTINUE_TRAVERSAL;
    }, bvhA, query);
    counts[tid] = count;
}

__global__ void fillKernel(IntersectionPair *pairs, const int *offsets, const cuBQL::Triangle *triA, cuBQL::bvh3f bvhA, const cuBQL::Triangle *triB, int N) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= N) return;
    int writePos = offsets[tid];
    cuBQL::Triangle b = triB[tid];
    cuBQL::box3f query = b.bounds();
    cuBQL::fixedBoxQuery::forEachLeaf([&](const uint32_t *ids, uint32_t num) {
        for (int i = 0; i < num; i++) {
            if (intersects(triA[ids[i]], b)) {
                pairs[writePos++] = { (int)ids[i], tid };
            }
        }
        return CUBQL_CONTINUE_TRAVERSAL;
    }, bvhA, query);
}

// ------------------------------------------------------------
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

    // 2. GPU ALLOCATION & UPLOAD
    t0 = cuBQL::getCurrentTime();
    cuBQL::Triangle *dA, *dB;
    CUBQL_CUDA_CALL(Malloc(&dA, NA * sizeof(cuBQL::Triangle)));
    CUBQL_CUDA_CALL(Malloc(&dB, NB * sizeof(cuBQL::Triangle)));
    CUBQL_CUDA_CALL(Memcpy(dA, hA.data(), NA * sizeof(cuBQL::Triangle), cudaMemcpyDefault));
    CUBQL_CUDA_CALL(Memcpy(dB, hB.data(), NB * sizeof(cuBQL::Triangle), cudaMemcpyDefault));
    CUBQL_CUDA_SYNC_CHECK();
    std::cout << "[Step 2] GPU Alloc & Upload:    " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - t0) << "s\n";

    // 3. BVH BUILD
    t0 = cuBQL::getCurrentTime();
    cuBQL::box3f *dBoxes;
    CUBQL_CUDA_CALL(Malloc(&dBoxes, NA * sizeof(cuBQL::box3f)));
    generateBoxes<<<cuBQL::divRoundUp(NA, 256), 256>>>(dBoxes, dA, NA);
    
    cuBQL::bvh3f bvh;
    cuBQL::gpuBuilder(bvh, dBoxes, NA, cuBQL::BuildConfig());
    CUBQL_CUDA_SYNC_CHECK();
    CUBQL_CUDA_CALL(Free(dBoxes));
    std::cout << "[Step 3] BVH Generation:        " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - t0) << "s\n";

    // 4. QUERY PASS 1 (COUNTING)
    t0 = cuBQL::getCurrentTime();
    int *dCounts, *dOffsets;
    CUBQL_CUDA_CALL(Malloc(&dCounts, NB * sizeof(int)));
    CUBQL_CUDA_CALL(Malloc(&dOffsets, NB * sizeof(int)));
    
    countKernel<<<cuBQL::divRoundUp(NB, 128), 128>>>(dCounts, dA, bvh, dB, NB);
    CUBQL_CUDA_SYNC_CHECK();
    std::cout << "[Step 4] Count Kernel (Pass 1): " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - t0) << "s\n";

    // 5. THRUST REDUCE & SCAN (STAYS ON GPU)
    t0 = cuBQL::getCurrentTime();
    thrust::device_ptr<int> dev_counts(dCounts);
    thrust::device_ptr<int> dev_offsets(dOffsets);
    long long totalIntersections = thrust::reduce(dev_counts, dev_counts + NB, 0LL, thrust::plus<long long>());
    thrust::exclusive_scan(dev_counts, dev_counts + NB, dev_offsets);
    CUBQL_CUDA_SYNC_CHECK();
    std::cout << "[Step 5] Thrust Reduce & Scan:  " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - t0) << "s\n";

    // 6. QUERY PASS 2 (FILLING)
    t0 = cuBQL::getCurrentTime();
    IntersectionPair *dPairs = nullptr;
    if (totalIntersections > 0) {
        CUBQL_CUDA_CALL(Malloc(&dPairs, totalIntersections * sizeof(IntersectionPair)));
        fillKernel<<<cuBQL::divRoundUp(NB, 128), 128>>>(dPairs, dOffsets, dA, bvh, dB, NB);
        CUBQL_CUDA_SYNC_CHECK();
    }
    std::cout << "[Step 6] Fill Kernel (Pass 2):  " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - t0) << "s\n";

    // 7. DOWNLOAD FROM GPU
    t0 = cuBQL::getCurrentTime();
    std::vector<IntersectionPair> hPairs;
    if (totalIntersections > 0) {
        hPairs.resize(totalIntersections);
        CUBQL_CUDA_CALL(Memcpy(hPairs.data(), dPairs, totalIntersections * sizeof(IntersectionPair), cudaMemcpyDefault));
    }
    std::cout << "[Step 7] Download Results:      " << cuBQL::prettyDouble(cuBQL::getCurrentTime() - t0) << "s\n";

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

    CUBQL_CUDA_CALL(Free(dA)); CUBQL_CUDA_CALL(Free(dB));
    CUBQL_CUDA_CALL(Free(dCounts)); CUBQL_CUDA_CALL(Free(dOffsets));
    if (dPairs) CUBQL_CUDA_CALL(Free(dPairs));
    cuBQL::cuda::free(bvh);

    return 0;
}
