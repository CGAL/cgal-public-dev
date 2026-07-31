#include "GPUPredicatesCheckDoubleAssisted.h"
#include "GPUPredicatesCommon.h"

// --------------------------------------------------------------------
// ANCHOR-SHIFTED CLASSIFICATION FUNCTION
// --------------------------------------------------------------------
__device__ inline PairStatus classifyPair(
    const double3& Aa, const double3& Ab, const double3& Ac,
    const double3& Ba, const double3& Bb, const double3& Bc) 
{    
    // 1. Pick Vertex A.a as the Local Anchor Origin
    const double3 anchor = Aa;

    // 2. Perform 15 exact DOUBLE subtractions in thread registers
    double3 relA1 = { Ab.x - anchor.x, Ab.y - anchor.y, Ab.z - anchor.z };
    double3 relA2 = { Ac.x - anchor.x, Ac.y - anchor.y, Ac.z - anchor.z };

    double3 relB0 = { Ba.x - anchor.x, Ba.y - anchor.y, Ba.z - anchor.z };
    double3 relB1 = { Bb.x - anchor.x, Bb.y - anchor.y, Bb.z - anchor.z };
    double3 relB2 = { Bc.x - anchor.x, Bc.y - anchor.y, Bc.z - anchor.z };

    // 3. Cast relative coordinates down to float3 (cuBQL::vec3f)
    const cuBQL::vec3f A0 = { 0.0f, 0.0f, 0.0f };
    const cuBQL::vec3f A1 = { (float)relA1.x, (float)relA1.y, (float)relA1.z };
    const cuBQL::vec3f A2 = { (float)relA2.x, (float)relA2.y, (float)relA2.z };

    const cuBQL::vec3f B0 = { (float)relB0.x, (float)relB0.y, (float)relB0.z };
    const cuBQL::vec3f B1 = { (float)relB1.x, (float)relB1.y, (float)relB1.z };
    const cuBQL::vec3f B2 = { (float)relB2.x, (float)relB2.y, (float)relB2.z };

    // 4. Dynamic Local Noise Bound Computation
    float max_x = fmaxf(fabsf(A1.x), fmaxf(fabsf(A2.x), fmaxf(fabsf(B0.x), fmaxf(fabsf(B1.x), fabsf(B2.x)))));
    float max_y = fmaxf(fabsf(A1.y), fmaxf(fabsf(A2.y), fmaxf(fabsf(B0.y), fmaxf(fabsf(B1.y), fabsf(B2.y)))));
    float max_z = fmaxf(fabsf(A1.z), fmaxf(fabsf(A2.z), fmaxf(fabsf(B0.z), fmaxf(fabsf(B1.z), fabsf(B2.z)))));
    float M_local = fmaxf(max_x, fmaxf(max_y, max_z));

    const float float_machine_epsilon = 1.1920929e-7f; 
    float eps = 16.0f * (M_local * M_local * M_local) * float_machine_epsilon;

    // 5. Geometrical Predicate Evaluation Passes
    int ob0 = orient3d_interval(A0, A1, A2, B0, eps);
    int ob1 = orient3d_interval(A0, A1, A2, B1, eps);
    int ob2 = orient3d_interval(A0, A1, A2, B2, eps);
    if (ob0 != 0 && ob1 != 0 && ob2 != 0) {
        if ((ob0 > 0 && ob1 > 0 && ob2 > 0) || (ob0 < 0 && ob1 < 0 && ob2 < 0)) return PAIR_NO;
    }

    int oa0 = orient3d_interval(B0, B1, B2, A0, eps);
    int oa1 = orient3d_interval(B0, B1, B2, A1, eps);
    int oa2 = orient3d_interval(B0, B1, B2, A2, eps);
    if (oa0 != 0 && oa1 != 0 && oa2 != 0) {
        if ((oa0 > 0 && oa1 > 0 && oa2 > 0) || (oa0 < 0 && oa1 < 0 && oa2 < 0)) return PAIR_NO;
    }

    int r;
    r = edgeTriInterval(A0, A1, B0, B1, B2, eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;
    r = edgeTriInterval(A1, A2, B0, B1, B2, eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;
    r = edgeTriInterval(A2, A0, B0, B1, B2, eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;

    r = edgeTriInterval(B0, B1, A0, A1, A2, eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;
    r = edgeTriInterval(B1, B2, A0, A1, A2, eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;
    r = edgeTriInterval(B2, B0, A0, A1, A2, eps); if (r == PAIR_GREEN) return PAIR_GREEN; if (r == PAIR_YELLOW) return PAIR_YELLOW;

    return PAIR_NO;
}

// --------------------------------------------------------------------
// CORE INTERSECTION KERNEL
// --------------------------------------------------------------------
__global__ void evaluateGeometricPairsKernelDoubleAssisted(
    int* outStatuses, 
    const int2* candidatePairs, 
    const double3* dVertsA,
    const uint3* dIndicesA,
    const double3* dVertsB,
    const uint3* dIndicesB,
    int numPairs) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= numPairs) return;

    int2 pair = candidatePairs[tid];
    
    // Read triangle indices (uint3)
    uint3 idxA = dIndicesA[pair.x];
    uint3 idxB = dIndicesB[pair.y];

    // Fetch double3 vertices directly
    double3 Aa = dVertsA[idxA.x];
    double3 Ab = dVertsA[idxA.y];
    double3 Ac = dVertsA[idxA.z];

    double3 Ba = dVertsB[idxB.x];
    double3 Bb = dVertsB[idxB.y];
    double3 Bc = dVertsB[idxB.z];

    PairStatus status = classifyPair(Aa, Ab, Ac, Ba, Bb, Bc);

    outStatuses[tid] = (int)status;
}

// --------------------------------------------------------------------
// EXPORTED PIPELINE WRAPPER STAGE
// --------------------------------------------------------------------
void evaluateAndCompactPairsDoubleAssisted(
    const int2* dCandidatePairs,
    int* dPairStatuses,
    const double3* dVertsA,
    const uint3* dIndicesA,
    const double3* dVertsB,
    const uint3* dIndicesB,
    int totalBatchPairs,
    cudaStream_t stream)
{
    if (totalBatchPairs <= 0) return;

    int threadsPerBlock = 256;
    int blocksPerGrid = (totalBatchPairs + threadsPerBlock - 1) / threadsPerBlock;
    
    evaluateGeometricPairsKernelDoubleAssisted<<<blocksPerGrid, threadsPerBlock, 0, stream>>>(
         dPairStatuses, dCandidatePairs, 
         dVertsA, dIndicesA, 
         dVertsB, dIndicesB, 
         totalBatchPairs
    );
}