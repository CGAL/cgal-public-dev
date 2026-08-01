#include "GPUPredicatesCheckDoubleAssisted.h"
#include "GPUPredicatesCommon.h"

// --------------------------------------------------------------------
// ANCHOR-SHIFTED + FP32 NORMALIZED CLASSIFICATION FUNCTION
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

    // 3. Cast relative coordinates down to float3
    const cuBQL::vec3f A0 = { 0.0f, 0.0f, 0.0f };
    cuBQL::vec3f A1 = { (float)relA1.x, (float)relA1.y, (float)relA1.z };
    cuBQL::vec3f A2 = { (float)relA2.x, (float)relA2.y, (float)relA2.z };

    cuBQL::vec3f B0 = { (float)relB0.x, (float)relB0.y, (float)relB0.z };
    cuBQL::vec3f B1 = { (float)relB1.x, (float)relB1.y, (float)relB1.z };
    cuBQL::vec3f B2 = { (float)relB2.x, (float)relB2.y, (float)relB2.z };

    // 4. Compute Local Maximum Coordinate Bound (Tree Reduction for ILP)
    float mx1 = fmaxf(fabsf(A1.x), fabsf(A1.y));
    float mx2 = fmaxf(fabsf(A1.z), fabsf(A2.x));
    float mx3 = fmaxf(fabsf(A2.y), fabsf(A2.z));
    float mx4 = fmaxf(fabsf(B0.x), fabsf(B0.y));
    float mx5 = fmaxf(fabsf(B0.z), fabsf(B1.x));
    float mx6 = fmaxf(fabsf(B1.y), fabsf(B1.z));
    float mx7 = fmaxf(fabsf(B2.x), fabsf(B2.y));
    float mx8 = fabsf(B2.z);

    float m12 = fmaxf(mx1, mx2);
    float m34 = fmaxf(mx3, mx4);
    float m56 = fmaxf(mx5, mx6);
    float m78 = fmaxf(mx7, mx8);

    float m1234 = fmaxf(m12, m34);
    float m5678 = fmaxf(m56, m78);
    float M_local = fmaxf(m1234, m5678);

    // 5. Early exit for extreme degeneracies (Avoids pointless math)
    if (M_local <= 1e-12f) {
        return PAIR_YELLOW;
    }

    // 6. Apply FP32 Normalization Scale using Fast Hardware Reciprocal
    float inv_scale = __frcp_rn(M_local);
    A1 = { A1.x * inv_scale, A1.y * inv_scale, A1.z * inv_scale };
    A2 = { A2.x * inv_scale, A2.y * inv_scale, A2.z * inv_scale };

    B0 = { B0.x * inv_scale, B0.y * inv_scale, B0.z * inv_scale };
    B1 = { B1.x * inv_scale, B1.y * inv_scale, B1.z * inv_scale };
    B2 = { B2.x * inv_scale, B2.y * inv_scale, B2.z * inv_scale };

    // 7. Constrained Noise Bound (Since M_local = 1.0 after normalization)
    const float float_machine_epsilon = 1.1920929e-7f; 
    const float eps = 16.0f * float_machine_epsilon;

    // 8. Geometrical Predicate Evaluation Passes
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

    // Fetch components explicitly via __ldg
    double3 Aa = { __ldg(&dVertsA[idxA.x].x), __ldg(&dVertsA[idxA.x].y), __ldg(&dVertsA[idxA.x].z) };
    double3 Ab = { __ldg(&dVertsA[idxA.y].x), __ldg(&dVertsA[idxA.y].y), __ldg(&dVertsA[idxA.y].z) };
    double3 Ac = { __ldg(&dVertsA[idxA.z].x), __ldg(&dVertsA[idxA.z].y), __ldg(&dVertsA[idxA.z].z) };

    double3 Ba = { __ldg(&dVertsB[idxB.x].x), __ldg(&dVertsB[idxB.x].y), __ldg(&dVertsB[idxB.x].z) };
    double3 Bb = { __ldg(&dVertsB[idxB.y].x), __ldg(&dVertsB[idxB.y].y), __ldg(&dVertsB[idxB.y].z) };
    double3 Bc = { __ldg(&dVertsB[idxB.z].x), __ldg(&dVertsB[idxB.z].y), __ldg(&dVertsB[idxB.z].z) };

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