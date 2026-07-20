#include "include/third-party/cubql/fixedBoxQueryv2.h"
#include <thrust/device_ptr.h>
#include <thrust/scan.h>
#include <thrust/execution_policy.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <atomic>
#include <thread>
#include <tbb/task_group.h>
#include <tbb/concurrent_vector.h>

#include "batchedCrossIntersectionV2.h"
#include "../custom_pipeline/GPUPredicatesCheckV2.h"
#include "TargetStatus.h"
#include "../src/CPU/YellowFilter.h"

#ifndef CUBQL_CUDA_CALL
#define CUBQL_CUDA_CALL(call)                                                                                          \
  {                                                                                                                    \
    cudaError_t err = cuda##call;                                                                                      \
    if(cudaSuccess != err) {                                                                                           \
      fprintf(stderr, "CUDA error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(err));      \
      exit(EXIT_FAILURE);                                                                                              \
    }                                                                                                                  \
  }
#endif




// --------------------------------------------------------------------
// WARP-AGGREGATED COMPACTION KERNEL
// --------------------------------------------------------------------

__global__ void compactGreenYellowKernel(const int2* __restrict__ candidatePairs,
                                         const int* __restrict__ pairStatuses,
                                         int totalPairs,
                                         int2* __restrict__ greenOut,
                                         int* __restrict__ d_greenCount,
                                         int2* __restrict__ yellowOut,
                                         int* __restrict__ d_yellowCount) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int status = (tid < totalPairs) ? pairStatuses[tid] : -1;

  // 1. Green Compaction via Warp Aggregation
  bool isGreen = (status == PAIR_GREEN);
  unsigned maskGreen = __ballot_sync(0xFFFFFFFF, isGreen);
  int warpGreenCount = __popc(maskGreen);
  if(warpGreenCount > 0) {
    int laneId = threadIdx.x & 31;
    int warpGreenOffset = 0;
    if(laneId == 0) {
      warpGreenOffset = atomicAdd(d_greenCount, warpGreenCount);
    }
    warpGreenOffset = __shfl_sync(0xFFFFFFFF, warpGreenOffset, 0);
    if(isGreen) {
      int greenLaneIdx = __popc(maskGreen & ((1u << laneId) - 1));
      greenOut[warpGreenOffset + greenLaneIdx] = candidatePairs[tid];
    }
  }

  // 2. Yellow Compaction via Warp Aggregation
  bool isYellow = (status == PAIR_YELLOW);
  unsigned maskYellow = __ballot_sync(0xFFFFFFFF, isYellow);
  int warpYellowCount = __popc(maskYellow);
  if(warpYellowCount > 0) {
    int laneId = threadIdx.x & 31;
    int warpYellowOffset = 0;
    if(laneId == 0) {
      warpYellowOffset = atomicAdd(d_yellowCount, warpYellowCount);
    }
    warpYellowOffset = __shfl_sync(0xFFFFFFFF, warpYellowOffset, 0);
    if(isYellow) {
      int yellowLaneIdx = __popc(maskYellow & ((1u << laneId) - 1));
      yellowOut[warpYellowOffset + yellowLaneIdx] = candidatePairs[tid];
    }
  }
}

// --------------------------------------------------------------------
// PRE-ALLOCATED MICROBATCH PIPELINE SLOT
// --------------------------------------------------------------------

constexpr int MICROBATCH_SIZE = 65536; // Size of each internal microbatch

struct MicrobatchSlot
{
  cudaStream_t stream = nullptr;

  // Fixed-size Device Buffers for MICROBATCH_SIZE
  int* d_pairStatuses = nullptr;
  int2* d_green = nullptr;
  int2* d_yellow = nullptr;
  int* d_greenCount = nullptr;
  int* d_yellowCount = nullptr;

  // Fixed-size Pinned Host Buffers
  int2* hGreenPinned = nullptr;
  int2* hYellowPinned = nullptr;
  int* hGreenCountPinned = nullptr;
  int* hYellowCountPinned = nullptr;

  std::atomic<bool> busy{false};

  void allocate() {
    CUBQL_CUDA_CALL(StreamCreateWithFlags(&stream, cudaStreamNonBlocking));

    CUBQL_CUDA_CALL(Malloc(&d_pairStatuses, MICROBATCH_SIZE * sizeof(int)));
    CUBQL_CUDA_CALL(Malloc(&d_green, MICROBATCH_SIZE * sizeof(int2)));
    CUBQL_CUDA_CALL(Malloc(&d_yellow, MICROBATCH_SIZE * sizeof(int2)));
    CUBQL_CUDA_CALL(Malloc(&d_greenCount, sizeof(int)));
    CUBQL_CUDA_CALL(Malloc(&d_yellowCount, sizeof(int)));

    CUBQL_CUDA_CALL(MallocHost((void**)&hGreenPinned, MICROBATCH_SIZE * sizeof(int2)));
    CUBQL_CUDA_CALL(MallocHost((void**)&hYellowPinned, MICROBATCH_SIZE * sizeof(int2)));
    CUBQL_CUDA_CALL(MallocHost((void**)&hGreenCountPinned, sizeof(int)));
    CUBQL_CUDA_CALL(MallocHost((void**)&hYellowCountPinned, sizeof(int)));
  }

  void destroy() {
    CUBQL_CUDA_CALL(Free(d_pairStatuses));
    CUBQL_CUDA_CALL(Free(d_green));
    CUBQL_CUDA_CALL(Free(d_yellow));
    CUBQL_CUDA_CALL(Free(d_greenCount));
    CUBQL_CUDA_CALL(Free(d_yellowCount));

    CUBQL_CUDA_CALL(FreeHost(hGreenPinned));
    CUBQL_CUDA_CALL(FreeHost(hYellowPinned));
    CUBQL_CUDA_CALL(FreeHost(hGreenCountPinned));
    CUBQL_CUDA_CALL(FreeHost(hYellowCountPinned));

    CUBQL_CUDA_CALL(StreamDestroy(stream));
  }
};

// --------------------------------------------------------------------
// EXISTING ASSEMBLY & TRAVERSAL KERNELS
// --------------------------------------------------------------------

__global__ void assembleChunkBuffersByBatchKernel_V3(uint32_t* d_BIter,
                                                  uint64_t* d_AIter,
                                                  const uint32_t* outPairsA,
                                                  const uint32_t* outPairsB,
                                                  const uint32_t* reverseMapB,
                                                  const uint32_t* outOffsetsB,
                                                  const uint32_t* outPrimsFlatB,
                                                  uint32_t outMarkedCountB,
                                                  uint32_t totalPrimsB,
                                                  int chunkStartBatchIdx,
                                                  int activeBatchesInChunk,
                                                  const int* d_chunkBatchOffsets) {
  int b = threadIdx.x + blockIdx.x * blockDim.x;
  if(b >= activeBatchesInChunk)
    return;

  int currentBatchArrayIdx = chunkStartBatchIdx + b;
  uint32_t bIdx = outPairsB[currentBatchArrayIdx];
  uint32_t markedIdxB = reverseMapB[bIdx];

  uint32_t startOffsetB = outOffsetsB[markedIdxB];
  uint32_t endOffsetB = (markedIdxB + 1 < outMarkedCountB) ? outOffsetsB[markedIdxB + 1] : totalPrimsB;
  int numPrims = (int)(endOffsetB - startOffsetB);

  if(numPrims <= 0)
    return;

  int writeSandboxOffset = d_chunkBatchOffsets[b];
  uint64_t startNodeIdxA = (uint64_t)(outPairsA[currentBatchArrayIdx]);

  for(int p = 0; p < numPrims; ++p) {
    int targetWriteIdx = writeSandboxOffset + p;
    d_BIter[targetWriteIdx] = outPrimsFlatB[startOffsetB + p];
    d_AIter[targetWriteIdx] = startNodeIdxA;
  }
}

__global__ void countAABBOverlapsKernel_Indirected_V3(int* pairCounts,
                                                   cuBQL::bvh3f bvhA,
                                                   const cuBQL::Triangle* triA,
                                                   const cuBQL::Triangle* triB,
                                                   const uint32_t* d_BIter,
                                                   uint32_t startOffsetB,
                                                   int numPrimsB,
                                                   const uint64_t* d_AIter) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= numPrimsB)
    return;

  uint32_t actualPrimIdB = d_BIter[startOffsetB + tid];
  uint64_t startNodeIdxA = d_AIter[startOffsetB + tid];
  cuBQL::Triangle b = triB[actualPrimIdB];
  cuBQL::box3f query = b.bounds();

  int count = 0;
  cuBQL::fixedBoxQueryv2::forEachLeaf(
      [&](const uint32_t* ids, uint32_t num) {
        for(uint32_t i = 0; i < num; i++) {
          if(triA[ids[i]].bounds().overlaps(query))
            count++;
        }
        return CUBQL_CONTINUE_TRAVERSAL;
      },
      bvhA, query, startNodeIdxA);

  pairCounts[tid] = count;
}

__global__ void fillAABBOverlapsKernel_Indirected_V3(int2* candidatePairs,
                                                  const int* offsets,
                                                  cuBQL::bvh3f bvhA,
                                                  const cuBQL::Triangle* triA,
                                                  const cuBQL::Triangle* triB,
                                                  const uint32_t* d_BIter,
                                                  uint32_t startOffsetB,
                                                  int numPrimsB,
                                                  const uint64_t* d_AIter) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= numPrimsB)
    return;

  int wPos = offsets[tid];
  uint32_t actualPrimIdB = d_BIter[startOffsetB + tid];
  uint64_t startNodeIdxA = d_AIter[startOffsetB + tid];
  cuBQL::Triangle b = triB[actualPrimIdB];
  cuBQL::box3f query = b.bounds();

  cuBQL::fixedBoxQueryv2::forEachLeaf(
      [&](const uint32_t* ids, uint32_t num) {
        for(uint32_t i = 0; i < num; i++) {
          if(triA[ids[i]].bounds().overlaps(query)) {
            candidatePairs[wPos++] = make_int2((int)ids[i], (int)actualPrimIdB);
          }
        }
        return CUBQL_CONTINUE_TRAVERSAL;
      },
      bvhA, query, startNodeIdxA);
}

// --------------------------------------------------------------------
// MAIN EXECUTABLE LOOP WITH INTERNAL MICROBATCHING
// ----------------------------------------------------------------fillAABBOverlapsKernel----

uint64_t executeBatchedCrossIntersectionLoopV3(Mesh& meshAcpu,
                                               Mesh& meshBcpu,
                                               int batchMultiplier,
                                               int totalBatches,
                                               const thrust::device_vector<uint32_t>& d_outPairsA,
                                               const thrust::device_vector<uint32_t>& d_outPairsB,
                                               const thrust::device_vector<uint32_t>& d_reverseMapB,
                                               const thrust::device_vector<uint32_t>& d_markedNodeIndicesB,
                                               const thrust::device_vector<uint32_t>& d_outOffsetsB,
                                               const thrust::device_vector<uint32_t>& d_outPrimsFlatB,
                                               const thrust::device_vector<uint32_t>& d_nodeDescendantCountsB,
                                               uint32_t h_outMarkedCountB,
                                               const cuBQL::bvh3f& bvhA,
                                               const cuBQL::Triangle* dMeshA,
                                               const cuBQL::Triangle* dMeshB,
                                               const float2* triAMetrics,
                                               const float2* triBMetrics,
                                               tbb::concurrent_vector<int2>& finalExactPairs,
                                               IntersectionTimeTracker& tracker,
                                               cudaStream_t mainStream) {
  double tTotalStart = cuBQL::getCurrentTime();
  double tAllocStart = tTotalStart;

  if(totalBatches <= 0)
    return 0;

  uint64_t finalCandidatePairs = 0;
  batchMultiplier = std::min(batchMultiplier, totalBatches);

  // 1. Download Lookup Tables to Host ONCE
  std::vector<uint32_t> h_outPairsA(d_outPairsA.size());
  std::vector<uint32_t> h_outPairsB(d_outPairsB.size());
  std::vector<uint32_t> h_reverseMapB(d_reverseMapB.size());
  std::vector<uint32_t> h_outOffsetsB(d_outOffsetsB.size());

  CUBQL_CUDA_CALL(MemcpyAsync(h_outPairsA.data(), thrust::raw_pointer_cast(d_outPairsA.data()),
                              d_outPairsA.size() * sizeof(uint32_t), cudaMemcpyDeviceToHost, mainStream));
  CUBQL_CUDA_CALL(MemcpyAsync(h_outPairsB.data(), thrust::raw_pointer_cast(d_outPairsB.data()),
                              d_outPairsB.size() * sizeof(uint32_t), cudaMemcpyDeviceToHost, mainStream));
  CUBQL_CUDA_CALL(MemcpyAsync(h_reverseMapB.data(), thrust::raw_pointer_cast(d_reverseMapB.data()),
                              d_reverseMapB.size() * sizeof(uint32_t), cudaMemcpyDeviceToHost, mainStream));
  CUBQL_CUDA_CALL(MemcpyAsync(h_outOffsetsB.data(), thrust::raw_pointer_cast(d_outOffsetsB.data()),
                              d_outOffsetsB.size() * sizeof(uint32_t), cudaMemcpyDeviceToHost, mainStream));
  CUBQL_CUDA_CALL(StreamSynchronize(mainStream));

  // 2. Pre-Calculate Maximum Primitives Bounds
  int maxChunkPrims = 0;
  int numChunks = (totalBatches + batchMultiplier - 1) / batchMultiplier;
  uint32_t totalPrimsB = (uint32_t)d_outPrimsFlatB.size();

  for(int c = 0; c < numChunks; ++c) {
    int startBatch = c * batchMultiplier;
    int activeInChunk = std::min(batchMultiplier, totalBatches - startBatch);
    int chunkPrims = 0;
    for(int b = 0; b < activeInChunk; ++b) {
      uint32_t bIdx = h_outPairsB[startBatch + b];
      uint32_t markedIdxB = h_reverseMapB[bIdx];
      uint32_t startOffset = h_outOffsetsB[markedIdxB];
      uint32_t endOffset = (markedIdxB + 1 < h_outMarkedCountB) ? h_outOffsetsB[markedIdxB + 1] : totalPrimsB;
      chunkPrims += (int)(endOffset - startOffset);
    }
    maxChunkPrims = std::max(maxChunkPrims, chunkPrims);
  }

  if(maxChunkPrims == 0)
    return 0;

  // 3. Pre-Allocate Static Sandboxes
  uint32_t* d_BIter = nullptr;
  uint64_t* d_AIter = nullptr;
  int* d_pairCounts = nullptr;
  int* d_offsets = nullptr;
  int* d_chunkBatchOffsets = nullptr;

  CUBQL_CUDA_CALL(Malloc(&d_BIter, maxChunkPrims * sizeof(uint32_t)));
  CUBQL_CUDA_CALL(Malloc(&d_AIter, maxChunkPrims * sizeof(uint64_t)));
  CUBQL_CUDA_CALL(Malloc(&d_pairCounts, maxChunkPrims * sizeof(int)));
  CUBQL_CUDA_CALL(Malloc(&d_offsets, maxChunkPrims * sizeof(int)));
  CUBQL_CUDA_CALL(Malloc(&d_chunkBatchOffsets, batchMultiplier * sizeof(int)));

  // Pre-allocate 3 Microbatch Slots
  constexpr int NUM_SLOTS = 3;
  MicrobatchSlot slots[NUM_SLOTS];

  struct SlotEvents
  {
    cudaEvent_t evEvalStart, evEvalEnd;
    cudaEvent_t evD2HStart, evD2HEnd;
    bool recorded = false;
  };
  SlotEvents slotEvents[NUM_SLOTS];

  for(int s = 0; s < NUM_SLOTS; ++s) {
    slots[s].allocate();
    CUBQL_CUDA_CALL(EventCreate(&slotEvents[s].evEvalStart));
    CUBQL_CUDA_CALL(EventCreate(&slotEvents[s].evEvalEnd));
    CUBQL_CUDA_CALL(EventCreate(&slotEvents[s].evD2HStart));
    CUBQL_CUDA_CALL(EventCreate(&slotEvents[s].evD2HEnd));
  }

  // mainStream Broad Phase Events
  cudaEvent_t evAssemblyStart, evAssemblyEnd;
  cudaEvent_t evCountStart, evCountEnd;
  cudaEvent_t evFillStart, evFillEnd;
  CUBQL_CUDA_CALL(EventCreate(&evAssemblyStart));
  CUBQL_CUDA_CALL(EventCreate(&evAssemblyEnd));
  CUBQL_CUDA_CALL(EventCreate(&evCountStart));
  CUBQL_CUDA_CALL(EventCreate(&evCountEnd));
  CUBQL_CUDA_CALL(EventCreate(&evFillStart));
  CUBQL_CUDA_CALL(EventCreate(&evFillEnd));
  bool hasPendingFillEvent = false;

  tbb::task_group cpuTaskGroup;

  std::atomic<uint64_t> atomicCpuTimeUs{0};
  std::atomic<int> atomicGreenPairs{0};
  std::atomic<int> atomicYellowPairs{0};

  const uint32_t* ptr_outPairsA = thrust::raw_pointer_cast(d_outPairsA.data());
  const uint32_t* ptr_outPairsB = thrust::raw_pointer_cast(d_outPairsB.data());
  const uint32_t* ptr_reverseMapB = thrust::raw_pointer_cast(d_reverseMapB.data());
  const uint32_t* ptr_outOffsetsB = thrust::raw_pointer_cast(d_outOffsetsB.data());
  const uint32_t* ptr_outPrimsFlatB = thrust::raw_pointer_cast(d_outPrimsFlatB.data());

  std::vector<int> h_chunkBatchOffsets(batchMultiplier);

  int2* d_candidatePairs = nullptr;
  size_t candidateCapacity = 0;

  cudaEvent_t candidateReadyEvent;
  CUBQL_CUDA_CALL(EventCreateWithFlags(&candidateReadyEvent, cudaEventDisableTiming));

  tracker.preallocateTimeMs = (cuBQL::getCurrentTime() - tAllocStart) * 1000.0;

  // 4. Main Chunk Loop
  for(int i = 0; i < totalBatches; i += batchMultiplier) {
    tracker.numberOfBatchLoops++;
    int activeBatchesInChunk = std::min(batchMultiplier, totalBatches - i);

    int totalChunkPrims = 0;
    for(int b = 0; b < activeBatchesInChunk; ++b) {
      h_chunkBatchOffsets[b] = totalChunkPrims;
      uint32_t bIdx = h_outPairsB[i + b];
      uint32_t markedIdxB = h_reverseMapB[bIdx];
      uint32_t startOffset = h_outOffsetsB[markedIdxB];
      uint32_t endOffset = (markedIdxB + 1 < h_outMarkedCountB) ? h_outOffsetsB[markedIdxB + 1] : totalPrimsB;
      totalChunkPrims += (int)(endOffset - startOffset);
    }

    if(totalChunkPrims == 0)
      continue;

    // --- ASSEMBLY PHASE (mainStream) ---
    CUBQL_CUDA_CALL(EventRecord(evAssemblyStart, mainStream));

    CUBQL_CUDA_CALL(MemcpyAsync(d_chunkBatchOffsets, h_chunkBatchOffsets.data(), activeBatchesInChunk * sizeof(int),
                                cudaMemcpyHostToDevice, mainStream));
    int assembleGridSize = (activeBatchesInChunk + 63) / 64;
    assembleChunkBuffersByBatchKernel_V3<<<assembleGridSize, 64, 0, mainStream>>>(
        d_BIter, d_AIter, ptr_outPairsA, ptr_outPairsB, ptr_reverseMapB, ptr_outOffsetsB, ptr_outPrimsFlatB,
        h_outMarkedCountB, totalPrimsB, i, activeBatchesInChunk, d_chunkBatchOffsets);

    CUBQL_CUDA_CALL(EventRecord(evAssemblyEnd, mainStream));

    // --- OVERLAP COUNTING PHASE (mainStream) ---
    CUBQL_CUDA_CALL(EventRecord(evCountStart, mainStream));

    int kernelGridSize = (totalChunkPrims + 127) / 128;
    countAABBOverlapsKernel_Indirected_V3<<<kernelGridSize, 128, 0, mainStream>>>(
        d_pairCounts, bvhA, const_cast<cuBQL::Triangle*>(dMeshA), const_cast<cuBQL::Triangle*>(dMeshB), d_BIter, 0,
        totalChunkPrims, d_AIter);

    thrust::device_ptr<int> dev_counts(d_pairCounts);
    thrust::device_ptr<int> dev_offsets(d_offsets);
    thrust::exclusive_scan(thrust::cuda::par.on(mainStream), dev_counts, dev_counts + totalChunkPrims, dev_offsets);

    int lastCount = 0, lastOffset = 0;
    CUBQL_CUDA_CALL(
        MemcpyAsync(&lastCount, d_pairCounts + totalChunkPrims - 1, sizeof(int), cudaMemcpyDeviceToHost, mainStream));
    CUBQL_CUDA_CALL(
        MemcpyAsync(&lastOffset, d_offsets + totalChunkPrims - 1, sizeof(int), cudaMemcpyDeviceToHost, mainStream));

    CUBQL_CUDA_CALL(EventRecord(evCountEnd, mainStream));

    CUBQL_CUDA_CALL(StreamSynchronize(mainStream));

    float assemblyMs = 0.0f, countMs = 0.0f;
    CUBQL_CUDA_CALL(EventElapsedTime(&assemblyMs, evAssemblyStart, evAssemblyEnd));
    CUBQL_CUDA_CALL(EventElapsedTime(&countMs, evCountStart, evCountEnd));
    tracker.assemblyPhaseMs += assemblyMs;
    tracker.executionPhaseMs += countMs;

    if(hasPendingFillEvent) {
      float fillMs = 0.0f;
      CUBQL_CUDA_CALL(EventElapsedTime(&fillMs, evFillStart, evFillEnd));
      tracker.executionPhaseMs += fillMs;
      hasPendingFillEvent = false;
    }

    int totalChunkPairs = lastCount + lastOffset;
    if(totalChunkPairs == 0)
      continue;

    finalCandidatePairs += totalChunkPairs;

    if((size_t)totalChunkPairs > candidateCapacity) {
      if(d_candidatePairs)
        CUBQL_CUDA_CALL(Free(d_candidatePairs));
      candidateCapacity = std::max((size_t)totalChunkPairs, candidateCapacity * 2);
      CUBQL_CUDA_CALL(Malloc(&d_candidatePairs, candidateCapacity * sizeof(int2)));
    }

    // --- OVERLAP FILLING PHASE (mainStream) ---
    CUBQL_CUDA_CALL(EventRecord(evFillStart, mainStream));

    fillAABBOverlapsKernel_Indirected_V3<<<kernelGridSize, 128, 0, mainStream>>>(
        d_candidatePairs, d_offsets, bvhA, const_cast<cuBQL::Triangle*>(dMeshA), const_cast<cuBQL::Triangle*>(dMeshB),
        d_BIter, 0, totalChunkPrims, d_AIter);

    CUBQL_CUDA_CALL(EventRecord(evFillEnd, mainStream));
    hasPendingFillEvent = true;

    CUBQL_CUDA_CALL(EventRecord(candidateReadyEvent, mainStream));

    int numMicrobatches = (totalChunkPairs + MICROBATCH_SIZE - 1) / MICROBATCH_SIZE;

    for(int m = 0; m < numMicrobatches; ++m) {
      int uOffset = m * MICROBATCH_SIZE;
      int uSize = std::min(MICROBATCH_SIZE, totalChunkPairs - uOffset);

      int slotIdx = m % NUM_SLOTS;
      MicrobatchSlot& slot = slots[slotIdx];
      SlotEvents& sEv = slotEvents[slotIdx];

      while(slot.busy.load(std::memory_order_acquire)) {
        std::this_thread::yield();
      }

      if(sEv.recorded) {
        float evalMs = 0.0f, d2hMs = 0.0f;
        CUBQL_CUDA_CALL(EventElapsedTime(&evalMs, sEv.evEvalStart, sEv.evEvalEnd));
        CUBQL_CUDA_CALL(EventElapsedTime(&d2hMs, sEv.evD2HStart, sEv.evD2HEnd));
        tracker.fineEvaluationPhaseMs += evalMs;
        tracker.DownloadAndClean += d2hMs;
        sEv.recorded = false;
      }

      slot.busy.store(true, std::memory_order_release);

      CUBQL_CUDA_CALL(StreamWaitEvent(slot.stream, candidateReadyEvent, 0));

      // 1. Predicate Check & Compaction
      CUBQL_CUDA_CALL(EventRecord(sEv.evEvalStart, slot.stream));

      evaluateAndCompactPairsV2(d_candidatePairs + uOffset, slot.d_pairStatuses, dMeshA, dMeshB, triAMetrics,
                                triBMetrics, uSize, slot.stream);

      CUBQL_CUDA_CALL(MemsetAsync(slot.d_greenCount, 0, sizeof(int), slot.stream));
      CUBQL_CUDA_CALL(MemsetAsync(slot.d_yellowCount, 0, sizeof(int), slot.stream));

      int compactGridSize = (uSize + 255) / 256;
      compactGreenYellowKernel<<<compactGridSize, 256, 0, slot.stream>>>(
          d_candidatePairs + uOffset, slot.d_pairStatuses, uSize, slot.d_green, slot.d_greenCount, slot.d_yellow,
          slot.d_yellowCount);

      CUBQL_CUDA_CALL(EventRecord(sEv.evEvalEnd, slot.stream));

      // 2. Exact Stream Transfers
      CUBQL_CUDA_CALL(EventRecord(sEv.evD2HStart, slot.stream));

      // Transfer counts first
      CUBQL_CUDA_CALL(
          MemcpyAsync(slot.hGreenCountPinned, slot.d_greenCount, sizeof(int), cudaMemcpyDeviceToHost, slot.stream));
      CUBQL_CUDA_CALL(
          MemcpyAsync(slot.hYellowCountPinned, slot.d_yellowCount, sizeof(int), cudaMemcpyDeviceToHost, slot.stream));

      // Transfer pair data asynchronously
      CUBQL_CUDA_CALL(
          MemcpyAsync(slot.hGreenPinned, slot.d_green, uSize * sizeof(int2), cudaMemcpyDeviceToHost, slot.stream));
      CUBQL_CUDA_CALL(
          MemcpyAsync(slot.hYellowPinned, slot.d_yellow, uSize * sizeof(int2), cudaMemcpyDeviceToHost, slot.stream));

      CUBQL_CUDA_CALL(EventRecord(sEv.evD2HEnd, slot.stream));
      sEv.recorded = true;

      // 3. Offload to TBB with precise host wall-time measurement
      // Inside the microbatch loop:
      struct CallbackCtx
      {
        MicrobatchSlot* slot;
        Mesh* meshA;
        Mesh* meshB;
        tbb::concurrent_vector<int2>* out;
        tbb::task_group* tg;
        std::atomic<uint64_t>* cpuTimeUs;
        std::atomic<int>* greenPairs;
        std::atomic<int>* yellowPairs;
      };

      auto* ctx = new CallbackCtx{&slot,         &meshAcpu,        &meshBcpu,         &finalExactPairs,
                                  &cpuTaskGroup, &atomicCpuTimeUs, &atomicGreenPairs, &atomicYellowPairs};

      CUBQL_CUDA_CALL(LaunchHostFunc(
          slot.stream,
          [](void* userData) {
            auto* c = static_cast<CallbackCtx*>(userData);
            MicrobatchSlot* s = c->slot;

            // Extract values while inside stream order (guaranteed stream sync)
            int greenCount = *s->hGreenCountPinned;
            int yellowCount = *s->hYellowCountPinned;

            // Fast-path: If no items, unlock slot immediately & exit
            if(greenCount == 0 && yellowCount == 0) {
              s->busy.store(false, std::memory_order_release);
              delete c;
              return;
            }

            // Offload actual computational work to TBB
            c->tg->run([c, s, greenCount, yellowCount]() {
              if(greenCount > 0) {
                c->out->grow_by(s->hGreenPinned, s->hGreenPinned + greenCount);
                c->greenPairs->fetch_add(greenCount, std::memory_order_relaxed);
              }

              if(yellowCount > 0) {
                auto tbbStart = std::chrono::high_resolution_clock::now();

                filterYellowPairsTBB(*c->meshA, *c->meshB, s->hYellowPinned, (size_t)yellowCount, *c->out);

                auto tbbEnd = std::chrono::high_resolution_clock::now();
                uint64_t durUs = std::chrono::duration_cast<std::chrono::microseconds>(tbbEnd - tbbStart).count();
                c->cpuTimeUs->fetch_add(durUs, std::memory_order_relaxed);
                c->yellowPairs->fetch_add(yellowCount, std::memory_order_relaxed);
              }

              // CRITICAL FIX: Only mark slot free AFTER TBB finishes processing memory!
              s->busy.store(false, std::memory_order_release);
              delete c;
            });
          },
          ctx));
    }
  }

  // --- CLEANUP & FINAL SYNCHRONIZATION ---
  double tCleanupStart = cuBQL::getCurrentTime();

  // 1. Ensure the main stream is done with all assembly/overlap tasks
  CUBQL_CUDA_CALL(StreamSynchronize(mainStream));
  if(hasPendingFillEvent) {
    float fillMs = 0.0f;
    CUBQL_CUDA_CALL(EventElapsedTime(&fillMs, evFillStart, evFillEnd));
    tracker.executionPhaseMs += fillMs;
  }

  // 2. Wait for GPU microbatch streams to finish enqueuing/executing their work
  for(int s = 0; s < NUM_SLOTS; ++s) {
    CUBQL_CUDA_CALL(StreamSynchronize(slots[s].stream));
    if(slotEvents[s].recorded) {
      float evalMs = 0.0f, d2hMs = 0.0f;
      CUBQL_CUDA_CALL(EventElapsedTime(&evalMs, slotEvents[s].evEvalStart, slotEvents[s].evEvalEnd));
      CUBQL_CUDA_CALL(EventElapsedTime(&d2hMs, slotEvents[s].evD2HStart, slotEvents[s].evD2HEnd));
      tracker.fineEvaluationPhaseMs += evalMs;
      tracker.DownloadAndClean += d2hMs;
      slotEvents[s].recorded = false;
    }
  }

  // 3. Wait for ALL TBB tasks to finish processing host pinned buffers!
  // This is the critical gate that prevents segmentation faults.
  cpuTaskGroup.wait();

  // 4. Now it is 100% safe to destroy events and free the slot resources
  for(int s = 0; s < NUM_SLOTS; ++s) {
    CUBQL_CUDA_CALL(EventDestroy(slotEvents[s].evEvalStart));
    CUBQL_CUDA_CALL(EventDestroy(slotEvents[s].evEvalEnd));
    CUBQL_CUDA_CALL(EventDestroy(slotEvents[s].evD2HStart));
    CUBQL_CUDA_CALL(EventDestroy(slotEvents[s].evD2HEnd));
    slots[s].destroy(); // Free pinned memory buffers
  }

  // Store total CPU predicate evaluation time in milliseconds
  tracker.CPUPredicates = atomicCpuTimeUs.load() / 1000.0;
  tracker.confirmedGreenPairs = atomicGreenPairs.load();
  tracker.confirmedYellowPairs = atomicYellowPairs.load();

  CUBQL_CUDA_CALL(EventDestroy(candidateReadyEvent));
  CUBQL_CUDA_CALL(EventDestroy(evAssemblyStart));
  CUBQL_CUDA_CALL(EventDestroy(evAssemblyEnd));
  CUBQL_CUDA_CALL(EventDestroy(evCountStart));
  CUBQL_CUDA_CALL(EventDestroy(evCountEnd));
  CUBQL_CUDA_CALL(EventDestroy(evFillStart));
  CUBQL_CUDA_CALL(EventDestroy(evFillEnd));

  if(d_candidatePairs)
    CUBQL_CUDA_CALL(Free(d_candidatePairs));

  CUBQL_CUDA_CALL(Free(d_BIter));
  CUBQL_CUDA_CALL(Free(d_AIter));
  CUBQL_CUDA_CALL(Free(d_pairCounts));
  CUBQL_CUDA_CALL(Free(d_offsets));
  CUBQL_CUDA_CALL(Free(d_chunkBatchOffsets));

  tracker.cleanupTimeMs = (cuBQL::getCurrentTime() - tCleanupStart) * 1000.0;
  tracker.totalLoopTimeMs = (cuBQL::getCurrentTime() - tTotalStart) * 1000.0;

  return finalCandidatePairs;
}