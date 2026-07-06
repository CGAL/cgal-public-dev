// // // --------------------------------------------------------------------
// //   // PRIMITIVE BOX CROSS INTERSECTION CHECK (BRUTEFORCE CHECK PROFILE)
// //   // --------------------------------------------------------------------
// //   std::cout << "==================================================" << std::endl;
// //   std::cout << "    EXECUTING PRIMITIVE CROSS-INTERSECTION CHECK  " << std::endl;
// //   std::cout << "==================================================" << std::endl;

// //   // Ensure any previous work on the stream is done before capturing start time
// //   CUBQL_CUDA_CALL(StreamSynchronize(s));
// //   double bruteForceStart = cuBQL::getCurrentTime();

// //   uint64_t totalPrimitiveIntersections = executePrimitiveBoxCrossIntersectionCheck<T, D>(
// //       batchMultiplier,
// //       d_intersectPairsA,
// //       d_intersectPairsB,
// //       outSortedPrimIDsA, outNodeOffsetsA, d_boxesA,
// //       outSortedPrimIDsB, outNodeOffsetsB, d_boxesB,
// //       s
// //   );
  
// //   // Synchronize to ensure the GPU finishes the check before capturing end time
// //   CUBQL_CUDA_CALL(StreamSynchronize(s));
// //   double bruteForceEnd = cuBQL::getCurrentTime();
// //   double bruteForceMs = (bruteForceEnd - bruteForceStart) * 1000.0;

// //   std::cout << " -> Bruteforce Check Execution Time : " << bruteForceMs << " ms" << std::endl;
// //   std::cout << " -> Total Primitive Overlaps Found  : " << totalPrimitiveIntersections << std::endl;
// //   std::cout << " -> Cell-Level Overlaps Baseline    : " << totalIntersections << std::endl;

// //   if (totalPrimitiveIntersections > 0) {
// //       std::cout << " -> [OK] Primitive intersection test produced active pairings." << std::endl;
// //   } else {
// //       std::cout << " -> [WARNING] Zero primitive-level overlaps found. Verify if dataset is actually disjoint!" << std::endl;
// //   }
// //   std::cout << "==================================================\n" << std::endl;

// template<typename T, int D>
// void testFastBVH(const cuBQL::box_t<T,D>* d_boxesA, int numPrimsA, int maxCellSizeA,
//                  const cuBQL::box_t<T,D>* d_boxesB, int numPrimsB, int maxCellSizeB,
//                  cuBQL::box_t<T,D> globalBoxA, cuBQL::box_t<T,D> globalBoxB, // Maintained signature
//                  cudaStream_t s, cuBQL::DeviceMemoryResource& memResource, int batchMultiplier, int userleafThreshold)
// {
//   std::cout << "\n==================================================" << std::endl;
//   std::cout << " [DEEP DEBUG] => Entering ZIMBLIFIED testFastBVH..." << std::endl;
//   std::cout << "==================================================" << std::endl;

//   // Check raw pointers and parameters right away
//   std::cout << " -> Device Pointer A : " << d_boxesA << " | Prims: " << numPrimsA << " | MaxCell: " << maxCellSizeA << std::endl;
//   std::cout << " -> Device Pointer B : " << d_boxesB << " | Prims: " << numPrimsB << " | MaxCell: " << maxCellSizeB << std::endl;
  
//   if (d_boxesA == d_boxesB) {
//       std::cout << " [CRITICAL WARNING] d_boxesA and d_boxesB point to the EXACT SAME memory address!" << std::endl;
//   } else {
//       std::cout << " [OK] Input pointers A and B are distinct memory blocks." << std::endl;
//   }

//   if (maxCellSizeA > 20 || maxCellSizeB > 20) {
//       std::cerr << "[CRITICAL ERROR] numIterations is too high! Aborting testFastBVH." << std::endl;
//       return;
//   }

//   // --------------------------------------------------------------------
//   // DYNAMIC GLOBAL BOX REDUCTION PHASE (FIXED)
//   // --------------------------------------------------------------------
//   std::cout << " [REDUCTION] => Launching computeGlobalBoxParallel for Mesh A and B..." << std::endl;
  
//   CUBQL_CUDA_CALL(StreamSynchronize(s));
//   double tReductionStart = cuBQL::getCurrentTime();

//   // Call the parallelized reduction function which safely implements its own internal tracking
//   cuBQL::box_t<T, D> calculatedGlobalBoxA = cuBQL::utils::computeGlobalBoxParallel<T, D>(d_boxesA, numPrimsA, s, memResource);
//   cuBQL::box_t<T, D> calculatedGlobalBoxB = cuBQL::utils::computeGlobalBoxParallel<T, D>(d_boxesB, numPrimsB, s, memResource);

//   CUBQL_CUDA_CALL(StreamSynchronize(s));
//   double tReductionMs = (cuBQL::getCurrentTime() - tReductionStart) * 1000.0;

//   // --- START ORIGINAL BENCHMARK TIMER ---
//   double tStart = cuBQL::getCurrentTime();

//   // --- MESH A RUN (Zimblified Edition) ---
//   cuBQL::box_t<T, D>* outNodeBoxesA = nullptr;
//   uint32_t* outSortedPrimIDsA = nullptr;
//   uint32_t* outNodeOffsetsA = nullptr;
//   uint32_t outTotalActiveCellsA = 0;

//   // Profile Mesh A Initialization
//   CUBQL_CUDA_CALL(StreamSynchronize(s));
//   double tInitAStart = cuBQL::getCurrentTime();

//   // Using calculatedGlobalBoxA here instead of the input parameter globalBoxA
//   cuBQL::gpuBuilder_v3::test_speedrun_initialization_linear(
//       d_boxesA, numPrimsA, (uint32_t)maxCellSizeA, globalBoxA,
//       outNodeBoxesA, outSortedPrimIDsA, outNodeOffsetsA, outTotalActiveCellsA,
//       s, memResource
//   );
//   CUBQL_CUDA_CALL(StreamSynchronize(s));
//   double tInitAMs = (cuBQL::getCurrentTime() - tInitAStart) * 1000.0;

//   // --- MESH B RUN (Zimblified Edition) ---
//   cuBQL::box_t<T, D>* outNodeBoxesB = nullptr;
//   uint32_t* outSortedPrimIDsB = nullptr;
//   uint32_t* outNodeOffsetsB = nullptr;
//   uint32_t outTotalActiveCellsB = 0;

//   // Profile Mesh B Initialization
//   CUBQL_CUDA_CALL(StreamSynchronize(s));
//   double tInitBStart = cuBQL::getCurrentTime();

//   // Using calculatedGlobalBoxB here instead of the input parameter globalBoxB
//   cuBQL::gpuBuilder_v3::test_speedrun_initialization_linear(
//       d_boxesB, numPrimsB, (uint32_t)maxCellSizeB, globalBoxB,
//       outNodeBoxesB, outSortedPrimIDsB, outNodeOffsetsB, outTotalActiveCellsB,
//       s, memResource
//   );
//   CUBQL_CUDA_CALL(StreamSynchronize(s));
//   double tInitBMs = (cuBQL::getCurrentTime() - tInitBStart) * 1000.0;

//   double tEnd = cuBQL::getCurrentTime();
//   double elapsedMs = (tEnd - tStart) * 1000.0;

//   // --------------------------------------------------------------------
//   // VERIFY BUILDING METRICS READBACK (No structural junk, 100% Active)
//   // --------------------------------------------------------------------
//   std::cout << "\n--------------------------------------------------" << std::endl;
//   std::cout << " [LINEAR INITIALIZATION METRICS OVERVIEW]" << std::endl;
//   std::cout << " -> Tree A: Total Compounded Active Cells = " << outTotalActiveCellsA << std::endl;
//   std::cout << " -> Tree B: Total Compounded Active Cells = " << outTotalActiveCellsB << std::endl;
//   std::cout << "--------------------------------------------------" << std::endl;

//   // Read back only the active, valid bounding boxes
//   std::vector<cuBQL::box_t<T,D>> h_boxesA(outTotalActiveCellsA);
//   std::vector<cuBQL::box_t<T,D>> h_boxesB(outTotalActiveCellsB);
//   std::vector<uint32_t> h_offsetsA(outTotalActiveCellsA + 1);
//   std::vector<uint32_t> h_offsetsB(outTotalActiveCellsB + 1);

//   CUBQL_CUDA_CALL(Memcpy(h_boxesA.data(), outNodeBoxesA, outTotalActiveCellsA * sizeof(cuBQL::box_t<T,D>), cudaMemcpyDeviceToHost));
//   CUBQL_CUDA_CALL(Memcpy(h_boxesB.data(), outNodeBoxesB, outTotalActiveCellsB * sizeof(cuBQL::box_t<T,D>), cudaMemcpyDeviceToHost));
//   CUBQL_CUDA_CALL(Memcpy(h_offsetsA.data(), outNodeOffsetsA, (outTotalActiveCellsA + 1) * sizeof(uint32_t), cudaMemcpyDeviceToHost));
//   CUBQL_CUDA_CALL(Memcpy(h_offsetsB.data(), outNodeOffsetsB, (outTotalActiveCellsB + 1) * sizeof(uint32_t), cudaMemcpyDeviceToHost));

//   // Calculate true primitive allocations via offsets
//   uint32_t sumPrimsA = h_offsetsA[outTotalActiveCellsA] - h_offsetsA[0];
//   uint32_t sumPrimsB = h_offsetsB[outTotalActiveCellsB] - h_offsetsB[0];

//   std::cout << "--------------------------------------------------" << std::endl;
//   std::cout << " [POPULATED CELLS CHECK SUMMARY]" << std::endl;
//   std::cout << " -> Tree A Active Nodes: " << outTotalActiveCellsA << " | Total Verified Prims = " << sumPrimsA << std::endl;
//   std::cout << " -> Tree B Active Nodes: " << outTotalActiveCellsB << " | Total Verified Prims = " << sumPrimsB << std::endl;
//   std::cout << "--------------------------------------------------\n" << std::endl;

//   // --------------------------------------------------------------------
//   // PARALLEL TEMP NODE CROSS CHECK
//   // --------------------------------------------------------------------
//   std::cout << "\n==================================================" << std::endl;
//   std::cout << "          LAUNCHING CROSS-CHECK KERNELS           " << std::endl;
//   std::cout << "==================================================" << std::endl;

//   thrust::device_vector<uint32_t> d_intersectPairsA;
//   thrust::device_vector<uint32_t> d_intersectPairsB;

//   double crossCheckStart = cuBQL::getCurrentTime();

//   uint32_t totalIntersections = executeBoxCrossCheck<T,D>(
//       outNodeBoxesA, outTotalActiveCellsA,
//       outNodeBoxesB, outTotalActiveCellsB,
//       d_intersectPairsA, d_intersectPairsB, s
//   );
//   CUBQL_CUDA_CALL(StreamSynchronize(s));

//   double crossCheckEnd = cuBQL::getCurrentTime();
//   double crossCheckMs = (crossCheckEnd - crossCheckStart) * 1000.0;

//   // Save historical metrics prior to post-build pruning actions
//   uint32_t initialCellsA = outTotalActiveCellsA;
//   uint32_t initialPrimsA = sumPrimsA;
//   uint32_t initialCellsB = outTotalActiveCellsB;
//   uint32_t initialPrimsB = sumPrimsB;

//   int currentPrimsNumA = (int)sumPrimsA;
//   int currentPrimsNumB = (int)sumPrimsB;

//   // --------------------------------------------------------------------
//   // PARALLEL STREAM COMPACTION & PRUNING ALGORITHM
//   // --------------------------------------------------------------------
//   std::cout << " [PIPELINE PRUNING] => Dispatched grid streaming parallel re-index compaction (POST-BUILD)..." << std::endl;

//   CUBQL_CUDA_CALL(StreamSynchronize(s));
//   double pruneStartA = cuBQL::getCurrentTime();
//   parallelPruneAndReindexAll(
//       thrust::raw_pointer_cast(d_intersectPairsA.data()), totalIntersections,
//       outSortedPrimIDsA, outNodeOffsetsA, outTotalActiveCellsA, currentPrimsNumA,
//       s, memResource
//   );
//   CUBQL_CUDA_CALL(StreamSynchronize(s));
//   double pruneMsA = (cuBQL::getCurrentTime() - pruneStartA) * 1000.0;

//   double pruneStartB = cuBQL::getCurrentTime();
//   parallelPruneAndReindexAll(
//       thrust::raw_pointer_cast(d_intersectPairsB.data()), totalIntersections,
//       outSortedPrimIDsB, outNodeOffsetsB, outTotalActiveCellsB, currentPrimsNumB,
//       s, memResource
//   );
//   CUBQL_CUDA_CALL(StreamSynchronize(s));
//   double pruneMsB = (cuBQL::getCurrentTime() - pruneStartB) * 1000.0;

//   // --------------------------------------------------------------------
//   // GPU BUILDER V4 FOREST EXPANSION
//   // --------------------------------------------------------------------
//   std::cout << " [FOREST EXPANSION] => Launching Level-by-Level Parallel Sub-Tree Compilation..." << std::endl;

//   const uint32_t numInitA = outTotalActiveCellsA + (outTotalActiveCellsA & 1u);
//   const uint32_t maxNodesA = 2u * (uint32_t)numPrimsA + numInitA + 2u;

//   const uint32_t numInitB = outTotalActiveCellsB + (outTotalActiveCellsB & 1u);
//   const uint32_t maxNodesB = 2u * (uint32_t)numPrimsB + numInitB + 2u;

//   thrust::device_vector<uint32_t> d_nodeDescendantCountsA(maxNodesA, 0);
//   thrust::device_vector<uint32_t> d_nodeDescendantCountsB(maxNodesB, 0);

//   cuBQL::BuildConfig buildConfig; 
//   buildConfig.makeLeafThreshold = userleafThreshold;

//   // --- Compile Forest A ---
//   cuBQL::BinaryBVH<T, D> bvhA;
//   CUBQL_CUDA_CALL(StreamSynchronize(s));
//   double forestAStart = cuBQL::getCurrentTime();

//   cuBQL::gpuBuilder_v4::build_forest<T, D>(
//       bvhA, d_boxesA, numPrimsA, outTotalActiveCellsA,
//       outSortedPrimIDsA, outNodeOffsetsA, buildConfig, thrust::raw_pointer_cast(d_nodeDescendantCountsA.data()), s, memResource
//   );
//   CUBQL_CUDA_CALL(StreamSynchronize(s));
//   double forestAMs = (cuBQL::getCurrentTime() - forestAStart) * 1000.0;

//   // --- Compile Forest B ---
//   cuBQL::BinaryBVH<T, D> bvhB;
//   CUBQL_CUDA_CALL(StreamSynchronize(s));
//   double forestBStart = cuBQL::getCurrentTime();

//   cuBQL::gpuBuilder_v4::build_forest<T, D>(
//       bvhB, d_boxesB, numPrimsB, outTotalActiveCellsB,
//       outSortedPrimIDsB, outNodeOffsetsB, buildConfig, thrust::raw_pointer_cast(d_nodeDescendantCountsB.data()), s, memResource
//   );
//   CUBQL_CUDA_CALL(StreamSynchronize(s));
//   double forestBMs = (cuBQL::getCurrentTime() - forestBStart) * 1000.0;

//   // --------------------------------------------------------------------
//   // FINAL PERFORMANCE & RESULTS READOUT
//   // --------------------------------------------------------------------
//   uint64_t totalValidActivePairs = (uint64_t)initialCellsA * initialCellsB;
//   double overlapPercentage = 0.0;
//   if (totalValidActivePairs > 0) {
//       overlapPercentage = (double)totalIntersections / (double)totalValidActivePairs * 100.0;
//   }

//   std::cout << "\n==================================================" << std::endl;
//   std::cout << "               FINAL BENCHMARK SUMMARY            " << std::endl;
//   std::cout << "==================================================" << std::endl;
//   std::cout << " -> Dynamic Global Box Reduction     : " << tReductionMs << " ms" << std::endl;
//   std::cout << " -> Mesh A Speedrun Init Speed      : " << tInitAMs << " ms" << std::endl;
//   std::cout << " -> Mesh B Speedrun Init Speed      : " << tInitBMs << " ms" << std::endl;
//   std::cout << " -> Cell Cross-Check Execution      : " << crossCheckMs << " ms" << std::endl;
//   std::cout << " ------------------------------------------------" << std::endl;
//   std::cout << " -> Mesh A Forest BVH Expansion     : " << forestAMs << " ms | Nodes: " << bvhA.numNodes << std::endl;
//   std::cout << " -> Mesh B Forest BVH Expansion     : " << forestBMs << " ms | Nodes: " << bvhB.numNodes << std::endl;
//   std::cout << " ------------------------------------------------" << std::endl;
//   std::cout << " -> Mesh A Structural Pruning (Post)  : " << pruneMsA << " ms" << std::endl;
//   std::cout << "    * Active Cells Surviving: " << initialCellsA << " -> " << outTotalActiveCellsA << " (Dropped: " << (initialCellsA - outTotalActiveCellsA) << ")" << std::endl;
//   std::cout << "    * Pack Prim Elements Size: " << initialPrimsA << " -> " << currentPrimsNumA << std::endl;
//   std::cout << " ------------------------------------------------" << std::endl;
//   std::cout << " -> Mesh B Structural Pruning (Post)  : " << pruneMsB << " ms" << std::endl;
//   std::cout << "    * Active Cells Surviving: " << initialCellsB << " -> " << outTotalActiveCellsB << " (Dropped: " << (initialCellsB - outTotalActiveCellsB) << ")" << std::endl;
//   std::cout << "    * Pack Prim Elements Size: " << initialPrimsB << " -> " << currentPrimsNumB << std::endl;
//   std::cout << " ------------------------------------------------" << std::endl;
//   std::cout << " -> Intersecting Pairs Detected     : " << totalIntersections << " pairs" << std::endl;
//   std::cout << " -> True Active Matrix Evaluated    : " << totalValidActivePairs << " pairings" << std::endl;
//   std::cout << " -> Percentage of Total Overlaps    : " << overlapPercentage << "%" << std::endl;
//   std::cout << "==================================================\n" << std::endl;

//   // Memory Cleanup Phase
//   if (outNodeBoxesA)     _FREE(outNodeBoxesA, s, memResource);
//   if (outSortedPrimIDsA) _FREE(outSortedPrimIDsA, s, memResource);
//   if (outNodeOffsetsA)   _FREE(outNodeOffsetsA, s, memResource);
//   if (outNodeBoxesB)     _FREE(outNodeBoxesB, s, memResource);
//   if (outSortedPrimIDsB) _FREE(outSortedPrimIDsB, s, memResource);
//   if (outNodeOffsetsB)   _FREE(outNodeOffsetsB, s, memResource);
// }