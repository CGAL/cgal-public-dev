#ifndef IntersectionTimeTracker_H
#define IntersectionTimeTracker_H

#include <iostream>

struct IntersectionTimeTracker {
    double totalLoopTimeMs       = 0.0; 
    double preallocateTimeMs     = 0.0; 
    double assemblyPhaseMs       = 0.0; 
    double executionPhaseMs      = 0.0; 
    double fineEvaluationPhaseMs = 0.0; // GPU Float predicates + compaction
    double gpuDoublePredicatesMs = 0.0; // GPU Double predicates pass on Yellow list
    double cleanupTimeMs         = 0.0; 
    double DownloadAndClean      = 0.0;
    double CPUPredicates         = 0.0; // CPU CGAL Exact Predicates (Orange List)
    
    int numberOfBatchLoops  = 0;
    int confirmedGreenPairs  = 0;
    int confirmedYellowPairs = 0; // Float yellow candidates
    int confirmedOrangePairs = 0; // Remaining unresolved pairs sent to CPU Exact

    void print() const {
        std::cout << "\n==================================================\n";
        std::cout << "          CROSS INTERSECTION TIME PROFILING       \n";
        std::cout << "==================================================\n";
        std::cout << " Preallocation Space Discovery : " << preallocateTimeMs     << " ms\n";
        std::cout << " Sandbox Assembly & Host Reads : " << assemblyPhaseMs       << " ms\n";
        std::cout << " CUDA Kernels & Scan Compute   : " << executionPhaseMs      << " ms\n";
        std::cout << " Fine GPU Float Evaluation     : " << fineEvaluationPhaseMs << " ms\n";
        std::cout << " Fine GPU Double Evaluation    : " << gpuDoublePredicatesMs << " ms\n";
        std::cout << " Upload and clean up           : " << DownloadAndClean      << " ms\n";
        std::cout << " Sandbox Cleanup & Free Cycles : " << cleanupTimeMs         << " ms\n";
        std::cout << " CPU Exact Predicates (CGAL)   : " << CPUPredicates         << " ms\n";
        std::cout << "--------------------------------------------------\n";
        std::cout << " Number of loops in batch      : " << numberOfBatchLoops    << "\n";
        std::cout << " Confirmed Green Pairs         : " << confirmedGreenPairs   << "\n";
        std::cout << " Confirmed Float Yellow Pairs  : " << confirmedYellowPairs  << "\n";
        std::cout << " Unresolved Orange Pairs (CPU) : " << confirmedOrangePairs  << "\n";
        std::cout << "--------------------------------------------------\n";
        std::cout << " Total Tracked Pipeline Time   : " << totalLoopTimeMs       << " ms\n";
        std::cout << "==================================================\n\n";
    }
};

#endif // IntersectionTimeTracker_H