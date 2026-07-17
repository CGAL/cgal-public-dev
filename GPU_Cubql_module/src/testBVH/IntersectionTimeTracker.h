#ifndef IntersectionTimeTracker_H
#define IntersectionTimeTracker_H

#include <iostream>

struct IntersectionTimeTracker {
    double totalLoopTimeMs      = 0.0; 
    double preallocateTimeMs    = 0.0; 
    double assemblyPhaseMs      = 0.0; 
    double executionPhaseMs     = 0.0; 
    double fineEvaluationPhaseMs = 0.0; // Added tracker for the exact geometric predicates + compaction
    double cleanupTimeMs        = 0.0; 
    double DownloadAndClean = 0.0;
    double CPUPredicates = 0.0;
    int numberOfBatchLoops = 0;
    int confirmedGreenPairs = 0;
    int confirmedYellowPairs = 0;

    void print() const {
        std::cout << "\n==================================================\n";
        std::cout << "          CROSS INTERSECTION TIME PROFILING       \n";
        std::cout << "==================================================\n";
        std::cout << " Preallocation Space Discovery : " << preallocateTimeMs << " ms\n";
        std::cout << " Sandbox Assembly & Host Reads : " << assemblyPhaseMs   << " ms\n";
        std::cout << " CUDA Kernels & Scan Compute   : " << executionPhaseMs  << " ms\n";
        std::cout << " Fine Geometric Evaluation     : " << fineEvaluationPhaseMs << " ms\n";
        std::cout << " Upload and clean up           : " << DownloadAndClean << " ms\n";
        std::cout << " Sandbox Cleanup & Free Cycles : " << cleanupTimeMs     << " ms\n";
        std::cout << "Number of loops in batch: " << numberOfBatchLoops  << "\n";
        std::cout << "--------------------------------------------------\n";
        std::cout << " Total Tracked Pipeline Time   : " << totalLoopTimeMs   << " ms\n";
        std::cout << "==================================================\n\n";
    }
};


#endif 