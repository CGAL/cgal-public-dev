#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>

// --- cuBQL INCLUDES ---
#include "include/loadOFFCGAL.h" 
#include "include/loadOFF.h"

// --- CGAL INCLUDES ---
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/intersections.h>
#include <CGAL/Real_timer.h>

// Types for Pure CGAL Pipeline
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;

// Types for Hybrid cuBQL + CGAL Pipeline
typedef Kernel::Point_3                             Point3;
typedef Kernel::Triangle_3                          Triangle3;

struct IntersectionPair {
    int idA;
    int idB;
};

struct CGALTrianglePair {
    Triangle3 triA;
    Triangle3 triB;
    IntersectionPair originalPair; 
};

// Mirror structural declaration from your kernels.cu file
struct GPUTimingBreakdown {
    double uploadTime = 0.0;
    double executionTime = 0.0; 
    double downloadTime = 0.0;
};

// External assembly mapping linked against updated kernels.cu
extern "C" long long runMeshIntersection(
    const cuBQL::Triangle* hA, int NA, 
    const cuBQL::Triangle* hB, int NB, 
    std::vector<IntersectionPair>& hPairs,
    GPUTimingBreakdown& outTimings); // <--- Updated to match your multi-out requirements

// Helper to convert cuBQL Triangle vectors to CGAL geometry objects
Triangle3 convertToCGAL(const cuBQL::Triangle& t) {
    return Triangle3(
        Point3(t.a.x, t.a.y, t.a.z),
        Point3(t.b.x, t.b.y, t.b.z),
        Point3(t.c.x, t.c.y, t.c.z)
    );
}

// Structure to capture fully granular metrics for both runs
struct CompleteMetrics {
    bool valid = false;
    double triA_count = 0;
    double triB_count = 0;
    
    // Pure CGAL Granular Metrics
    double cgal_read_time = 0;
    double cgal_query_time = 0;      
    double cgal_total_routine_time = 0; 
    double cgal_overlaps = 0;
    
    // Hybrid cuBQL + CGAL Metrics (Now tracking 5 custom metrics + Total)
    double hybrid_read_time = 0;          // Metric 1: Disk -> Host CPU
    double hybrid_upload_time = 0;        // Metric 2: Host -> Device Allocation & Copy
    double hybrid_bvh_gpu_time = 0;       // Metric 3: Pure Device GPU Kernels
    double hybrid_download_time = 0;      // Metric 4: Device -> Host Output Copy
    double hybrid_conv_time = 0;          // Metric 5: Part A Conversion
    double hybrid_filter_time = 0;        // Metric 6: Part B Filtering
    double hybrid_total_routine_time = 0; // Total Summary Time
    
    double hybrid_aabb_hits = 0;
    double hybrid_overlaps = 0;
};

// Executive execution runner processing a singular file pairing through both pipelines
CompleteMetrics run_comparative_test(const std::string& fileA, const std::string& fileB) {
    CompleteMetrics metric;
    
    // ==========================================
    // 1. PURE CGAL METHOD (Granular Timing)
    // ==========================================
    double cgal_routine_start = cuBQL::getCurrentTime();
    
    double cgal_read_start = cuBQL::getCurrentTime();
    Mesh m1, m2;
    std::ifstream in1(fileA);
    std::ifstream in2(fileB);
    if (!(in1 >> m1) || !(in2 >> m2)) {
        return metric; // Return invalid if reading fails
    }
    in1.close();
    in2.close();
    metric.cgal_read_time = cuBQL::getCurrentTime() - cgal_read_start;

    metric.triA_count = (double)num_faces(m1);
    metric.triB_count = (double)num_faces(m2);

    CGAL::Real_timer cgal_timer;
    cgal_timer.start();
    std::vector<std::pair<face_descriptor, face_descriptor>> cgal_intersected_tris;
    
    //Note: Implicit BVH AABB-tree construction happens right inside this execution routine
    CGAL::Polygon_mesh_processing::internal::compute_face_face_intersection(
        m1, m2, 
        std::back_inserter(cgal_intersected_tris),
        CGAL::parameters::all_default(), 
        CGAL::parameters::all_default()
    );

    cgal_timer.stop();
    
    metric.cgal_query_time = cgal_timer.time();
    metric.cgal_overlaps = (double)cgal_intersected_tris.size();
    metric.cgal_total_routine_time = cuBQL::getCurrentTime() - cgal_routine_start;

    // ==========================================
    // 2. HYBRID CUBQL + CGAL METHOD
    // ==========================================
    double hybrid_routine_start = cuBQL::getCurrentTime();
    
    // 1. Read Time
    double hybrid_read_start = cuBQL::getCurrentTime();
    std::vector<cuBQL::Triangle> hA = cuBQL::samples::loadOFF(fileA);
    std::vector<cuBQL::Triangle> hB = cuBQL::samples::loadOFF(fileB);
    metric.hybrid_read_time = cuBQL::getCurrentTime() - hybrid_read_start;
    
    int NA = hA.size();
    int NB = hB.size();

    // 2, 3, 4. GPU Engine Extraction Steps (Upload, Compute, Download)
    std::vector<IntersectionPair> hPairs;
    GPUTimingBreakdown gpuTimings;
    long long totalAABBHits = runMeshIntersection(hA.data(), NA, hB.data(), NB, hPairs, gpuTimings);
    
    metric.hybrid_upload_time   = gpuTimings.uploadTime;
    metric.hybrid_bvh_gpu_time  = gpuTimings.executionTime;
    metric.hybrid_download_time = gpuTimings.downloadTime;
    metric.hybrid_aabb_hits     = (double)totalAABBHits;

    // 5. Phase A: Conversion Loop
    double conv_t0 = cuBQL::getCurrentTime();
    std::vector<CGALTrianglePair> cgalPairs;
    cgalPairs.reserve(hPairs.size());
    for (const auto& pair : hPairs) {
        cgalPairs.push_back({
            convertToCGAL(hA[pair.idA]),
            convertToCGAL(hB[pair.idB]),
            pair
        });
    }
    double conv_t1 = cuBQL::getCurrentTime();
    metric.hybrid_conv_time = conv_t1 - conv_t0;

    // 6. Phase B: Exact Intersections Filtering Loop
    double inter_t0 = cuBQL::getCurrentTime();
    std::vector<IntersectionPair> exactPairs;
    exactPairs.reserve(hPairs.size());
    for (const auto& cgalPair : cgalPairs) {
        if (CGAL::do_intersect(cgalPair.triA, cgalPair.triB)) {
            exactPairs.push_back(cgalPair.originalPair);
        }
    }
    double inter_t1 = cuBQL::getCurrentTime();
    
    metric.hybrid_filter_time = inter_t1 - inter_t0;
    metric.hybrid_overlaps = (double)exactPairs.size();
    
    // 7. Complete Macro Routine Time
    metric.hybrid_total_routine_time = cuBQL::getCurrentTime() - hybrid_routine_start;
    
    metric.valid = true;
    return metric;
}

int main(int argc, char** argv) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <basePath/> <output.csv> <K_categories> <pairCount>\n";
        std::cerr << "Example: " << argv[0] << " /path/to/Data/BigData/ results.csv 11 20\n";
        return 1;
    }

    std::string basePath = argv[1];
    std::string outCsvName = argv[2];
    int maxCategories = std::stoi(argv[3]);
    int pairCount = std::stoi(argv[4]);

    if (!basePath.empty() && basePath.back() != '/') {
        basePath += '/';
    }

    std::ofstream csv(outCsvName);
    if (!csv.is_open()) {
        std::cerr << "Critical Error: Could not generate file output at " << outCsvName << "\n";
        return 1;
    }
    
    // Updated data schema headers to reflect the newly partitioned timeline profile
    csv << "Category,"
        << "Avg_TriA,Avg_TriB,"
        << "CGAL_Avg_Overlaps,CGAL_Avg_ReadTime_S,CGAL_Avg_QueryTime_S,CGAL_Avg_TotalRoutineTime_S,"
        << "Hybrid_Avg_AABBHits,Hybrid_Avg_ExactOverlaps,"
        << "Hybrid_Avg_ReadTime_S,Hybrid_Avg_UploadTime_S,Hybrid_Avg_GPU_BVH_Time_S,Hybrid_Avg_DownloadTime_S,Hybrid_Avg_ConversionTime_S,Hybrid_Avg_FilteringTime_S,Hybrid_Avg_TotalRoutineTime_S\n";

    std::cout << "======================================================================\n";
    std::cout << "Starting Automated Batch Evaluation Task\n";
    std::cout << " Target Directory : " << basePath << "\n";
    std::cout << " Category Limit   : [0 to " << maxCategories - 1 << "]\n";
    std::cout << " File Pairs/Cat   : " << pairCount << "\n";
    std::cout << "======================================================================\n";

    for (int cat = 0; cat < maxCategories; ++cat) {
        std::cout << "Benchmarking Category " << cat << "..." << std::endl;
        
        CompleteMetrics catSum; 
        int validPairs = 0;

        for (int i = 0; i < pairCount; ++i) {
            std::string fileA = basePath + std::to_string(cat) + "/A_" + std::to_string(i) + ".off";
            std::string fileB = basePath + std::to_string(cat) + "/B_" + std::to_string(i) + ".off";

            CompleteMetrics res = run_comparative_test(fileA, fileB);
            
            if (res.valid) {
                // Per-iteration validation logic
                if (res.cgal_overlaps != res.hybrid_overlaps) {
                    std::cerr << "\n######################################################################\n"
                              << " CRITICAL ERROR: PIPELINE INTERSECTION MISMATCH DETECTED!\n"
                              << " Category : " << cat << " | Pair Index : " << i << "\n"
                              << " File A   : " << fileA << "\n"
                              << " File B   : " << fileB << "\n"
                              << " CGAL Exact Intersections: " << (long long)res.cgal_overlaps << "\n"
                              << " Hybrid Exact Intersections: " << (long long)res.hybrid_overlaps << "\n"
                              << " Delta    : " << std::abs((long long)res.cgal_overlaps - (long long)res.hybrid_overlaps) << "\n"
                              << "######################################################################\n\n";
                } else {
                    std::cout << "   -> Pair " << i << " Verified: CGAL (" << (long long)res.cgal_overlaps 
                              << ") == Hybrid (" << (long long)res.hybrid_overlaps << ")\n";
                }

                catSum.triA_count               += res.triA_count;
                catSum.triB_count               += res.triB_count;
                
                // CGAL accumulation
                catSum.cgal_overlaps            += res.cgal_overlaps;
                catSum.cgal_read_time           += res.cgal_read_time;
                catSum.cgal_query_time          += res.cgal_query_time;
                catSum.cgal_total_routine_time  += res.cgal_total_routine_time;
                
                // Hybrid accumulation
                catSum.hybrid_aabb_hits          += res.hybrid_aabb_hits;
                catSum.hybrid_overlaps           += res.hybrid_overlaps;
                catSum.hybrid_read_time          += res.hybrid_read_time;
                catSum.hybrid_upload_time        += res.hybrid_upload_time;
                catSum.hybrid_bvh_gpu_time       += res.hybrid_bvh_gpu_time;
                catSum.hybrid_download_time      += res.hybrid_download_time;
                catSum.hybrid_conv_time          += res.hybrid_conv_time;
                catSum.hybrid_filter_time        += res.hybrid_filter_time;
                catSum.hybrid_total_routine_time += res.hybrid_total_routine_time;
                
                validPairs++;
            }
        }

        if (validPairs > 0) {
            csv << cat << ","
                << catSum.triA_count / validPairs << ","
                << catSum.triB_count / validPairs << ","
                << catSum.cgal_overlaps / validPairs << ","
                << catSum.cgal_read_time / validPairs << ","
                << catSum.cgal_query_time / validPairs << ","
                << catSum.cgal_total_routine_time / validPairs << ","
                << catSum.hybrid_aabb_hits / validPairs << ","
                << catSum.hybrid_overlaps / validPairs << ","
                << catSum.hybrid_read_time / validPairs << ","
                << catSum.hybrid_upload_time / validPairs << ","
                << catSum.hybrid_bvh_gpu_time / validPairs << ","
                << catSum.hybrid_download_time / validPairs << ","
                << catSum.hybrid_conv_time / validPairs << ","
                << catSum.hybrid_filter_time / validPairs << ","
                << catSum.hybrid_total_routine_time / validPairs << "\n";
                
            std::cout << " > Category " << cat << " complete. (" << validPairs << " valid pairs)\n"
                      << "   |-- CGAL Routine Total:   " << (catSum.cgal_total_routine_time / validPairs) << "s  (Query: " << (catSum.cgal_query_time / validPairs) << "s)\n"
                      << "   |-- Hybrid Routine Total: " << (catSum.hybrid_total_routine_time / validPairs) << "s  (Pure GPU: " << (catSum.hybrid_bvh_gpu_time / validPairs) << "s)\n";
        } else {
            std::cout << " > Category " << cat << " skipped (No valid files found).\n";
        }
    }

    csv.close();
    std::cout << "\nUnified benchmark completed. Output saved to '" << outCsvName << "'\n";
    return 0;
}