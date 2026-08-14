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

// --- TBB INCLUDES ----
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/concurrent_vector.h>

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

struct GPUTimingBreakdown {
    double uploadTime = 0.0;
    double executionTime = 0.0; 
    double downloadTime = 0.0;
};

// External assembly mapping linked against updated kernels.cu
// Added extra bool parameter to control GPU workload characteristics directly
extern "C" long long runMeshIntersection(
    const cuBQL::Triangle* hA, int NA, 
    const cuBQL::Triangle* hB, int NB, 
    std::vector<IntersectionPair>& hPairs,
    GPUTimingBreakdown& outTimings,
    bool doMoreWork); 

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
    
    // Hybrid cuBQL + CGAL Metrics
    double hybrid_read_time = 0;          
    double hybrid_upload_time = 0;        
    double hybrid_bvh_gpu_time = 0;       
    double hybrid_download_time = 0;      
    double hybrid_conv_time = 0;          // Step A: CGAL Mesh -> cuBQL Conversion Time
    double hybrid_filter_time = 0;        // Steps C + D: Combined Extraction & Filtering Time
    double hybrid_total_routine_time = 0; 
    
    double hybrid_aabb_hits = 0;
    double hybrid_overlaps = 0;
};

// Helper to convert a CGAL face to a cuBQL Triangle format safely
cuBQL::Triangle convertCgalFaceToCuBQL(const Mesh& mesh, face_descriptor face) {
    auto conn = mesh.vertices_around_face(mesh.halfedge(face));
    auto it = conn.begin();
    
    Point3 pA = mesh.point(*it++);
    Point3 pB = mesh.point(*it++);
    Point3 pC = mesh.point(*it);

    cuBQL::Triangle tri;
    tri.a = { (float)pA.x(), (float)pA.y(), (float)pA.z() };
    tri.b = { (float)pB.x(), (float)pB.y(), (float)pB.z() };
    tri.c = { (float)pC.x(), (float)pC.y(), (float)pC.z() };
    return tri;
}

// Executive execution runner processing a singular file pairing through both pipelines
CompleteMetrics run_comparative_test(const std::string& fileA, const std::string& fileB, bool doMoreWork) {
    CompleteMetrics metric;
    
    // ==========================================
    // 1. UNIFIED DISK READ (Using CGAL Engine)
    // ==========================================
    double read_start = cuBQL::getCurrentTime();
    Mesh m1, m2;
    std::ifstream in1(fileA);
    std::ifstream in2(fileB);
    if (!(in1 >> m1) || !(in2 >> m2)) {
        return metric; // Return invalid if reading fails
    }
    in1.close();
    in2.close();
    double total_read_time = cuBQL::getCurrentTime() - read_start;

    metric.cgal_read_time = total_read_time;
    metric.hybrid_read_time = total_read_time;

    metric.triA_count = (double)num_faces(m1);
    metric.triB_count = (double)num_faces(m2);

    // ==========================================
    // 2. PURE CGAL METHOD Execution
    // ==========================================
    double cgal_routine_start = cuBQL::getCurrentTime();
    
    CGAL::Real_timer cgal_timer;
    cgal_timer.start();
    std::vector<std::pair<face_descriptor, face_descriptor>> cgal_intersected_tris;
    
    CGAL::Polygon_mesh_processing::internal::compute_face_face_intersection(
        m1, m2, 
        std::back_inserter(cgal_intersected_tris),
        CGAL::parameters::all_default(), 
        CGAL::parameters::all_default()
    );

    cgal_timer.stop();
    
    metric.cgal_query_time = cgal_timer.time();
    metric.cgal_overlaps = (double)cgal_intersected_tris.size();
    metric.cgal_total_routine_time = (cuBQL::getCurrentTime() - cgal_routine_start) + metric.cgal_read_time;

    // ==========================================
    // 3. HYBRID CUBQL + CGAL METHOD 
    // ==========================================
    double hybrid_routine_start = cuBQL::getCurrentTime();
    
    // ------------------------------------------
    // Step A: Extract / Convert CGAL Mesh into cuBQL Triangles
    // ------------------------------------------
    double conv_t0 = cuBQL::getCurrentTime();
    
    std::vector<cuBQL::Triangle> hA; 
    hA.reserve(num_faces(m1));
    std::vector<face_descriptor> facesA; 
    facesA.reserve(num_faces(m1));
    for (auto f : m1.faces()) {
        hA.push_back(convertCgalFaceToCuBQL(m1, f));
        facesA.push_back(f);
    }

    std::vector<cuBQL::Triangle> hB; 
    hB.reserve(num_faces(m2));
    std::vector<face_descriptor> facesB; 
    facesB.reserve(num_faces(m2));
    for (auto f : m2.faces()) {
        hB.push_back(convertCgalFaceToCuBQL(m2, f));
        facesB.push_back(f);
    }
    
    double conv_t1 = cuBQL::getCurrentTime();
    metric.hybrid_conv_time = conv_t1 - conv_t0; // Isolated Timer for Data Conversion
    
    int NA = hA.size();
    int NB = hB.size();

    // ------------------------------------------
    // Step B: Execute Broad-Phase Intersection via GPU BVH
    // ------------------------------------------
    std::vector<IntersectionPair> hPairs;
    GPUTimingBreakdown gpuTimings;
    
    // Propagated control flag directly into the backend GPU hook
    long long totalAABBHits = runMeshIntersection(hA.data(), NA, hB.data(), NB, hPairs, gpuTimings, doMoreWork);
    
    metric.hybrid_upload_time   = gpuTimings.uploadTime;
    metric.hybrid_bvh_gpu_time  = gpuTimings.executionTime;
    metric.hybrid_download_time = gpuTimings.downloadTime;
    metric.hybrid_aabb_hits     = (double)totalAABBHits;

    // ------------------------------------------
    // Combined Steps C & D: Narrow-Phase Extraction & Intersection Test On-The-Fly
    // ------------------------------------------
    double filter_start = cuBQL::getCurrentTime();
    
    auto& coords1 = m1.points();
    auto& coords2 = m2.points();

    // Use a concurrent vector to allow safe parallel push_backs
    tbb::concurrent_vector<IntersectionPair> exactPairs_concurrent;

    tbb::parallel_for(tbb::blocked_range<size_t>(0, hPairs.size()), [&](const tbb::blocked_range<size_t>& r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
            const auto& pair = hPairs[i];
            
            face_descriptor fA = facesA[pair.idA];
            face_descriptor fB = facesB[pair.idB];

            auto halfedgeA = m1.halfedge(fA);
            auto vA0 = m1.target(halfedgeA);
            auto vA1 = m1.target(m1.next(halfedgeA));
            auto vA2 = m1.target(m1.next(m1.next(halfedgeA)));
            Triangle3 triA(coords1[vA0], coords1[vA1], coords1[vA2]);

            auto halfedgeB = m2.halfedge(fB);
            auto vB0 = m2.target(halfedgeB);
            auto vB1 = m2.target(m2.next(halfedgeB));
            auto vB2 = m2.target(m2.next(m2.next(halfedgeB)));
            Triangle3 triB(coords2[vB0], coords2[vB1], coords2[vB2]);

            if (CGAL::do_intersect(triA, triB)) {
                exactPairs_concurrent.push_back(pair);
            }
        }
    });

    //  std::vector<IntersectionPair> exactPairs(exactPairs_concurrent.begin(), exactPairs_concurrent.end());
    double filter_end = cuBQL::getCurrentTime();
    
    metric.hybrid_filter_time = filter_end - filter_start;
    // Read the size directly from the concurrent vector (this is thread-safe and fast)
    metric.hybrid_overlaps    = (double)exactPairs_concurrent.size();
    
    metric.hybrid_total_routine_time = (cuBQL::getCurrentTime() - hybrid_routine_start) + metric.hybrid_read_time;
    
    metric.valid = true;
    return metric;
   
}

int main(int argc, char** argv) {
    if (argc < 7) { // Bumped from 6 to 7
        std::cerr << "Usage: " << argv[0] << " <basePath/> <output.csv> <K_categories> <pairCount> <speedUp (1/0)> <numThreads>\n";
        std::cerr << "Example: " << argv[0] << " /path/Data/ results.csv 11 20 1 30\n";
        return 1;
    }

    std::string basePath = argv[1];
    std::string outCsvName = argv[2];
    int maxCategories = std::stoi(argv[3]);
    int pairCount = std::stoi(argv[4]);
    bool doMoreWork = (std::stoi(argv[5]) != 0);
    int numThreads = std::stoi(argv[6]); // Parse thread count

    // Limit TBB's global thread pool to your specific allocation (e.g., 30)
    tbb::global_control global_limit(
        tbb::global_control::max_allowed_parallelism, numThreads
    );

    std::cout << "======================================================================\n";
    std::cout << "Starting Automated Batch Evaluation Task (Unified & Streamlined)\n";
    std::cout << " GPU Extra Work Parameter : " << (doMoreWork ? "TRUE" : "FALSE") << "\n";
    std::cout << " TBB Thread Worker Limit  : " << numThreads << "\n";
    std::cout << "======================================================================\n";
    
    // ... rest of your main loop remains exactly the same ...
    if (!basePath.empty() && basePath.back() != '/') {
        basePath += '/';
    }

    std::ofstream csv(outCsvName);
    if (!csv.is_open()) {
        std::cerr << "Critical Error: Could not generate file output at " << outCsvName << "\n";
        return 1;
    }
    
    csv << "Category,"
        << "Avg_TriA,Avg_TriB,"
        << "CGAL_Avg_Overlaps,CGAL_Avg_ReadTime_S,CGAL_Avg_QueryTime_S,CGAL_Avg_TotalRoutineTime_S,"
        << "Hybrid_Avg_AABBHits,Hybrid_Avg_ExactOverlaps,"
        << "Hybrid_Avg_ReadTime_S,Hybrid_Avg_UploadTime_S,Hybrid_Avg_GPU_BVH_Time_S,Hybrid_Avg_DownloadTime_S,Hybrid_Avg_ConversionTime_S,Hybrid_Avg_FilteringTime_S,Hybrid_Avg_TotalRoutineTime_S\n";

    std::cout << "======================================================================\n";
    std::cout << "Starting Automated Batch Evaluation Task (Unified & Streamlined)\n";
    std::cout << " GPU Extra Work Parameter : " << (doMoreWork ? "TRUE" : "FALSE") << "\n";
    std::cout << "======================================================================\n";

    for (int cat = 0; cat < maxCategories; ++cat) {
        std::cout << "Benchmarking Category " << cat << "..." << std::endl;
        
        CompleteMetrics catSum; 
        int validPairs = 0;

        for (int i = 0; i < pairCount; ++i) {
            std::string fileA = basePath + std::to_string(cat) + "/A_" + std::to_string(i) + ".off";
            std::string fileB = basePath + std::to_string(cat) + "/B_" + std::to_string(i) + ".off";

            CompleteMetrics res = run_comparative_test(fileA, fileB, doMoreWork);
            
            if (res.valid) {
                if (res.cgal_overlaps != res.hybrid_overlaps) {
                    std::cerr << "\n######################################################################\n"
                              << " CRITICAL ERROR: PIPELINE INTERSECTION MISMATCH DETECTED!\n"
                              << " Category : " << cat << " | Pair Index : " << i << "\n"
                              << " CGAL Exact Intersections: " << (long long)res.cgal_overlaps << "\n"
                              << " Hybrid Exact Intersections: " << (long long)res.hybrid_overlaps << "\n"
                              << "######################################################################\n\n";
                } else {
                    std::cout << "   -> Pair " << i << " Verified: CGAL (" << (long long)res.cgal_overlaps 
                              << ") == Hybrid (" << (long long)res.hybrid_overlaps << ")\n";
                }

                catSum.triA_count               += res.triA_count;
                catSum.triB_count               += res.triB_count;
                catSum.cgal_overlaps            += res.cgal_overlaps;
                catSum.cgal_read_time           += res.cgal_read_time;
                catSum.cgal_query_time          += res.cgal_query_time;
                catSum.cgal_total_routine_time  += res.cgal_total_routine_time;
                
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
                      << "   |-- CGAL Routine Total:   " << (catSum.cgal_total_routine_time / validPairs) << "s\n"
                      << "   |-- Hybrid Routine Total: " << (catSum.hybrid_total_routine_time / validPairs) << "s\n";
        }
    }

    csv.close();
    return 0;
}