#include "TestSuite.h"
#include "CPU/RotationTools.h"

#include "CPU/PolyscopeBridge.h"

#include "GPUIntersector/StandaloneBVHPipeline.h"

#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>



// Assuming ApplicationState is accessible via include context in main.cpp or its header
#include "ApplicationState.h"

TestSuite::TestSuite(ApplicationState& appState) : app_(appState) {}

std::set<std::pair<int, int>> TestSuite::convertToCanonicalSet(const std::vector<int2>& pairs) {
    std::set<std::pair<int, int>> canonicalSet;
    for (const auto& p : pairs) {
        canonicalSet.insert({p.x, p.y});
    }
    return canonicalSet;
}

std::vector<int2> TestSuite::runCGALClassicGroundTruth(
    const double3& rotA, const double3& transA,
    const double3& rotB, const double3& transB,
    double& outTimeMs) 
{
    Mesh meshA_transformed, meshB_transformed;
    Point3 centerA(app_.normCenterA.x, app_.normCenterA.y, app_.normCenterA.z);
    Point3 centerB(app_.normCenterB.x, app_.normCenterB.y, app_.normCenterB.z);

    transformCgalMesh(app_.meshA, meshA_transformed, centerA, rotA, transA);
    transformCgalMesh(app_.meshB, meshB_transformed, centerB, rotB, transB);

    using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;
    std::vector<std::pair<face_descriptor, face_descriptor>> cgal_intersected_tris;

    auto tStart = std::chrono::high_resolution_clock::now();
    CGAL::Polygon_mesh_processing::internal::compute_face_face_intersection(
        meshA_transformed, meshB_transformed,
        std::back_inserter(cgal_intersected_tris),
        CGAL::parameters::all_default(),
        CGAL::parameters::all_default()
    );
    auto tEnd = std::chrono::high_resolution_clock::now();
    outTimeMs = std::chrono::duration<double, std::milli>(tEnd - tStart).count();

    std::vector<int2> results;
    results.reserve(cgal_intersected_tris.size());
    for (const auto& pair : cgal_intersected_tris) {
        results.push_back(make_int2(static_cast<int>(pair.first), static_cast<int>(pair.second)));
    }
    return results;
}

std::vector<int2> TestSuite::runMainGPUPipeline(
    const double3& rotA, const double3& transA,
    const double3& rotB, const double3& transB,
    const TestConfig& config,
    double& outTimeMs) 
{
    if (!app_.controller.isGPUAllocated()) {
        app_.controller.reconstructGPU(app_.stats);
    }

    float tGPU, tCPU, tTV, tAT, tGB;
    app_.controller.setTransformBoth(rotA, transA, rotB, transB, tGPU, tCPU, tTV, tAT, tGB);

    int2* outPairs = nullptr;
    size_t outCount = 0;

    auto tStart = std::chrono::high_resolution_clock::now();
    app_.controller.runIntersectionPipeline(
        config.batchMultiplier, config.dualTreeSteps, 0, outPairs, outCount, 
        app_.stats, config.enableGpuPrecision // Passed as integer mode
    );
    auto tEnd = std::chrono::high_resolution_clock::now();
    outTimeMs = std::chrono::duration<double, std::milli>(tEnd - tStart).count();

    std::vector<int2> results;
    if (outPairs && outCount > 0) {
        results.assign(outPairs, outPairs + outCount);
        std::free(outPairs);
    }
    return results;
}

std::vector<int2> TestSuite::runStandalonePipeline(
    const double3& rotA, const double3& transA,
    const double3& rotB, const double3& transB,
    const TestConfig& config,
    double& outTimeMs) 
{
    Point3 centerA(app_.normCenterA.x, app_.normCenterA.y, app_.normCenterA.z);
    Point3 centerB(app_.normCenterB.x, app_.normCenterB.y, app_.normCenterB.z);

    ExecutionStats testStats;
    int2* outPairs = nullptr;
    size_t outCount = 0;

    auto tStart = std::chrono::high_resolution_clock::now();
    runGridForestIntersectionPipeline(
        app_.meshA, app_.meshB, app_.hVertsA.data(), static_cast<int>(app_.hVertsA.size()), 
        app_.hIndicesA.data(), static_cast<int>(app_.hIndicesA.size()), config.queryDescentLevel, 
        app_.hVertsB.data(), static_cast<int>(app_.hVertsB.size()), app_.hIndicesB.data(), 
        static_cast<int>(app_.hIndicesB.size()), config.referenceDescentLevel, config.batchMultiplier, 
        config.dualTreeSteps, config.enableGpuPrecision, config.leafThreshold, testStats, 
        outPairs, outCount, centerA, centerB, rotA, transA, rotB, transB
    );
    auto tEnd = std::chrono::high_resolution_clock::now();
    outTimeMs = std::chrono::duration<double, std::milli>(tEnd - tStart).count();

    std::vector<int2> results;
    if (outPairs && outCount > 0) {
        results.assign(outPairs, outPairs + outCount);
        std::free(outPairs);
    }
    return results;
}

VerificationResult TestSuite::evaluateStep(
    int stepIdx,
    const double3& rotA, const double3& transA,
    const double3& rotB, const double3& transB,
    const TestConfig& config) 
{
    VerificationResult res;
    res.stepIndex = stepIdx;
    res.rotA = rotA; res.transA = transA;
    res.rotB = rotB; res.transB = transB;

    auto gtPairs = runCGALClassicGroundTruth(rotA, transA, rotB, transB, res.cgalTimeMs);
    auto gtSet = convertToCanonicalSet(gtPairs);
    res.groundTruthCount = gtSet.size();

    if (config.testMainGpuPipeline) {
        auto mainPairs = runMainGPUPipeline(rotA, transA, rotB, transB, config, res.gpuMainTimeMs);
        auto mainSet = convertToCanonicalSet(mainPairs);
        res.gpuMainCount = mainSet.size();
        res.gpuMainExactMatch = (gtSet == mainSet);

        for (const auto& pair : mainSet) {
            if (gtSet.find(pair) == gtSet.end()) res.mainFalsePositives++;
        }
        for (const auto& pair : gtSet) {
            if (mainSet.find(pair) == mainSet.end()) res.mainFalseNegatives++;
        }
    }

    if (config.testStandalonePipeline) {
        auto standalonePairs = runStandalonePipeline(rotA, transA, rotB, transB, config, res.standaloneTimeMs);
        auto standaloneSet = convertToCanonicalSet(standalonePairs);
        res.standaloneCount = standaloneSet.size();
        res.standaloneExactMatch = (gtSet == standaloneSet);

        for (const auto& pair : standaloneSet) {
            if (gtSet.find(pair) == gtSet.end()) res.standaloneFalsePositives++;
        }
        for (const auto& pair : gtSet) {
            if (standaloneSet.find(pair) == standaloneSet.end()) res.standaloneFalseNegatives++;
        }
    }

    return res;
}

void TestSuite::runSuite(const TestConfig& config) {
    if (!app_.isLoaded) {
        std::cout << "Error: No meshes loaded. Run 'load' before invoking TestSuite.\n";
        return;
    }

    // 1. Fetch active viewport base positions from Polyscope
    float3 baseRotA_f{0.0f, 0.0f, 0.0f}, baseTransA_f{0.0f, 0.0f, 0.0f};
    float3 baseRotB_f{0.0f, 0.0f, 0.0f}, baseTransB_f{0.0f, 0.0f, 0.0f};

    if (!PolyscopeBridge::getCurrentTransforms(baseRotA_f, baseTransA_f, baseRotB_f, baseTransB_f)) {
        std::cout << "[TestSuite] Warning: Could not fetch Polyscope transforms. Falling back to origin.\n";
    }

    double3 baseRotA   = make_double3(baseRotA_f.x, baseRotA_f.y, baseRotA_f.z);
    double3 baseTransA = make_double3(baseTransA_f.x, baseTransA_f.y, baseTransA_f.z);
    double3 baseRotB   = make_double3(baseRotB_f.x, baseRotB_f.y, baseRotB_f.z);
    double3 baseTransB = make_double3(baseTransB_f.x, baseTransB_f.y, baseTransB_f.z);

    std::cout << "\n=========================================================================\n";
    std::cout << "               STARTING AUTOMATED PIPELINE TEST SUITE                    \n";
    std::cout << "=========================================================================\n";
    std::cout << "  Base Pos Mesh A : Trans(" << baseTransA.x << ", " << baseTransA.y << ", " << baseTransA.z << ")\n";
    std::cout << "  Base Pos Mesh B : Trans(" << baseTransB.x << ", " << baseTransB.y << ", " << baseTransB.z << ")\n";
    std::cout << "  Steps: " << config.numSteps 
              << " | Trans Perturb: [-" << config.maxTranslation << ", +" << config.maxTranslation 
              << "] | Rot Perturb: [-" << config.maxRotationDeg << "°, +" << config.maxRotationDeg << "°]"
              << " | Seed: " << config.seed << "\n\n";

    size_t totalMainPassed = 0, totalStandalonePassed = 0;

    std::mt19937 rng(config.seed);
    const double maxRotRad = config.maxRotationDeg * (M_PI / 180.0);

    std::uniform_real_distribution<double> rotDist(-maxRotRad, maxRotRad);
    std::uniform_real_distribution<double> transDist(-config.maxTranslation, config.maxTranslation);

    for (int i = 0; i < config.numSteps; ++i) {
        // Mesh A remains centered at its base transform
        double3 rotA   = baseRotA;
        double3 transA = baseTransA;

        // Apply random perturbations relative to Mesh B's current starting position
        double3 rotB = make_double3(
            baseRotB.x + rotDist(rng),
            baseRotB.y + rotDist(rng),
            baseRotB.z + rotDist(rng)
        );
        double3 transB = make_double3(
            baseTransB.x + transDist(rng),
            baseTransB.y + transDist(rng),
            baseTransB.z + transDist(rng)
        );

        VerificationResult res = evaluateStep(i, rotA, transA, rotB, transB, config);

        if (res.gpuMainExactMatch) totalMainPassed++;
        if (res.standaloneExactMatch) totalStandalonePassed++;

        std::cout << "Step [" << std::setw(2) << i + 1 << "/" << config.numSteps << "] "
                  << "GT: " << std::setw(4) << res.groundTruthCount << " tris | "
                  << "Main GPU: " << std::setw(4) << res.gpuMainCount 
                  << " (" << (res.gpuMainExactMatch ? "MATCH" : "MISMATCH") << ", FP:" 
                  << res.mainFalsePositives << " FN:" << res.mainFalseNegatives << ") " 
                  << std::fixed << std::setprecision(1) << res.gpuMainTimeMs << "ms | "
                  << "Standalone: " << std::setw(4) << res.standaloneCount 
                  << " (" << (res.standaloneExactMatch ? "MATCH" : "MISMATCH") << ") " 
                  << res.standaloneTimeMs << "ms | "
                  << "CGAL: " << res.cgalTimeMs << "ms\n";
    }

    std::cout << "-------------------------------------------------------------------------\n";
    std::cout << "  Main GPU Pipeline Accuracy      : " << totalMainPassed << "/" << config.numSteps 
              << " (" << (100.0 * totalMainPassed / config.numSteps) << "%)\n";
    std::cout << "  Standalone BVH Pipeline Accuracy : " << totalStandalonePassed << "/" << config.numSteps 
              << " (" << (100.0 * totalStandalonePassed / config.numSteps) << "%)\n";
    std::cout << "=========================================================================\n\n";
}