#ifndef TEST_SUITE_H
#define TEST_SUITE_H

#include <vector>
#include <set>
#include <string>
#include <cuda_runtime.h>
#include <vector_types.h>

#include "CPU/CgalDefinitions.h"
#include "GPUIntersector/KernelBVHController.h"

// Forward declaration of ApplicationState from main context
struct ApplicationState;

struct TestConfig {
    int numSteps = 8;
    double maxTranslation = 1.0;     // Maximum offset in scene units
    double maxRotationDeg = 180.0;   // Rotation range [-180, 180] degrees
    unsigned int seed = 1337;        // Random seed for reproducibility
    bool testMainGpuPipeline = true;
    bool testStandalonePipeline = true;
    int enableGpuPrecision = 1;      // Integer precision mode
    int queryDescentLevel = 12;
    int referenceDescentLevel = 12;
    int batchMultiplier = 2147483647;
    int dualTreeSteps = 8;
    int leafThreshold = 20;
};

struct VerificationResult {
    int stepIndex = 0;
    double3 rotA{0,0,0}, transA{0,0,0};
    double3 rotB{0,0,0}, transB{0,0,0};

    size_t groundTruthCount = 0;
    size_t gpuMainCount = 0;
    size_t standaloneCount = 0;

    bool gpuMainExactMatch = false;
    bool standaloneExactMatch = false;

    size_t mainFalsePositives = 0;
    size_t mainFalseNegatives = 0;
    size_t standaloneFalsePositives = 0;
    size_t standaloneFalseNegatives = 0;

    double cgalTimeMs = 0.0;
    double gpuMainTimeMs = 0.0;
    double standaloneTimeMs = 0.0;
};

class TestSuite {
public:
    explicit TestSuite(ApplicationState& appState);

    // Runs full automated sweep across translation and rotation ranges
    void runSuite(const TestConfig& config);

    // Evaluates a single iteration with explicit mesh transforms
    VerificationResult evaluateStep(
        int stepIdx,
        const double3& rotA, const double3& transA,
        const double3& rotB, const double3& transB,
        const TestConfig& config
    );

private:
    ApplicationState& app_;

    std::set<std::pair<int, int>> convertToCanonicalSet(const std::vector<int2>& pairs);

    std::vector<int2> runCGALClassicGroundTruth(
        const double3& rotA, const double3& transA,
        const double3& rotB, const double3& transB,
        double& outTimeMs
    );

    std::vector<int2> runMainGPUPipeline(
        const double3& rotA, const double3& transA,
        const double3& rotB, const double3& transB,
        const TestConfig& config,
        double& outTimeMs
    );

    std::vector<int2> runStandalonePipeline(
        const double3& rotA, const double3& transA,
        const double3& rotB, const double3& transB,
        const TestConfig& config,
        double& outTimeMs
    );
};

#endif // TEST_SUITE_H