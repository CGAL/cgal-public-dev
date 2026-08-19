#ifndef APPLICATION_STATE_H
#define APPLICATION_STATE_H

#include <vector>
#include <memory>
#include <cuda_runtime.h>
#include <vector_types.h>
#include <tbb/global_control.h>

#include "GPUIntersector/KernelBVHController.h"
#include "CPU/CgalDefinitions.h"
#include "CPU/MeshTriangleDegeneracyVisualizer.h"
#include "CPU/ExecutionTimingVisualizer.h"

struct ApplicationState
{
  KernelBVHController controller;
  ExecutionStats stats;
  Mesh meshA, meshB;
  std::vector<double3> hVertsA, hVertsB;
  std::vector<uint3> hIndicesA, hIndicesB;
  std::unique_ptr<tbb::global_control> tbbControl;
  MeshTriangleDegeneracyVisualizer edgeVisualizer;
  ExecutionTimingVisualizer timingVisualizer;

  bool isLoaded = false;
  double currentScaleFactor = 1.0;
  double currentCenterX = 0.0, currentCenterY = 0.0, currentCenterZ = 0.0;
  double currentMaxSpan = 0.0;
  bool currentScaledToUnit = false;

  double3 origCenterA{0.0, 0.0, 0.0}, origCenterB{0.0, 0.0, 0.0};
  float3 normCenterA{0.0f, 0.0f, 0.0f}, normCenterB{0.0f, 0.0f, 0.0f};
};

#endif // APPLICATION_STATE_H