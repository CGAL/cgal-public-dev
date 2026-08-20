#pragma once

#include <vector>
#include <vector_types.h>

// -------------------------------------------------------------
// POD Structs shared across NVCC and host C++ callers
// -------------------------------------------------------------
struct ExactPoint3POD {
  double x, y, z;
};

struct IntersectionSegment {
  ExactPoint3POD p0;
  ExactPoint3POD p1;
};

// Forward declaration of mesh types to prevent including CGAL in pure header calls
#include "CPU/CgalDefinitions.h"

// -------------------------------------------------------------
// Host API Interface
// -------------------------------------------------------------
std::vector<IntersectionSegment> computeExactIntersectionCurves(
    const Mesh& meshA,
    const Mesh& meshB,
    const std::vector<int2>& candidatePairs,
    const Point3& centerA, const double3& rotDegA, const double3& transA,
    const Point3& centerB, const double3& rotDegB, const double3& transB);