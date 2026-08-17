#ifndef TRIANGLE_DOUBLE_H
#define TRIANGLE_DOUBLE_H

#include <cuda_runtime.h>
#include <vector_types.h>
#include <algorithm>

// Optional: Include cuBQL box types if you want built-in bounding box support
//#include "cuBQL/bvh.h" 

struct TriangleDouble {
    double3 a;
    double3 b;
    double3 c;

    // Default constructors for device & host initialization
    __host__ __device__ TriangleDouble() = default;

    __host__ __device__ TriangleDouble(double3 p0, double3 p1, double3 p2)
        : a(p0), b(p1), c(p2) {}

 
};

#endif // TRIANGLE_DOUBLE_H