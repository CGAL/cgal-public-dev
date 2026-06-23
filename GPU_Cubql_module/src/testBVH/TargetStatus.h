#ifndef TARGET_STATUS_H
#define TARGET_STATUS_H

#include <thrust/device_vector.h>

struct IsTargetPairStatus {
    int target;
    __host__ __device__ bool operator()(const int status) const {
        return status == target;
    }
};


#endif // BATCHED_CROSS_INTERSECTION_H