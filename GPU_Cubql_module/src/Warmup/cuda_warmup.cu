// cuda_warmup.cu
#include "cuda_warmup.h"
#include <cuda_runtime.h>
#include <thrust/device_vector.h>

extern "C" void warmupCUDA() {
    // 1. Force CUDA context creation & driver initialization
    cudaFree(0);

    // 2. Warm up Thrust allocator and internal kernel templates
    thrust::device_vector<int> dummy(1024);

    // 3. Ensure all setup operations finish and clocks settle
    cudaDeviceSynchronize();
}