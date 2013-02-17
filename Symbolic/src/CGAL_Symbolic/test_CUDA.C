
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <cuda_runtime.h>
#include <math.h>

int main( int argc, char** argv) {

    int deviceCount;
    
    if(cudaGetDeviceCount(&deviceCount) != cudaSuccess || 
            deviceCount == 0 ) {
        printf("There is no device supporting CUDA");
        return 0;
    }

    int dev;
    for (dev = 0; dev < deviceCount; ++dev) {
        cudaDeviceProp deviceProp;
        
        if(cudaGetDeviceProperties(&deviceProp, dev) != cudaSuccess) {
            printf("There is no device supporting CUDA");
            return 0;
        }

        if (dev == 0) {
            if (deviceProp.major < 1) {
                printf("There is no device supporting CUDA.");
                return 0;
            } else if (deviceCount == 1)
                printf("There is 1 device supporting CUDA");
            else
                printf("There are %d devices supporting CUDA", deviceCount);
        }


        return ((deviceProp.major * 10) + deviceProp.minor);

//         printf("Device %d: %s", dev, deviceProp.name);
//         printf("  Major revision number:                         %d",
//                deviceProp.major);
//         printf("  Minor revision number:                         %d",
//                deviceProp.minor);
//         printf("  Total amount of global memory:                 %d bytes",
//                deviceProp.totalGlobalMem);
//         printf("  Total amount of constant memory:               %d bytes",
//                deviceProp.totalConstMem); 
//         printf("  Total amount of shared memory per block:       %d bytes",
//                deviceProp.sharedMemPerBlock);
//         printf("  Total number of registers available per block: %d",
//                deviceProp.regsPerBlock);
//         printf("  Warp size:                                     %d",
//                deviceProp.warpSize);
//         printf("  Maximum number of threads per block:           %d",
//                deviceProp.maxThreadsPerBlock);
//         printf("  Maximum sizes of each dimension of a block:    %d x %d x %d",
//                deviceProp.maxThreadsDim[0],
//                deviceProp.maxThreadsDim[1],
//                deviceProp.maxThreadsDim[2]);
//         printf("  Maximum sizes of each dimension of a grid:     %d x %d x %d",
//                deviceProp.maxGridSize[0],
//                deviceProp.maxGridSize[1],
//                deviceProp.maxGridSize[2]);
//         printf("  Maximum memory pitch:                          %d bytes",
//                deviceProp.memPitch);
//         printf("  Texture alignment:                             %d bytes",
//                deviceProp.textureAlignment);
//         printf("  Clock rate:                                    %d kilohertz",
//                deviceProp.clockRate);
    }
}
