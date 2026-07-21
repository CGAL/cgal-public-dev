// cuda_warmup.h
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// Function to initialize CUDA context, ramp up GPU clocks, and warm up Thrust
void warmupCUDA();

#ifdef __cplusplus
}
#endif