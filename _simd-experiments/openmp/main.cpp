#include <random>
#include <chrono>
#include <iostream>
#include <algorithm>

#include <CGAL/Real_timer.h>

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

// Naive function
float func_naive(float a, float b) {
    return std::max(a, b);
}

// Openmp-vectorized function
#pragma omp declare simd

float func_openmp(float a, float b) {
    return std::max(a, b);
}

// Test harness
int main() {

    std::size_t N = 1000000000;

    // Generate test data
    std::random_device rd;
    std::mt19937 re{rd()};
    std::uniform_real_distribution<float> distribution(0, 10);
    std::vector<float> A(N), B(N);
    std::generate(A.begin(), A.end(), [&]() { return distribution(re); });
    std::generate(B.begin(), B.end(), [&]() { return distribution(re); });

    // Use CGAL's timing utility
    CGAL::Real_timer timer_naive, timer_openmp;

    // Repeat test many times, average times

    std::vector<float> C_naive(N);
    timer_naive.start();
    for (std::size_t i = 0; i < N; ++i)
        C_naive[i] = func_naive(A[i], B[i]);
    timer_naive.stop();

    std::vector<float> C_openmp(N);
    timer_openmp.start();
#pragma omp simd
    for (std::size_t i = 0; i < N; ++i)
        C_openmp[i] = func_openmp(A[i], B[i]);
    timer_openmp.stop();

    // Report results
    std::cout << "Naive: " << timer_naive.time() << "\n"
              << "Openmp:  " << timer_openmp.time();

}
