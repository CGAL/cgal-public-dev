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
float func_naive(const float &a, const float &b) {
  return a + b;
}

void __attribute__((noinline))
func_naive(const std::vector<float> &A, const std::vector<float> &B, std::vector<float> &C) {
  for (std::size_t i = 0; i < C.capacity(); ++i)
    C[i] = func_naive(A[i], B[i]);
}

// Openmp-vectorized function
#pragma omp declare simd
#pragma omp declare simd linear(ref(a)) linear(ref(b))

float func_openmp(const float &a, const float &b) {
  return a + b;
}

void __attribute__((noinline))
func_openmp(const std::vector<float> &A, const std::vector<float> &B, std::vector<float> &C) {
  #pragma omp for simd
  for (std::size_t i = 0; i < C.capacity(); ++i)
    C[i] = func_openmp(A[i], B[i]);
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

  // Repeat test many times, average results

  std::vector<float> C_naive{};
  C_naive.reserve(N);
  std::vector<float> C_openmp{};
  C_openmp.reserve(N);

  timer_naive.start();
  func_naive(A, B, C_naive);
  timer_naive.stop();

  timer_openmp.start();
  func_openmp(A, B, C_openmp);
  timer_openmp.stop();

  // Report results
  std::cout << "Naive:   " << timer_naive.time() << "\n"
            << "Openmp:  " << timer_openmp.time();

}
