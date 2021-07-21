
#include <iostream>

int main() {
  std::size_t N = 100'000;
  std::size_t R = 1'000;

  // All test data will be randomly generated

  // Shapes used for testing will be confined to a cubic region

  // Generate rays to cast

  // Generate boxes to cast

  // Generate boxes to hit

  // Generate primitives (triangles) to hit

  // Benchmark R times, so that tests are interleaved
  for (int r = 0; r < R; ++r) {

    // Time bbox-bbox intersection

    // Time ray-bbox intersection

    // Time ray-primitive intersection

  }

  // Divide times to produce averages

  // Display results

}