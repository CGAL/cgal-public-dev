
#include "../util.h"

#include "vbbox.h"

#include "intersection_strategies/xsimd.h"
#include "intersection_strategies/implicit.h"

#include <iostream>
#include <numeric>

template<typename T>
struct VQuery {
  Ray<T> ray;
  VBBox<T> vbox;
};

template<typename T>
auto pack_queries(std::vector<Query<T>> queries) {
  std::vector<VQuery<T>> vqueries;

  for (const Query<T> &query : queries)
    vqueries.push_back({query.ray, VBBox<T>(query.boxes)});

  return vqueries;
}

int main() {


  long N = 3742217;
  long R = 100;

  // Load test data
  auto file = std::ifstream("../data/remeshing_intersections_3742217.txt");
  if (!file.is_open()) return EXIT_FAILURE;
  auto queries = load_queries<double>(file, N);
  std::cout << "Loaded "
            << std::accumulate(queries.begin(), queries.end(), 0,
                               [](const auto &a, const Query<double> &b) { return a + b.boxes.size(); })
            << " scenarios "
            << "divided into " << queries.size() << " queries." << std::endl;

  // Convert test data to vector format
  auto vqueries = pack_queries(queries);

  std::vector<double> explicit_times, implicit_times;

  for (int i = 0; i < R; ++i) {
    std::cout << i + 1 << "/" << R << std::endl;

    std::vector<bool> explicit_results, implicit_results;

    for (const auto &query : vqueries) {

      explicit_times.push_back(time([&] {
        xsimd::intersect(query.vbox, query.ray, explicit_results);
      }));

      implicit_times.push_back(time([&] {
        implicit::intersect(query.vbox, query.ray, implicit_results);
      }));

    }

    // Check results for correctness
    if (implicit_results != explicit_results)
      throw std::logic_error("Incorrect results");

    std::cout << "\tHit-rate: "
              << std::accumulate(implicit_results.begin(), implicit_results.end(), 0.0) * 100.0 /
                 (double) implicit_results.size()
              << "%" << std::endl;
  }

  std::cout << std::endl;
  std::cout << "{| class=\"wikitable\"" << std::endl;
  std::cout << "|+ Time to Complete " << N << " Intersection Tests" << std::endl;
  std::cout << "|-" << std::endl;
  std::cout << "! Implementation !! Time" << std::endl;
  std::cout << "|-" << std::endl;
  std::cout << "| Explicit SIMD || "
            << std::accumulate(explicit_times.begin(), explicit_times.end(), 0.0) / (double) R
            << " ms"
            << std::endl;
  std::cout << "|-" << std::endl;
  std::cout << "| Implicit SIMD || "
            << std::accumulate(implicit_times.begin(), implicit_times.end(), 0.0) / (double) R
            << " ms"
            << std::endl;
  std::cout << "|-" << std::endl;
  std::cout << "|}";
}