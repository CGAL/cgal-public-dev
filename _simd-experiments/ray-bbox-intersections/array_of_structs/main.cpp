#include <iostream>
#include <random>
#include <chrono>
#include <algorithm>

#include "../util.h"

#include "../vector3.h"
#include "../ray.h"
#include "../bbox.h"

#include "intersection_strategies/smits_method.h"
#include "intersection_strategies/improved.h"
#include "intersection_strategies/clarified.h"
#include "intersection_strategies/branchless.h"
#include "intersection_strategies/xsimd.h"

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

int main() {


  long N = 3742217;
  long R = 100;

  // Load test data
  auto file = std::ifstream("../data/remeshing_intersections_3742217.txt");
  if (!file.is_open()) return EXIT_FAILURE;
  auto queries = load_queries<double>(file, N);
  std::cout << "Loaded " << std::accumulate(queries.begin(), queries.end(), 0,
                               [](const auto &a, const Query<double> &b) { return a + b.boxes.size(); })
            << " scenarios "
            << "divided into " << queries.size() << " queries." << std::endl;


  std::vector<double> smits_method_times, improved_times, clarified_times, branchless_times, xsimd_times;

  for (int i = 0; i < R; ++i) {
    std::cout << i + 1 << "/" << R  << std::endl;

    std::vector<bool> smits_method_results, improved_results, clarified_results, branchless_results, xsimd_results;

    for (const auto &query : queries) {

      smits_method_times.push_back(time([&] {
        for (const auto &bbox : query.boxes)
          smits_method_results.push_back(smits_method::intersect(bbox, query.ray));
      }));

      improved_times.push_back(time([&] {
        for (const auto &bbox : query.boxes)
          improved_results.push_back(improved::intersect(bbox, query.ray));
      }));
//
//      clarified_times.push_back(time([&] {
//        const auto &ray = query.first;
//        for (const auto &bbox : query.second)
//          sum += clarified::intersect(bbox, ray);
//      }));
//
//      branchless_times.push_back(time([&] {
//        const auto &ray = query.first;
//        for (const auto &bbox : query.second)
//          sum += branchless::intersect(bbox, ray);
//      }));
//
//      xsimd_times.push_back(time([&] {
//        const auto &ray = query.first;
//        for (const auto &bbox : query.second)
//          sum += xsimd::intersect(bbox, ray);
//      }));
//

    }

    if (smits_method_results != improved_results) throw std::logic_error("Incorrect results");

    std::cout << "\tHitrate: "
              << std::accumulate(smits_method_results.begin(), smits_method_results.end(), 0.0) * 100.0 /
                 (double) smits_method_results.size()
              << "%" << std::endl;
  }

  std::cout << std::endl;
  std::cout << "{| class=\"wikitable\"" << std::endl;
  std::cout << "|+ Time to Complete " << N << " Intersection Tests" << std::endl;
  std::cout << "| Smits' Method || "
            << std::accumulate(smits_method_times.begin(), smits_method_times.end(), 0.0) / (double) R
            << " ms"
            << std::endl;
  std::cout << "|-" << std::endl;
  std::cout << "| Improved || "
            << std::accumulate(improved_times.begin(), improved_times.end(), 0.0) / (double) R
            << " ms"
            << std::endl;
  std::cout << "|-" << std::endl;
  std::cout << "| Clarified || "
            << std::accumulate(clarified_times.begin(), clarified_times.end(), 0.0) / (double) R
            << " ms"
            << std::endl;
  std::cout << "|-" << std::endl;
  std::cout << "| Branchless || "
            << std::accumulate(branchless_times.begin(), branchless_times.end(), 0.0) / (double) R
            << " ms"
            << std::endl;
  std::cout << "|-" << std::endl;
  std::cout << "| XSimd || "
            << std::accumulate(xsimd_times.begin(), xsimd_times.end(), 0.0) / (double) R
            << " ms"
            << std::endl;
  std::cout << "|-" << std::endl;
  std::cout << "|}";
}