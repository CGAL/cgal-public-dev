#include <iostream>
#include <random>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <vector>

#include <xsimd/xsimd.hpp>

#include "../ray.h"
#include "../bbox.h"

#include "../util.h"

#include "intersection_strategies/implicit.h"
#include "intersection_strategies/xsimd.h"

template<typename T, std::size_t N>
struct XQuery {
  Ray<T> ray;
  std::vector<BBox<xsimd::batch<T, N>>> boxes;
  std::vector<BBox<T>> leftover_boxes;
};

template<typename T, std::size_t N>
auto pack_queries(std::vector<Query<T>> queries) {
  std::vector<XQuery<T, N>> xqueries;

  for (Query<T> &query : queries) {
    std::vector<BBox<xsimd::batch<T, N>>> xboxes;

    // Take chunkes of size N until there are fewer than N boxes remaining
    while (query.boxes.size() > N) {

      // Consume N boxes to create a xsimd-bbox
      xsimd::batch<T, N> xmin, ymin, zmin, xmax, ymax, zmax;
      for (std::size_t i = 0; i < N; ++i) {
        const BBox<T> &box = query.boxes.back();
        query.boxes.pop_back();

        xmin[i] = box.min().x();
        ymin[i] = box.min().y();
        zmin[i] = box.min().z();

        xmax[i] = box.max().x();
        ymax[i] = box.max().y();
        zmax[i] = box.max().z();
      }

      // Assemble the bbox and add it to this query
      xboxes.emplace_back(Vector3{xmin, ymin, zmin}, Vector3{xmax, ymax, zmax});
    }

    // The remainder of the boxes make up the leftover list
    const auto &leftover_boxes = query.boxes;


    // Create a new query in the list
    xqueries.push_back({query.ray, xboxes, leftover_boxes});
  }

  return xqueries;
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
  auto vqueries = pack_queries<double, 4>(queries);

  std::vector<double> explicit_times, implicit_times;

  for (int i = 0; i < R; ++i) {
    std::cout << i + 1 << "/" << R << std::endl;

    std::vector<bool> explicit_results, implicit_results;

    for (const auto &query : vqueries) {

      explicit_times.push_back(time([&] {
        xsimd::intersect(query.boxes, query.leftover_boxes, query.ray, explicit_results);
      }));

      implicit_times.push_back(time([&] {
        implicit::intersect(query.boxes, query.leftover_boxes, query.ray, implicit_results);
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
  std::cout << "|+ Time to Complete " << N << " Intersection Tests, AoSoA Data" << std::endl;
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

