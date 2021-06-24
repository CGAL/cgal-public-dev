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
      for (int i = 0; i < N; ++i) {
        const BBox<T> &box = query.boxes.back();
        query.boxes.pop_back();

        xmin[i] = box.min().x()[i];
        ymin[i] = box.min().y()[i];
        zmin[i] = box.min().z()[i];

        xmax[i] = box.max().x()[i];
        ymax[i] = box.max().y()[i];
        zmax[i] = box.max().z()[i];
      }

      // Assemble the bbox and add it to this query
      xboxes.emplace_back({xmin, ymin, zmin}, {xmax, ymax, zmax});
    }

    // The remainder of the boxes make up the leftover list
    const auto &leftover_boxes = query.boxes;

    // Create a new query in the list
    xqueries.emplace_back(query.ray, xboxes, leftover_boxes);
  }
}

int main() {

}

