#ifndef RAY_BBOX_INTERSECTIONS_ARRAY_OF_STRUCTS_LOAD_SCENARIOS_H
#define RAY_BBOX_INTERSECTIONS_ARRAY_OF_STRUCTS_LOAD_SCENARIOS_H

#include <fstream>
#include <unordered_map>

#include "ray.h"
#include "bbox.h"

template<typename T>
struct Query {
  Ray<T> ray;
  std::vector<BBox<T>> boxes;
};

template<typename T>
auto load_queries(std::ifstream &file, std::size_t N) {

  // We're going to pair each ray with the bounding boxes it's queried against
  std::vector<Query<T>> queries;

  // Read in each pair of Ray & BBox
  Ray<T> ray{{0, 0, 0}, {0, 0, 0}};
  BBox<T> box{{0, 0, 0}, {0, 0, 0}};
  for (std::size_t i = 0; i < N; ++i) {
    file >> ray >> box;

    // If this ray is different from the last, begin a new query
    if (queries.empty() || queries.back().ray != ray)
      queries.push_back({ray, std::vector<BBox<T>>()});

    // Add the latest box to the existing query
    queries.back().boxes.push_back(box);
  }

  return queries;
}

#endif //RAY_BBOX_INTERSECTIONS_ARRAY_OF_STRUCTS_LOAD_SCENARIOS_H
