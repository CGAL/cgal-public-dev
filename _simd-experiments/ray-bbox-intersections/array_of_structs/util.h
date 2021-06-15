#ifndef RAY_BBOX_INTERSECTIONS_LOAD_H
#define RAY_BBOX_INTERSECTIONS_LOAD_H

#include <iostream>
#include <random>
#include <chrono>
#include <algorithm>
#include <fstream>
#include "ray.h"
#include "bbox.h"

std::vector<std::pair<Ray<double>, std::vector<BBox<double>>>> load_scenarios(std::ifstream &file, long N) {
  std::vector<std::pair<Ray<double>, std::vector<BBox<double>>>> scenarios;

  T px, py, pz, qx, qy, qz,
          bxmin, bymin, bzmin, bxmax, bymax, bzmax;

  std::vector<T> px, py, pz, qx, qy, qz,
          bxmin, bymin, bzmin, bxmax, bymax, bzmax;

  for (int i = 0; i < N; ++i) {

    // Read ray data
    file >> px >> py >> pz
         >> qx >> qy >> qz;
    Vector3 origin = {px, py, pz};
    Vector3 direction = {qx, qy, qz};
    Ray ray = {origin, direction};

    // Read box data
    file >> bxmin >> bymin >> bzmin
         >> bxmax >> bymax >> bzmax;
    Vector3 min = {bxmin, bymin, bzmin};
    Vector3 max = {bxmax, bymax, bzmax};
    BBox box = {min, max};

    // Only create a new scenario when the query ray has changed
    if (scenarios.empty() || !(ray == scenarios.back().first))

    scenarios.back().second.push_back(box);
  }

  long count = 0;
  for (const auto &scenario : scenarios) count += scenario.second.size();

  return scenarios;
}

double time(const std::function<void(void)> &f) {

  auto start = std::chrono::_V2::high_resolution_clock::now();
  {
    f();
  }
  auto end = std::chrono::_V2::high_resolution_clock::now();

  return (end - start).count();
}

#endif //RAY_BBOX_INTERSECTIONS_LOAD_H
