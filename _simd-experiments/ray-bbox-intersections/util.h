//
// Created by jackcamp on 6/15/21.
//

#ifndef RAY_BBOX_INTERSECTIONS_LOAD_H
#define RAY_BBOX_INTERSECTIONS_LOAD_H

#include "ray.h"
#include "bbox.h"

#include "xray.h"
#include "xbbox.h"

std::vector<std::pair<Ray<double>, std::vector<BBox<double>>>> load_scenarios(std::ifstream &file, long N) {
  std::vector<std::pair<Ray<double>, std::vector<BBox<double>>>> scenarios;

  double px, py, pz, qx, qy, qz,
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
      scenarios.emplace_back(ray, std::vector<BBox<double>>());

    scenarios.back().second.push_back(box);
  }

  long count = 0;
  for (const auto &scenario : scenarios) count += scenario.second.size();

  return scenarios;
}

std::vector<std::pair<XRay<double, 4>, std::vector<XBBox<double, 4>>>>
load_vector_scenarios(std::ifstream &file, long N) {
  auto scalar_scenarios = load_scenarios(file, N);
  std::vector<std::pair<XRay<double, 4>, std::vector<XBBox<double, 4>>>> vector_scenarios;

  for (const auto &scenario : scalar_scenarios) {
    vector_scenarios.emplace_back(XRay<double, 4>(scenario.first), std::vector<XBBox<double, 4>>());
    const auto &boxes = scenario.second;

    for (int b = 3; b < boxes.size(); b += 4) {

      vector_scenarios.back().second.emplace_back(std::array<BBox<double>, 4>{
              boxes[b - 3], boxes[b - 2], boxes[b - 1], boxes[b]
      });
    }
  }

  return vector_scenarios;
}

#endif //RAY_BBOX_INTERSECTIONS_LOAD_H
