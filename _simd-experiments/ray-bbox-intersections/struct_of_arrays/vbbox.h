
#ifndef RAY_BBOX_INTERSECTIONS_ARRAY_OF_BATCHED_STRUCTS_VBBOX_H
#define RAY_BBOX_INTERSECTIONS_ARRAY_OF_BATCHED_STRUCTS_VBBOX_H

#include "../bbox.h"

#include <xsimd/xsimd.hpp>
#include <iostream>

template<typename T>
struct VBBox : BBox<std::vector<T>> {

  explicit VBBox(const std::vector<BBox<T>> &boxes) :
  BBox<std::vector<T>>({{}, {}, {}}, {{}, {}, {}}) {

    std::vector<Vector3<T>> mins, maxs;

    for (const BBox<T> &b : boxes) {
      mins.push_back(b.min);
      maxs.push_back(b.max);
    }

    for (const auto &v : mins) {
      this->min.x.push_back(v.x);
      this->min.y.push_back(v.y);
      this->min.z.push_back(v.z);
    }

    for (const auto &v : maxs) {
      this->max.x.push_back(v.x);
      this->max.y.push_back(v.y);
      this->max.z.push_back(v.z);
    }
  }
};

#endif //RAY_BBOX_INTERSECTIONS_ARRAY_OF_BATCHED_STRUCTS_VBBOX_H
