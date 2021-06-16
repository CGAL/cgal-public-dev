
#ifndef RAY_BBOX_INTERSECTIONS_ARRAY_OF_BATCHED_STRUCTS_VBBOX_H
#define RAY_BBOX_INTERSECTIONS_ARRAY_OF_BATCHED_STRUCTS_VBBOX_H

#include "../bbox.h"

template<typename T>
struct VBBox : BBox<std::vector<T>> {

  VBBox(const std::vector<BBox<T>> &boxes) {

    std::vector<Vector3<T>> mins, maxs;

    for (const BBox<T> &b : boxes) {
      mins.push_back(b.min);
      maxs.push_back(b.max);
    }
  }
};

#endif //RAY_BBOX_INTERSECTIONS_ARRAY_OF_BATCHED_STRUCTS_VBBOX_H
