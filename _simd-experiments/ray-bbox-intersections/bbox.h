#ifndef RAY_BBOX_INTERSECTIONS_BBOX_H
#define RAY_BBOX_INTERSECTIONS_BBOX_H

#include "vector3.h"

#include <cassert>

template<typename T>
struct BBox {

  Vector3<T> min, max;

  BBox(const Vector3<T> &min, const Vector3<T> &max) : min(min), max(max) {
    assert(min.x <= max.x);
    assert(min.y <= max.y);
    assert(min.z <= max.z);
  }

  std::array<std::reference_wrapper<const Vector3<T>>, 2> bounds() const {
    return {std::cref(min), std::cref(max)};
  }

  friend std::istream &operator>>(std::istream &input, BBox &bbox) {
    input >> bbox.min >> bbox.max;
    assert(bbox.min.x <= bbox.max.x);
    assert(bbox.min.y <= bbox.max.y);
    assert(bbox.min.z <= bbox.max.z);
    return input;
  }
};

#endif //RAY_BBOX_INTERSECTIONS_BBOX_H
