#ifndef RAY_BBOX_INTERSECTIONS_BBOX_H
#define RAY_BBOX_INTERSECTIONS_BBOX_H

#include "vector3.h"

#include <cassert>
#include <functional>

template<typename T>
class BBox {
protected:

  std::array<Vector3<T>, 2> _bounds;

public:

  BBox(const Vector3<T> &min, const Vector3<T> &max) : _bounds{min, max} {}

  const std::array<Vector3<T>, 2> &bounds() const { return _bounds; }

  const Vector3<T> &min() const { return _bounds[0]; }

  const Vector3<T> &max() const { return _bounds[1]; }

  inline friend std::istream &operator>>(std::istream &input, BBox &bbox) {
    input >> bbox._bounds[0] >> bbox._bounds[1];

    for (int i = 0; i < 3; ++i)
      assert(bbox.min().arr()[i] <= bbox.max().arr()[i]);

    return input;
  }
};

#endif //RAY_BBOX_INTERSECTIONS_BBOX_H
