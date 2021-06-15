#ifndef RAY_BBOX_INTERSECTIONS_RAY_H
#define RAY_BBOX_INTERSECTIONS_RAY_H

#include <array>
#include "vector3.h"

template <typename T>
class Ray {
private:
  Vector3<T> o, d;

  Vector3<T> inv{0, 0, 0};
  std::array<int, 3> s;

public:

  Ray(const Vector3<T> &origin, const Vector3<T> &direction) : o(origin), d(direction) {

    inv = {
            1 / direction.x(),
            1 / direction.y(),
            1 / direction.z()
    };

    s = {
            inv.x() < 0,
            inv.y() < 0,
            inv.z() < 0
    };

  }

  [[nodiscard]] const Vector3<T> &origin() const { return o; };

  [[nodiscard]] const Vector3<T> &direction() const { return d; };

  [[nodiscard]] const Vector3<T> &inv_direction() const { return inv; };

  [[nodiscard]] const decltype(s) &sign() const { return s; }

  bool operator==(const Ray &other) const {
    return other.direction() == this->direction() &&
           other.origin() == this->origin();
  }
};

#endif //RAY_BBOX_INTERSECTIONS_RAY_H
