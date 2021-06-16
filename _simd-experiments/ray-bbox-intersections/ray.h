#ifndef RAY_BBOX_INTERSECTIONS_RAY_H
#define RAY_BBOX_INTERSECTIONS_RAY_H

#include <array>
#include "vector3.h"

template<typename T>
struct Ray {

  Vector3<T> origin, direction;

  Vector3<T> inv_direction{0, 0, 0};
  Vector3<short> sign{0, 0, 0};

  Ray(const Vector3<T> &origin, const Vector3<T> &direction) : origin(origin), direction(direction) {

    inv_direction = {
            1 / direction.x,
            1 / direction.y,
            1 / direction.z
    };

    sign = {
            inv_direction.x < 0,
            inv_direction.y < 0,
            inv_direction.z < 0
    };

  }

  bool operator==(const Ray &other) const {
    return other.direction == this->direction &&
           other.origin == this->origin;
  }

  bool operator!=(const Ray &other) const {
    return !operator==(other);
  }

  friend std::istream &operator>>(std::istream &input, Ray &ray) {
    input >> ray.origin >> ray.direction;

    ray.inv_direction = {
            1 / ray.direction.x,
            1 / ray.direction.y,
            1 / ray.direction.z
    };

    ray.sign = {
            ray.inv_direction.x < 0,
            ray.inv_direction.y < 0,
            ray.inv_direction.z < 0
    };

    return input;
  }
};

#endif //RAY_BBOX_INTERSECTIONS_RAY_H
