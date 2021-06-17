#ifndef RAY_BBOX_INTERSECTIONS_RAY_H
#define RAY_BBOX_INTERSECTIONS_RAY_H

#include <array>
#include "vector3.h"

template<typename T>
class Ray {
protected:

  Vector3<T> _origin, _direction;

  Vector3<T> _inv_direction{0, 0, 0};
  Vector3<int> _sign{0, 0, 0};

public:

  Ray(const Vector3<T> &origin, const Vector3<T> &direction) : _origin(origin), _direction(direction) {

    _inv_direction = {
            1 / direction.x(),
            1 / direction.y(),
            1 / direction.z()
    };

    _sign = {
            _inv_direction.x() < 0,
            _inv_direction.y() < 0,
            _inv_direction.z() < 0
    };
  }

  const Vector3<T> &origin() const { return _origin; }
  const Vector3<T> &direction() const { return _direction; }
  const Vector3<T> &inv_direction() const { return _inv_direction; }
  const Vector3<int> &sign() const { return _sign; }

  bool operator==(const Ray &other) const {
    return other.direction() == this->direction() &&
           other.origin() == this->origin();
  }

  bool operator!=(const Ray &other) const {
    return !operator==(other);
  }

  friend std::istream &operator>>(std::istream &input, Ray &ray) {
    input >> ray._origin >> ray._direction;

    ray._inv_direction = {
            1 / ray.direction().x(),
            1 / ray.direction().y(),
            1 / ray.direction().z()
    };

    ray._sign = {
            ray.inv_direction().x() < 0,
            ray.inv_direction().y() < 0,
            ray.inv_direction().z() < 0
    };

    return input;
  }
};

#endif //RAY_BBOX_INTERSECTIONS_RAY_H
