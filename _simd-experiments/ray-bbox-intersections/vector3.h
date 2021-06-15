#ifndef RAY_BBOX_INTERSECTIONS_VECTOR3_H
#define RAY_BBOX_INTERSECTIONS_VECTOR3_H

template <typename T>
class Vector3 {
private:

  T val[3];

public:

  Vector3(T x, T y, T z) : val{x, y, z} {}

  [[nodiscard]] const decltype(val) &arr() const { return val; }

  [[nodiscard]] const T &x() const { return val[0]; }

  [[nodiscard]] const T &y() const { return val[1]; }

  [[nodiscard]] const T &z() const { return val[2]; }

  bool operator==(const Vector3 &other) const {
    return other.x() == this->x() &&
           other.y() == this->y() &&
           other.z() == this->z();
  }
};

#endif //RAY_BBOX_INTERSECTIONS_VECTOR3_H
