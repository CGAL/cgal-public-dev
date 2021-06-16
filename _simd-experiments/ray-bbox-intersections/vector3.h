#ifndef RAY_BBOX_INTERSECTIONS_VECTOR3_H
#define RAY_BBOX_INTERSECTIONS_VECTOR3_H

template<typename T>
class Vector3 {
protected:

  std::array<T, 3> _arr;

public:

  Vector3(T x, T y, T z) : _arr{x, y, z} {}

  const T &x() const { return _arr[0]; };
  T &x() { return _arr[0]; };

  const T &y() const { return _arr[1]; };
  T &y() { return _arr[1]; };

  const T &z() const { return _arr[2]; };
  T &z() { return _arr[2]; };

  const std::array<T, 3> &arr() const { return _arr; };
  std::array<T, 3> &arr() { return _arr; };

  bool operator==(const Vector3 &other) const {
    return other.x() == this->x() &&
           other.y() == this->y() &&
           other.z() == this->z();
  }

  inline friend std::istream &operator>>(std::istream &input, Vector3 &vector3) {
    input >> vector3._arr[0] >> vector3._arr[1] >> vector3._arr[2];
    return input;
  }
};

#endif //RAY_BBOX_INTERSECTIONS_VECTOR3_H
