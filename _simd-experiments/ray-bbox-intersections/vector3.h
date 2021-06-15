#ifndef RAY_BBOX_INTERSECTIONS_VECTOR3_H
#define RAY_BBOX_INTERSECTIONS_VECTOR3_H

template <typename T>
struct Vector3 {

  T x, y, z;

  Vector3(T x, T y, T z) : x(x), y(y), z(z) {}

  bool operator==(const Vector3 &other) const {
    return other.x == this->x &&
           other.y == this->y &&
           other.z == this->z;
  }

  friend std::istream &operator>>(std::istream &input, Vector3 &vector3) {
    input >> vector3.x >> vector3.y >> vector3.z;
    return input;
  }
};

#endif //RAY_BBOX_INTERSECTIONS_VECTOR3_H
