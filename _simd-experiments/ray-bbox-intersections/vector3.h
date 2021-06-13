#ifndef RAY_BBOX_INTERSECTIONS_VECTOR3_H
#define RAY_BBOX_INTERSECTIONS_VECTOR3_H

class Vector3 {
private:

  double val[3];

public:

  Vector3(double x, double y, double z) : val{x, y, z} {}

  [[nodiscard]] const decltype(val) &arr() const { return val; }

  [[nodiscard]] const double &x() const { return val[0]; }

  [[nodiscard]] const double &y() const { return val[1]; }

  [[nodiscard]] const double &z() const { return val[2]; }

  bool operator==(const Vector3 &other) const {
    return other.x() == this->x() &&
           other.y() == this->y() &&
           other.z() == this->z();
  }
};

#endif //RAY_BBOX_INTERSECTIONS_VECTOR3_H
