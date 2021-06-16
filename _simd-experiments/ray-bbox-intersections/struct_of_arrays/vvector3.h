#ifndef RAY_BBOX_INTERSECTIONS_ARRAY_OF_BATCHED_STRUCTS_VVECTOR3_H
#define RAY_BBOX_INTERSECTIONS_ARRAY_OF_BATCHED_STRUCTS_VVECTOR3_H

#include "../vector3.h"

#include <vector>

template<typename T>
struct VVector3 : Vector3<std::vector<T>> {

  VVector3(const std::vector<Vector3<T>> &vectors) {
    for (const auto &v : vectors) {
      this->x.push_back(v.x);
      this->y.push_back(v.y);
      this->z.push_back(v.z);
    }
  }
};

#endif //RAY_BBOX_INTERSECTIONS_ARRAY_OF_BATCHED_STRUCTS_VVECTOR3_H
