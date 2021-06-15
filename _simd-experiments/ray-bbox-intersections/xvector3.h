
#ifndef RAY_BBOX_INTERSECTIONS_XVECTOR_H
#define RAY_BBOX_INTERSECTIONS_XVECTOR_H

#include "vector3.h"

#include <xsimd/xsimd.hpp>

template<typename T, int N>
class XVector3 : Vector3<xsimd::batch<T, N>> {
public:

  XVector3(std::array<Vector3<T> N> vectors) {

  }
};

#endif //RAY_BBOX_INTERSECTIONS_XVECTOR_H
