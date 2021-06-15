#ifndef RAY_BBOX_INTERSECTIONS_XRAY_H
#define RAY_BBOX_INTERSECTIONS_XRAY_H

#include "ray.h"

#include <xsimd/xsimd.hpp>

template<typename T, int N>
class XRay : Ray<xsimd::batch<T, N>> {
public:

  XRay(Ray<T> ray) {

  }
};

#endif //RAY_BBOX_INTERSECTIONS_XRAY_H
