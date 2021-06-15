
#ifndef RAY_BBOX_INTERSECTIONS_XVECTOR_H
#define RAY_BBOX_INTERSECTIONS_XVECTOR_H

#include <xsimd/xsimd.hpp>

template<typename T, int N>
struct XVector3 {

  xsimd::batch<T, N> x, y, z;

  XVector3(std::array<T, N> x, std::array<T, N> y, std::array<T, N> z) {
    x = xsimd::batch<T, N>(x);
    y = xsimd::batch<T, N>(y);
    z = xsimd::batch<T, N>(z);
  }
};

#endif //RAY_BBOX_INTERSECTIONS_XVECTOR_H
