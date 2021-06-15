#ifndef RAY_BBOX_INTERSECTIONS_XRAY_H
#define RAY_BBOX_INTERSECTIONS_XRAY_H

#include <xsimd/xsimd.hpp>

template<typename T, int N>
struct XRay {

  XVector3<T, N> o, d;

  XRay(XVector3<T, N> origin, XVector3<T, N> direction) : o(origin), d(direction) {

  }
};

#endif //RAY_BBOX_INTERSECTIONS_XRAY_H
