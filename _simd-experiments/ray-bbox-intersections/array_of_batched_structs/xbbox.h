

#ifndef RAY_BBOX_INTERSECTIONS_XBBOX_H
#define RAY_BBOX_INTERSECTIONS_XBBOX_H

#include <xsimd/xsimd.hpp>

template<typename T, int N>
struct XBBox {

  std::array<XVector3<T, N>, 2> arr;

  explicit XBBox(XVector3<T, N> min, XVector3<T, N> max) : arr(min, max) {}

};

#endif //RAY_BBOX_INTERSECTIONS_XBBOX_H
