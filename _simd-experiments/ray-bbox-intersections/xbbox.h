

#ifndef RAY_BBOX_INTERSECTIONS_XBBOX_H
#define RAY_BBOX_INTERSECTIONS_XBBOX_H

#include "bbox.h"

#include <xsimd/xsimd.hpp>

template<typename T, int N>
class XBBox : BBox<xsimd::batch<T, N>> {
public:

  explicit XBBox(std::array<BBox<T>, N> boxes) {

  }
};

#endif //RAY_BBOX_INTERSECTIONS_XBBOX_H
