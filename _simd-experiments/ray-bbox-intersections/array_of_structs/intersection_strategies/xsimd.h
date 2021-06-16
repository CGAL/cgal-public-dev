#ifndef RAY_BBOX_INTERSECTIONS_XSIMD_H
#define RAY_BBOX_INTERSECTIONS_XSIMD_H

#include "../../bbox.h"
#include "../../ray.h"
#include "../../vector3.h"

#include <xsimd/xsimd.hpp>

namespace xsimd {

  template<typename T>
  inline bool intersect(const BBox<T> &bbox, const Ray<T> &ray) {

    // Determine intermediate value for minimum bounds()
    xsimd::batch<double, 4> min(
            bbox.bounds()[ray.sign.x].get().x,
            bbox.bounds()[ray.sign.y].get().y,
            bbox.bounds()[ray.sign.z].get().z,
            std::numeric_limits<double>::max()
    );

    // Determine intermediate value for maximum bounds()
    xsimd::batch<double, 4> max(
            bbox.bounds()[1 - ray.sign.x].get().x,
            bbox.bounds()[1 - ray.sign.y].get().y,
            bbox.bounds()[1 - ray.sign.z].get().z,
            std::numeric_limits<double>::min()
    );

    // Apply transform to all bounds()
    xsimd::batch<double, 4> origin(ray.origin.x, ray.origin.y, ray.origin.z, 1);
    xsimd::batch<double, 4> inv_direction(ray.inv_direction.x, ray.inv_direction.y, ray.inv_direction.z, 1);
    min = (min - origin) * inv_direction;
    max = (max - origin) * inv_direction;

    // Determine the minimum bounds() of the region where x, y, and z overlap
    double overall_min = std::max({min[0], min[1], min[2]});

    // Determine the maximum bounds() of the region where x, y, and z overlap
    double overall_max = std::min({max[0], max[1], max[2]});

    // The ray intercepts if this region exists and overlaps with the bounds() provided
    return (overall_max >= overall_min);
  }

}

#endif //RAY_BBOX_INTERSECTIONS_XSIMD_H
