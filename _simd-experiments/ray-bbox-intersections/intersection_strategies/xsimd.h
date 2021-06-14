#ifndef RAY_BBOX_INTERSECTIONS_XSIMD_H
#define RAY_BBOX_INTERSECTIONS_XSIMD_H

#include "../bbox.h"
#include "../ray.h"
#include "../vector3.h"

#include <xsimd/xsimd.hpp>

namespace xs = xsimd;

bool intersect_xsimd(const BBox &bbox, const Ray &ray, float rmin, float rmax) {

  // Determine intermediate value for minimum bounds
  xs::batch<double, 4> min(
          bbox.bounds()[ray.sign()[0]].x(),
          bbox.bounds()[ray.sign()[1]].y(),
          bbox.bounds()[ray.sign()[2]].z(),
          std::numeric_limits<double>::max()
  );

  // Determine intermediate value for maximum bounds
  xs::batch<double, 4> max(
          bbox.bounds()[1 - ray.sign()[0]].x(),
          bbox.bounds()[1 - ray.sign()[1]].y(),
          bbox.bounds()[1 - ray.sign()[2]].z(),
          std::numeric_limits<double>::min()
  );

  // Apply transform to all bounds
  xs::batch<double, 4> origin(ray.origin().x(), ray.origin().y(), ray.origin().z(), 1);
  xs::batch<double, 4> inv_direction(ray.inv_direction().x(), ray.inv_direction().y(), ray.inv_direction().z(), 1);
  min = (min - origin) * inv_direction;
  max = (max - origin) * inv_direction;

  // Determine the minimum bounds of the region where x, y, and z overlap
  double overall_min = std::max({min[0], min[1], min[2]});

  // Determine the maximum bounds of the region where x, y, and z overlap
  double overall_max = std::min({max[0], max[1], max[2]});

  // The ray intercepts if this region exists and overlaps with the bounds provided
  return (overall_max > overall_min) && (overall_min < rmax) && (overall_max > rmin);
}

#endif //RAY_BBOX_INTERSECTIONS_XSIMD_H
