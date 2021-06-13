#ifndef RAY_BBOX_INTERSECTIONS_STD_SIMD_H
#define RAY_BBOX_INTERSECTIONS_STD_SIMD_H

#include "../bbox.h"
#include "../ray.h"
#include "../vector3.h"

#include <experimental/simd>

bool intersect_std_simd(const BBox &bbox, const Ray &ray, float rmin, float rmax) {

  // TODO: This is temporarily identical to the clarified approach

  std::experimental::fixed_size_simd<double, 3> min;

  // Determine intermediate value for minimum bounds
  float xmin = box.bounds[ray.sign[0]].x;
  float ymin = box.bounds[ray.sign[1]].y;
  float zmin = box.bounds[ray.sign[2]].z;
  /* pad_min = FLOAT_MAX; */

  // Determine intermediate value for maximum bounds
  float xmax = box.bounds[1 - ray.sign[0]].x;
  float ymax = box.bounds[1 - ray.sign[1]].y;
  float zmax = box.bounds[1 - ray.sign[2]].z;
  /* pad_max = FLOAT_MIN; */

  // Apply transform to all bounds
  float xmin = (xmin - r.origin.x) * r.inv_direction.x;
  float ymin = (ymin - r.origin.y) * r.inv_direction.y;
  float zmin = (zmin - r.origin.z) * r.inv_direction.z;
  /* pad_min = (pading - 0.0)      * 1.0; /*
  float xmax = (zmax - r.origin.x) * r.inv_direction.x;
  float ymax = (ymax - r.origin.y) * r.inv_direction.y;
  float zmax = (zmax - r.origin.z) * r.inv_direction.z;
  /* pad_max = (pading - 0.0)      * 1.0; /*

  // Determine the minimum bounds of the region where x, y, and z overlap
  min = MAX(xmin, ymin, zmin /*, pad_min */);

  // Determine the maximum bounds of the region where x, y, and z overlap
  max = MIN(xmax, ymax, zmax /*, pad_max */);

  // The ray intercepts if this region exists and overlaps with the bounds provided
  return (max > min) && (min < rmax) && (max > rmin);

}

#endif //RAY_BBOX_INTERSECTIONS_STD_SIMD_H
