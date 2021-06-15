#ifndef RAY_BBOX_INTERSECTIONS_SMITS_METHOD_H
#define RAY_BBOX_INTERSECTIONS_SMITS_METHOD_H

#include "../../bbox.h"
#include "../../ray.h"
#include "../../vector3.h"

// As explained [here](http://people.csail.mit.edu/amy/papers/box-jgt.ps)
// This implementation generally sticks close to the one described in the paper

namespace smits_method {

  template<typename T>
  inline bool intersect(const BBox<T> &bbox, const Ray<T> &ray) {
    double tmin, tmax, tymin, tymax, tzmin, tzmax;

    if (ray.direction.x >= 0) {
      tmin = (bbox.bounds()[0].get().x - ray.origin.x) / ray.direction.x;
      tmax = (bbox.bounds()[1].get().x - ray.origin.x) / ray.direction.x;
    } else {
      tmin = (bbox.bounds()[1].get().x - ray.origin.x) / ray.direction.x;
      tmax = (bbox.bounds()[0].get().x - ray.origin.x) / ray.direction.x;
    }

    if (ray.direction.y >= 0) {
      tymin = (bbox.bounds()[0].get().y - ray.origin.y) / ray.direction.y;
      tymax = (bbox.bounds()[1].get().y - ray.origin.y) / ray.direction.y;
    } else {
      tymin = (bbox.bounds()[1].get().y - ray.origin.y) / ray.direction.y;
      tymax = (bbox.bounds()[0].get().y - ray.origin.y) / ray.direction.y;
    }

    if ((tmin > tymax) || (tymin > tmax)) return false;
    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;

    if (ray.direction.z >= 0) {
      tzmin = (bbox.bounds()[0].get().z - ray.origin.z) / ray.direction.z;
      tzmax = (bbox.bounds()[1].get().z - ray.origin.z) / ray.direction.z;
    } else {
      tzmin = (bbox.bounds()[1].get().z - ray.origin.z) / ray.direction.z;
      tzmax = (bbox.bounds()[0].get().z - ray.origin.z) / ray.direction.z;
    }

    if ((tmin > tzmax) || (tzmin > tmax)) return false;

    return true;
  }

}

#endif //RAY_BBOX_INTERSECTIONS_SMITS_METHOD_H
