#ifndef RAY_BBOX_INTERSECTIONS_IMPROVED_H
#define RAY_BBOX_INTERSECTIONS_IMPROVED_H

#include "../../bbox.h"
#include "../../ray.h"
#include "../../vector3.h"

// As explained [here](http://people.csail.mit.edu/amy/papers/box-jgt.ps)
// This implementation generally sticks close to the one described in the paper
namespace improved {

  template<typename T>
  inline bool intersect(const BBox<T> &bbox, const Ray<T> &ray) {

    double tmin, tmax, tymin, tymax, tzmin, tzmax;

    tmin = (bbox.bounds()[ray.sign[0]].get().x - ray.origin.x) * ray.inv_direction.x;
    tmax = (bbox.bounds()[1 - ray.sign[0]].get().x - ray.origin.x) * ray.inv_direction.x;

    tymin = (bbox.bounds()[ray.sign[1]].get().y - ray.origin.y) * ray.inv_direction.y;
    tymax = (bbox.bounds()[1 - ray.sign[1]].get().y - ray.origin.y) * ray.inv_direction.y;

    if ((tmin > tymax) || (tymin > tmax)) return false;
    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;

    tzmin = (bbox.bounds()[ray.sign[2]].get().z - ray.origin.z) * ray.inv_direction.z;
    tzmax = (bbox.bounds()[1 - ray.sign[2]].get().z - ray.origin.z) * ray.inv_direction.z;

    return !((tmin > tzmax) || (tzmin > tmax));

  }

}

#endif //RAY_BBOX_INTERSECTIONS_IMPROVED_H
