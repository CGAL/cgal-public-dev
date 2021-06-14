#ifndef RAY_BBOX_INTERSECTIONS_IMPROVED_H
#define RAY_BBOX_INTERSECTIONS_IMPROVED_H

#include "../bbox.h"
#include "../ray.h"
#include "../vector3.h"

// As explained [here](http://people.csail.mit.edu/amy/papers/box-jgt.ps)
// This implementation generally sticks close to the one described in the paper
namespace improved {

  inline bool intersect(const BBox &bbox, const Ray &ray, float t0, float t1) {

    double tmin, tmax, tymin, tymax, tzmin, tzmax;

    tmin = (bbox.bounds()[ray.sign()[0]].x() - ray.origin().x()) * ray.inv_direction().x();
    tmax = (bbox.bounds()[1 - ray.sign()[0]].x() - ray.origin().x()) * ray.inv_direction().x();

    tymin = (bbox.bounds()[ray.sign()[1]].y() - ray.origin().y()) * ray.inv_direction().y();
    tymax = (bbox.bounds()[1 - ray.sign()[1]].y() - ray.origin().y()) * ray.inv_direction().y();

    if ((tmin > tymax) || (tymin > tmax)) return false;
    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;

    tzmin = (bbox.bounds()[ray.sign()[2]].z() - ray.origin().z()) * ray.inv_direction().z();
    tzmax = (bbox.bounds()[1 - ray.sign()[2]].z() - ray.origin().z()) * ray.inv_direction().z();

    if ((tmin > tzmax) || (tzmin > tmax)) return false;
    if (tzmin > tmin) tmin = tzmin;
    if (tzmax < tmax) tmax = tzmax;

    return (tmin < t1) && (tmax > t0);

  }

}

#endif //RAY_BBOX_INTERSECTIONS_IMPROVED_H
