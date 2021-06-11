#ifndef RAY_BBOX_INTERSECTIONS_SMITS_METHOD_H
#define RAY_BBOX_INTERSECTIONS_SMITS_METHOD_H

#include "../bbox.h"
#include "../ray.h"
#include "../vector3.h"

// As explained [here](http://people.csail.mit.edu/amy/papers/box-jgt.ps)
// This implementation generally sticks close to the one described in the paper

bool intersect_smits_method(const BBox &bbox, const Ray &ray, float t0, float t1) {
    double tmin, tmax, tymin, tymax, tzmin, tzmax;

    if (ray.direction().x() >= 0) {
        tmin = (bbox.bounds()[0].x() - ray.origin().x()) / ray.direction().x();
        tmax = (bbox.bounds()[1].x() - ray.origin().x()) / ray.direction().x();
    } else {
        tmin = (bbox.bounds()[1].x() - ray.origin().x()) / ray.direction().x();
        tmax = (bbox.bounds()[0].x() - ray.origin().x()) / ray.direction().x();
    }

    if (ray.direction().y() >= 0) {
        tymin = (bbox.bounds()[0].y() - ray.origin().y()) / ray.direction().y();
        tymax = (bbox.bounds()[1].y() - ray.origin().y()) / ray.direction().y();
    } else {
        tymin = (bbox.bounds()[1].y() - ray.origin().y()) / ray.direction().y();
        tymax = (bbox.bounds()[0].y() - ray.origin().y()) / ray.direction().y();
    }

    if ((tmin > tymax) || (tymin > tmax)) return false;
    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;

    if (ray.direction().z() >= 0) {
        tzmin = (bbox.bounds()[0].z() - ray.origin().z()) / ray.direction().z();
        tzmax = (bbox.bounds()[1].z() - ray.origin().z()) / ray.direction().z();
    } else {
        tzmin = (bbox.bounds()[1].z() - ray.origin().z()) / ray.direction().z();
        tzmax = (bbox.bounds()[0].z() - ray.origin().z()) / ray.direction().z();
    }

    if ((tmin > tzmax) || (tzmin > tmax)) return false;
    if (tzmin > tmin) tmin = tzmin;
    if (tzmax < tmax) tmax = tzmax;

    return (tmin < t1) && (tmax > t0);
}

#endif //RAY_BBOX_INTERSECTIONS_SMITS_METHOD_H
