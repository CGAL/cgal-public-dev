#ifndef RAY_BBOX_INTERSECTIONS_EXPLICIT_H
#define RAY_BBOX_INTERSECTIONS_EXPLICIT_H

#include "../bbox.h"
#include "../ray.h"
#include "../vector3.h"

bool intersect_explicit(const BBox &bbox, const Ray &ray, float rmin, float rmax) {

  // Determine bounds for x and y
  double xmin = (bbox.bounds()[ray.sign()[0]].x() - ray.origin().x()) * ray.inv_direction().x();
  double xmax = (bbox.bounds()[1 - ray.sign()[0]].x() - ray.origin().x()) * ray.inv_direction().x();

  double ymin = (bbox.bounds()[ray.sign()[1]].y() - ray.origin().y()) * ray.inv_direction().y();
  double ymax = (bbox.bounds()[1 - ray.sign()[1]].y() - ray.origin().y()) * ray.inv_direction().y();

  // If the x and y bounds don't overlap, the ray doesn't intersect with the box
  if (xmin > ymax || ymin > xmax) return false;

  // Determine the bounds of the overlapping region
  double min = std::max(xmin, ymin);
  double max = std::min(xmax, ymax);

  // Determine bounds for z
  double zmin = (bbox.bounds()[ray.sign()[2]].z() - ray.origin().z()) * ray.inv_direction().z();
  double zmax = (bbox.bounds()[1 - ray.sign()[2]].z() - ray.origin().z()) * ray.inv_direction().z();

  // If the z bounds don't overlap with the existing region, the ray doesn't intercept
  if (min > zmax || zmin > max) return false;

  // Update the bounds to find the region where x, y, and z overlap
  min = std::max(min, zmin);
  max = std::min(max, zmax);

  // The ray intercepts if this region overlaps with the bounds provided
  return (min < rmax) && (max > rmin);
}


#endif //RAY_BBOX_INTERSECTIONS_EXPLICIT_H
