#ifndef RAY_BBOX_INTERSECTIONS_BRANCHLESS_H
#define RAY_BBOX_INTERSECTIONS_BRANCHLESS_H

// As explained [here](https://cgal.geometryfactory.com/CGAL/Members/wiki/GSoC2021/AABB_tree#Examining_Ray-BBox_Intersection)
// This implementation is a modification of the version described in the paper,
// with the goal of improving auto-vectorization

namespace branchless {

  template<typename T>
  inline bool intersect(const BBox<T> &bbox, const Ray<T> &ray) {

    // Determine bounds x, y, and z
    double xmin = (bbox.bounds()[ray.sign()[0]].x() - ray.origin().x()) * ray.inv_direction().x();
    double xmax = (bbox.bounds()[1 - ray.sign()[0]].x() - ray.origin().x()) * ray.inv_direction().x();
    double ymin = (bbox.bounds()[ray.sign()[1]].y() - ray.origin().y()) * ray.inv_direction().y();
    double ymax = (bbox.bounds()[1 - ray.sign()[1]].y() - ray.origin().y()) * ray.inv_direction().y();
    double zmin = (bbox.bounds()[ray.sign()[2]].z() - ray.origin().z()) * ray.inv_direction().z();
    double zmax = (bbox.bounds()[1 - ray.sign()[2]].z() - ray.origin().z()) * ray.inv_direction().z();

    // Determine the bounds of the overlapping region
    double min = std::max({xmin, ymin, zmin});
    double max = std::min({xmax, ymax, zmax});

    // The ray intercepts if this region overlaps with the interval provided
    return (max > min);
  }

}


#endif //RAY_BBOX_INTERSECTIONS_BRANCHLESS_H
