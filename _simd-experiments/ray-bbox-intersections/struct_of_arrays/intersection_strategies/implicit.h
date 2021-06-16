#ifndef RAY_BBOX_INTERSECTIONS_ARRAY_OF_STRUCTS_IMPLICIT_H
#define RAY_BBOX_INTERSECTIONS_ARRAY_OF_STRUCTS_IMPLICIT_H

#include <iostream>

namespace implicit {

  template<typename T>
  inline bool intersect(const VBBox<T> &bbox, const Ray<T> &ray, std::size_t i) {

    // Determine bounds x, y, and z
    double xmin = (bbox.bounds()[ray.sign().x()].x()[i] - ray.origin().x()) * ray.inv_direction().x();
    double xmax = (bbox.bounds()[1 - ray.sign().x()].x()[i] - ray.origin().x()) * ray.inv_direction().x();
    double ymin = (bbox.bounds()[ray.sign().y()].y()[i] - ray.origin().y()) * ray.inv_direction().y();
    double ymax = (bbox.bounds()[1 - ray.sign().y()].y()[i] - ray.origin().y()) * ray.inv_direction().y();
    double zmin = (bbox.bounds()[ray.sign().z()].z()[i] - ray.origin().z()) * ray.inv_direction().z();
    double zmax = (bbox.bounds()[1 - ray.sign().z()].z()[i] - ray.origin().z()) * ray.inv_direction().z();

    // Determine the bounds of the overlapping region
    double min = std::max({xmin, ymin, zmin});
    double max = std::min({xmax, ymax, zmax});

    // The ray intercepts if this region overlaps with the interval provided
    return (max >= min);
  }

  template<typename T>
  void intersect(const VBBox<T> &vbbox, const Ray<T> &ray, std::vector<bool> &results) {
    for (std::size_t i = 0; i < vbbox.min().x().size(); ++i)
      results.push_back(intersect(vbbox, ray, i));
  }

}

#endif //RAY_BBOX_INTERSECTIONS_ARRAY_OF_STRUCTS_IMPLICIT_H
