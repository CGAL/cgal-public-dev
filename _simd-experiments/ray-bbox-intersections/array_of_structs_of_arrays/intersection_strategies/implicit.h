#ifndef RAY_BBOX_INTERSECTIONS_ARRAY_OF_STRUCTS_OF_ARRAYS_IMPLICIT_H
#define RAY_BBOX_INTERSECTIONS_ARRAY_OF_STRUCTS_OF_ARRAYS_IMPLICIT_H

#include "../../bbox.h"
#include "../../ray.h"

#include "../../struct_of_arrays/intersection_strategies/implicit.h"
#include "../../array_of_structs/intersection_strategies/branchless.h"

namespace implicit {

  template<typename T, std::size_t N>
  inline bool intersect(const BBox<xsimd::batch<T, N>> &bbox, const Ray<T> &ray, std::size_t i) {

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

  template<typename T, std::size_t N>
  void intersect(const std::vector<BBox<xsimd::batch<T, N>>> &boxes, const std::vector<BBox<T>> &leftover_boxes,
                 const Ray<T> &ray, std::vector<bool> &results) {

    // Check for intersections with all batched boxes
    for (const auto &xbox : boxes) {
      for (std::size_t j = 0; j < N; ++j)
        results.push_back(intersect(xbox, ray, j));
    }

    // Check for intersections with leftovers
    for (const auto &box : leftover_boxes) {
      results.push_back(branchless::intersect(box, ray));
    }
  }
}


#endif //RAY_BBOX_INTERSECTIONS_ARRAY_OF_STRUCTS_OF_ARRAYS_IMPLICIT_H
