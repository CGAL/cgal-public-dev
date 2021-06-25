#ifndef RAY_BBOX_INTERSECTIONS_ARRAY_OF_STRUCTS_OF_ARRAYS_XSIMD_H
#define RAY_BBOX_INTERSECTIONS_ARRAY_OF_STRUCTS_OF_ARRAYS_XSIMD_H

#include "../../bbox.h"
#include "../../ray.h"

#include "../../struct_of_arrays/intersection_strategies/xsimd.h"
#include "../../array_of_structs/intersection_strategies/branchless.h"

namespace xsimd {

  template<typename T, std::size_t N>
  void intersect(const std::vector<BBox<xsimd::batch<T, N>>> &boxes, const std::vector<BBox<T>> &leftover_boxes,
                 const Ray<T> &ray, std::vector<bool> &results) {

    // First, broadcast the ray
    Vector3<xsimd::batch<T, N>> origin{
            xsimd::batch<T, N>(ray.origin().x()),
            xsimd::batch<T, N>(ray.origin().y()),
            xsimd::batch<T, N>(ray.origin().z())
    };
    Vector3<xsimd::batch<T, N>> inv_direction{
            xsimd::batch<T, N>(ray.inv_direction().x()),
            xsimd::batch<T, N>(ray.inv_direction().y()),
            xsimd::batch<T, N>(ray.inv_direction().z())
    };
    Vector3<int> sign = ray.sign();

    // Check for intersections with all batched boxes
    for (const auto &xbox : boxes) {
      const auto &r = intersect(xbox, origin, inv_direction, sign);
      for (std::size_t j = 0; j < N; ++j) {
        results.push_back(r[j]);
      }
    }

    // Check for intersections with leftovers
    for (const auto &box : leftover_boxes) {
      results.push_back(branchless::intersect(box, ray));
    }
  }
}

#endif //RAY_BBOX_INTERSECTIONS_ARRAY_OF_STRUCTS_OF_ARRAYS_XSIMD_H
