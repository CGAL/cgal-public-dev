#ifndef RAY_BBOX_INTERSECTIONS_ARRAY_OF_STRUCTS_XSIMD_H
#define RAY_BBOX_INTERSECTIONS_ARRAY_OF_STRUCTS_XSIMD_H

#include "../../ray.h"
#include "../vbbox.h"

#include <xsimd/xsimd.hpp>

// Created with the help of xsimd's documentation
// https://xsimd.readthedocs.io/en/latest/basic_usage.html

namespace xsimd {

  template<typename T>
  void intersect(const VBBox<T> &vbbox, const Ray<T> &ray, std::vector<bool> &results) {

    std::size_t data_size = vbbox.min.x.size();
    constexpr std::size_t batch_size = xsimd::simd_type<double>::size;
    std::size_t aligned_size = data_size - (data_size % batch_size);

    // Load each batch from the vectorized box type
    for (std::size_t i = 0; i < aligned_size; i += batch_size) {

      auto xmin = xsimd::load_unaligned(&vbbox.min.x[i]);
      auto ymin = xsimd::load_unaligned(&vbbox.min.y[i]);
      auto zmin = xsimd::load_unaligned(&vbbox.min.z[i]);

      auto xmax = xsimd::load_unaligned(&vbbox.max.x[i]);
      auto ymax = xsimd::load_unaligned(&vbbox.max.y[i]);
      auto zmax = xsimd::load_unaligned(&vbbox.max.z[i]);

      auto min = Vector3(xmin, ymin, zmin);
      auto max = Vector3(xmax, ymax, zmax);

      auto box = BBox(min, max);

      // Perform intersection tests for that batch, save the results to the output array
      auto result = intersect(box, ray);

      // Add results to output
      // TODO: Is there a faster way to save the results?
      for (std::size_t j = 0; j < batch_size; ++j) {
        results.push_back(result[j]);
      }
    }

    // Perform scalar operations on leftover values
    for (std::size_t i = aligned_size; i < data_size; ++i) {
      // TODO This is just a placeholder!
      results.push_back(false);
    }
  }

  template<typename T, std::size_t N>
  inline xsimd::batch_bool<T, N> intersect(const BBox<xsimd::batch<T, N>> &xbbox, const Ray<T> &ray) {

    // TODO
    return xbbox.min.x > 0;

  }
}

#endif //RAY_BBOX_INTERSECTIONS_ARRAY_OF_STRUCTS_XSIMD_H
