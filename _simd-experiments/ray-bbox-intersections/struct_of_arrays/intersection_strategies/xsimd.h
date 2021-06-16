#ifndef RAY_BBOX_INTERSECTIONS_ARRAY_OF_STRUCTS_XSIMD_H
#define RAY_BBOX_INTERSECTIONS_ARRAY_OF_STRUCTS_XSIMD_H

#include "../../ray.h"
#include "../vbbox.h"

#include <xsimd/xsimd.hpp>

// Created with the help of xsimd's documentation
// https://xsimd.readthedocs.io/en/latest/basic_usage.html

namespace xsimd {

  template<typename T>
  inline std::vector<bool, xsimd::aligned_allocator<T, 16>> intersect(const VBBox<T> &bbox, const Ray<T> &ray) {

    std::size_t data_size = bbox.min.x.size();
    constexpr std::size_t batch_size = xsimd::simd_type<double>::size;
    std::size_t aligned_size = data_size - (data_size % batch_size);

    for (int i = 0; i < aligned_size; i += batch_size) {

      // TODO Build an xsimd bbox
    }

    std::vector<bool> results;
  }
}

#endif //RAY_BBOX_INTERSECTIONS_ARRAY_OF_STRUCTS_XSIMD_H
