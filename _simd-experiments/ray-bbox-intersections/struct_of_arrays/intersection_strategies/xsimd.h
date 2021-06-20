#ifndef RAY_BBOX_INTERSECTIONS_ARRAY_OF_STRUCTS_XSIMD_H
#define RAY_BBOX_INTERSECTIONS_ARRAY_OF_STRUCTS_XSIMD_H

#include "../../ray.h"
#include "../vbbox.h"

#include "implicit.h"

#include <xsimd/xsimd.hpp>

// Created with the help of xsimd's documentation
// https://xsimd.readthedocs.io/en/latest/basic_usage.html

namespace xsimd {

  template<typename T, std::size_t N>
  inline xsimd::batch_bool<T, N> intersect(const BBox<xsimd::batch<T, N>> &xbbox,
                                           const Vector3<xsimd::batch<T, N>> &origin,
                                           const Vector3<xsimd::batch<T, N>> &inv_direction, const Vector3<int> sign) {

    // When the ray is negative, flip the box's bounds

    xsimd::batch<T, N> min_bound_x = (sign.x()) ? xbbox.max().x() : xbbox.min().x();
    xsimd::batch<T, N> max_bound_x = (!sign.x()) ? xbbox.max().x() : xbbox.min().x();

    xsimd::batch<T, N> min_bound_y = (sign.y()) ? xbbox.max().y() : xbbox.min().y();
    xsimd::batch<T, N> max_bound_y = (!sign.y()) ? xbbox.max().y() : xbbox.min().y();

    xsimd::batch<T, N> min_bound_z = (sign.z()) ? xbbox.max().z() : xbbox.min().z();
    xsimd::batch<T, N> max_bound_z = (!sign.z()) ? xbbox.max().z() : xbbox.min().z();

    // Calculate bounds for each axis

    xsimd::batch<T, N> min_x = (min_bound_x - origin.x()) * inv_direction.x();
    xsimd::batch<T, N> max_x = (max_bound_x - origin.x()) * inv_direction.x();

    xsimd::batch<T, N> min_y = (min_bound_y - origin.y()) * inv_direction.y();
    xsimd::batch<T, N> max_y = (max_bound_y - origin.y()) * inv_direction.y();

    xsimd::batch<T, N> min_z = (min_bound_z - origin.z()) * inv_direction.z();
    xsimd::batch<T, N> max_z = (max_bound_z - origin.z()) * inv_direction.z();

    // Consolidate bounds into the segment of intersection

    xsimd::batch<T, N> max = xsimd::min(max_x, xsimd::min(max_y, max_z));
    xsimd::batch<T, N> min = xsimd::max(min_x, xsimd::max(min_y, min_z));

    // Intersection exists if the segment has non-negative length

    return max >= min;
  }

  template<typename T, std::size_t N>
  inline xsimd::batch_bool<T, N> intersect(const BBox<xsimd::batch<T, N>> &xbbox, const Ray<T> &ray) {

    // When the ray is negative, flip the box's bounds

    xsimd::batch<T, N> min_bound_x = (ray.sign().x()) ? xbbox.max().x() : xbbox.min().x();
    xsimd::batch<T, N> max_bound_x = (!ray.sign().x()) ? xbbox.max().x() : xbbox.min().x();

    xsimd::batch<T, N> min_bound_y = (ray.sign().y()) ? xbbox.max().y() : xbbox.min().y();
    xsimd::batch<T, N> max_bound_y = (!ray.sign().y()) ? xbbox.max().y() : xbbox.min().y();

    xsimd::batch<T, N> min_bound_z = (ray.sign().z()) ? xbbox.max().z() : xbbox.min().z();
    xsimd::batch<T, N> max_bound_z = (!ray.sign().z()) ? xbbox.max().z() : xbbox.min().z();

    // Calculate bounds for each axis

    xsimd::batch<T, N> min_x = (min_bound_x - ray.origin().x()) * ray.inv_direction().x();
    xsimd::batch<T, N> max_x = (max_bound_x - ray.origin().x()) * ray.inv_direction().x();

    xsimd::batch<T, N> min_y = (min_bound_y - ray.origin().y()) * ray.inv_direction().y();
    xsimd::batch<T, N> max_y = (max_bound_y - ray.origin().y()) * ray.inv_direction().y();

    xsimd::batch<T, N> min_z = (min_bound_z - ray.origin().z()) * ray.inv_direction().z();
    xsimd::batch<T, N> max_z = (max_bound_z - ray.origin().z()) * ray.inv_direction().z();

    // Consolidate bounds into the segment of intersection

    xsimd::batch<T, N> max = xsimd::min(max_x, xsimd::min(max_y, max_z));
    xsimd::batch<T, N> min = xsimd::max(min_x, xsimd::max(min_y, min_z));

    // Intersection exists if the segment has non-negative length

    return max >= min;
  }

  template<typename T>
  void intersect(const VBBox<T> &vbbox, const Ray<T> &ray, std::vector<bool> &results) {

    std::size_t data_size = vbbox.min().x().size();
    constexpr std::size_t batch_size = xsimd::simd_type<double>::size;
    std::size_t aligned_size = data_size - (data_size % batch_size);

    // Create a broadcasted-batch ray
    Vector3<xsimd::batch<T, batch_size>> origin{
            xsimd::batch<T, batch_size>(ray.origin().x()),
            xsimd::batch<T, batch_size>(ray.origin().y()),
            xsimd::batch<T, batch_size>(ray.origin().z())
    };
    Vector3<xsimd::batch<T, batch_size>> inv_direction{
            xsimd::batch<T, batch_size>(ray.inv_direction().x()),
            xsimd::batch<T, batch_size>(ray.inv_direction().y()),
            xsimd::batch<T, batch_size>(ray.inv_direction().z())
    };
    Vector3<int> sign = ray.sign();

    // Load each batch from the vectorized box type
    for (std::size_t i = 0; i < aligned_size; i += batch_size) {

      auto xmin = xsimd::load_unaligned(&vbbox.min().x()[i]);
      auto ymin = xsimd::load_unaligned(&vbbox.min().y()[i]);
      auto zmin = xsimd::load_unaligned(&vbbox.min().z()[i]);

      auto xmax = xsimd::load_unaligned(&vbbox.max().x()[i]);
      auto ymax = xsimd::load_unaligned(&vbbox.max().y()[i]);
      auto zmax = xsimd::load_unaligned(&vbbox.max().z()[i]);

      auto min = Vector3(xmin, ymin, zmin);
      auto max = Vector3(xmax, ymax, zmax);

      auto box = BBox(min, max);

      // Perform intersection tests for that batch, save the results to the output array
      auto result = intersect(box, origin, inv_direction, sign);

      // Add results to output
      // TODO: Is there a faster way to save the results?
      for (std::size_t j = 0; j < batch_size; ++j) {
        results.push_back(result[j]);
      }
    }

    // Perform scalar operations on leftover values
    for (std::size_t i = aligned_size; i < data_size; ++i)
      results.push_back(implicit::intersect(vbbox.getr(i), ray));
  }
}

#endif //RAY_BBOX_INTERSECTIONS_ARRAY_OF_STRUCTS_XSIMD_H
