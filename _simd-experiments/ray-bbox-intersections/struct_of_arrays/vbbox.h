
#ifndef RAY_BBOX_INTERSECTIONS_ARRAY_OF_BATCHED_STRUCTS_VBBOX_H
#define RAY_BBOX_INTERSECTIONS_ARRAY_OF_BATCHED_STRUCTS_VBBOX_H

#include "../bbox.h"

#include <xsimd/xsimd.hpp>
#include <iostream>

template<typename T>
struct VBBox : public BBox<std::vector<T>> {
public:

  explicit VBBox(const std::vector<BBox<T>> &boxes) :
          BBox<std::vector<T>>({{},
                                {},
                                {}},
                               {{},
                                {},
                                {}}) {

    std::vector<Vector3<T>> mins, maxs;

    for (const BBox<T> &b : boxes) {
      mins.push_back(b.min());
      maxs.push_back(b.max());
    }

    for (const auto &v : mins) {
      this->_bounds[0].x().push_back(v.x());
      this->_bounds[0].y().push_back(v.y());
      this->_bounds[0].z().push_back(v.z());
    }

    for (const auto &v : maxs) {
      this->_bounds[1].x().push_back(v.x());
      this->_bounds[1].y().push_back(v.y());
      this->_bounds[1].z().push_back(v.z());
    }
  }

  BBox<const T *> getp(std::size_t i) const {
    return {Vector3<const T *>{&(this->min().x()[i]), &(this->min().y()[i]), &(this->min().z()[i])},
            Vector3<const T *>{&(this->max().x()[i]), &(this->max().y()[i]), &(this->max().z()[i])}};
  }

  BBox<std::reference_wrapper<const T>> getr(std::size_t i) const {
    return {Vector3<std::reference_wrapper<const T>>{std::cref(this->min().x()[i]),
                                                     std::cref(this->min().y()[i]),
                                                     std::cref(this->min().z()[i])},
            Vector3<std::reference_wrapper<const T>>{std::cref(this->max().x()[i]),
                                                     std::cref(this->max().y()[i]),
                                                     std::cref(this->max().z()[i])}
    };
  }
};

#endif //RAY_BBOX_INTERSECTIONS_ARRAY_OF_BATCHED_STRUCTS_VBBOX_H
