#ifndef CGAL_INTERNAL_SPLAT_SURFACE_RECONSTRUCTION_3_SURFACE_APPROXIMATION_H
#define CGAL_INTERNAL_SPLAT_SURFACE_RECONSTRUCTION_3_SURFACE_APPROXIMATION_H

#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Kernel/global_functions_3.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

#include "box_grid.h"
#include "geometry_utils.h"

namespace CGAL {
namespace internal {
namespace Splat_surface_reconstruction_3 {

template <typename Kernel_>
class Surface_approximation {
  public:
    using Kernel   = Kernel_;
    using FT       = typename Kernel::FT;
    using Point_3  = typename Kernel::Point_3;
    using Point_2  = typename Kernel::Point_2;
    using Vector_3 = typename Kernel::Vector_3;
    using Index    = std::size_t;

    struct Initial_seed {
      Point_3 first;
      Point_3 second;
      Vector_3 normal;
    };

    struct Feature_segment {
      Point_3 point_1;
      Point_3 point_2;
      Vector_3 normal;

      Feature_segment(const Point_3& point_1, const Point_3& point_2, const Vector_3& normal)
        : point_1(point_1)
        , point_2(point_2)
        , normal(normal)
        {}
    };

    virtual ~Surface_approximation() = default;

    virtual Initial_seed get_initial_seed() const = 0;

    virtual Vector_3 box_normal(const Point_3& p) const = 0;

    virtual std::vector<std::pair<Point_3, Vector_3>> intersect_circle(const Point_3& center, const Vector_3& normal, const FT radius) const = 0;

    virtual std::vector<Feature_segment> find_feature_segments(FT minimum_angle) = 0;

    virtual std::vector<std::pair<Point_3, Vector_3>> find_initial_feature_vertex_candidates(FT minimum_angle) = 0;

    virtual void write_to_obj(const char * filename) const = 0;
};

} // namespace Splat_surface_reconstruction_3
} // namespace internal
} // namespace CGAL

#endif