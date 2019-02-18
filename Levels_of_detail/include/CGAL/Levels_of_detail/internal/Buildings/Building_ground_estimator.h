#ifndef CGAL_LEVELS_OF_DETAIL_BUILDING_GROUND_ESTIMATOR_H
#define CGAL_LEVELS_OF_DETAIL_BUILDING_GROUND_ESTIMATOR_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/number_utils.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utilities.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Building_ground_estimator {

  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Triangle_2 = typename Traits::Triangle_2;
    using Plane_3 = typename Traits::Plane_3;

    Building_ground_estimator(
      const std::vector<Triangle_2>& footprint,
      const Plane_3& ground_plane) :
    m_footprint(footprint),
    m_ground_plane(ground_plane)
    { }

    void estimate(std::vector<Point_3>& ground) const {
      
      FT minx = internal::max_value<FT>();
      FT miny = internal::max_value<FT>();
      FT maxx = -internal::max_value<FT>();
      FT maxy = -internal::max_value<FT>();

      ground.clear();
      for (std::size_t i = 0; i < m_footprint.size(); ++i) {

        for (std::size_t j = 0; j < 3; ++j) {
          const Point_2& p = m_footprint[i][j];

          minx = CGAL::min(minx, p.x());
          miny = CGAL::min(miny, p.y());

          maxx = CGAL::max(minx, p.x());
          maxy = CGAL::max(maxy, p.y());
        }
      }

      ground.resize(4);

      ground[0] = internal::position_on_plane_3(Point_2(minx, miny), m_ground_plane);
      ground[1] = internal::position_on_plane_3(Point_2(maxx, miny), m_ground_plane);
      ground[2] = internal::position_on_plane_3(Point_2(maxx, maxy), m_ground_plane);
      ground[3] = internal::position_on_plane_3(Point_2(minx, maxy), m_ground_plane);
    }

  private:
    const std::vector<Triangle_2>& m_footprint;
    const Plane_3& m_ground_plane;

  }; // Building_ground_estimator

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_BUILDING_GROUND_ESTIMATOR_H
