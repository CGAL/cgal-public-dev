#ifndef CGAL_LEVELS_OF_DETAIL_WALL_ESTIMATOR_H
#define CGAL_LEVELS_OF_DETAIL_WALL_ESTIMATOR_H

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
  class Wall_estimator {

  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Segment_2 = typename Traits::Segment_2;
    using Plane_3 = typename Traits::Plane_3;

    Wall_estimator(
      const std::vector<Segment_2>& boundaries,
      const Plane_3& ground_plane,
      const FT building_height) :
    m_boundaries(boundaries),
    m_ground_plane(ground_plane),
    m_building_height(building_height)
    { }

    void estimate(std::vector< std::vector<Point_3> >& walls) const {
      
      walls.clear();
      walls.reserve(m_boundaries.size());

      for (std::size_t i = 0; i < m_boundaries.size(); ++i)
        estimate_wall(m_boundaries[i], walls);
    }

  private:
    const std::vector<Segment_2>& m_boundaries;
    const Plane_3& m_ground_plane;
    const FT m_building_height;

    void estimate_wall(
      const Segment_2& boundary,
      std::vector< std::vector<Point_3> >& walls) const {

      const Point_2& source = boundary.source();
      const Point_2& target = boundary.target();
      
      std::vector<Point_3> wall(4);

      wall[0] = internal::position_on_plane_3(source, m_ground_plane);
      wall[1] = internal::position_on_plane_3(target, m_ground_plane); 

      const FT height1 = m_building_height - wall[1].z();
      const FT height0 = m_building_height - wall[0].z();

      wall[2] = Point_3(wall[1].x(), wall[1].y(), wall[1].z() + height1);
      wall[3] = Point_3(wall[0].x(), wall[0].y(), wall[0].z() + height0);

      walls.push_back(wall);
    }

  }; // Wall_estimator

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_WALL_ESTIMATOR_H
