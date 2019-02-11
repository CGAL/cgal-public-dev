#ifndef CGAL_LEVEL_OF_DETAIL_BUILDINGS_CLUSTERING_H
#define CGAL_LEVEL_OF_DETAIL_BUILDINGS_CLUSTERING_H

// STL includes.
#include <vector>
#include <algorithm>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/structures.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits, 
  typename InputRange, 
  typename PointMap>
  class Buildings_clustering {
  
  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Triangle_2 = typename Traits::Triangle_2;
    
    using Iterator = typename Input_range::const_iterator;
    using Iterators = std::vector<Iterator>;

    using Building = Building<Traits>;

  private:
    const Input_range& m_input_range;
    const Point_map m_point_map;

  public:
    Buildings_clustering(
      const Input_range& input_range, 
      Point_map point_map) : 
    m_input_range(input_range), 
    m_point_map(point_map)
    { }

    void create(
      const std::vector<Building>& buildings,
      std::vector<Iterators>& clusters) const {
      
      clusters.clear();
      clusters.resize(buildings.size());

      for (auto it = m_input_range.begin(); it != m_input_range.end(); ++it) {
        const Point_3& point = get(m_point_map, *it);

        const int building_index = 
        find_building_index(point, buildings);

        if (building_index >= 0)
          clusters[static_cast<std::size_t>(building_index)].push_back(it);
      }
    }

  private:
    int find_building_index(
      const Point_3& query, 
      const std::vector<Building>& buildings) const {

      const Point_2 point = Point_2(query.x(), query.y());

      for (std::size_t i = 0; i < buildings.size(); ++i) {
        const auto& triangles = buildings[i].footprint;

        for (std::size_t j = 0; j < triangles.size(); ++j) {
          const Triangle_2& triangle = triangles[j];

          const auto type = triangle.bounded_side(point);
          if (type == CGAL::ON_BOUNDARY || type == CGAL::ON_BOUNDED_SIDE)
            return static_cast<int>(i);
        }
      }
      return -1;
    }

  }; // Buildings_clustering
  
} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDINGS_CLUSTERING_H
