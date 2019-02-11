#ifndef CGAL_LEVELS_OF_DETAIL_BUILDING_HEIGHT_ESTIMATOR_H
#define CGAL_LEVELS_OF_DETAIL_BUILDING_HEIGHT_ESTIMATOR_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/number_utils.h>

// Internal includes.
#include <CGAL/Levels_of_detail/enumerations.h>
#include <CGAL/Levels_of_detail/internal/structures.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits, 
  typename InputRange, 
  typename PointMap>
  class Building_height_estimator {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;

    using Iterator = typename Input_range::const_iterator;
    using Iterators = std::vector<Iterator>;
    using Building = Building<Traits>;

    Building_height_estimator(
      const std::vector<Iterators>& clusters, 
      Point_map point_map) : 
    m_clusters(clusters), 
    m_point_map(point_map)
    { }

    void compute_heights(
      const Extrusion_type extrusion_type,
      std::vector<Building>& buildings) const {

      switch (extrusion_type) {

        case Extrusion_type::MIN:
        compute_min_heights(buildings);
        return;

        case Extrusion_type::AVERAGE:
        compute_average_heights(buildings);
        return;

        case Extrusion_type::MAX:
        compute_max_heights(buildings);
        return;

        default:
        compute_max_heights(buildings);
        return;
      }
    }

  private:
    const std::vector<Iterators>& m_clusters;
    const Point_map m_point_map;
    
    void compute_min_heights(std::vector<Building>& buildings) const {

      for (std::size_t i = 0; i < buildings.size(); ++i) {
        const Iterators& cluster = m_clusters[i];

        FT height = internal::max_value<FT>();
        for (std::size_t j = 0; j < cluster.size(); ++j)
          height = CGAL::min(
            height, get(m_point_map, *(cluster[j])).z());
        buildings[i].height = height;
      }
    }

    void compute_average_heights(std::vector<Building>& buildings) const {

      for (std::size_t i = 0; i < buildings.size(); ++i) {
        const Iterators& cluster = m_clusters[i];

        FT height = FT(0);
        for (std::size_t j = 0; j < cluster.size(); ++j)
          height += get(m_point_map, *(cluster[j])).z();
        buildings[i].height = height / static_cast<FT>(cluster.size());
      }
    }

    void compute_max_heights(std::vector<Building>& buildings) const {

      for (std::size_t i = 0; i < buildings.size(); ++i) {
        const Iterators& cluster = m_clusters[i];

        FT height = -internal::max_value<FT>();
        for (std::size_t j = 0; j < cluster.size(); ++j)
          height = CGAL::max(
            height, get(m_point_map, *(cluster[j])).z());
        buildings[i].height = height;
      }
    }

  }; // Building_height_estimator
    
} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_BUILDING_HEIGHT_ESTIMATOR_H
