#ifndef CGAL_LEVELS_OF_DETAIL_GROUND_H
#define CGAL_LEVELS_OF_DETAIL_GROUND_H

// LOD includes.
#include <CGAL/Levels_of_detail/internal/Utilities.h>

namespace CGAL {

namespace Levels_of_detail {

namespace internal {

  template<class DataStructure>
  class Ground {

  public:
    using Data_structure = DataStructure;

    using Traits = typename DataStructure::Traits;
    using Point_3 = typename Traits::Point_3;

    Ground(Data_structure &data_structure) :
    m_data(data_structure)
    { }

    void make_planar() {

      if (m_data.verbose) 
        std::cout << "* computing planar ground" << std::endl;

      internal::fit_plane_to_points_3(
        m_data.ground_points(), m_data.point_map, 
        m_data.ground_plane);

      internal::compute_bounding_box_3(
        m_data.ground_points(), m_data.point_map, m_data.ground_plane, 
        m_data.planar_ground);
    }

    void make_smooth() {

    }

    template<typename VerticesOutputIterator>
    void return_as_polygon(VerticesOutputIterator vertices) const {

      CGAL_precondition(!m_data.planar_ground.empty());
      for (size_t i = 0; i < m_data.planar_ground.size(); ++i)
        *(vertices++) = m_data.planar_ground[i];
    }

    template<typename VerticesOutputIterator, typename FacesOutputIterator>
    void return_as_triangle_soup(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

    }

  private:
    Data_structure &m_data;
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_GROUND_H
