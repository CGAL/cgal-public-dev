#ifndef CGAL_LEVELS_OF_DETAIL_BUILDINGS_H
#define CGAL_LEVELS_OF_DETAIL_BUILDINGS_H

// LOD includes.
#include <CGAL/Levels_of_detail/internal/Filtering/Alpha_shapes_filtering.h>

namespace CGAL {

namespace Levels_of_detail {

namespace internal {

  template<class DataStructure>
  class Buildings {

  public:
    using Data_structure = DataStructure;

    using Traits = typename DataStructure::Traits;
    using FT = typename Traits::FT;

    using Alpha_shapes_filtering = Alpha_shapes_filtering<Traits>;

    Buildings(Data_structure &data_structure) :
    m_data(data_structure)
    { }

    void detect_building_boundaries(
      const FT alpha_shape_size,
      const FT grid_cell_width,
      const FT region_growing_scale,
      const FT region_growing_noise_level,
      const FT region_growing_normal_threshold,
      const FT region_growing_minimum_length) {

      if (m_data.verbose) 
        std::cout << "- Computing building boundaries" 
        << std::endl;

      extract_building_boundary_points(alpha_shape_size, grid_cell_width);
    }

    template<typename OutputIterator>
    void return_building_boundary_points(OutputIterator output) const {

      CGAL_precondition(!m_data.building_boundary_points_2.empty());
      std::copy(
        boost::make_transform_iterator(
          m_data.building_boundary_points_2.begin(),
          internal::Point_3_from_point_2_and_plane<Traits>(m_data.ground_plane)),
        boost::make_transform_iterator(
          m_data.building_boundary_points_2.end(),
          internal::Point_3_from_point_2_and_plane<Traits>(m_data.ground_plane)),
        output);
    }

  private:
    Data_structure &m_data;

    void extract_building_boundary_points(
      const FT alpha_shape_size, 
      const FT grid_cell_width) {

      m_data.building_boundary_points_2.clear();

      apply_alpha_shapes_filtering(alpha_shape_size, grid_cell_width);
      apply_grid_based_filtering(grid_cell_width);

      if (m_data.verbose)
        std::cout << "-> " << m_data.building_boundary_points_2.size()
        << " boundary point(s) extracted" << std::endl;
    }

    void apply_alpha_shapes_filtering(
      const FT alpha_shape_size, 
      const FT sampling) {

      if (m_data.verbose) 
        std::cout << "* alpha shapes filtering" 
        << std::endl;

      const std::size_t numb_pts = m_data.building_boundary_points().size();
      const std::size_t numi_pts = m_data.building_interior_points().size();

      CGAL_precondition(numb_pts >= 3 || numi_pts >= 3);
      Alpha_shapes_filtering filtering(alpha_shape_size);

      if (numb_pts >= 3)
        filtering.add_points(
          m_data.building_boundary_points(),
          m_data.point_map);

      if (numi_pts >= 3)
        filtering.add_points(
          m_data.building_interior_points(),
          m_data.point_map);

      filtering.get_filtered_points(
        m_data.building_boundary_points_2, sampling);
    }

    void apply_grid_based_filtering(const FT grid_cell_width) {

        if (m_data.verbose) 
        std::cout << "* grid-based filtering" 
        << std::endl;


    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_BUILDINGS_H
