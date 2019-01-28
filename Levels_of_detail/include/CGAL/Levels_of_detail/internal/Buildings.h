#ifndef CGAL_LEVELS_OF_DETAIL_BUILDINGS_H
#define CGAL_LEVELS_OF_DETAIL_BUILDINGS_H

// STL includes.
#include <vector>
#include <algorithm>

// Boost includes.
#include <boost/iterator/transform_iterator.hpp>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utilities.h>

#include <CGAL/Levels_of_detail/internal/Simplification/Thinning_2.h>
#include <CGAL/Levels_of_detail/internal/Simplification/Grid_based_filtering_2.h>
#include <CGAL/Levels_of_detail/internal/Simplification/Alpha_shapes_filtering_2.h>

#include <CGAL/Levels_of_detail/internal/Shape_detection/Region_growing.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Estimate_normals_2.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Points_k_nearest_neighbors_connectivity.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Points_2_least_squares_line_fit_conditions.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename DataStructure>
  class Buildings {

  public:
    using Data_structure = DataStructure;

    using Traits = typename DataStructure::Traits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;

    using Grid_based_filtering_2 = Grid_based_filtering_2<Traits>;
    using Alpha_shapes_filtering_2 = Alpha_shapes_filtering_2<Traits>;

    using Points_2 = std::vector<Point_2>;
    using Connectivity_2 = Points_k_nearest_neighbors_connectivity<Traits, Point_2>;
    using Normals_estimator_2 = Estimate_normals_2<Traits, Connectivity_2>;
    using Conditions_2 = Points_2_least_squares_line_fit_conditions<Traits>;
    using Region_growing_2 = Region_growing<Points_2, Connectivity_2, Conditions_2>;
    
    Buildings(Data_structure& data_structure) :
    m_data(data_structure)
    { }
    
    void detect_building_boundaries(
      const FT alpha_shape_size,
      const FT grid_cell_width,
      const FT region_growing_search_size,
      const FT region_growing_noise_level,
      const FT region_growing_angle,
      const FT region_growing_minimum_length) {

      if (m_data.verbose) 
        std::cout << "- Computing building boundaries" 
        << std::endl;

      extract_boundary_points_2(
        alpha_shape_size, 
        grid_cell_width);

      extract_wall_points_2(
        region_growing_search_size,
        region_growing_noise_level,
        region_growing_angle,
        region_growing_minimum_length);
    }

    template<typename OutputIterator>
    void return_boundary_points(OutputIterator output) const {

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

    template<typename OutputIterator>
    void return_wall_points(OutputIterator output) const {

      std::vector<std::pair<Point_3, int> > data(
        m_data.building_boundary_points_2.size());

      for (std::size_t i = 0; i < m_data.building_boundary_points_2.size(); ++i)
        data[i] = std::make_pair(
                  internal::position_on_plane_3(
                    m_data.building_boundary_points_2[i], 
                    m_data.ground_plane), -1);

        for(std::size_t i = 0; i < m_data.building_wall_points_2.size(); ++i)
          for(std::size_t j = 0; j < m_data.building_wall_points_2[i].size(); ++j)
            data[m_data.building_wall_points_2[i][j]].second = i;

        std::copy(data.begin(), data.end(), output);
    }

    template<typename OutputIterator>
    void return_boundary_edges(OutputIterator output) const {

      std::copy(
        boost::make_transform_iterator(
          m_data.building_wall_points_2.begin(),
          internal::Segment_3_from_points_and_plane<Traits>(
            m_data.building_boundary_points_2, m_data.ground_plane)),
        boost::make_transform_iterator(
          m_data.building_wall_points_2.end(),
          internal::Segment_3_from_points_and_plane<Traits>(
            m_data.building_boundary_points_2, m_data.ground_plane)),
        output);
    }

  private:
    Data_structure& m_data;

    // Building boundaries.
    void extract_boundary_points_2(
      const FT alpha_shape_size, 
      const FT grid_cell_width) {

      m_data.building_boundary_points_2.clear();
      const FT sampling = grid_cell_width;

      apply_alpha_shapes_filtering_2(alpha_shape_size, sampling);
      apply_grid_based_filtering_2(grid_cell_width);

      if (m_data.verbose)
        std::cout << "-> " << m_data.building_boundary_points_2.size()
        << " boundary point(s) extracted" << std::endl;
    }

    void apply_alpha_shapes_filtering_2(
      const FT alpha_shape_size, 
      const FT sampling) {

      if (m_data.verbose) 
        std::cout << "* alpha shapes filtering" 
        << std::endl;

      const std::size_t numb_pts = m_data.building_boundary_points().size();
      const std::size_t numi_pts = m_data.building_interior_points().size();

      CGAL_precondition(numb_pts >= 3 || numi_pts >= 3);
      Alpha_shapes_filtering_2 filtering(alpha_shape_size);

      if (numb_pts >= 3)
        filtering.add_points(
          m_data.building_boundary_points(),
          m_data.point_map);

      if (numi_pts >= 3)
        filtering.add_points(
          m_data.building_interior_points(),
          m_data.point_map);

      filtering.get_filtered_points(
        sampling, m_data.building_boundary_points_2);
    }

    void apply_grid_based_filtering_2(const FT grid_cell_width) {

        if (m_data.verbose) 
        std::cout << "* grid-based filtering" 
        << std::endl;

        const Grid_based_filtering_2 filtering(grid_cell_width);
				filtering.apply(m_data.building_boundary_points_2);
    }

    // Building walls.
    void extract_wall_points_2(
      const FT region_growing_search_size,
      const FT region_growing_noise_level,
      const FT region_growing_angle,
      const FT region_growing_minimum_length) {

      if (m_data.verbose) 
        std::cout << "* region growing" 
        << std::endl;

      m_data.building_wall_points_2.clear();
      const auto& points = m_data.building_boundary_points_2;

      Connectivity_2 connectivity(
        points, 
        region_growing_search_size);

      Normals_estimator_2 estimator(
        points, 
        connectivity);

      Conditions_2 conditions(
        points, 
        estimator.normals(),
        region_growing_noise_level,
        region_growing_angle,
        region_growing_minimum_length);

      std::vector<std::size_t> indices(points.size());
      for (std::size_t i = 0; i < points.size(); ++i)
        indices[i] = i;
      std::stable_sort(
        indices.begin(), indices.end(), estimator.sorter());

      Region_growing_2 region_growing(
        indices,
        connectivity,
        conditions);

      region_growing.detect(
        m_data.building_wall_points_2);

      if (m_data.verbose) 
        std::cout << "-> " << m_data.building_wall_points_2.size()
        << " wall(s) extracted" 
        << std::endl;
    }

  }; // Buildings

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_BUILDINGS_H
