#ifndef CGAL_LEVELS_OF_DETAIL_BUILDINGS_H
#define CGAL_LEVELS_OF_DETAIL_BUILDINGS_H

// STL includes.
#include <algorithm>

// Boost includes.
#include <boost/iterator/transform_iterator.hpp>

// LOD includes.
#include <CGAL/Levels_of_detail/internal/Utilities.h>

#include <CGAL/Levels_of_detail/internal/Filtering/Grid_based_filtering.h>
#include <CGAL/Levels_of_detail/internal/Filtering/Alpha_shapes_filtering.h>

#include <CGAL/Levels_of_detail/internal/Tools/Kd_tree_creator.h>
#include <CGAL/Levels_of_detail/internal/Estimations/Tree_based_lines_estimator.h>

namespace CGAL {

namespace Levels_of_detail {

namespace internal {

  template<class DataStructure>
  class Buildings {

  public:
    using Data_structure = DataStructure;

    using Traits = typename DataStructure::Traits;
    
    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;

    using Grid_based_filtering = Grid_based_filtering<Traits>;
    using Alpha_shapes_filtering = Alpha_shapes_filtering<Traits>;

    using Tree = Kd_tree_creator<Traits>;
    using Estimator = Tree_based_lines_estimator<Traits, Tree>;

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

      extract_building_boundary_points_2(
        alpha_shape_size, 
        grid_cell_width);

      extract_building_wall_points_2(
        region_growing_scale,
        region_growing_noise_level,
        region_growing_normal_threshold,
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
                  internal::position_on_plane(
                    m_data.building_boundary_points_2[i], 
                    m_data.ground_plane), -1);

        for(std::size_t i = 0; i < m_data.building_wall_points_2.size(); ++i)
          for(std::size_t j = 0; j < m_data.building_wall_points_2[i].size(); ++j)
            data[m_data.building_wall_points_2[i][j]].second = i;

        std::copy(data.begin(), data.end(), output);
    }

  private:
    Data_structure &m_data;

    // Building boundaries.
    void extract_building_boundary_points_2(
      const FT alpha_shape_size, 
      const FT grid_cell_width) {

      m_data.building_boundary_points_2.clear();
      const FT sampling = grid_cell_width;

      apply_alpha_shapes_filtering(alpha_shape_size, sampling);
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
        sampling, m_data.building_boundary_points_2);
    }

    void apply_grid_based_filtering(const FT grid_cell_width) {

        if (m_data.verbose) 
        std::cout << "* grid-based filtering" 
        << std::endl;

        const Grid_based_filtering filtering(grid_cell_width);
				filtering.apply(m_data.building_boundary_points_2);
    }

    // Building walls.
    void extract_building_wall_points_2(
      const FT region_growing_scale,
      const FT region_growing_noise_level,
      const FT region_growing_normal_threshold,
      const FT region_growing_minimum_length) {

      if (m_data.verbose) 
        std::cout << "* region growing" 
        << std::endl;

      m_data.building_wall_points_2.clear();
      
      Tree tree(
        m_data.building_boundary_points_2, 
        region_growing_scale);

      Estimator estimator(
        m_data.building_boundary_points_2, tree);

      
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_BUILDINGS_H
