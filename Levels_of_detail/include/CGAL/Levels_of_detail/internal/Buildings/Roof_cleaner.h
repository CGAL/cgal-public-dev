#ifndef CGAL_LEVELS_OF_DETAIL_ROOF_CLEANER_H
#define CGAL_LEVELS_OF_DETAIL_ROOF_CLEANER_H

// STL includes.
#include <vector>
#include <algorithm>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/number_utils.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utilities.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits, 
  typename InputRange,
  typename PointMap>
  class Roof_cleaner {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;
    using Vector_3 = typename Traits::Vector_3;
    using Line_3 = typename Traits::Line_3;

    using Search_traits = CGAL::Search_traits_3<Traits>;
		using Neighbor_search = CGAL::Orthogonal_k_neighbor_search<Search_traits>;
		using Fuzzy_sphere = CGAL::Fuzzy_sphere<Search_traits>;
		using Tree = typename Neighbor_search::Tree;

    Roof_cleaner(
      const Input_range& input_range, 
      const Point_map point_map,
      const std::vector<Vector_3>& normals,
      const FT scale) :
    m_input_range(input_range),
    m_point_map(point_map),
    m_normals(normals),
    m_max_percentage(FT(90)),
    m_angle_threshold(FT(15)),
    m_scale(scale),
    m_distance_threshold(m_scale / FT(2))
    { }

    void clean(std::vector< std::vector<std::size_t> >& roofs) const {

      std::vector<std::size_t> indices;
      set_default_indices(roofs.size(), indices);

      apply_size_criteria(roofs, indices);
      apply_scale_based_criteria(roofs, indices);
      apply_vertical_criteria(roofs, indices);
      apply_thin_criteria(roofs, indices);

      update(indices, roofs);
    }

  private:
    const Input_range& m_input_range;
    const Point_map m_point_map;
    const std::vector<Vector_3>& m_normals;

    const FT m_scale;
    const FT m_max_percentage;
    const FT m_angle_threshold;
    const FT m_distance_threshold;

    class Compare_size {
                
    public:
      Compare_size(
        const std::vector< std::vector<std::size_t> >& roofs) : 
      m_roofs(roofs) 
      { }
                
      bool operator()(const std::size_t i, const size_t j) const { 

        CGAL_precondition(i >= 0 && i < m_roofs.size());
        CGAL_precondition(j >= 0 && j < m_roofs.size());

        return m_roofs[i].size() > m_roofs[j].size();
      }

    private:
      const std::vector< std::vector<std::size_t> >& m_roofs;
    };

    void set_default_indices(
      const std::size_t num_indices,
      std::vector<std::size_t>& indices) const {
                
      indices.clear();
      indices.resize(num_indices);

      for (std::size_t i = 0; i < num_indices; ++i)
        indices[i] = i;
    }

    void apply_size_criteria(
      const std::vector< std::vector<std::size_t> >& roofs, 
      std::vector<std::size_t>& indices) const {

      CGAL_precondition(
        m_max_percentage >= FT(0) && m_max_percentage <= FT(100));

      sort_indices_by_size(roofs, indices);
      remove_outliers(roofs, indices);
    }

    void sort_indices_by_size(
      const std::vector< std::vector<std::size_t> >& roofs, 
      std::vector<std::size_t>& indices) const {

      Compare_size cmp(roofs);
      std::sort(indices.begin(), indices.end(), cmp);
    }

    void remove_outliers(
      const std::vector< std::vector<std::size_t> >& roofs, 
      std::vector<std::size_t>& indices) const {
                
      const std::size_t num_total_points = 
      get_total_number_of_points(roofs, indices);

      const FT scale = m_max_percentage / FT(100);

      const std::size_t num_points_to_keep = 
        static_cast<std::size_t>(
          std::ceil(
            CGAL::to_double(
              scale * static_cast<FT>(num_total_points))));
                
      std::size_t curr_num_points = 0;
      for (std::size_t i = 0; i < indices.size(); ++i) {
                    
        const std::size_t index = indices[i];
        curr_num_points += roofs[index].size();

        if (curr_num_points >= num_points_to_keep) {
                        
          indices.erase(indices.begin() + i + 1, indices.end());
          break;
        }
      }
    }

    std::size_t get_total_number_of_points(
      const std::vector< std::vector<std::size_t> >& roofs, 
      const std::vector<std::size_t>& indices) const {

      std::size_t num_total_points = 0;
      for (std::size_t i = 0; i < indices.size(); ++i)
        num_total_points += roofs[indices[i]].size();

      return num_total_points;
    }

    void apply_vertical_criteria(
      const std::vector< std::vector<std::size_t> >& roofs, 
      std::vector<std::size_t>& indices) const {

      std::vector<std::size_t> new_indices;
      for (std::size_t i = 0; i < indices.size(); ++i) {

        const std::size_t index = indices[i];
        if (!is_vertical_shape(roofs[index])) 
            new_indices.push_back(index);
      }
      indices = new_indices;
    }

    bool is_vertical_shape(const std::vector<std::size_t>& indices) const {

      Vector_3 m; internal::average_vector_3(m_normals, indices, m);
		  const Vector_3 n = Vector_3(FT(0), FT(0), FT(1));

      const FT angle_deg = internal::angle_deg_3(m, n);
      const FT angle_dif = CGAL::abs(FT(90) - CGAL::abs(angle_deg));

      if (angle_dif < m_angle_threshold) 
        return true;

      return false;
    }

    void apply_scale_based_criteria(
      const std::vector< std::vector<std::size_t> >& roofs, 
      std::vector<std::size_t>& indices) const {

      std::vector<std::size_t> new_indices;
      for (std::size_t i = 0; i < indices.size(); ++i) {

        const std::size_t index = indices[i];
        if (!is_within_scale_bounds(roofs[index])) 
          new_indices.push_back(index);
      }
      indices = new_indices;
    }

    bool is_within_scale_bounds(
      const std::vector<std::size_t>& shape_indices) const {
                
      CGAL_precondition(m_scale > FT(0));

      Point_3 barycentre;
      internal::barycenter_3(
        m_input_range, m_point_map, shape_indices, barycentre);

      std::vector<Point_3> points(shape_indices.size());
      for (std::size_t i = 0; i < shape_indices.size(); ++i) 
        points[i] = 
          get(m_point_map, *(m_input_range.begin() + shape_indices[i]));

      Tree tree(points.begin(), points.end());
      const Fuzzy_sphere sphere(barycentre, m_scale);

      std::vector<Point_3> result;
      tree.search(std::back_inserter(result), sphere);

      if (result.size() == points.size()) 
        return true;
      
      return false;
    }

    void apply_thin_criteria(
      const std::vector< std::vector<std::size_t> >& roofs, 
      std::vector<std::size_t>& indices) const {

      std::vector<std::size_t> new_indices;
      for (std::size_t i = 0; i < indices.size(); ++i) {

        const std::size_t index = indices[i];
        if (!is_thin(roofs[index])) 
          new_indices.push_back(index);
      }
      indices = new_indices;
    }

    bool is_thin(const std::vector<std::size_t>& shape_indices) const {

      Line_3 line; internal::line_from_points_3(
        m_input_range, m_point_map, shape_indices, line);

      const FT average_distance = internal::average_distance_to_line_3(
        m_input_range, m_point_map, shape_indices, line);

      return average_distance < m_distance_threshold;
    }

    void update(
      const std::vector<std::size_t>& indices, 
      std::vector< std::vector<std::size_t> >& roofs) const {

      std::vector< std::vector<std::size_t> > new_roofs;
      for (std::size_t i = 0; i < indices.size(); ++i)
        if (roofs[indices[i]].size() > 2)
          new_roofs.push_back(roofs[indices[i]]);

      roofs = new_roofs;
    }

  }; // Roof_cleaner

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_ROOF_CLEANER_H
