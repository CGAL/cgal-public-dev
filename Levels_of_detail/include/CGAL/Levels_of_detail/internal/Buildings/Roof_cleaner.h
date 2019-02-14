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

	  typename Traits::Compute_squared_length_3 squared_length;
		typename Traits::Compute_scalar_product_3 dot_product;
		typename Traits::Construct_cross_product_vector_3 cross_product;

    using Search_traits = CGAL::Search_traits_3<Traits>;
		using Neighbor_search = CGAL::Orthogonal_k_neighbor_search<Search_traits>;
		using Fuzzy_sphere = CGAL::Fuzzy_sphere<Search_traits>;
		using Tree = typename Neighbor_search::Tree;

    Roof_cleaner(
      const Input_range& input_range, 
      const Point_map point_map,
      const std::vector<Vector_3>& normals,
      const FT min_size) :
    m_input_range(input_range),
    m_point_map(point_map),
    m_normals(normals),
    m_max_percentage(FT(80)),
    m_angle_threshold(FT(25)),
    m_min_size(min_size)
    { }

    void clean(std::vector< std::vector<std::size_t> >& roofs) const {

      std::vector<std::size_t> indices;
      set_default_indices(roofs.size(), indices);

      apply_size_criteria(roofs, indices);
      apply_scale_based_criteria(roofs, indices);
      apply_vertical_criteria(roofs, indices);
      
      // apply_thin_criteria(roofs, indices);

      update(indices, roofs);
    }

  private:
    const Input_range& m_input_range;
    const Point_map m_point_map;
    const std::vector<Vector_3>& m_normals;

    const FT m_min_size;
    const FT m_max_percentage;
    const FT m_angle_threshold;

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

    bool is_vertical_shape(
      const std::vector<std::size_t>& shape_indices) const {

			Vector_3 shape_normal;
		  set_shape_normal(shape_indices, shape_normal);

			Vector_3 ground_normal;
			set_ground_normal(ground_normal);

      const FT angle = compute_angle(shape_normal, ground_normal);
      const FT angle_diff = CGAL::abs(FT(90) - CGAL::abs(angle));

      if (angle_diff < m_angle_threshold) 
        return true;

      return false;
    }

    void set_shape_normal(
      const std::vector<std::size_t>& shape_indices, 
      Vector_3& m) const {
				
      FT x = FT(0), y = FT(0), z = FT(0);
      for (std::size_t i = 0; i < shape_indices.size(); ++i) {
                    
        const std::size_t index = shape_indices[i];
        const Vector_3& normal = m_normals[index];

        x += normal.x();
        y += normal.y();
        z += normal.z();
      }

      x /= static_cast<FT>(shape_indices.size());
      y /= static_cast<FT>(shape_indices.size());
      z /= static_cast<FT>(shape_indices.size());

      m = Vector_3(x, y, z);
		}

		void set_ground_normal(Vector_3& n) const {
			n = Vector_3(FT(0), FT(0), FT(1));
		}

    FT compute_angle(
      const Vector_3& m, 
      const Vector_3& n) const {
				
			const auto cross = cross_product(m, n);
			
      const FT length = 
      static_cast<FT>(
        CGAL::sqrt(
          CGAL::to_double(
            squared_length(cross))));

			const FT dot = dot_product(m, n);

			FT angle_rad = 
      static_cast<FT>(
        std::atan2(
          CGAL::to_double(length), 
          CGAL::to_double(dot)));
                
      const FT half_pi = 
      static_cast<FT>(CGAL_PI) / FT(2);
      
      if (angle_rad > half_pi) 
        angle_rad = static_cast<FT>(CGAL_PI) - angle_rad;

			const FT angle_deg = 
        angle_rad * FT(180) / static_cast<FT>(CGAL_PI);
      
      return angle_deg;
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
                
      CGAL_precondition(m_min_size > FT(0));

      Point_3 barycentre;
      compute_barycentre(shape_indices, barycentre);

      std::vector<Point_3> points;
      set_points(shape_indices, points);

      Tree tree(points.begin(), points.end());
      const Fuzzy_sphere sphere(barycentre, m_min_size);

      std::vector<Point_3> result;
      tree.search(std::back_inserter(result), sphere);

      if (result.size() == points.size()) 
        return true;
      
      return false;
    }

    void compute_barycentre(
      const std::vector<std::size_t>& shape_indices, 
      Point_3& barycentre) const {
                
      FT x = FT(0), y = FT(0), z = FT(0);
      for (std::size_t i = 0; i < shape_indices.size(); ++i) {
        
        const Point_3& point = 
          get(m_point_map, *(m_input_range.begin() + shape_indices[i]));

        x += point.x();
        y += point.y();
        z += point.z();
      }

      x /= static_cast<FT>(shape_indices.size());
      y /= static_cast<FT>(shape_indices.size());
      z /= static_cast<FT>(shape_indices.size());

      barycentre = Point_3(x, y, z);
    }

    void set_points(
      const std::vector<std::size_t>& shape_indices, 
      std::vector<Point_3>& points) const {

      points.clear();
      points.resize(shape_indices.size());

      for (std::size_t i = 0; i < shape_indices.size(); ++i) 
        points[i] = 
          get(m_point_map, *(m_input_range.begin() + shape_indices[i]));
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
