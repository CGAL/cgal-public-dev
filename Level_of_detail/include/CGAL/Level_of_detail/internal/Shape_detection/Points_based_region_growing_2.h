#ifndef CGAL_LEVEL_OF_DETAIL_POINTS_BASED_REGION_GROWING_2_H
#define CGAL_LEVEL_OF_DETAIL_POINTS_BASED_REGION_GROWING_2_H

// STL includes.
#include <map>
#include <list>

// CGAL includes.
#include <CGAL/number_utils.h>

// LOD includes.
#include <CGAL/Level_of_detail/internal/Fitting/Line_to_points_fitter.h>

namespace CGAL {

namespace Level_of_detail {

template<class GeomTraits, class Tree>
class Points_based_region_growing_2 {

public:
  using Kernel = GeomTraits;

  using FT 	   = typename Kernel::FT;
  using Line_2   = typename Kernel::Line_2;
  using Point_2  = typename Kernel::Point_2;
  using Vector_2 = typename Kernel::Vector_2;

  using Shape_index 			= std::vector<int>;
  using Line_to_points_fitter = Line_to_points_fitter<Kernel>;

  typename Kernel::Compute_squared_distance_2 squared_distance_2;

  Points_based_region_growing_2(const FT epsilon, const FT cluster_epsilon,
                                const FT normal_threshold, const size_t min_points, const Tree &tree) :
    m_epsilon(epsilon),
    m_cluster_epsilon(cluster_epsilon),
    m_normal_threshold(normal_threshold),
    m_min_points(min_points), 
    m_tree(tree)
  { }

  void detect(const std::vector<std::size_t>& indices,
              const std::vector<Point_2>& points,
              const std::vector<Vector_2>& normals,
              std::vector<std::vector<std::size_t> > &output) const {

    output.clear();
    CGAL_precondition(indices.size() > 1);

    // Main structures.
    std::vector<std::size_t> index_container;
    size_t available_points = points.size();

    Shape_index shape_index (points.size(), -1);

    // Main loop.
    int class_index = -1;
    for (std::size_t i = 0; i < indices.size(); ++ i) {

      std::size_t idx = indices[i];
      if (shape_index[idx] >= 0) continue;

      // Get query point and its normal.
      const Point_2  &position = points[idx];
      const Vector_2 &normal   = normals[idx];

      // Get optimal line and its normal.
      Line_2   optimal_line; 
      Vector_2 optimal_normal;
      get_line_with_normal(position, normal, optimal_line, optimal_normal);

      // Initialize containers.
      shape_index[idx] = ++class_index;
      index_container.clear();
      index_container.push_back(idx);

      // Propagate through all neighbouring points.
      propagate(idx, points, normals, class_index, 
                shape_index, index_container, optimal_line, optimal_normal);

      // Try to create a new region.
      if (index_container.size() >= m_min_points) { // if the global condition is satisfied

        output.push_back(index_container);
        available_points -= index_container.size();

      } else { // otherwise, cancel region

        --class_index;
        shape_index[idx] = -1;
        for (std::size_t j = 0; j < index_container.size(); ++ j)
          shape_index[index_container[j]] = -1;
      }
    }
  }

private:
  const FT 	 m_epsilon;
  const FT 	 m_cluster_epsilon;
  const FT 	 m_normal_threshold;
  const size_t m_min_points;

  const Tree &m_tree;

  void get_line_with_normal(const Point_2 &point, const Vector_2 &normal, Line_2 &output_line, Vector_2 &output_normal) const {
				
    const Line_2 tmp_line(point, normal);
    output_line = tmp_line.perpendicular(point);

    const FT length = static_cast<FT>(CGAL::sqrt(CGAL::to_double(normal * normal)));
    output_normal = normal / length;
  }

  void propagate(const std::size_t point_id,
                 const std::vector<Point_2>& points,
                 const std::vector<Vector_2>& normals,
                 const int class_index,
                 Shape_index &shape_index,
                 std::vector<std::size_t>& index_container,
                 Line_2 &optimal_line,
                 Vector_2 &optimal_normal) const {

    // Initialize containers.
    std::vector<std::size_t> neighbors;
    std::vector<std::size_t> current_ring;
        		
    std::vector<std::size_t> former_ring;
    former_ring.push_back(point_id);

    // Main loop.
    bool propagation = true;
    do {
					
      propagation = false;
      for (std::size_t i = 0; i < former_ring.size(); ++ i) {

        // Find neighbours of the current point.
        const Point_2 &query = points[former_ring[i]];
        m_tree.search_2(query, neighbors);

        // Check neighbours if they are valid for propagation.
        handle_neighbours(neighbors, points, normals, class_index, optimal_line, optimal_normal, 
                          current_ring, shape_index, propagation);
      }

      // Update containers with current data.
      update_containers(current_ring, former_ring, index_container);

      // Update conditions: refit the line.
      if (index_container.size() < 2) continue;
      refit_the_line(points, index_container, optimal_line, optimal_normal);

    } while (propagation);
  }

  void handle_neighbours(
    const std::vector<std::size_t>& neighbors,
    const std::vector<Point_2>& points,
    const std::vector<Vector_2>& normals, 
    const int class_index,
    const Line_2 &optimal_line,
    const Vector_2 &optimal_normal,
    std::vector<std::size_t>& current_ring,
    Shape_index &shape_index,
    bool &propagation) const {

    CGAL_precondition(neighbors.size() > 0);
    for (std::size_t i = 0; i < neighbors.size(); ++ i) {
				 	
      if (shape_index[neighbors[i]] >= 0) continue;

      // Get neighbour's position and normal.
      const Point_2  &position = points[neighbors[i]];
      const Vector_2 &normal   = normals[neighbors[i]];

      // Compute distance to the line and cos.
      const FT squared_distance = squared_distance_2(position, optimal_line);
      const FT cos_angle 		  = CGAL::abs(normal * optimal_normal);

      // Check local conditions.
      CGAL_precondition(m_epsilon > FT(0) && m_normal_threshold > FT(0) && m_normal_threshold < FT(1));
      if (squared_distance > m_epsilon * m_epsilon || cos_angle < m_normal_threshold) continue; // not satisfied

      // If satisfied, add point to the region.
      shape_index[neighbors[i]] = class_index;
      current_ring.push_back(neighbors[i]);
      propagation = true;
    }
  }

  void update_containers(std::vector<std::size_t>& current_ring,
                         std::vector<std::size_t>& former_ring,
                         std::vector<std::size_t>& index_container) const {

    former_ring.clear();
    for (std::size_t i = 0; i < current_ring.size(); ++ i) {
      former_ring.push_back(current_ring[i]);
      index_container.push_back(current_ring[i]);
    }
    current_ring.clear();
  }

  void refit_the_line(const std::vector<Point_2>& points,
                      const std::vector<std::size_t>& index_container,
                      Line_2 &optimal_line, Vector_2 &optimal_normal) const {

    const Line_to_points_fitter line_to_points_fitter;
    line_to_points_fitter.fit_line_2(index_container, points, optimal_line);
				
    optimal_normal  = (optimal_line.perpendicular(optimal_line.point(0))).to_vector();
    optimal_normal /= static_cast<FT>(CGAL::sqrt(CGAL::to_double(optimal_normal * optimal_normal)));
  }
};

} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_POINTS_BASED_REGION_GROWING_2_H
