#ifndef CGAL_LEVEL_OF_DETAIL_POINTS_BASED_REGION_GROWING_2_H
#define CGAL_LEVEL_OF_DETAIL_POINTS_BASED_REGION_GROWING_2_H

// STL includes.
#include <map>
#include <list>

// CGAL includes.
#include <CGAL/number_utils.h>

// LOD includes.
#include <CGAL/Level_of_detail/Tools/Fitting/Line_to_points_fitter.h>

namespace CGAL {

	namespace Level_of_detail {

		template<class InputKernel, class TreeWrapper>
		class Points_based_region_growing_2 {

		public:
		    using Kernel = InputKernel;
            using Tree   = TreeWrapper;

			using FT 	   = typename Kernel::FT;
			using Line_2   = typename Kernel::Line_2;
			using Point_2  = typename Kernel::Point_2;
			using Vector_2 = typename Kernel::Vector_2;

			using Point_identifier          = typename Tree::Element_identifier;
			using Neighbours 	            = typename Tree::Neighbours;
			using Const_neighbours_iterator = typename Neighbours::const_iterator;

			using Shape_index 			= std::map<Point_identifier, int>;
			using Line_to_points_fitter = LOD::Line_to_points_fitter<Kernel>;

			typename Kernel::Compute_squared_distance_2 squared_distance_2;

			Points_based_region_growing_2(const FT epsilon, const FT cluster_epsilon, const FT normal_threshold, const size_t min_points, const Tree &tree) :
			m_epsilon(epsilon),
			m_cluster_epsilon(cluster_epsilon),
			m_normal_threshold(normal_threshold),
			m_min_points(min_points), 
			m_tree(tree)
			{ }

			template<class Elements, class Point_map, class Normal_map, class Index_container>
			void detect(const Elements &elements, const Point_map &point_map, const Normal_map &normal_map, std::list<Index_container> &output) const {
				
				using Const_elements_iterator = typename Elements::const_iterator;
				using Const_index_iterator 	  = typename Index_container::const_iterator;

				output.clear();
				CGAL_precondition(elements.size() > 1);

				// Main structures.
				Index_container index_container;
				size_t available_points = elements.size();

				Shape_index shape_index;
				for (Const_elements_iterator ce_it = elements.begin(); ce_it != elements.end(); ++ce_it) shape_index[*ce_it] = -1;

				// Main loop.
				int class_index = -1;
      			for (Const_elements_iterator ce_it = elements.begin(); ce_it != elements.end(); ++ce_it) {
					
					const Point_identifier &point_id = *ce_it;
					if (shape_index.at(point_id) >= 0) continue;

					// Get query point and its normal.
					const Point_2  &position = get(point_map, point_id);
					const Vector_2 &normal   = get(normal_map, point_id);

					// Get optimal line and its normal.
					Line_2   optimal_line; 
					Vector_2 optimal_normal;
					get_line_with_normal(position, normal, optimal_line, optimal_normal);

					// Initialize containers.
					shape_index[point_id] = ++class_index;
					index_container.clear();
					index_container.push_back(point_id);

					// Propagate through all neighbouring points.
					propagate(point_id, point_map, normal_map, class_index, 
					shape_index, index_container, optimal_line, optimal_normal);

					// Try to create a new region.
					if (index_container.size() >= m_min_points) { // if the global condition is satisfied

						output.push_back(index_container);
      			 		available_points -= index_container.size();

      			 	} else { // otherwise, cancel region

						--class_index;
          		 		shape_index[point_id] = -1;
          		 		for (Const_index_iterator ic_it = index_container.begin(); ic_it != index_container.end(); ++ic_it) shape_index[*ic_it] = -1;
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

			template<class Point_map, class Normal_map, class Index_container>
			void propagate(const Point_identifier &point_id, const Point_map &point_map, const Normal_map &normal_map, const int class_index,
			Shape_index &shape_index, Index_container &index_container, Line_2 &optimal_line, Vector_2 &optimal_normal) const {

				using Const_index_iterator = typename Index_container::const_iterator;

				// Initialize containers.
				Neighbours 		neighbours;
				Index_container current_ring;
        		
				Index_container former_ring;
				former_ring.push_back(point_id);

				// Main loop.
				bool propagation = true;
				do {
					
					propagation = false;
					for (Const_index_iterator fr_it = former_ring.begin(); fr_it != former_ring.end(); ++fr_it) {

						// Find neighbours of the current point.
						const Point_2 &query = get(point_map, *fr_it);
						m_tree.search_2(query, neighbours);

						// Check neighbours if they are valid for propagation.
						handle_neighbours(neighbours, point_map, normal_map, class_index, optimal_line, optimal_normal, 
						current_ring, shape_index, propagation);
					}

					// Update containers with current data.
					update_containers(current_ring, former_ring, index_container);

					// Update conditions: refit the line.
					if (index_container.size() < 2) continue;
					refit_the_line(point_map, index_container, optimal_line, optimal_normal);

				} while (propagation);
			}

			template<class Point_map, class Normal_map, class Index_container>
			void handle_neighbours(
			const Neighbours &neighbours, const Point_map &point_map, const Normal_map &normal_map, 
			const int class_index, const Line_2 &optimal_line, const Vector_2 &optimal_normal,
			Index_container &current_ring, Shape_index &shape_index, bool &propagation) const {

				CGAL_precondition(neighbours.size() > 0);
				for (Const_neighbours_iterator cn_it = neighbours.begin(); cn_it != neighbours.end(); ++cn_it) {
				 	
					const Point_identifier &neighbour_id = get(m_tree.point_identifier_map(), *cn_it);
					if (shape_index.at(neighbour_id) >= 0) continue;

					// Get neighbour's position and normal.
					const Point_2  &position = get(point_map, neighbour_id);
					const Vector_2 &normal   = get(normal_map, neighbour_id);

					// Compute distance to the line and cos.
					const FT squared_distance = squared_distance_2(position, optimal_line);
					const FT cos_angle 		  = CGAL::abs(normal * optimal_normal);

					// Check local conditions.
					CGAL_precondition(m_epsilon > FT(0) && m_normal_threshold > FT(0) && m_normal_threshold < FT(1));
					if (squared_distance > m_epsilon * m_epsilon || cos_angle < m_normal_threshold) continue; // not satisfied

					// If satisfied, add point to the region.
					shape_index[neighbour_id] = class_index;
              		current_ring.push_back(neighbour_id);
              		propagation = true;
				}
			}

			template<class Index_container>
			void update_containers(Index_container &current_ring, Index_container &former_ring, Index_container &index_container) const {
				using Const_index_iterator = typename Index_container::const_iterator;

				former_ring.clear();
          		for (Const_index_iterator cr_it = current_ring.begin(); cr_it != current_ring.end(); ++cr_it) {
	            	
	            	former_ring.push_back(*cr_it);
	            	index_container.push_back(*cr_it);
	          	}
          		current_ring.clear();
			}

			template<class Point_map, class Index_container>
			void refit_the_line(const Point_map &point_map, const Index_container &index_container, Line_2 &optimal_line, Vector_2 &optimal_normal) const {

				const Line_to_points_fitter line_to_points_fitter;
				line_to_points_fitter.fit_line_2(index_container, point_map, optimal_line);
				
				optimal_normal  = (optimal_line.perpendicular(optimal_line.point(0))).to_vector();
				optimal_normal /= static_cast<FT>(CGAL::sqrt(CGAL::to_double(optimal_normal * optimal_normal)));
			}
		};

	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_POINTS_BASED_REGION_GROWING_2_H