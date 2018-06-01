#ifndef CGAL_LEVEL_OF_DETAIL_POINTS_BASED_REGION_GROWING_2_H
#define CGAL_LEVEL_OF_DETAIL_POINTS_BASED_REGION_GROWING_2_H

// STL includes.
#include <map>
#include <list>
#include <cmath>
#include <algorithm>
#include <unordered_set>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_2.h>

namespace CGAL {

	namespace Level_of_detail {

		template<class InputKernel, class PointIdentifier>
		class Points_based_region_growing_2 {

		public:
		    using Kernel           = InputKernel;
            using Point_identifier = PointIdentifier;

			using FT 	   = typename Kernel::FT;
			using Point_2  = typename Kernel::Point_2;
			using Vector_2 = typename Kernel::Vector_2;

			using Search_traits_2 = CGAL::Search_traits_2<Kernel>;
			using Search_circle   = CGAL::Fuzzy_sphere<Search_traits_2>;
			using Search_tree     = CGAL::Kd_tree<Search_traits_2>;

			Points_based_region_growing_2(const FT epsilon, const FT cluster_epsilon, const FT normal_threshold, const FT min_points) :
			m_epsilon(epsilon),
			m_cluster_epsilon(cluster_epsilon),
			m_normal_threshold(normal_threshold),
			m_min_points(min_points)
			{ }

			template<class Elements, class Point_map, class Normal_map, class Range>
			void detect(const Elements &elements, const Point_map &point_map, const Normal_map &normal_map, std::list<Range> &output) const {
				CGAL_precondition(elements.size() > 1);

				// // Tree.
				// Search_tree tree;
				// create_tree(elements, point_map, tree);
				
				// // Main structures.
				// size_t available_points = elements.size();
				// std::map<Point_identifier, int> shape_index;
				// for (typename Elements::const_iterator element = elements.begin(); element != elements.end(); ++element) shape_index[*element] = -1;
				// std::vector<Point_identifier> index_container;

				// // Main loop.
				// int class_index = -1;
      			// for (typename Elements::const_iterator element = elements.begin(); element != elements.end(); ++element, ++i) {
      				
      			// 	const Point_identifier &point_index = *element;
      			// 	if (shape_index.at(point_index) >= 0) continue;

      			// 	// Get query point and its normal.
      			// 	const Point_2 &query_pos = points.at(point_index);
      			// 	const Normal &query_norm = normals.at(point_index);

      			// 	// Get optimal line and its normal.
      			// 	Line_2 optimal_line; Normal line_normal;
      			// 	get_line_with_normal(query_pos, query_norm, optimal_line, line_normal);

      			// 	// Initialize containers.
				// 	shape_index[point_index] = ++class_index;

      			// 	index_container.clear();
        		// 	index_container.push_back(point_index);

        		// 	// Propagate region through all neighbouring points.
        		// 	propagate(shape_index, index_container, optimal_line, line_normal, internal, 
        		// 		class_index, points, normals, point_index, tree);

        		// 	// Update input data.
      			// 	if (index_container.size() >= m_min_points) {

      			// 		update_input_data(boundaries, boundaries_projected, class_index, index_container, points);
      			// 		available_points -= index_container.size();

      			// 	} else {

      			// 		--class_index;
          		// 		shape_index[point_index] = -1;

          		// 		for (Index_iterator it = index_container.begin(); it != index_container.end(); ++it)
            	// 			shape_index[*it] = -1;
      			// 	}
      			// }

      			// // Update all unassigned points.
      			// update_unassigned_points(points, shape_index);
      			// assert(points.size() == available_points);

      			// // Save internal log.
				// // if (m_save_info) internal.save("tmp" + std::string(PSR) + "internal_rg");

				// // Save found shape indices.
				// if (m_save_info) {

				// 	log.out << "" + std::string(PN) + "Shape indices (" << shape_index.size() << ") :";
				// 	for (size_t i = 0; i < shape_index.size(); ++i) log.out << shape_index[i] << " ";
				// 	log.out << std::endl;
				// }
			}

		private:
			const FT m_epsilon;
			const FT m_cluster_epsilon;
			const FT m_normal_threshold;
			const FT m_min_points;

			// void create_tree(const Elements &elements, const Point_map &point_map, Search_tree &tree) const {
                
            //     tree.clear();
            //     for (typename Elements::const_iterator element = elements.begin(); element != elements.end(); ++element) {
                    
            //         const Point_2 &point = get(point_map, *element);
            //         tree.insert(point);
            //     }
            // }
		};

	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_POINTS_BASED_REGION_GROWING_2_H