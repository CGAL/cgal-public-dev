#ifndef CGAL_LEVEL_OF_DETAIL_SEMANTIC_DATA_SPLITTER_H
#define CGAL_LEVEL_OF_DETAIL_SEMANTIC_DATA_SPLITTER_H

// LOD includes.
#include <CGAL/Level_of_detail/Enumerations.h>

namespace CGAL {

	namespace Level_of_detail {

		class Semantic_data_splitter {

		public:
			template<class Input_range, class Semantic_element_map, class Points>
			void split_semantics(const Input_range &input_range, const Semantic_element_map &semantic_element_map, 
			Points &ground_points, Points &building_boundary_points, Points &building_interior_points) const {
				
				using Const_points_iterator = typename Input_range::const_iterator;

				ground_points.clear();
				building_boundary_points.clear();
				building_interior_points.clear();

				for (Const_points_iterator point = input_range.begin(); point != input_range.end(); ++point) {
					
					const Semantic_label label = get(semantic_element_map, *point);
					switch (label) {

						case Semantic_label::GROUND:
							ground_points.push_back(point);
							break;

						case Semantic_label::BUILDING_BOUNDARY:
							building_boundary_points.push_back(point);
							break;

						case Semantic_label::BUILDING_INTERIOR:
							building_interior_points.push_back(point);
							break;

						default:
							break;
					}
				}
			}
		};
	
	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_SEMANTIC_DATA_SPLITTER_H