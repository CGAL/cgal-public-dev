#ifndef CGAL_LEVEL_OF_DETAIL_SEMANTIC_DATA_SPLITTER_H
#define CGAL_LEVEL_OF_DETAIL_SEMANTIC_DATA_SPLITTER_H

// LOD includes.
#include <CGAL/Level_of_detail/Enumerations.h>

namespace CGAL {

	namespace Level_of_detail {

		template<class LodDataStructure>
		class Semantic_data_splitter {

		public:
			using LOD_data_structure   = LodDataStructure;
			using Input_range 		   = typename LodDataStructure::Input_range;
			using Input_range_iterator = typename Input_range::const_iterator;

			template<class Semantic_element_map>
			void split_semantics(LOD_data_structure &lod_data_structure, const Semantic_element_map &semantic_element_map) const {

				lod_data_structure.ground_points().clear();
				lod_data_structure.building_boundary_points().clear();
				lod_data_structure.building_interior_points().clear();

				const Input_range &input_range = lod_data_structure.input_range();
				for (Input_range_iterator point = input_range.begin(); point != input_range.end(); ++point) {
					
					const Semantic_label label = get(semantic_element_map, *point);
					switch (label) {

						case Semantic_label::GROUND:
							lod_data_structure.ground_points().push_back(point);
							break;

						case Semantic_label::BUILDING_BOUNDARY:
							lod_data_structure.building_boundary_points().push_back(point);
							break;

						case Semantic_label::BUILDING_INTERIOR:
							lod_data_structure.building_interior_points().push_back(point);
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