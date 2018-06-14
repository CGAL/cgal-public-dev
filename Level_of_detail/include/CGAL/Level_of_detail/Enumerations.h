#ifndef CGAL_LEVEL_OF_DETAIL_ENUMERATIONS_H
#define CGAL_LEVEL_OF_DETAIL_ENUMERATIONS_H

namespace CGAL {

	namespace Level_of_detail {

		enum class Semantic_label { 
			
			UNASSIGNED = 0,
			GROUND = 1, 
			BUILDING_INTERIOR = 2, 
			BUILDING_BOUNDARY = 3
		};

		enum class Colour_map_type {
			RANDOM = 0,
			WHITE = 1,
			BLACK = 2,
			GROUND_DEFAULT = 3,
			WALL_DEFAULT = 4,
			ROOF_DEFAULT = 5,
		};

		enum class Visibility_label {
			OUTSIDE = 0,
			INSIDE = 1
		};

		enum class Flat_roof_type {
			AVERAGE = 0,
			FIXED_10_METERS = 1
		};

	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_ENUMERATIONS_H