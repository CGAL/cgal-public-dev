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

	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_ENUMERATIONS_H