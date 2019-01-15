#ifndef CGAL_LEVELS_OF_DETAIL_ENUMERATIONS_H
#define CGAL_LEVELS_OF_DETAIL_ENUMERATIONS_H

namespace CGAL {

  namespace Levels_of_detail {

    enum class Semantic_label { 
			
			UNASSIGNED = 0,
			GROUND = 1,
			BUILDING_BOUNDARY = 2,
      BUILDING_INTERIOR = 3, 
      VEGETATION = 4
		};

  } // Levels_of_detail

} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_ENUMERATIONS_H
