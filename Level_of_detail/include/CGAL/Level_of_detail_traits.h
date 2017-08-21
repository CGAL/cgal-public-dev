#ifndef CGAL_LEVEL_OF_DETAIL_TRAITS_H
#define CGAL_LEVEL_OF_DETAIL_TRAITS_H

// New CGAL includes.
#include <CGAL/Loader/Level_of_detail_loader_stub.h>
#include <CGAL/Preprocessor/Level_of_detail_preprocessor.h>
#include <CGAL/Selector/Level_of_detail_selector.h>
#include <CGAL/Selector/Level_of_detail_selection_strategy.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class OutputContainer>
		struct Level_of_detail_traits {

			typedef KernelTraits 							       Kernel;
			typedef OutputContainer                                Container;
			typedef Level_of_detail_loader_stub<Kernel, Container> Loader;
			typedef Level_of_detail_preprocessor<Kernel>           Preprocessor;

			typedef Level_of_detail_clutter<Kernel, Container> 			 ClutterStrategy;
			typedef Level_of_detail_ground<Kernel, Container> 			 GroundStrategy;
			typedef Level_of_detail_building_boundary<Kernel, Container> BuildingBoundaryStrategy;
			typedef Level_of_detail_building_interior<Kernel, Container> BuildingInteriorStrategy;

			typedef Level_of_detail_selector<Kernel, ClutterStrategy> 		   ClutterSelector;
			typedef Level_of_detail_selector<Kernel, GroundStrategy> 		   GroundSelector;
			typedef Level_of_detail_selector<Kernel, BuildingBoundaryStrategy> BuildingBoundarySelector;
			typedef Level_of_detail_selector<Kernel, BuildingInteriorStrategy> BuildingInteriorSelector;
			
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_TRAITS_H