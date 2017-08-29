#ifndef CGAL_LEVEL_OF_DETAIL_TRAITS_H
#define CGAL_LEVEL_OF_DETAIL_TRAITS_H

// STL includes.
#include <map>
#include <vector>

// New CGAL includes.
#include <CGAL/Loader/Level_of_detail_loader_stub.h>
#include <CGAL/Preprocessor/Level_of_detail_preprocessor.h>
#include <CGAL/Selector/Level_of_detail_selector.h>
#include <CGAL/Selector/Level_of_detail_selection_strategy.h>
#include <CGAL/Regularizer/Level_of_detail_regularizer.h>
#include <CGAL/Projector/Level_of_detail_projector.h>
#include <CGAL/Structuring_2/Level_of_detail_structuring_2.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class OutputContainer>
		struct Level_of_detail_traits {

			typedef KernelTraits 							        Kernel;
			typedef OutputContainer                                 Container;
			typedef Level_of_detail_loader_stub<Kernel, Container>  Loader;
			typedef Level_of_detail_preprocessor<Kernel, Container> Preprocessor;

			typedef Level_of_detail_clutter<Kernel, Container> 			 Clutter_strategy;
			typedef Level_of_detail_ground<Kernel, Container> 			 Ground_strategy;
			typedef Level_of_detail_building_boundary<Kernel, Container> Building_boundary_strategy;
			typedef Level_of_detail_building_interior<Kernel, Container> Building_interior_strategy;

			typedef Level_of_detail_selector<Kernel, Clutter_strategy> 		     Clutter_selector;
			typedef Level_of_detail_selector<Kernel, Ground_strategy> 		     Ground_selector;
			typedef Level_of_detail_selector<Kernel, Building_boundary_strategy> Building_boundary_selector;
			typedef Level_of_detail_selector<Kernel, Building_interior_strategy> Building_interior_selector;

			typedef std::map<int, std::vector<int> >          						Planes;
			typedef Level_of_detail_vertical_regularizer<Kernel, Container, Planes> Vertical_regularizer;
			
			typedef std::map<int, typename Kernel::Point_2>           	   				   Projected;
			typedef Level_of_detail_simple_projector<Kernel, Container, Planes, Projected> Ground_projector;

			typedef Level_of_detail_structuring_2<Kernel> Structuring_2;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_TRAITS_H