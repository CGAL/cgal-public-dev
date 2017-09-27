#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_ROOF_FITTER_2_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_ROOF_FITTER_2_H

// STL includes.
#include <map>
#include <vector>
#include <cassert>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class CDTInput, class ContainerInput>
		class Level_of_detail_building_roof_fitter_2 {

		public:
			typedef KernelTraits   Kernel;
			typedef CDTInput       CDT;
			typedef ContainerInput Container;

			typedef typename CDT::Face_handle Face_handle;

			// Extra.
			using Buildings = std::map<int, std::vector<Face_handle> >;

			Level_of_detail_building_roof_fitter_2() { }

			void fit_roof_heights(const CDT &, const Buildings &, const Container &) const {
				
				// to be implemented
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_ROOF_FITTER_2_H