#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_OUTLINER_2_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_OUTLINER_2_H

// STL includes.
#include <map>
#include <vector>
#include <cassert>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class CDTInput>
		class Level_of_detail_building_outliner_2 {

		public:
			typedef KernelTraits Kernel;
			typedef CDTInput     CDT;

			typedef typename CDT::Face_handle Face_handle;

			// Extra.
			using Buildings = std::map<int, std::vector<Face_handle> >;

			Level_of_detail_building_outliner_2() { }

			void order_walls(const CDT &, const Buildings &) const {
				
				// to be implemented
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_OUTLINER_2_H