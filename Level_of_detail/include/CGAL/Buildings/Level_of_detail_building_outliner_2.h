#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_OUTLINER_2_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_OUTLINER_2_H

// STL includes.
#include <map>
#include <vector>
#include <cassert>

// New CGAL includes.
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class CDTInput>
		class Level_of_detail_building_outliner_2 {

		public:
			typedef KernelTraits Kernel;
			typedef CDTInput     CDT;

			typedef typename Kernel::FT FT;

			typedef typename CDT::Vertex_handle Vertex_handle;
			typedef typename CDT::Face_handle   Face_handle;

			// Extra.
			using Building  = CGAL::LOD::Building<FT, Vertex_handle, Face_handle>;
			using Buildings = std::map<int, Building>;

			Level_of_detail_building_outliner_2() { }

			void find_boundaries(const CDT &, Buildings &) const {
				
				// to be implemented
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_OUTLINER_2_H