#ifndef CGAL_LEVEL_OF_DETAIL_RECONSTRUCTION_0_H
#define CGAL_LEVEL_OF_DETAIL_RECONSTRUCTION_0_H

// STL includes.
#include <vector>
#include <iostream>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class CDTInput, class Visibility>
		class Level_of_detail_reconstruction_0 {

		public:
			typedef KernelTraits Kernel;
			typedef CDTInput     CDT;
			typedef Visibility   Visibility_result;

			typedef typename Kernel::Segment_2 Segment;
			typedef std::vector<Segment> Lod_0_result;

			void reconstruct(const CDT &, const Visibility_result &, Lod_0_result &) const {

				// to be implemented
				// std::cout << "lod_0" << std::endl;
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_RECONSTRUCTION_0_H