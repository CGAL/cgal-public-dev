#ifndef CGAL_LEVEL_OF_DETAIL_RECONSTRUCTION_H
#define CGAL_LEVEL_OF_DETAIL_RECONSTRUCTION_H

// STL includes.
#include <iostream>

namespace CGAL {

	namespace Level_of_detail {

		template<class InputTraits>
		class Level_of_detail_reconstruction {

		public:
			Level_of_detail_reconstruction() { }

			void build_lod0() const {

				std::cout << "building lod0" << std::endl;
			}

			void build_lod1() const {

				std::cout << "building lod1" << std::endl;
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_RECONSTRUCTION_H