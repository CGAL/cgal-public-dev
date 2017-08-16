#ifndef CGAL_LEVEL_OF_DETAIL_LOADER_STUB_H
#define CGAL_LEVEL_OF_DETAIL_LOADER_STUB_H

#include <CGAL/Loader/Level_of_detail_loader.h>

namespace CGAL {

	namespace LOD {

		template<class Traits, class PointSet_3>
		class Level_of_detail_loader_stub : public Level_of_detail_loader<Traits, PointSet_3> {
		
		public:
			typedef Traits 					 Kernel;
			typedef PointSet_3               Container;
			typedef typename Kernel::Point_3 Point;

			void load(const std::string &, Container &container) const override {

				// Her I assume that we use Point_set_3 package.
				container.insert(Point(0, 0, 0));
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_LOADER_H